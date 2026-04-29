#' @importFrom rlang %||%
NULL

# ── Logit helpers ──────────────────────────────────────────────────────────────

.logit <- function(p) {
  p <- pmax(pmin(p, 1 - .Machine$double.eps), .Machine$double.eps)
  log(p / (1 - p))
}

.inv_logit <- function(x) 1 / (1 + exp(-x))

# ── Mixture of normals ─────────────────────────────────────────────────────────

#' Fit a mixture of normal distributions to a vector of draws
#'
#' @param draws Numeric vector of posterior draws.
#' @param K_max Maximum number of mixture components to try (selected by BIC).
#' @return Named list with `weights`, `means`, `sds`, `K`.
#' @noRd
.fit_mixture_of_normals <- function(draws, K_max = 6L) {
  K_max <- min(K_max, length(unique(round(draws, 4))))
  K_max <- max(K_max, 1L)

  if (K_max == 1L) {
    return(list(
      weights = 1.0,
      means   = mean(draws),
      sds     = max(sd(draws), 1e-6),
      K       = 1L
    ))
  }

  best_bic <- Inf
  best_fit <- NULL
  best_K   <- 1L

  for (k in seq_len(K_max)) {
    fit <- tryCatch(
      mclust::Mclust(draws, G = k, modelNames = "V", verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    if (fit$bic < best_bic || is.null(best_fit)) {
      # mclust uses -2*logLik + k*log(n), lower is better — but mclust
      # reports BIC as 2*logLik - k*log(n) (higher = better).  We maximise.
      # Re-use the maximisation direction mclust uses:
      if (is.null(best_fit) || fit$bic > best_bic) {
        best_bic <- fit$bic
        best_fit <- fit
        best_K   <- k
      }
    }
  }

  if (is.null(best_fit)) {
    return(list(
      weights = 1.0,
      means   = mean(draws),
      sds     = max(sd(draws), 1e-6),
      K       = 1L
    ))
  }

  params  <- best_fit$parameters
  weights <- as.numeric(params$pro)
  means   <- as.numeric(params$mean)
  sds     <- as.numeric(sqrt(params$variance$sigmasq))

  # Ensure lengths are consistent (mclust may drop components)
  K_actual <- length(weights)
  if (length(means) < K_actual) means <- rep(means, K_actual)[seq_len(K_actual)]
  if (length(sds)   < K_actual) sds   <- rep(sds,   K_actual)[seq_len(K_actual)]

  list(
    weights = weights / sum(weights),  # re-normalise for safety
    means   = means,
    sds     = pmax(sds, 1e-6),
    K       = K_actual
  )
}

#' Evaluate mixture of normals PDF at points x
#' @noRd
.eval_mixture_density <- function(x, mix) {
  dens <- numeric(length(x))
  for (k in seq_len(mix$K)) {
    dens <- dens + mix$weights[k] * dnorm(x, mix$means[k], mix$sds[k])
  }
  dens
}

#' Draw samples from a mixture of normals
#' @noRd
.sample_mixture <- function(mix, n) {
  comp <- sample(seq_len(mix$K), size = n, replace = TRUE, prob = mix$weights)
  rnorm(n, mean = mix$means[comp], sd = mix$sds[comp])
}

# ── Stan compilation and caching ──────────────────────────────────────────────

#' Check that CmdStan is available
#' @noRd
.check_cmdstan <- function() {
  tryCatch(
    cmdstanr::cmdstan_version(),
    error = function(e) {
      stop(
        "CmdStan is not installed or not found. ",
        "Please run `cmdstanr::install_cmdstan()` before fitting models.",
        call. = FALSE
      )
    }
  )
}

#' Hash a Stan code string (for cache keying)
#' @noRd
.hash_stan_code <- function(code) {
  digest::digest(code, algo = "sha256", serialize = FALSE)
}

#' Return path to the model cache directory
#' @noRd
.cache_dir <- function() {
  cache <- tools::R_user_dir("historicalBorrow", which = "cache")
  if (!dir.exists(cache)) dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  cache
}

#' Compile a Stan model, using the file-based cache if available
#'
#' @param code Character string of Stan code.
#' @param model_name Short name used in the cached filename.
#' @return A `CmdStanModel` object.
#' @noRd
.compile_stan <- function(code, model_name = "model") {
  .check_cmdstan()

  hash       <- substr(.hash_stan_code(code), 1, 12)
  stan_file  <- file.path(.cache_dir(), paste0(model_name, "_", hash, ".stan"))

  if (!file.exists(stan_file)) {
    writeLines(code, stan_file)
  }

  cmdstanr::cmdstan_model(stan_file, compile = TRUE)
}

# ── ESS helpers ───────────────────────────────────────────────────────────────

#' Morita ESS for a binary MAP prior (approximation)
#'
#' Uses the Fisher information of the binomial at the MAP mean probability.
#' Reference: Morita et al. (2008) Biometrics.
#' @noRd
.ess_binary_map <- function(mix, n_ref = 1L) {
  p_mean <- .inv_logit(sum(mix$weights * mix$means))
  fisher_prior <- 1 / sum(mix$weights * mix$sds^2)  # 1 / Var(logit(p))
  fisher_binom <- p_mean * (1 - p_mean)              # Var(X/n) * n = p*(1-p)
  ess <- fisher_prior / fisher_binom
  round(max(ess, 0), 2)
}

#' Variance-based ESS for continuous MAP prior
#' @noRd
.ess_continuous_map <- function(mix, sigma_ref = 1.0) {
  prior_var <- sum(mix$weights * mix$sds^2 +
                     mix$weights * (mix$means - sum(mix$weights * mix$means))^2)
  ess <- sigma_ref^2 / prior_var
  round(max(ess, 0), 2)
}

#' ESS for count MAP prior (log-rate scale)
#' @noRd
.ess_count_map <- function(mix) {
  lambda_mean <- exp(sum(mix$weights * mix$means))
  prior_var   <- sum(mix$weights * mix$sds^2)
  ess <- 1 / (prior_var * lambda_mean)
  round(max(ess, 0), 2)
}

#' Compute ESS appropriate for outcome type
#' @noRd
.compute_ess <- function(mix, outcome_type) {
  switch(outcome_type,
    binary     = .ess_binary_map(mix),
    continuous = .ess_continuous_map(mix),
    count      = .ess_count_map(mix),
    NA_real_
  )
}
