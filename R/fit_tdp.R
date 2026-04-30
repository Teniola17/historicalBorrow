#' Fit a Time-Dependent Prior (TDP) from historical data with temporal decay
#'
#' Implements an Ornstein-Uhlenbeck continuous-time process to model the
#' evolution of treatment effect parameters over time. This allows borrowing
#' from historical studies while automatically down-weighting older data via
#' a correlation decay parameter `rho`.
#'
#' @param data A `bb_data` object from [prepare_data()]. The historical data
#'   must include a `time` column (numeric, e.g., study year, days since epoch).
#' @param hist_times Optional numeric vector of times for historical studies.
#'   If provided, overrides the `time` column in `data$hist_data`.
#' @param current_time Numeric; the time at which to predict (current trial date).
#'   If `NULL`, uses `max(hist_times) + 1`.
#' @param hyperpriors Named list of prior specifications. Keys:
#'   \describe{
#'     \item{`mu_mean`, `mu_sd`}{Mean and SD of the long-term mean `mu`.}
#'     \item{`mu1_mean`, `mu1_sd`}{Mean and SD of the first study's parameter `mu1`.}
#'     \item{`sigma1_mean`, `sigma1_sd`}{Mean and SD of initial SD `sigma1`.}
#'     \item{`tau_mean`, `tau_sd`}{Mean and SD of temporal noise `tau`.}
#'     \item{`rho_alpha`, `rho_beta`}{Beta hyperprior on correlation decay `rho`.}
#'   }
#' @param chains, iter_warmup, iter_sampling, seed MCMC control parameters.
#' @param show_messages Logical; show Stan messages.
#' @param ... Additional arguments to `cmdstanr`'s `$sample()`.
#'
#' @return An S3 object of class `tdp_prior`:
#'   \describe{
#'     \item{`theta_star_draws`}{Posterior draws of `theta_star` (predictive parameter).}
#'     \item{`outcome_type`}{Outcome type string.}
#'     \item{`current_time`}{The prediction time.}
#'     \item{`hyperpriors`}{Hyperprior specifications used.}
#'     \item{`ess_prior`}{Equivalent historical sample size (via Morita method).}
#'     \item{`draws`}{Full posterior draws.}
#'     \item{`stan_fit`}{CmdStanFit object.}
#'     \item{`data`}{The bb_data input.}
#'     \item{`mixture`}{`list(weights, means, sds, K)` — mixture approximation to theta_star.}
#'   }
#'
#' @seealso [robustify_tdp()], [fit_borrowing_model()]
#'
#' @export
fit_tdp <- function(data,
                    hist_times        = NULL,
                    current_time      = NULL,
                    hyperpriors       = list(),
                    chains            = 4L,
                    iter_warmup       = 1000L,
                    iter_sampling     = 2000L,
                    seed              = 42L,
                    show_messages     = FALSE,
                    ...) {

  if (!inherits(data, "bb_data"))
    stop("`data` must be a `bb_data` object from prepare_data().", call. = FALSE)

  .check_cmdstan()

  # ── Extract times ──────────────────────────────────────────────────────────
  times_hist <- .extract_times(data, hist_times)
  t_star     <- current_time %||% (max(times_hist) + 1)

  if (t_star <= max(times_hist))
    warning("current_time is not after the latest historical study. ",
            "Prediction may be extrapolation in reverse.", call. = FALSE)

  # ── Fill hyperpriors ──────────────────────────────────────────────────────
  hp <- .fill_tdp_hyperpriors(hyperpriors, data$outcome)

  # ── Generate Stan code ────────────────────────────────────────────────────
  code <- .generate_stan_code_tdp(
    outcome_type = data$outcome,
    hyperpriors  = hp
  )

  # ── Build Stan data ───────────────────────────────────────────────────────
  stan_d <- .build_stan_data_tdp(data, times_hist, t_star)

  # ── Compile & sample ──────────────────────────────────────────────────────
  model <- .compile_stan(code, model_name = paste0("tdp_", data$outcome))

  fit <- model$sample(
    data            = stan_d,
    chains          = chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    seed            = seed,
    show_messages   = show_messages,
    show_exceptions = show_messages,
    ...
  )

  draws <- posterior::as_draws_df(fit$draws())
  theta_star_draws <- as.numeric(posterior::subset_draws(draws, variable = "theta_star"))

  # ── Fit mixture of normals ────────────────────────────────────────────────
  mix <- .fit_mixture_of_normals(theta_star_draws, K_max = 6L)
  ess <- .compute_ess(mix, data$outcome)

  structure(
    list(
      theta_star_draws = theta_star_draws,
      outcome_type     = data$outcome,
      current_time     = t_star,
      hyperpriors      = hp,
      ess_prior        = ess,
      draws            = draws,
      stan_fit         = fit,
      data             = data,
      mixture          = mix
    ),
    class = "tdp_prior"
  )
}

# ── Helper: extract times ─────────────────────────────────────────────────────

.extract_times <- function(data, hist_times) {
  if (!is.null(hist_times)) {
    checkmate::assert_numeric(hist_times, len = data$n_hist,
      .var.name = "hist_times")
    return(as.numeric(hist_times))
  }

  if ("time" %in% names(data$hist_data)) {
    return(as.numeric(data$hist_data$time))
  }

  stop(
    "No times provided. Either:\n",
    "  1. Add a 'time' column to hist_data in prepare_data(), or\n",
    "  2. Pass hist_times = c(...) to fit_tdp().",
    call. = FALSE
  )
}

# ── Helper: fill TDP hyperpriors ──────────────────────────────────────────────

.fill_tdp_hyperpriors <- function(hp, outcome_type) {
  defaults <- switch(outcome_type,
    binary = list(
      mu_mean      = 0,
      mu_sd        = 2.5,
      mu1_mean     = 0,
      mu1_sd       = 2.5,
      sigma1_mean  = 0,
      sigma1_sd    = 1,
      tau_mean     = 0,
      tau_sd       = 0.5,
      rho_alpha    = 2,
      rho_beta     = 2
    ),
    continuous = list(
      mu_mean      = 0,
      mu_sd        = 100,
      mu1_mean     = 0,
      mu1_sd       = 50,
      sigma1_mean  = 0,
      sigma1_sd    = 10,
      tau_mean     = 0,
      tau_sd       = 10,
      rho_alpha    = 2,
      rho_beta     = 2
    ),
    count = list(
      mu_mean      = 0,
      mu_sd        = 2.5,
      mu1_mean     = 0,
      mu1_sd       = 2.5,
      sigma1_mean  = 0,
      sigma1_sd    = 1,
      tau_mean     = 0,
      tau_sd       = 0.5,
      rho_alpha    = 2,
      rho_beta     = 2
    )
  )
  modifyList(defaults, hp)
}

# ── Helper: build Stan data for TDP ───────────────────────────────────────────

.build_stan_data_tdp <- function(data, times_hist, t_star) {
  c(
    list(
      H      = data$n_hist,
      time   = times_hist,
      t_star = t_star
    ),
    data$stan_data_map
  )
}
