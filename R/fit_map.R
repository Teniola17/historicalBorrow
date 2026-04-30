#' Fit a Meta-Analytic-Predictive (MAP) prior from historical data
#'
#' Estimates the MAP prior by fitting a Bayesian hierarchical model to the
#' historical studies in `data`. Stan code is generated automatically from the
#' outcome type and hyperprior specifications. The posterior predictive
#' distribution for a new (exchangeable) study is approximated as a
#' K-component mixture of normal distributions, which is then used in
#' [fit_borrowing_model()].
#'
#' @param data A `bb_data` object created by [prepare_data()].
#' @param model_type Currently only `"standard"` (exchangeable hierarchical
#'   model) is supported. Future versions will add `"robust"` directly here.
#' @param mixture_components Maximum number of normal components to fit to the
#'   MAP predictive posterior. `NULL` lets the BIC choose up to 6 components.
#' @param hyperpriors Named list of hyperprior overrides. Supported keys:
#'   `mu_mean`, `mu_sd`, `tau_scale`. See Details.
#' @param chains Number of MCMC chains (default `4`).
#' @param iter_warmup Warm-up iterations per chain (default `1000`).
#' @param iter_sampling Sampling iterations per chain (default `2000`).
#' @param seed Integer random seed for reproducibility.
#' @param show_messages Logical; whether to show Stan compilation messages.
#' @param ... Additional arguments passed to `cmdstanr`'s `$sample()` method.
#'
#' @details
#' ## Hyperpriors
#' Default hyperpriors:
#' * `mu_mean = 0`, `mu_sd = 2.5` (weakly informative on logit/log scale)
#' * `tau_scale = 0.5` (half-normal on between-study SD)
#'
#' For continuous outcomes, wider defaults (`mu_sd = 100`, `tau_scale = 10`)
#' are used automatically.
#'
#' ## MAP transfer
#' After MCMC, the predictive parameter `mu_pred` is extracted from all chains.
#' A mixture of normal distributions is fitted to these draws via the EM
#' algorithm (`mclust`), selecting the number of components by BIC. This
#' mixture is stored in the returned object and passed as Stan data when
#' calling [fit_borrowing_model()].
#'
#' @return An S3 object of class `map_prior` with components:
#'   \describe{
#'     \item{`mixture`}{`list(weights, means, sds, K)` — mixture of normals approximation.}
#'     \item{`mu_pred_draws`}{Full numeric vector of `mu_pred` posterior draws.}
#'     \item{`outcome_type`}{Outcome type string.}
#'     \item{`model_type`}{`"standard"`.}
#'     \item{`hyperpriors`}{List of hyperprior values used.}
#'     \item{`ess_prior`}{Approximate equivalent historical sample size.}
#'     \item{`stan_fit`}{The `CmdStanFit` object.}
#'     \item{`draws`}{Full posterior draws as a `draws_df`.}
#'     \item{`data`}{The `bb_data` input.}
#'     \item{`mixture_components`}{Number of mixture components selected.}
#'   }
#'
#' @seealso [robustify_map()], [fit_borrowing_model()]
#'
#' @examples
#' \dontrun{
#' hist <- data.frame(
#'   study = paste0("S", 1:4),
#'   n = c(100L, 120L, 90L, 110L),
#'   r = c(30L, 40L, 25L, 38L)
#' )
#' curr <- data.frame(study = "C", n = 80L, r = 22L)
#' pd  <- prepare_data(hist, curr, outcome = "binary")
#' mp  <- fit_map(pd, seed = 42)
#' print(mp)
#' plot(mp)
#' }
#'
#' @export
fit_map <- function(data,
                    model_type         = "standard",
                    mixture_components = NULL,
                    hyperpriors        = list(),
                    chains             = 4L,
                    iter_warmup        = 1000L,
                    iter_sampling      = 2000L,
                    seed               = 42L,
                    show_messages      = FALSE,
                    ...) {

  if (!inherits(data, "bb_data"))
    stop("`data` must be a `bb_data` object from prepare_data().", call. = FALSE)

  model_type <- match.arg(model_type, "standard")
  .check_cmdstan()

  # Apply sensible outcome-specific defaults before generating Stan code
  hyperpriors <- .fill_map_hyperpriors(hyperpriors, data$outcome)

  code <- .generate_stan_code(
    outcome_type    = data$outcome,
    model_structure = paste0("map_", model_type),
    prior_spec      = hyperpriors
  )

  model  <- .compile_stan(code, model_name = paste0("map_", data$outcome))
  stan_d <- data$stan_data_map

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

  draws    <- posterior::as_draws_df(fit$draws())
  mu_draws <- posterior::extract_variable(draws, variable = "mu_pred")

  K_max <- mixture_components %||% 6L
  mix   <- .fit_mixture_of_normals(mu_draws, K_max = K_max)
  ess   <- .compute_ess(mix, data$outcome)

  structure(
    list(
      mixture            = mix,
      mu_pred_draws      = mu_draws,
      outcome_type       = data$outcome,
      model_type         = model_type,
      hyperpriors        = hyperpriors,
      ess_prior          = ess,
      stan_fit           = fit,
      draws              = draws,
      data               = data,
      mixture_components = mix$K
    ),
    class = "map_prior"
  )
}

# ── Hyperprior defaults ────────────────────────────────────────────────────────

.fill_map_hyperpriors <- function(hp, outcome_type) {
  defaults <- switch(outcome_type,
    binary = list(mu_mean = 0, mu_sd = 2.5,  tau_scale = 0.5),
    continuous = list(mu_mean = 0, mu_sd = 100, tau_scale = 10),
    count  = list(mu_mean = 0, mu_sd = 2.5,  tau_scale = 0.5)
  )
  modifyList(defaults, hp)
}
