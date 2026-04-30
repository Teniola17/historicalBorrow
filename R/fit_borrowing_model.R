#' Fit a Bayesian borrowing model using a MAP, RMAP, power, or commensurate prior
#'
#' Generates Stan code from the specified model type, compiles it (using a
#' file-based cache to avoid recompilation on repeated calls), runs MCMC via
#' `cmdstanr`, and returns a fitted model object.
#'
#' @param current_data A `bb_data` object from [prepare_data()].
#' @param prior A `map_prior` or `rmap_prior` object from [fit_map()] /
#'   [robustify_map()]. Required for `model = "map"` or `"rmap"`.
#'   Ignored (and may be `NULL`) for `model = "power"` or `"commensurate"`.
#' @param model One of:
#'   \describe{
#'     \item{`"map"`}{Uses the MAP mixture prior from `prior$mixture`.}
#'     \item{`"rmap"`}{Uses the robust MAP mixture from `prior$mixture_rmap`.}
#'     \item{`"tdp"`}{Uses the TDP mixture prior from `prior$mixture`.}
#'     \item{`"rtdp"`}{Uses the robust TDP mixture from `prior$mixture_rtdp`.}
#'     \item{`"power"`}{Power prior — scales historical log-likelihood by `a0`.}
#'     \item{`"commensurate"`}{Commensurate prior — `tau_comm` controls borrowing.}
#'   }
#' @param hyperpriors Named list of model hyperpriors. Relevant keys:
#'   * MAP/RMAP: none (mixture passed as Stan data).
#'   * Power prior: `a0_fixed` (scalar `[0,1]`, or `-1` to estimate),
#'     `a0_alpha`, `a0_beta`.
#'   * Commensurate: `mu_hist`, `sig_hist`.
#' @param chains Number of MCMC chains (default `4`).
#' @param iter_warmup Warm-up iterations per chain (default `1000`).
#' @param iter_sampling Sampling iterations per chain (default `2000`).
#' @param seed Integer random seed.
#' @param show_messages Logical; show Stan compilation/sampling messages.
#' @param ... Additional arguments forwarded to `cmdstanr`'s `$sample()`.
#'
#' @return An S3 object of class `borrowing_fit`:
#'   \describe{
#'     \item{`fit`}{The `CmdStanFit` object.}
#'     \item{`draws`}{Posterior draws as a `draws_df`.}
#'     \item{`prior`}{The prior object supplied.}
#'     \item{`data`}{The `bb_data` object supplied.}
#'     \item{`model`}{Model type string.}
#'     \item{`outcome_type`}{Outcome type string.}
#'     \item{`stan_code`}{Generated Stan code string.}
#'     \item{`diagnostics`}{`NULL` until [diagnostics()] is called.}
#'   }
#'
#' @seealso [fit_map()], [robustify_map()], [diagnostics()],
#'   [plot_prior_posterior()], [decision_criteria()]
#'
#' @examples
#' \dontrun{
#' pd  <- prepare_data(hist_data, curr_data, "binary")
#' mp  <- fit_map(pd)
#' rmp <- robustify_map(mp, weight = 0.8)
#' bf  <- fit_borrowing_model(pd, prior = rmp, model = "rmap")
#' summary(bf)
#' }
#'
#' @export
fit_borrowing_model <- function(current_data,
                                prior         = NULL,
                                model         = c("map", "rmap", "power", "commensurate"),
                                hyperpriors   = list(),
                                chains        = 4L,
                                iter_warmup   = 1000L,
                                iter_sampling = 2000L,
                                seed          = 42L,
                                show_messages = FALSE,
                                ...) {

  if (!inherits(current_data, "bb_data"))
    stop("`current_data` must be a `bb_data` object.", call. = FALSE)

  model <- match.arg(model)
  .check_cmdstan()

  model_structure <- switch(model,
    map          = "borrowing_map",
    rmap         = "borrowing_rmap",
    tdp          = "borrowing_tdp",
    rtdp         = "borrowing_rtdp",
    power        = "borrowing_power",
    commensurate = "borrowing_commensurate"
  )

  # ── Build Stan data ──────────────────────────────────────────────────────────
  stan_d <- .build_stan_data(current_data, prior, model, hyperpriors)

  # ── Generate & compile Stan code ────────────────────────────────────────────
  hp_for_gen <- .build_hp_for_gen(current_data, prior, model, hyperpriors)
  code <- .generate_stan_code(
    outcome_type    = current_data$outcome,
    model_structure = model_structure,
    prior_spec      = hp_for_gen
  )

  compiled <- .compile_stan(
    code,
    model_name = paste0(model, "_", current_data$outcome)
  )

  # ── Sample ──────────────────────────────────────────────────────────────────
  fit <- compiled$sample(
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

  structure(
    list(
      fit          = fit,
      draws        = draws,
      prior        = prior,
      data         = current_data,
      model        = model,
      outcome_type = current_data$outcome,
      stan_code    = code,
      diagnostics  = NULL
    ),
    class = "borrowing_fit"
  )
}

# ── Stan data builders ────────────────────────────────────────────────────────

.build_stan_data <- function(current_data, prior, model, hyperpriors) {
  curr <- current_data$stan_data_current

  if (model %in% c("map", "rmap", "tdp", "rtdp")) {
    mix <- if (model == "rmap") {
      prior$mixture_rmap
    } else if (model == "rtdp") {
      prior$mixture_rtdp
    } else {
      prior$mixture
    }
    if (is.null(mix))
      stop(
        if (model == "rmap")
          "prior$mixture_rmap is NULL. Did you call robustify_map()?"
        else if (model == "rtdp")
          "prior$mixture_rtdp is NULL. Did you call robustify_tdp()?"
        else if (model %in% c("tdp", "map"))
          "prior$mixture is NULL."
        else
          "prior object missing mixture.",
        call. = FALSE
      )
    c(
      list(
        K           = mix$K,
        mix_weights = mix$weights,
        mix_means   = mix$means,
        mix_sds     = mix$sds
      ),
      curr
    )

  } else if (model == "power") {
    hist       <- current_data$stan_data_map
    a0_fixed   <- hyperpriors$a0_fixed %||% -1
    a0_flag    <- as.integer(a0_fixed >= 0)
    a0_val     <- if (a0_fixed >= 0) a0_fixed else 0.5  # placeholder when estimated
    a0_alpha   <- hyperpriors$a0_alpha %||% 1
    a0_beta    <- hyperpriors$a0_beta  %||% 1

    base <- c(hist, curr)
    # Rename map data keys to _hist suffix for power prior Stan model
    base <- .rename_for_power(base, current_data$outcome)
    c(base, list(
      a0_fixed_flag = a0_flag,
      a0_fixed_val  = a0_val,
      a0_alpha      = a0_alpha,
      a0_beta       = a0_beta
    ))

  } else {  # commensurate
    curr
  }
}

.rename_for_power <- function(data_list, outcome_type) {
  renames <- switch(outcome_type,
    binary = c(J = "J_hist", n = "n_hist", r = "r_hist",
               n_curr = "n_curr", r_curr = "r_curr"),
    continuous = c(J = "J_hist", y = "y_hist", se = "se_hist",
                   n_curr = "n_curr", y_curr = "y_curr", se_curr = "se_curr"),
    count = c(J = "J_hist", events = "events_hist", log_n = "log_n_hist",
              events_curr = "events_curr", log_n_curr = "log_n_curr")
  )
  new_list <- data_list
  for (old in names(renames)) {
    new <- renames[[old]]
    if (old %in% names(new_list)) {
      new_list[[new]] <- new_list[[old]]
      if (old != new) new_list[[old]] <- NULL
    }
  }
  new_list
}

.build_hp_for_gen <- function(current_data, prior, model, hyperpriors) {
  if (model == "power") {
    a0_fixed <- hyperpriors$a0_fixed %||% -1
    list(
      a0_fixed = a0_fixed,
      a0_alpha = hyperpriors$a0_alpha %||% 1,
      a0_beta  = hyperpriors$a0_beta  %||% 1,
      mu_mean  = hyperpriors$mu_mean  %||% 0,
      mu_sd    = hyperpriors$mu_sd    %||% 2.5
    )
  } else if (model == "commensurate") {
    list(
      mu_hist  = hyperpriors$mu_hist  %||% 0,
      sig_hist = hyperpriors$sig_hist %||% 2.5,
      mu_mean  = hyperpriors$mu_mean  %||% 0,
      mu_sd    = hyperpriors$mu_sd    %||% 2.5
    )
  } else {
    list()  # MAP/RMAP: mixture passed as data, no hyperprior injection needed
  }
}
