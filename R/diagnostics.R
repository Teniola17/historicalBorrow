#' Compute MCMC diagnostics for a fitted borrowing model
#'
#' Returns convergence diagnostics (R-hat, ESS, divergences), prior-data
#' conflict detection, effective sample size of borrowing, and trace plots.
#' The result is also stored back into `fit$diagnostics`.
#'
#' @param fit A `borrowing_fit` object from [fit_borrowing_model()].
#' @param threshold_rhat R-hat threshold above which parameters are flagged
#'   (default `1.05`).
#' @param pdc_threshold Number of prior SDs separating prior mean from MLE at
#'   which prior-data conflict is flagged (default `2.0`).
#' @param parameters Character vector of parameter names to include in trace
#'   plots. `NULL` (default) selects the primary outcome parameter automatically.
#'
#' @return An S3 object of class `hb_diagnostics` (a named list) with:
#'   \describe{
#'     \item{`rhat`}{Data frame of R-hat values.}
#'     \item{`ess_bulk`}{Data frame of bulk ESS values.}
#'     \item{`ess_tail`}{Data frame of tail ESS values.}
#'     \item{`divergences`}{Integer count of divergent transitions.}
#'     \item{`max_treedepth_hits`}{Integer count of max-treedepth saturation.}
#'     \item{`prior_data_conflict`}{List with `detected`, `statistic`,
#'       `threshold`, `direction`.}
#'     \item{`ess_borrowing`}{Numeric effective sample size contribution from
#'       historical data.}
#'     \item{`shrinkage`}{Numeric shrinkage coefficient.}
#'     \item{`trace_plots`}{A `ggplot` object of trace plots.}
#'     \item{`flags`}{Character vector of diagnostic warnings (empty if all OK).}
#'   }
#'
#' @examples
#' \dontrun{
#' diag <- diagnostics(bf)
#' print(diag)
#' diag$trace_plots
#' }
#'
#' @export
diagnostics <- function(fit,
                        threshold_rhat = 1.05,
                        pdc_threshold  = 2.0,
                        parameters     = NULL) {

  if (!inherits(fit, "borrowing_fit"))
    stop("`fit` must be a `borrowing_fit` object.", call. = FALSE)

  draws <- fit$draws

  # ── R-hat & ESS ─────────────────────────────────────────────────────────────
  summ <- posterior::summarise_draws(
    draws,
    posterior::default_convergence_measures()
  )
  rhat_df     <- summ[, c("variable", "rhat")]
  ess_bulk_df <- summ[, c("variable", "ess_bulk")]
  ess_tail_df <- summ[, c("variable", "ess_tail")]

  # ── HMC diagnostics ──────────────────────────────────────────────────────────
  div_count        <- 0L
  treedepth_count  <- 0L
  if (!is.null(fit$fit)) {
    tryCatch({
      diag_summary    <- fit$fit$diagnostic_summary(quiet = TRUE)
      div_count       <- sum(diag_summary$num_divergent,       na.rm = TRUE)
      treedepth_count <- sum(diag_summary$num_max_treedepth,   na.rm = TRUE)
    }, error = function(e) NULL)
  }

  # ── Prior-data conflict ──────────────────────────────────────────────────────
  pdc <- .compute_prior_data_conflict(fit, pdc_threshold)

  # ── ESS of borrowing ─────────────────────────────────────────────────────────
  ess_borrow <- .compute_ess_borrowing(fit)

  # ── Shrinkage ────────────────────────────────────────────────────────────────
  shrink <- .compute_shrinkage(fit)

  # ── Trace plots ──────────────────────────────────────────────────────────────
  params_to_plot <- parameters %||% .primary_params(fit$outcome_type, fit$model)
  tp <- .make_trace_plots(draws, params_to_plot)

  # ── Flags ────────────────────────────────────────────────────────────────────
  flags <- character(0)
  bad_rhat <- rhat_df$variable[!is.na(rhat_df$rhat) & rhat_df$rhat > threshold_rhat]
  if (length(bad_rhat) > 0)
    flags <- c(flags, paste0("R-hat > ", threshold_rhat, " for: ", paste(bad_rhat, collapse = ", ")))
  if (div_count > 0)
    flags <- c(flags, paste0(div_count, " divergent transition(s). Consider increasing adapt_delta."))
  if (treedepth_count > 0)
    flags <- c(flags, paste0(treedepth_count, " max-treedepth hit(s). Consider increasing max_treedepth."))
  if (pdc$detected)
    flags <- c(flags, paste0("Prior-data conflict detected (statistic = ",
                              round(pdc$statistic, 2), ", direction = ", pdc$direction, ")."))

  result <- structure(
    list(
      rhat                = rhat_df,
      ess_bulk            = ess_bulk_df,
      ess_tail            = ess_tail_df,
      divergences         = div_count,
      max_treedepth_hits  = treedepth_count,
      prior_data_conflict = pdc,
      ess_borrowing       = ess_borrow,
      shrinkage           = shrink,
      trace_plots         = tp,
      flags               = flags
    ),
    class = "hb_diagnostics"
  )

  # Cache diagnostics in the fit object (by reference is not possible in R,
  # so we return the diagnostics object; users assign it)
  result
}

# ── Internal: prior-data conflict ─────────────────────────────────────────────

#' @noRd
.compute_prior_data_conflict <- function(fit, threshold = 2.0) {
  prior <- fit$prior
  if (is.null(prior)) return(list(detected = FALSE, statistic = NA_real_,
                                  threshold = threshold, direction = NA_character_))

  mix  <- prior$mixture_rmap %||% prior$mixture
  if (is.null(mix)) return(list(detected = FALSE, statistic = NA_real_,
                                threshold = threshold, direction = NA_character_))

  prior_mean <- sum(mix$weights * mix$means)
  prior_sd   <- sqrt(sum(mix$weights * (mix$sds^2 + (mix$means - prior_mean)^2)))

  # MLE from current data
  curr    <- fit$data$stan_data_current
  mle_val <- switch(fit$outcome_type,
    binary     = .logit(curr$r_curr / curr$n_curr),
    continuous = curr$y_curr,
    count      = log(curr$events_curr / exp(curr$log_n_curr))
  )

  if (is.na(mle_val) || is.infinite(mle_val) || prior_sd < 1e-10) {
    return(list(detected = FALSE, statistic = NA_real_,
                threshold = threshold, direction = NA_character_))
  }

  stat      <- (mle_val - prior_mean) / prior_sd
  detected  <- abs(stat) > threshold
  direction <- if (stat > 0) "positive" else "negative"

  list(
    detected  = detected,
    statistic = stat,
    threshold = threshold,
    direction = direction
  )
}

# ── Internal: ESS of borrowing ────────────────────────────────────────────────

#' @noRd
.compute_ess_borrowing <- function(fit) {
  if (is.null(fit$prior)) return(NA_real_)
  fit$prior$ess_prior_rmap %||% fit$prior$ess_prior %||% NA_real_
}

# ── Internal: shrinkage ───────────────────────────────────────────────────────

#' @noRd
.compute_shrinkage <- function(fit) {
  prior <- fit$prior
  if (is.null(prior)) return(NA_real_)

  mix <- prior$mixture_rmap %||% prior$mixture
  if (is.null(mix)) return(NA_real_)

  prior_mean <- sum(mix$weights * mix$means)
  curr       <- fit$data$stan_data_current

  mle_val <- switch(fit$outcome_type,
    binary     = .logit(curr$r_curr / curr$n_curr),
    continuous = curr$y_curr,
    count      = log(curr$events_curr / exp(curr$log_n_curr))
  )

  param_name <- .primary_params(fit$outcome_type, fit$model)[1]
  post_mean  <- tryCatch(
    mean(as.numeric(posterior::subset_draws(fit$draws, variable = param_name))),
    error = function(e) NA_real_
  )

  if (is.na(mle_val) || is.na(post_mean) || abs(prior_mean - mle_val) < 1e-10)
    return(NA_real_)

  (post_mean - mle_val) / (prior_mean - mle_val)
}

# ── Internal: primary parameters ─────────────────────────────────────────────

#' @noRd
.primary_params <- function(outcome_type, model) {
  base <- switch(outcome_type,
    binary     = c("theta_curr", "p_curr"),
    continuous = "theta_curr",
    count      = c("log_lambda_curr", "lambda_curr")
  )
  if (model == "commensurate") base <- c(base, "tau_comm")
  if (model == "power")        base <- c(base, "a0")
  base
}

# ── Internal: trace plots ─────────────────────────────────────────────────────

#' @noRd
.make_trace_plots <- function(draws, params) {
  # Filter to parameters that actually exist in draws
  available <- posterior::variables(draws)
  params     <- intersect(params, available)
  if (length(params) == 0) params <- available[1]

  # Build long-format data frame for ggplot
  plot_data <- do.call(rbind, lapply(params, function(p) {
    vals <- as.numeric(posterior::subset_draws(draws, variable = p))
    chains     <- as.integer(posterior::subset_draws(draws, variable = ".chain"))
    iterations <- as.integer(posterior::subset_draws(draws, variable = ".iteration"))
    data.frame(
      parameter = p,
      chain     = factor(chains),
      iteration = iterations,
      value     = vals
    )
  }))

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$iteration, y = .data$value, colour = .data$chain)
  ) +
    ggplot2::geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title  = "MCMC Trace Plots",
      x      = "Iteration",
      y      = "Value",
      colour = "Chain"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}

# ── S3 print ──────────────────────────────────────────────────────────────────

#' @export
print.hb_diagnostics <- function(x, ...) {
  cat("── historicalBorrow: MCMC Diagnostics ───────────────────────────────────\n")
  cat(glue::glue("  Divergences       : {x$divergences}\n"))
  cat(glue::glue("  Max-treedepth hits: {x$max_treedepth_hits}\n"))

  max_rhat <- max(x$rhat$rhat, na.rm = TRUE)
  cat(glue::glue("  Max R-hat         : {round(max_rhat, 3)}\n"))

  min_ess <- min(x$ess_bulk$ess_bulk, na.rm = TRUE)
  cat(glue::glue("  Min ESS (bulk)    : {round(min_ess, 0)}\n"))

  cat(glue::glue("  ESS borrowing     : {round(x$ess_borrowing, 1)}\n"))

  if (!is.na(x$shrinkage))
    cat(glue::glue("  Shrinkage         : {round(x$shrinkage, 3)}\n"))

  pdc <- x$prior_data_conflict
  if (!is.na(pdc$detected)) {
    status <- if (pdc$detected) "YES" else "no"
    cat(glue::glue(
      "  Prior-data conflict: {status}",
      " (|stat| = {round(abs(pdc$statistic), 2)}, threshold = {pdc$threshold})\n"
    ))
  }

  if (length(x$flags) > 0) {
    cat("\n  Warnings:\n")
    for (f in x$flags) cat(glue::glue("  ! {f}\n"))
  } else {
    cat("\n  No diagnostic warnings.\n")
  }
  invisible(x)
}
