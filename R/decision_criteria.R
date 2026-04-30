#' Compute Go/No-Go decision criteria from posterior draws
#'
#' Evaluates a one-sided posterior probability threshold rule common in Phase II
#' Bayesian clinical trials. Returns a Go decision if
#' `P(theta > theta_threshold | data) >= probability_threshold`
#' (or `<` when `direction = "less"`).
#'
#' @param fit A `borrowing_fit` object.
#' @param theta_threshold Numeric threshold for the primary parameter (on the
#'   **natural** scale: probability for binary, mean for continuous, rate for
#'   count). For example, `theta_threshold = 0.25` for a binary outcome tests
#'   whether the response probability exceeds 25%.
#' @param direction `"greater"` (default) or `"less"`.
#' @param probability_threshold Minimum posterior probability required for a Go
#'   decision (default `0.8`).
#' @param parameter Optional character string naming the draw variable to use.
#'   If `NULL`, the natural-scale outcome parameter is chosen automatically
#'   (`"p_curr"` for binary, `"theta_curr"` for continuous, `"lambda_curr"` for
#'   count).
#'
#' @return An S3 object of class `decision_result`:
#'   \describe{
#'     \item{`go`}{`TRUE` if the posterior probability meets the threshold.}
#'     \item{`probability`}{Posterior probability of the criterion.}
#'     \item{`theta_threshold`}{The value threshold used.}
#'     \item{`probability_threshold`}{The probability threshold used.}
#'     \item{`direction`}{`"greater"` or `"less"`.}
#'     \item{`parameter`}{Name of the parameter used.}
#'     \item{`n_draws`}{Number of posterior draws used.}
#'     \item{`posterior_mean`}{Posterior mean of the parameter.}
#'     \item{`posterior_ci`}{Named vector with `lo` and `hi` (95% CI).}
#'   }
#'
#' @examples
#' \dontrun{
#' dc <- decision_criteria(
#'   fit,
#'   theta_threshold      = 0.25,
#'   direction            = "greater",
#'   probability_threshold = 0.80
#' )
#' print(dc)
#' }
#'
#' @export
decision_criteria <- function(fit,
                               theta_threshold,
                               direction             = c("greater", "less"),
                               probability_threshold = 0.80,
                               parameter             = NULL) {

  if (!inherits(fit, "borrowing_fit"))
    stop("`fit` must be a `borrowing_fit` object.", call. = FALSE)

  checkmate::assert_number(theta_threshold)
  checkmate::assert_number(probability_threshold, lower = 0, upper = 1)
  direction <- match.arg(direction)

  # ── Resolve parameter name ────────────────────────────────────────────────
  param_name <- parameter %||% switch(fit$outcome_type,
    binary     = "p_curr",
    continuous = "theta_curr",
    count      = "lambda_curr"
  )

  draws_vec <- tryCatch(
    posterior::extract_variable(fit$draws, variable = param_name),
    error = function(e) {
      stop(glue::glue(
        "Parameter '{param_name}' not found in posterior draws. ",
        "Use `parameter = ` to specify a valid name."
      ), call. = FALSE)
    }
  )

  n_draws <- length(draws_vec)

  # ── Posterior probability ─────────────────────────────────────────────────
  post_prob <- if (direction == "greater") {
    mean(draws_vec > theta_threshold, na.rm = TRUE)
  } else {
    mean(draws_vec < theta_threshold, na.rm = TRUE)
  }

  go <- post_prob >= probability_threshold

  # ── Summary stats ─────────────────────────────────────────────────────────
  post_mean <- mean(draws_vec, na.rm = TRUE)
  post_ci   <- quantile(draws_vec, c(0.025, 0.975), na.rm = TRUE)
  names(post_ci) <- c("lo", "hi")

  structure(
    list(
      go                    = go,
      probability           = post_prob,
      theta_threshold       = theta_threshold,
      probability_threshold = probability_threshold,
      direction             = direction,
      parameter             = param_name,
      n_draws               = n_draws,
      posterior_mean        = post_mean,
      posterior_ci          = post_ci
    ),
    class = "decision_result"
  )
}

# ── S3 print ──────────────────────────────────────────────────────────────────

#' @export
print.decision_result <- function(x, ...) {
  verdict <- if (x$go) "GO  ✓" else "NO-GO  ✗"
  cat("── historicalBorrow: Decision Criteria ──────────────────────────────────\n")
  cat(glue::glue("  Verdict            : {verdict}\n"))
  cat(glue::glue("  Parameter          : {x$parameter}\n"))
  dir_sym  <- if (x$direction == "greater") ">" else "<"
  cat(glue::glue(
    "  Criterion          : P({x$parameter} {dir_sym} {x$theta_threshold}) >= {x$probability_threshold}\n"
  ))
  cat(glue::glue(
    "  Posterior prob     : {round(x$probability, 4)} ",
    "({if (x$go) 'meets' else 'does NOT meet'} threshold)\n"
  ))
  cat(glue::glue("  Posterior mean     : {round(x$posterior_mean, 4)}\n"))
  cat(glue::glue(
    "  95% CI             : [{round(x$posterior_ci['lo'], 4)}, {round(x$posterior_ci['hi'], 4)}]\n"
  ))
  cat(glue::glue("  Draws used         : {x$n_draws}\n"))
  invisible(x)
}
