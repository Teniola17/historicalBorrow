#' Prepare historical and current trial data for Bayesian borrowing
#'
#' Validates and structures data for use with [fit_map()] and
#' [fit_borrowing_model()]. No Stan compilation occurs at this step.
#'
#' @param hist_data A data frame of historical study summaries. Required columns
#'   depend on `outcome`:
#'   * **binary**: `study`, `n` (integer), `r` (integer events)
#'   * **continuous**: `study`, `n`, `mean`, `sd`
#'   * **count**: `study`, `n` (person-time, numeric), `events` (integer)
#' @param current_data A data frame for the current trial. Same column
#'   requirements as `hist_data`. The response/events column may be `NA` if
#'   the current trial has not yet observed an outcome (forward planning).
#' @param outcome One of `"binary"`, `"continuous"`, or `"count"`.
#' @param study_col Name of the column that identifies each study. Defaults to
#'   `"study"`.
#'
#' @return An object of class `bb_data` (a named list) containing:
#'   \describe{
#'     \item{`outcome`}{Outcome type string.}
#'     \item{`hist_data`}{Standardised historical data frame.}
#'     \item{`current_data`}{Standardised current data frame.}
#'     \item{`n_hist`}{Number of historical studies.}
#'     \item{`stan_data_map`}{List formatted for the MAP Stan model.}
#'     \item{`stan_data_current`}{List formatted for borrowing Stan models.}
#'   }
#'
#' @examples
#' hist <- data.frame(
#'   study = paste0("S", 1:4),
#'   n = c(100L, 120L, 90L, 110L),
#'   r = c(30L, 40L, 25L, 38L)
#' )
#' curr <- data.frame(study = "C", n = 80L, r = 22L)
#' pd <- prepare_data(hist, curr, outcome = "binary")
#' print(pd)
#'
#' @export
prepare_data <- function(hist_data,
                         current_data,
                         outcome   = c("binary", "continuous", "count"),
                         study_col = "study") {

  outcome <- match.arg(outcome)

  checkmate::assert_data_frame(hist_data,   min.rows = 2)
  checkmate::assert_data_frame(current_data, min.rows = 1)
  checkmate::assert_string(study_col)

  # Rename study column to "study" if needed
  hist_data    <- .rename_study_col(hist_data,    study_col)
  current_data <- .rename_study_col(current_data, study_col)

  switch(outcome,
    binary     = .prepare_binary(hist_data, current_data),
    continuous = .prepare_continuous(hist_data, current_data),
    count      = .prepare_count(hist_data, current_data)
  )
}

# ── Internal: rename study column ─────────────────────────────────────────────

.rename_study_col <- function(df, study_col) {
  if (study_col != "study") {
    if (!study_col %in% names(df)) {
      stop(glue::glue("Column '{study_col}' not found in data."), call. = FALSE)
    }
    names(df)[names(df) == study_col] <- "study"
  }
  df
}

# ── Internal: binary ──────────────────────────────────────────────────────────

.prepare_binary <- function(hist, curr) {
  req <- c("study", "n", "r")
  .check_cols(hist, req, "hist_data (binary)")

  hist$n <- as.integer(hist$n)
  hist$r <- as.integer(hist$r)

  if (any(hist$n <= 0L, na.rm = TRUE))
    stop("All `n` values in hist_data must be positive.", call. = FALSE)
  if (any(hist$r > hist$n, na.rm = TRUE))
    stop("Some r > n in hist_data: events cannot exceed sample size.", call. = FALSE)
  if (any(hist$r < 0L, na.rm = TRUE))
    stop("r values in hist_data must be non-negative.", call. = FALSE)

  curr$n <- as.integer(curr$n)
  if ("r" %in% names(curr)) curr$r <- as.integer(curr$r)

  stan_map <- list(
    J = nrow(hist),
    n = hist$n,
    r = hist$r
  )

  stan_curr <- list(
    n_curr = curr$n[1],
    r_curr = if ("r" %in% names(curr) && !is.na(curr$r[1])) curr$r[1] else 0L
  )

  .make_bb_data("binary", hist, curr, nrow(hist), stan_map, stan_curr)
}

# ── Internal: continuous ──────────────────────────────────────────────────────

.prepare_continuous <- function(hist, curr) {
  req <- c("study", "n", "mean", "sd")
  .check_cols(hist, req, "hist_data (continuous)")

  hist$n    <- as.integer(hist$n)
  hist$mean <- as.numeric(hist$mean)
  hist$sd   <- as.numeric(hist$sd)

  if (any(hist$sd <= 0, na.rm = TRUE))
    stop("All `sd` values in hist_data must be positive.", call. = FALSE)
  if (any(hist$n <= 0L, na.rm = TRUE))
    stop("All `n` values in hist_data must be positive.", call. = FALSE)

  hist$se <- hist$sd / sqrt(hist$n)

  curr$n    <- as.integer(curr$n)
  curr$mean <- if ("mean" %in% names(curr)) as.numeric(curr$mean) else NA_real_
  curr$sd   <- if ("sd"   %in% names(curr)) as.numeric(curr$sd)   else NA_real_
  curr$se   <- if (!is.na(curr$sd[1])) curr$sd / sqrt(curr$n) else NA_real_

  stan_map <- list(
    J  = nrow(hist),
    n  = hist$n,
    y  = hist$mean,
    se = hist$se
  )

  stan_curr <- list(
    n_curr  = curr$n[1],
    y_curr  = if (!is.na(curr$mean[1])) curr$mean[1] else 0.0,
    se_curr = if (!is.na(curr$se[1]))   curr$se[1]   else 1.0
  )

  .make_bb_data("continuous", hist, curr, nrow(hist), stan_map, stan_curr)
}

# ── Internal: count ────────────────────────────────────────────────────────────

.prepare_count <- function(hist, curr) {
  req <- c("study", "n", "events")
  .check_cols(hist, req, "hist_data (count)")

  hist$n      <- as.numeric(hist$n)
  hist$events <- as.integer(hist$events)

  if (any(hist$n <= 0, na.rm = TRUE))
    stop("All `n` (person-time) values in hist_data must be positive.", call. = FALSE)
  if (any(hist$events < 0L, na.rm = TRUE))
    stop("`events` values in hist_data must be non-negative.", call. = FALSE)

  hist$log_n <- log(hist$n)

  curr$n      <- as.numeric(curr$n)
  curr$events <- if ("events" %in% names(curr)) as.integer(curr$events) else NA_integer_

  stan_map <- list(
    J      = nrow(hist),
    events = hist$events,
    log_n  = hist$log_n
  )

  stan_curr <- list(
    events_curr = if (!is.na(curr$events[1])) curr$events[1] else 0L,
    log_n_curr  = log(curr$n[1])
  )

  .make_bb_data("count", hist, curr, nrow(hist), stan_map, stan_curr)
}

# ── Internal: helpers ─────────────────────────────────────────────────────────

.check_cols <- function(df, required, df_name) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      glue::glue(
        "{df_name} is missing required columns: {paste(missing, collapse = ', ')}."
      ),
      call. = FALSE
    )
  }
}

.make_bb_data <- function(outcome, hist, curr, n_hist, stan_map, stan_curr) {
  structure(
    list(
      outcome           = outcome,
      hist_data         = hist,
      current_data      = curr,
      n_hist            = n_hist,
      stan_data_map     = stan_map,
      stan_data_current = stan_curr
    ),
    class = "bb_data"
  )
}

# ── S3 methods ────────────────────────────────────────────────────────────────

#' @export
print.bb_data <- function(x, ...) {
  cat("── historicalBorrow: Prepared Data ──────────────────────────────────────\n")
  cat(glue::glue("  Outcome type : {x$outcome}\n"))
  cat(glue::glue("  Historical   : {x$n_hist} studies\n"))

  hist_n <- switch(x$outcome,
    binary     = sum(x$hist_data$n),
    continuous = sum(x$hist_data$n),
    count      = sum(x$hist_data$events)
  )
  label <- if (x$outcome == "count") "total events" else "total patients"
  cat(glue::glue("               ({hist_n} {label})\n"))

  curr_n <- switch(x$outcome,
    binary     = x$stan_data_current$n_curr,
    continuous = x$stan_data_current$n_curr,
    count      = exp(x$stan_data_current$log_n_curr)
  )
  cat(glue::glue("  Current trial: {curr_n} {if (x$outcome == 'count') 'person-time' else 'patients'}\n"))
  invisible(x)
}

#' @export
summary.bb_data <- function(object, ...) {
  cat("── historicalBorrow: Data Summary ───────────────────────────────────────\n")
  cat(glue::glue("Outcome: {object$outcome}\n\n"))
  cat("Historical studies:\n")
  print(object$hist_data)
  cat("\nCurrent trial:\n")
  print(object$current_data)
  invisible(object)
}
