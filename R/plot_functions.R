#' Plot prior and posterior density overlay
#'
#' Draws the MAP/RMAP prior density and the posterior density of the primary
#' treatment parameter on the same panel, with optional credible interval
#' ribbons.
#'
#' @param fit A `borrowing_fit` object.
#' @param transform `"natural"` (default) plots on the natural scale (probability
#'   for binary, mean for continuous, rate for count). `"logit"` / `"log"` plots
#'   on the link scale.
#' @param ci_width Width of the credible interval ribbon (default `0.95`).
#' @param n_grid Number of grid points for the prior density curve (default `500`).
#'
#' @return A `ggplot` object.
#' @export
plot_prior_posterior <- function(fit,
                                 transform = c("natural", "link"),
                                 ci_width  = 0.95,
                                 n_grid    = 500L) {

  if (!inherits(fit, "borrowing_fit"))
    stop("`fit` must be a `borrowing_fit` object.", call. = FALSE)

  transform <- match.arg(transform)
  alpha     <- (1 - ci_width) / 2

  outcome  <- fit$outcome_type
  prior    <- fit$prior
  mix      <- if (!is.null(prior)) {
    prior$mixture_rmap %||% prior$mixture
  } else NULL

  # Posterior draws of primary parameter (link scale)
  param_link <- switch(outcome,
    binary     = "theta_curr",
    continuous = "theta_curr",
    count      = "log_lambda_curr"
  )
  post_link <- tryCatch(
    posterior::extract_variable(fit$draws, variable = param_link),
    error = function(e) NULL
  )

  # Natural-scale parameter if available
  param_nat <- switch(outcome,
    binary     = "p_curr",
    continuous = "theta_curr",
    count      = "lambda_curr"
  )
  post_nat <- tryCatch(
    posterior::extract_variable(fit$draws, variable = param_nat),
    error = function(e) post_link
  )

  post_vals <- if (transform == "natural") post_nat else post_link

  if (is.null(post_vals))
    stop("Could not extract posterior draws.", call. = FALSE)

  # ── Prior density on appropriate scale ──────────────────────────────────────
  x_range <- range(post_vals, na.rm = TRUE)
  padding  <- diff(x_range) * 0.4
  x_grid   <- seq(x_range[1] - padding, x_range[2] + padding, length.out = n_grid)

  prior_df <- NULL
  if (!is.null(mix)) {
    mix_for_plot <- mix
    if (transform == "natural") {
      # Sample from mixture on link scale and transform
      mix_samples <- .sample_mixture(mix_for_plot, n = 20000L)
      nat_samples <- switch(outcome,
        binary     = .inv_logit(mix_samples),
        continuous = mix_samples,
        count      = exp(mix_samples)
      )
      dens_prior <- density(nat_samples, from = x_grid[1], to = x_grid[n_grid], n = n_grid)
      prior_df   <- data.frame(x = dens_prior$x, density = dens_prior$y, source = "Prior")
    } else {
      prior_dens <- .eval_mixture_density(x_grid, mix_for_plot)
      prior_df   <- data.frame(x = x_grid, density = prior_dens, source = "Prior")
    }
  }

  # ── Posterior kernel density ─────────────────────────────────────────────────
  post_dens  <- density(post_vals, from = x_grid[1], to = x_grid[n_grid], n = n_grid)
  post_df    <- data.frame(x = post_dens$x, density = post_dens$y, source = "Posterior")

  plot_data  <- if (!is.null(prior_df)) rbind(prior_df, post_df) else post_df

  # ── Labels ───────────────────────────────────────────────────────────────────
  x_label <- switch(outcome,
    binary     = if (transform == "natural") "Response probability" else "logit(p)",
    continuous = "Mean",
    count      = if (transform == "natural") "Event rate" else "log(rate)"
  )

  ci_lo <- quantile(post_vals, alpha,       na.rm = TRUE)
  ci_hi <- quantile(post_vals, 1 - alpha,   na.rm = TRUE)
  post_mode <- post_dens$x[which.max(post_dens$y)]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$x, y = .data$density,
    colour = .data$source, linetype = .data$source
  )) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_vline(
      xintercept = c(ci_lo, ci_hi),
      linetype   = "dashed",
      colour     = "steelblue",
      linewidth  = 0.4,
      alpha      = 0.7
    ) +
    ggplot2::scale_colour_manual(
      values = c(Prior = "#E66101", Posterior = "#5E81AC")
    ) +
    ggplot2::labs(
      title    = "Prior vs Posterior",
      subtitle = glue::glue("{round(ci_width * 100)}% credible interval: [{round(ci_lo, 3)}, {round(ci_hi, 3)}]"),
      x        = x_label,
      y        = "Density",
      colour   = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  p
}

#' Plot borrowing effect — credible interval comparison
#'
#' Shows side-by-side 95% credible intervals for the primary parameter under
#' the fitted borrowing model versus a reference no-borrowing analysis.
#'
#' @param fit A `borrowing_fit` object.
#' @param no_borrow_fit Optional `borrowing_fit` with `model = "map"` fitted
#'   using a flat / very vague prior, used as the "no borrowing" reference.
#'   If `NULL`, the MLE from the current data is used as the reference point.
#' @param ci_width Credible interval width (default `0.95`).
#'
#' @return A `ggplot` object.
#' @export
plot_borrowing_effect <- function(fit,
                                  no_borrow_fit = NULL,
                                  ci_width      = 0.95) {

  if (!inherits(fit, "borrowing_fit"))
    stop("`fit` must be a `borrowing_fit` object.", call. = FALSE)

  alpha    <- (1 - ci_width) / 2
  outcome  <- fit$outcome_type

  param_nat <- switch(outcome,
    binary     = "p_curr",
    continuous = "theta_curr",
    count      = "lambda_curr"
  )

  .get_ci <- function(f, label) {
    draws <- tryCatch(
      posterior::extract_variable(f$draws, variable = param_nat),
      error = function(e) NULL
    )
    if (is.null(draws)) return(NULL)
    data.frame(
      model = label,
      mean  = mean(draws, na.rm = TRUE),
      lo    = quantile(draws, alpha, na.rm = TRUE),
      hi    = quantile(draws, 1 - alpha, na.rm = TRUE)
    )
  }

  rows <- list(.get_ci(fit, paste0("With borrowing\n(", toupper(fit$model), ")")))

  if (!is.null(no_borrow_fit) && inherits(no_borrow_fit, "borrowing_fit")) {
    rows <- c(rows, list(.get_ci(no_borrow_fit, "No borrowing")))
  } else {
    # Use MLE as reference point
    curr <- fit$data$stan_data_current
    mle  <- switch(outcome,
      binary     = curr$r_curr / curr$n_curr,
      continuous = curr$y_curr,
      count      = exp(curr$events_curr - curr$log_n_curr)
    )
    rows <- c(rows, list(data.frame(
      model = "Current data\n(MLE)",
      mean  = mle,
      lo    = NA_real_,
      hi    = NA_real_
    )))
  }

  plot_data       <- do.call(rbind, rows[!sapply(rows, is.null)])
  plot_data$model <- factor(plot_data$model, levels = rev(plot_data$model))

  x_label <- switch(outcome,
    binary     = "Response probability",
    continuous = "Mean",
    count      = "Event rate"
  )

  ggplot2::ggplot(plot_data, ggplot2::aes(
    x    = .data$mean,
    y    = .data$model,
    xmin = .data$lo,
    xmax = .data$hi,
    colour = .data$model
  )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(width = 0.2, linewidth = 0.8, na.rm = TRUE) +
    ggplot2::labs(
      title  = "Borrowing Effect",
      subtitle = glue::glue("{round(ci_width * 100)}% credible intervals"),
      x      = x_label,
      y      = NULL,
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

#' Plot shrinkage of the posterior towards the MAP prior
#'
#' Compares the current-data MLE with the posterior mean incorporating
#' borrowing, and annotates with the shrinkage coefficient.
#'
#' @param fit A `borrowing_fit` object.
#'
#' @return A `ggplot` object.
#' @export
plot_shrinkage <- function(fit) {

  if (!inherits(fit, "borrowing_fit"))
    stop("`fit` must be a `borrowing_fit` object.", call. = FALSE)

  outcome <- fit$outcome_type
  prior   <- fit$prior
  mix     <- if (!is.null(prior)) prior$mixture_rmap %||% prior$mixture else NULL

  param_nat <- switch(outcome,
    binary     = "p_curr",
    continuous = "theta_curr",
    count      = "lambda_curr"
  )

  post_vals <- tryCatch(
    posterior::extract_variable(fit$draws, variable = param_nat),
    error = function(e) NULL
  )
  if (is.null(post_vals)) stop("Cannot extract posterior draws.", call. = FALSE)

  post_mean <- mean(post_vals, na.rm = TRUE)

  curr <- fit$data$stan_data_current
  mle  <- switch(outcome,
    binary     = curr$r_curr / curr$n_curr,
    continuous = curr$y_curr,
    count      = exp(curr$events_curr - curr$log_n_curr)
  )

  prior_mean_nat <- if (!is.null(mix)) {
    lm <- sum(mix$weights * mix$means)
    switch(outcome,
      binary     = .inv_logit(lm),
      continuous = lm,
      count      = exp(lm)
    )
  } else NA_real_

  shrink <- .compute_shrinkage(fit)

  pts <- data.frame(
    label = c("Current MLE", "Posterior mean"),
    value = c(mle, post_mean)
  )

  x_label <- switch(outcome,
    binary = "Response probability", continuous = "Mean", count = "Event rate")

  p <- ggplot2::ggplot(pts, ggplot2::aes(
    x = .data$value, y = .data$label, colour = .data$label
  )) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_segment(
      data = data.frame(x = mle, xend = post_mean, y = 1.5, yend = 1.5),
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      arrow     = ggplot2::arrow(length = ggplot2::unit(0.2, "cm"), ends = "last"),
      colour    = "grey40",
      linewidth = 0.6,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title    = "Shrinkage towards MAP Prior",
      subtitle = if (!is.na(shrink))
        glue::glue("Shrinkage coefficient: {round(shrink, 3)}")
      else
        "Shrinkage: N/A",
      x      = x_label,
      y      = NULL,
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  if (!is.na(prior_mean_nat)) {
    p <- p + ggplot2::geom_vline(
      xintercept = prior_mean_nat,
      linetype   = "dashed",
      colour     = "#E66101",
      linewidth  = 0.5
    ) +
      ggplot2::annotate(
        "text", x = prior_mean_nat, y = 2.4,
        label  = "Prior mean",
        colour = "#E66101",
        size   = 3,
        hjust  = -0.1
      )
  }

  p
}
