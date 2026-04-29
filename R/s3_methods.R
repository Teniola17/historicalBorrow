# S3 methods for historicalBorrow objects
# print / summary / plot dispatchers

# ── map_prior ─────────────────────────────────────────────────────────────────

#' @export
print.map_prior <- function(x, ...) {
  cat("── historicalBorrow: MAP Prior ──────────────────────────────────────────\n")
  cat(glue::glue("  Outcome type : {x$outcome_type}\n"))
  cat(glue::glue("  Model        : {x$model_type}\n"))
  cat(glue::glue("  Studies used : {x$data$n_hist}\n"))
  cat(glue::glue("  Mixture K    : {x$mixture$K} component(s)\n"))
  cat(glue::glue("  ESS (prior)  : ~{x$ess_prior} patients\n\n"))

  mix_tbl <- data.frame(
    Component = seq_len(x$mixture$K),
    Weight    = round(x$mixture$weights, 4),
    Mean      = round(x$mixture$means,   4),
    SD        = round(x$mixture$sds,     4)
  )
  cat("  Mixture components (MAP predictive posterior):\n")
  print(mix_tbl, row.names = FALSE)
  invisible(x)
}

#' @export
summary.map_prior <- function(object, ...) {
  print(object)
  if (!is.null(object$draws)) {
    cat("\n  Posterior summary (mu, tau):\n")
    summ <- posterior::summarise_draws(
      posterior::subset_draws(object$draws, variable = c("mu", "tau")),
      mean, sd, ~quantile(.x, c(0.025, 0.975))
    )
    print(as.data.frame(summ))
  }
  invisible(object)
}

#' @export
plot.map_prior <- function(x, ...) {
  mix      <- x$mixture
  draws    <- x$mu_pred_draws

  x_range  <- range(draws, na.rm = TRUE)
  padding  <- diff(x_range) * 0.4
  x_grid   <- seq(x_range[1] - padding, x_range[2] + padding, length.out = 500L)
  mix_dens <- .eval_mixture_density(x_grid, mix)
  mix_df   <- data.frame(x = x_grid, density = mix_dens)

  p <- ggplot2::ggplot() +
    ggplot2::geom_histogram(
      data    = data.frame(draws = draws),
      mapping = ggplot2::aes(x = .data$draws, y = ggplot2::after_stat(density)),
      bins    = 50, fill = "steelblue", alpha = 0.4, colour = "white"
    ) +
    ggplot2::geom_line(
      data    = mix_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$density),
      colour  = "#E66101",
      linewidth = 1.0
    ) +
    ggplot2::labs(
      title    = "MAP Prior",
      subtitle = glue::glue(
        "Outcome: {x$outcome_type} | ",
        "K = {mix$K} | ESS ≈ {x$ess_prior}"
      ),
      x = glue::glue("mu_pred ({x$outcome_type} link scale)"),
      y = "Density"
    ) +
    ggplot2::theme_bw()

  p
}

# ── rmap_prior ────────────────────────────────────────────────────────────────

#' @export
print.rmap_prior <- function(x, ...) {
  cat("── historicalBorrow: Robust MAP Prior ───────────────────────────────────\n")
  cat(glue::glue("  Outcome type  : {x$outcome_type}\n"))
  cat(glue::glue("  MAP weight    : {x$rmap_weight}\n"))
  cat(glue::glue("  Vague weight  : {round(1 - x$rmap_weight, 4)}\n"))
  cat(glue::glue("  MAP ESS       : ~{x$ess_prior} patients\n"))
  cat(glue::glue("  RMAP ESS      : ~{x$ess_prior_rmap} patients\n\n"))

  mix <- x$mixture_rmap
  n   <- mix$K
  labels <- c(
    paste0("MAP-", seq_len(n - 1)),
    "Vague"
  )
  mix_tbl <- data.frame(
    Component = labels,
    Weight    = round(mix$weights, 4),
    Mean      = round(mix$means,   4),
    SD        = round(mix$sds,     4)
  )
  cat("  Extended mixture (RMAP):\n")
  print(mix_tbl, row.names = FALSE)
  invisible(x)
}

#' @export
summary.rmap_prior <- function(object, ...) {
  print(object)
  invisible(object)
}

#' @export
plot.rmap_prior <- function(x, ...) {
  # Show MAP mixture and the vague component separately
  mix  <- x$mixture_rmap
  K    <- mix$K
  n_g  <- 500L

  # Extent of the plot range
  x_range <- c(
    min(mix$means - 3 * mix$sds),
    max(mix$means + 3 * mix$sds)
  )
  x_grid  <- seq(x_range[1], x_range[2], length.out = n_g)

  # Full RMAP density
  rmap_dens <- .eval_mixture_density(x_grid, mix)

  # MAP-only density (first K-1 components, re-normalised)
  map_mix <- list(
    weights = mix$weights[-K] / sum(mix$weights[-K]),
    means   = mix$means[-K],
    sds     = mix$sds[-K],
    K       = K - 1L
  )
  map_dens  <- .eval_mixture_density(x_grid, map_mix) * sum(mix$weights[-K])

  vague_dens <- mix$weights[K] *
    dnorm(x_grid, mix$means[K], mix$sds[K])

  plot_df <- rbind(
    data.frame(x = x_grid, density = rmap_dens, source = "RMAP (total)"),
    data.frame(x = x_grid, density = map_dens,  source = "MAP component"),
    data.frame(x = x_grid, density = vague_dens, source = "Vague component")
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$x, y = .data$density,
                 colour = .data$source, linetype = .data$source)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(
      values = c(
        "RMAP (total)"   = "black",
        "MAP component"  = "#E66101",
        "Vague component" = "steelblue"
      )
    ) +
    ggplot2::labs(
      title    = "Robust MAP Prior",
      subtitle = glue::glue(
        "MAP weight = {x$rmap_weight} | RMAP ESS ≈ {x$ess_prior_rmap}"
      ),
      x        = glue::glue("Parameter ({x$outcome_type} link scale)"),
      y        = "Density",
      colour   = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

# ── borrowing_fit ─────────────────────────────────────────────────────────────

#' @export
print.borrowing_fit <- function(x, ...) {
  cat("── historicalBorrow: Borrowing Model Fit ────────────────────────────────\n")
  cat(glue::glue("  Outcome : {x$outcome_type}\n"))
  cat(glue::glue("  Model   : {toupper(x$model)}\n\n"))

  params <- .primary_params(x$outcome_type, x$model)
  available <- posterior::variables(x$draws)
  params <- intersect(params, available)

  if (length(params) > 0) {
    summ <- posterior::summarise_draws(
      posterior::subset_draws(x$draws, variable = params),
      mean, sd,
      ~quantile(.x, 0.025, na.rm = TRUE),
      ~quantile(.x, 0.975, na.rm = TRUE)
    )
    names(summ)[3:4] <- c("q2.5", "q97.5")
    summ[, 2:5] <- round(summ[, 2:5], 4)
    cat("  Posterior summary:\n")
    print(as.data.frame(summ), row.names = FALSE)
  }
  invisible(x)
}

#' @export
summary.borrowing_fit <- function(object, ...) {
  print(object)
  cat("\n")
  if (!is.null(object$diagnostics)) {
    print(object$diagnostics)
  } else {
    cat("  Run diagnostics(fit) for convergence summaries.\n")
  }
  invisible(object)
}

#' @export
plot.borrowing_fit <- function(x,
                                type = c("prior_posterior", "effect",
                                         "shrinkage", "trace"),
                                ...) {
  type <- match.arg(type)
  switch(type,
    prior_posterior = plot_prior_posterior(x, ...),
    effect          = plot_borrowing_effect(x, ...),
    shrinkage       = plot_shrinkage(x, ...),
    trace           = {
      diag <- diagnostics(x, ...)
      diag$trace_plots
    }
  )
}
