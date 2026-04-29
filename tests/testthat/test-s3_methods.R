# Shared mock helpers ──────────────────────────────────────────────────────────

.mock_map_s3 <- function(outcome = "binary", K = 2L) {
  set.seed(10)
  mix <- list(
    weights = rep(1/K, K),
    means   = seq(-0.5, 0.5, length.out = K),
    sds     = rep(0.4, K),
    K       = K
  )
  curr_data <- structure(
    list(outcome = outcome, n_hist = 4L,
         hist_data = data.frame(study = paste0("S", 1:4)),
         stan_data_current = list(n_curr = 80L, r_curr = 22L)),
    class = "bb_data"
  )
  structure(
    list(
      mixture            = mix,
      mu_pred_draws      = rnorm(800, -0.3, 0.4),
      outcome_type       = outcome,
      model_type         = "standard",
      hyperpriors        = list(),
      ess_prior          = 22.0,
      stan_fit           = NULL,
      draws              = NULL,
      data               = curr_data,
      mixture_components = K
    ),
    class = "map_prior"
  )
}

.mock_rmap_s3 <- function() {
  mp  <- .mock_map_s3()
  rmp <- robustify_map(mp, weight = 0.8)
  rmp
}

.mock_bf_s3 <- function(outcome = "binary") {
  set.seed(20)
  n_draws <- 400L
  draws_df <- posterior::as_draws_df(posterior::draws_matrix(
    theta_curr = rnorm(n_draws, -0.2, 0.35),
    p_curr     = plogis(rnorm(n_draws, -0.2, 0.35)),
    .nchains   = 4L
  ))
  prior <- .mock_map_s3(outcome)
  curr  <- list(n_curr = 80L, r_curr = 22L)
  bd    <- structure(
    list(outcome = outcome, stan_data_current = curr),
    class = "bb_data"
  )
  structure(
    list(
      fit = NULL, draws = draws_df, prior = prior,
      data = bd, model = "map", outcome_type = outcome,
      stan_code = "", diagnostics = NULL
    ),
    class = "borrowing_fit"
  )
}

# ── print.map_prior ───────────────────────────────────────────────────────────

test_that("print.map_prior produces MAP header", {
  mp <- .mock_map_s3()
  expect_output(print(mp), "MAP Prior")
  expect_output(print(mp), "Mixture")
  expect_output(print(mp), "binary")
})

test_that("print.map_prior shows ESS", {
  mp <- .mock_map_s3()
  expect_output(print(mp), "ESS")
})

# ── plot.map_prior ────────────────────────────────────────────────────────────

test_that("plot.map_prior returns a ggplot", {
  mp <- .mock_map_s3()
  p  <- plot(mp)
  expect_s3_class(p, "ggplot")
})

# ── print.rmap_prior ──────────────────────────────────────────────────────────

test_that("print.rmap_prior produces Robust MAP header", {
  rmp <- .mock_rmap_s3()
  expect_output(print(rmp), "Robust MAP")
})

test_that("print.rmap_prior shows Vague component", {
  rmp <- .mock_rmap_s3()
  expect_output(print(rmp), "Vague")
})

# ── plot.rmap_prior ───────────────────────────────────────────────────────────

test_that("plot.rmap_prior returns a ggplot", {
  rmp <- .mock_rmap_s3()
  p   <- plot(rmp)
  expect_s3_class(p, "ggplot")
})

# ── print.borrowing_fit ───────────────────────────────────────────────────────

test_that("print.borrowing_fit produces header", {
  bf <- .mock_bf_s3()
  expect_output(print(bf), "Borrowing Model")
})

test_that("print.borrowing_fit shows posterior summary", {
  bf <- .mock_bf_s3()
  expect_output(print(bf), "Posterior summary")
})

# ── plot.borrowing_fit ────────────────────────────────────────────────────────

test_that("plot.borrowing_fit type prior_posterior returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot(bf, type = "prior_posterior")
  expect_s3_class(p, "ggplot")
})

test_that("plot.borrowing_fit type shrinkage returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot(bf, type = "shrinkage")
  expect_s3_class(p, "ggplot")
})

test_that("plot.borrowing_fit type effect returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot(bf, type = "effect")
  expect_s3_class(p, "ggplot")
})

# ── plot_prior_posterior standalone ──────────────────────────────────────────

test_that("plot_prior_posterior returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot_prior_posterior(bf)
  expect_s3_class(p, "ggplot")
})

test_that("plot_prior_posterior transform link returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot_prior_posterior(bf, transform = "link")
  expect_s3_class(p, "ggplot")
})

# ── plot_borrowing_effect standalone ─────────────────────────────────────────

test_that("plot_borrowing_effect returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot_borrowing_effect(bf)
  expect_s3_class(p, "ggplot")
})

# ── plot_shrinkage standalone ─────────────────────────────────────────────────

test_that("plot_shrinkage returns ggplot", {
  bf <- .mock_bf_s3()
  p  <- plot_shrinkage(bf)
  expect_s3_class(p, "ggplot")
})
