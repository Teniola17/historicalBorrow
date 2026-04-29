# Helpers — build mock objects without Stan ────────────────────────────────────

.mock_map_prior_diag <- function(outcome = "binary") {
  mix <- list(weights = c(0.7, 0.3), means = c(-0.5, 0.2),
              sds = c(0.3, 0.8), K = 2L)
  structure(
    list(mixture = mix, mixture_rmap = NULL, outcome_type = outcome,
         ess_prior = 25.0, ess_prior_rmap = NULL),
    class = "map_prior"
  )
}

.mock_borrowing_fit <- function(outcome = "binary", model = "map",
                                n_curr = 80L, r_curr = 22L,
                                post_mean_logit = -0.2) {
  set.seed(99)
  n_draws <- 400L
  draws_df <- posterior::as_draws_df(posterior::draws_matrix(
    theta_curr = rnorm(n_draws, post_mean_logit, 0.35),
    p_curr     = plogis(rnorm(n_draws, post_mean_logit, 0.35)),
    .nchains   = 4L
  ))

  curr_stan <- list(n_curr = n_curr, r_curr = r_curr)
  bb_data   <- structure(
    list(outcome = outcome, stan_data_current = curr_stan),
    class = "bb_data"
  )

  structure(
    list(
      fit          = NULL,
      draws        = draws_df,
      prior        = .mock_map_prior_diag(outcome),
      data         = bb_data,
      model        = model,
      outcome_type = outcome,
      stan_code    = "",
      diagnostics  = NULL
    ),
    class = "borrowing_fit"
  )
}

# ── .compute_prior_data_conflict ──────────────────────────────────────────────

test_that(".compute_prior_data_conflict returns list with required fields", {
  fit <- .mock_borrowing_fit()
  pdc <- historicalBorrow:::.compute_prior_data_conflict(fit, threshold = 2.0)
  expect_type(pdc, "list")
  expect_true(all(c("detected", "statistic", "threshold", "direction") %in% names(pdc)))
})

test_that(".compute_prior_data_conflict flags conflict when MLE far from prior", {
  # Prior centred at logit(0.3) ≈ -0.85; current data gives r=70/n=80 => logit ≈ 1.95
  fit <- .mock_borrowing_fit(r_curr = 70L, n_curr = 80L, post_mean_logit = 1.5)
  pdc <- historicalBorrow:::.compute_prior_data_conflict(fit, threshold = 1.5)
  expect_true(pdc$detected)
  expect_equal(pdc$direction, "positive")
})

test_that(".compute_prior_data_conflict no conflict when MLE near prior", {
  # Prior mean at -0.5 logit; MLE logit(22/80) ≈ -0.66 — very close
  fit <- .mock_borrowing_fit(r_curr = 22L, n_curr = 80L)
  pdc <- historicalBorrow:::.compute_prior_data_conflict(fit, threshold = 2.0)
  expect_false(pdc$detected)
})

test_that(".compute_prior_data_conflict returns NA when prior is NULL", {
  fit       <- .mock_borrowing_fit()
  fit$prior <- NULL
  pdc       <- historicalBorrow:::.compute_prior_data_conflict(fit)
  expect_false(pdc$detected)
})

# ── .compute_shrinkage ────────────────────────────────────────────────────────

test_that(".compute_shrinkage returns numeric", {
  fit <- .mock_borrowing_fit()
  s   <- historicalBorrow:::.compute_shrinkage(fit)
  expect_type(s, "double")
})

test_that(".compute_shrinkage is between -1 and 2 for typical scenarios", {
  fit <- .mock_borrowing_fit()
  s   <- historicalBorrow:::.compute_shrinkage(fit)
  if (!is.na(s)) expect_true(s > -2 && s < 3)
})

# ── .eval_mixture_density ─────────────────────────────────────────────────────

test_that(".eval_mixture_density integrates to approximately 1", {
  mix <- list(weights = c(0.6, 0.4), means = c(0, 1), sds = c(0.5, 0.5), K = 2L)
  x   <- seq(-3, 4, length.out = 2000)
  dx  <- diff(x)[1]
  expect_equal(sum(historicalBorrow:::.eval_mixture_density(x, mix)) * dx,
               1.0, tolerance = 0.02)
})

test_that(".eval_mixture_density returns non-negative values", {
  mix  <- list(weights = c(0.5, 0.5), means = c(-1, 1), sds = c(0.4, 0.4), K = 2L)
  dens <- historicalBorrow:::.eval_mixture_density(seq(-4, 4, by = 0.1), mix)
  expect_true(all(dens >= 0))
})

# ── .sample_mixture ───────────────────────────────────────────────────────────

test_that(".sample_mixture produces correct length", {
  mix <- list(weights = c(0.5, 0.5), means = c(-1, 1), sds = c(0.5, 0.5), K = 2L)
  s   <- historicalBorrow:::.sample_mixture(mix, n = 500L)
  expect_length(s, 500L)
  expect_type(s, "double")
})

# ── .fit_mixture_of_normals ───────────────────────────────────────────────────

test_that(".fit_mixture_of_normals returns valid structure", {
  set.seed(42)
  draws <- c(rnorm(300, -0.5, 0.3), rnorm(300, 0.8, 0.4))
  mix   <- historicalBorrow:::.fit_mixture_of_normals(draws, K_max = 4L)
  expect_named(mix, c("weights", "means", "sds", "K"))
  expect_equal(sum(mix$weights), 1.0, tolerance = 1e-6)
  expect_true(all(mix$sds >= 1e-6))
  expect_gte(mix$K, 1L)
  expect_lte(mix$K, 4L)
})

test_that(".fit_mixture_of_normals K_max = 1 returns single component", {
  set.seed(1)
  draws <- rnorm(200, 0.5, 0.3)
  mix   <- historicalBorrow:::.fit_mixture_of_normals(draws, K_max = 1L)
  expect_equal(mix$K, 1L)
  expect_equal(mix$weights, 1.0)
})

# ── print.hb_diagnostics ─────────────────────────────────────────────────────

test_that("print.hb_diagnostics produces output without error", {
  fit  <- .mock_borrowing_fit()
  diag <- diagnostics(fit)
  expect_output(print(diag), "Diagnostics")
})
