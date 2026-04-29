# Mock borrowing_fit for decision criteria tests ──────────────────────────────

.mock_bf_dc <- function(outcome = "binary", p_true = 0.3, n_draws = 400L) {
  set.seed(77)
  logit_p <- qlogis(p_true)
  draws_df <- posterior::as_draws_df(posterior::draws_matrix(
    theta_curr = rnorm(n_draws, logit_p, 0.3),
    p_curr     = plogis(rnorm(n_draws, logit_p, 0.3)),
    .nchains   = 4L
  ))

  curr <- list(n_curr = 80L, r_curr = round(p_true * 80))
  bd   <- structure(
    list(outcome = outcome, stan_data_current = curr),
    class = "bb_data"
  )

  structure(
    list(
      fit = NULL, draws = draws_df, prior = NULL,
      data = bd, model = "rmap", outcome_type = outcome,
      stan_code = "", diagnostics = NULL
    ),
    class = "borrowing_fit"
  )
}

test_that("decision_criteria returns decision_result", {
  fit <- .mock_bf_dc(p_true = 0.35)
  dc  <- decision_criteria(fit, theta_threshold = 0.20, direction = "greater")
  expect_s3_class(dc, "decision_result")
})

test_that("Go when posterior prob exceeds probability_threshold", {
  # p_true = 0.5 >> 0.20, so P(p > 0.20) should be very high
  fit <- .mock_bf_dc(p_true = 0.50)
  dc  <- decision_criteria(fit, theta_threshold = 0.20,
                            probability_threshold = 0.80)
  expect_true(dc$go)
  expect_gt(dc$probability, 0.80)
})

test_that("No-Go when posterior prob below probability_threshold", {
  # p_true = 0.10 << 0.40, so P(p > 0.40) should be very low
  fit <- .mock_bf_dc(p_true = 0.10)
  dc  <- decision_criteria(fit, theta_threshold = 0.40,
                            probability_threshold = 0.80)
  expect_false(dc$go)
  expect_lt(dc$probability, 0.80)
})

test_that("direction = 'less' computes correct probability", {
  fit <- .mock_bf_dc(p_true = 0.10)
  dc  <- decision_criteria(fit, theta_threshold = 0.30, direction = "less",
                            probability_threshold = 0.80)
  expect_true(dc$go)
  expect_gt(dc$probability, 0.80)
})

test_that("probability is in [0, 1]", {
  fit <- .mock_bf_dc()
  dc  <- decision_criteria(fit, theta_threshold = 0.25)
  expect_gte(dc$probability, 0)
  expect_lte(dc$probability, 1)
})

test_that("go is TRUE exactly when probability >= probability_threshold", {
  fit  <- .mock_bf_dc(p_true = 0.35)
  dc   <- decision_criteria(fit, theta_threshold = 0.20,
                             probability_threshold = 0.80)
  expect_equal(dc$go, dc$probability >= dc$probability_threshold)
})

test_that("custom parameter name is used", {
  fit <- .mock_bf_dc()
  dc  <- decision_criteria(fit, theta_threshold = 0.0,
                            parameter = "theta_curr")
  expect_equal(dc$parameter, "theta_curr")
})

test_that("error on invalid parameter name", {
  fit <- .mock_bf_dc()
  expect_error(
    decision_criteria(fit, theta_threshold = 0.25, parameter = "nonexistent_var"),
    regexp = "not found"
  )
})

test_that("error on non-borrowing_fit input", {
  expect_error(decision_criteria(list(), theta_threshold = 0.25), regexp = "borrowing_fit")
})

test_that("posterior_ci has lo and hi", {
  fit <- .mock_bf_dc()
  dc  <- decision_criteria(fit, theta_threshold = 0.25)
  expect_named(dc$posterior_ci, c("lo", "hi"))
  expect_lt(dc$posterior_ci["lo"], dc$posterior_ci["hi"])
})

test_that("n_draws matches actual draws", {
  fit <- .mock_bf_dc(n_draws = 400L)
  dc  <- decision_criteria(fit, theta_threshold = 0.25)
  expect_equal(dc$n_draws, 400L)
})

test_that("print.decision_result produces Go/No-Go output", {
  fit <- .mock_bf_dc(p_true = 0.50)
  dc  <- decision_criteria(fit, theta_threshold = 0.20)
  expect_output(print(dc), "GO|NO-GO")
  expect_output(print(dc), "Posterior prob")
})
