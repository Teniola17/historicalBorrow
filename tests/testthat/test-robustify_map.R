# Mock map_prior to avoid Stan compilation in CI
.make_mock_map <- function(outcome = "binary", K = 2L) {
  mix <- list(
    weights = rep(1 / K, K),
    means   = seq(-0.5, 0.5, length.out = K),
    sds     = rep(0.4, K),
    K       = K
  )
  structure(
    list(
      mixture            = mix,
      mu_pred_draws      = rnorm(400, -0.3, 0.4),
      outcome_type       = outcome,
      model_type         = "standard",
      hyperpriors        = list(),
      ess_prior          = 25.0,
      stan_fit           = NULL,
      draws              = NULL,
      data               = NULL,
      mixture_components = K
    ),
    class = "map_prior"
  )
}

test_that("robustify_map returns rmap_prior inheriting map_prior", {
  mp  <- .make_mock_map()
  rmp <- robustify_map(mp, weight = 0.8)
  expect_s3_class(rmp, "rmap_prior")
  expect_s3_class(rmp, "map_prior")
})

test_that("mixture_rmap has K + 1 components", {
  mp  <- .make_mock_map(K = 2L)
  rmp <- robustify_map(mp, weight = 0.8)
  expect_equal(rmp$mixture_rmap$K, 3L)
})

test_that("mixture_rmap weights sum to 1", {
  mp  <- .make_mock_map()
  rmp <- robustify_map(mp, weight = 0.7)
  expect_equal(sum(rmp$mixture_rmap$weights), 1.0, tolerance = 1e-10)
})

test_that("vague component weight equals 1 - weight", {
  mp  <- .make_mock_map(K = 2L)
  rmp <- robustify_map(mp, weight = 0.8)
  expect_equal(rmp$mixture_rmap$weights[3], 0.2, tolerance = 1e-10)
})

test_that("MAP component weights sum to weight", {
  mp  <- .make_mock_map(K = 3L)
  rmp <- robustify_map(mp, weight = 0.75)
  expect_equal(sum(rmp$mixture_rmap$weights[1:3]), 0.75, tolerance = 1e-10)
})

test_that("robustify_map errors on weight = 0", {
  mp <- .make_mock_map()
  expect_error(robustify_map(mp, weight = 0), regexp = "lower")
})

test_that("robustify_map errors on weight = 1", {
  mp <- .make_mock_map()
  expect_error(robustify_map(mp, weight = 1), regexp = "upper")
})

test_that("robustify_map errors on weight > 1", {
  mp <- .make_mock_map()
  expect_error(robustify_map(mp, weight = 1.5))
})

test_that("custom vague prior mean and sd are used", {
  mp  <- .make_mock_map(K = 1L)
  rmp <- robustify_map(mp, weight = 0.8,
                       vague_prior = list(mean = 1.0, sd = 5.0))
  expect_equal(rmp$mixture_rmap$means[2], 1.0)
  expect_equal(rmp$mixture_rmap$sds[2],   5.0)
})

test_that("ess_prior_rmap is weight * ess_prior", {
  mp  <- .make_mock_map()
  rmp <- robustify_map(mp, weight = 0.8)
  expect_equal(rmp$ess_prior_rmap, 0.8 * mp$ess_prior, tolerance = 0.01)
})

test_that("ess_prior_rmap is less than ess_prior", {
  mp  <- .make_mock_map()
  rmp <- robustify_map(mp, weight = 0.8)
  expect_lt(rmp$ess_prior_rmap, mp$ess_prior)
})

test_that("rmap_weight slot is stored correctly", {
  mp  <- .make_mock_map()
  rmp <- robustify_map(mp, weight = 0.65)
  expect_equal(rmp$rmap_weight, 0.65)
})

test_that("original MAP mixture is unchanged", {
  mp_orig <- .make_mock_map()
  rmp     <- robustify_map(mp_orig, weight = 0.8)
  expect_equal(rmp$mixture$K,       mp_orig$mixture$K)
  expect_equal(rmp$mixture$weights, mp_orig$mixture$weights)
})

test_that("robustify_map errors on non-map_prior input", {
  expect_error(robustify_map(list(a = 1), weight = 0.8), regexp = "map_prior")
})
