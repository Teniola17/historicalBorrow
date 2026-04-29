test_that("prepare_data returns bb_data for binary outcome", {
  hist <- data.frame(
    study = paste0("S", 1:4),
    n     = c(100L, 120L, 90L, 110L),
    r     = c(30L, 40L, 25L, 38L)
  )
  curr <- data.frame(study = "C", n = 80L, r = 22L)
  pd   <- prepare_data(hist, curr, outcome = "binary")

  expect_s3_class(pd, "bb_data")
  expect_equal(pd$outcome, "binary")
  expect_equal(pd$n_hist, 4L)
  expect_equal(pd$stan_data_map$J, 4L)
  expect_equal(pd$stan_data_map$n, c(100L, 120L, 90L, 110L))
  expect_equal(pd$stan_data_map$r, c(30L, 40L, 25L, 38L))
  expect_equal(pd$stan_data_current$n_curr, 80L)
  expect_equal(pd$stan_data_current$r_curr, 22L)
})

test_that("prepare_data returns bb_data for continuous outcome", {
  hist <- data.frame(
    study = paste0("S", 1:3),
    n     = c(50L, 60L, 55L),
    mean  = c(5.2, 4.8, 5.5),
    sd    = c(1.1, 0.9, 1.3)
  )
  curr <- data.frame(study = "C", n = 45L, mean = 5.1, sd = 1.0)
  pd   <- prepare_data(hist, curr, outcome = "continuous")

  expect_s3_class(pd, "bb_data")
  expect_equal(pd$outcome, "continuous")
  expect_true("se" %in% names(pd$hist_data))
  expect_equal(pd$hist_data$se, hist$sd / sqrt(hist$n), tolerance = 1e-10)
})

test_that("prepare_data returns bb_data for count outcome", {
  hist <- data.frame(
    study  = paste0("S", 1:3),
    n      = c(200.0, 180.0, 220.0),
    events = c(12L, 10L, 15L)
  )
  curr <- data.frame(study = "C", n = 150.0, events = NA_integer_)
  pd   <- prepare_data(hist, curr, outcome = "count")

  expect_s3_class(pd, "bb_data")
  expect_equal(pd$stan_data_map$log_n, log(c(200, 180, 220)), tolerance = 1e-10)
})

test_that("prepare_data errors when r > n", {
  hist <- data.frame(study = "S1", n = 10L, r = 15L)
  curr <- data.frame(study = "C",  n = 50L, r = NA_integer_)
  expect_error(prepare_data(hist, curr, outcome = "binary"), regexp = "r > n")
})

test_that("prepare_data errors on missing required columns", {
  hist <- data.frame(study = "S1", n = 100L)  # missing 'r'
  curr <- data.frame(study = "C",  n = 50L)
  expect_error(
    prepare_data(hist, curr, outcome = "binary"),
    regexp = "missing required columns"
  )
})

test_that("prepare_data accepts custom study column name", {
  hist <- data.frame(
    trial = paste0("T", 1:3),
    n     = c(100L, 110L, 90L),
    r     = c(20L, 25L, 18L)
  )
  curr <- data.frame(trial = "Current", n = 80L, r = 10L)
  pd   <- prepare_data(hist, curr, outcome = "binary", study_col = "trial")
  expect_equal(names(pd$hist_data)[1], "study")
})

test_that("prepare_data requires at least 2 historical studies", {
  hist <- data.frame(study = "S1", n = 100L, r = 30L)
  curr <- data.frame(study = "C",  n = 80L)
  expect_error(prepare_data(hist, curr, outcome = "binary"))
})

test_that("prepare_data errors on non-positive n for continuous", {
  hist <- data.frame(study = "S1", n = 0L, mean = 5.0, sd = 1.0)
  curr <- data.frame(study = "C",  n = 40L)
  expect_error(
    prepare_data(
      rbind(hist, data.frame(study = "S2", n = 50L, mean = 5.0, sd = 1.0)),
      curr, outcome = "continuous"
    ),
    regexp = "positive"
  )
})

test_that("print.bb_data produces output", {
  hist <- data.frame(study = paste0("S", 1:3), n = c(100L, 110L, 90L), r = c(20L, 25L, 18L))
  curr <- data.frame(study = "C", n = 80L, r = 10L)
  pd   <- prepare_data(hist, curr, "binary")
  expect_output(print(pd), "historicalBorrow")
  expect_output(print(pd), "binary")
})

test_that("summary.bb_data produces output", {
  hist <- data.frame(study = paste0("S", 1:3), n = c(100L, 110L, 90L), r = c(20L, 25L, 18L))
  curr <- data.frame(study = "C", n = 80L, r = 10L)
  pd   <- prepare_data(hist, curr, "binary")
  expect_output(summary(pd), "Historical")
})

test_that("prepare_data handles negative events in count data", {
  hist <- data.frame(
    study  = paste0("S", 1:2),
    n      = c(200.0, 180.0),
    events = c(-1L, 10L)
  )
  curr <- data.frame(study = "C", n = 150.0)
  expect_error(prepare_data(hist, curr, outcome = "count"), regexp = "non-negative")
})
