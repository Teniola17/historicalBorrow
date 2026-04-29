test_that(".generate_stan_code returns character for all outcome/structure combos", {
  outcomes   <- c("binary", "continuous", "count")
  structures <- c("map_standard", "borrowing_map", "borrowing_rmap",
                  "borrowing_power", "borrowing_commensurate")

  for (ot in outcomes) {
    for (ms in structures) {
      code <- historicalBorrow:::.generate_stan_code(ot, ms, list())
      expect_type(code, "character", info = paste(ot, ms))
      expect_gt(nchar(code), 50L, info = paste(ot, ms))
    }
  }
})

test_that(".generate_stan_code MAP binary uses non-centred parameterisation", {
  code <- historicalBorrow:::.generate_stan_code("binary", "map_standard", list())
  expect_match(code, "eta")
  expect_match(code, "mu \\+ tau \\* eta")
  expect_match(code, "std_normal")
})

test_that(".generate_stan_code MAP continuous uses non-centred parameterisation", {
  code <- historicalBorrow:::.generate_stan_code("continuous", "map_standard", list())
  expect_match(code, "mu \\+ tau \\* eta")
})

test_that(".generate_stan_code MAP contains mu_pred in generated quantities", {
  for (ot in c("binary", "continuous", "count")) {
    code <- historicalBorrow:::.generate_stan_code(ot, "map_standard", list())
    expect_match(code, "mu_pred", info = ot)
  }
})

test_that(".generate_stan_code injects hyperprior values into MAP model", {
  code <- historicalBorrow:::.generate_stan_code(
    "binary", "map_standard",
    list(mu_mean = 1.5, mu_sd = 5.0, tau_scale = 0.3)
  )
  expect_match(code, "1.5")
  expect_match(code, "5\\.0")
  expect_match(code, "0\\.3")
})

test_that(".generate_stan_code borrowing contains log_sum_exp", {
  for (ot in c("binary", "continuous", "count")) {
    code <- historicalBorrow:::.generate_stan_code(ot, "borrowing_map", list())
    expect_match(code, "log_sum_exp", info = ot)
  }
})

test_that(".generate_stan_code borrowing contains K==1 branch", {
  code <- historicalBorrow:::.generate_stan_code("binary", "borrowing_map", list())
  expect_match(code, "K == 1")
})

test_that(".generate_stan_code borrowing_rmap produces same code as borrowing_map", {
  code_map  <- historicalBorrow:::.generate_stan_code("binary", "borrowing_map",  list())
  code_rmap <- historicalBorrow:::.generate_stan_code("binary", "borrowing_rmap", list())
  expect_equal(code_map, code_rmap)
})

test_that(".generate_stan_code power prior with fixed a0 injects literal value", {
  code <- historicalBorrow:::.generate_stan_code(
    "binary", "borrowing_power",
    list(a0_fixed = 0.7)
  )
  expect_match(code, "0\\.7")
  expect_false(grepl("real<lower=0, upper=1> a0", code))
})

test_that(".generate_stan_code power prior with random a0 has a0 parameter", {
  code <- historicalBorrow:::.generate_stan_code(
    "binary", "borrowing_power",
    list(a0_fixed = -1, a0_alpha = 1, a0_beta = 1)
  )
  expect_match(code, "real<lower=0, upper=1> a0")
  expect_match(code, "beta\\(")
})

test_that(".generate_stan_code commensurate prior has tau_comm parameter", {
  for (ot in c("binary", "continuous", "count")) {
    code <- historicalBorrow:::.generate_stan_code(ot, "borrowing_commensurate", list())
    expect_match(code, "tau_comm", info = ot)
    expect_match(code, "cauchy", info = ot)
  }
})

test_that(".generate_stan_code binary MAP generates log_lik in GQ", {
  code <- historicalBorrow:::.generate_stan_code("binary", "map_standard", list())
  expect_match(code, "log_lik")
  expect_match(code, "binomial_logit_lpmf")
})

test_that(".generate_stan_code count MAP uses poisson_log", {
  code <- historicalBorrow:::.generate_stan_code("count", "map_standard", list())
  expect_match(code, "poisson_log")
  expect_match(code, "log_n")
})

test_that(".generate_stan_code continuous borrowing uses normal_lpdf", {
  code <- historicalBorrow:::.generate_stan_code("continuous", "borrowing_map", list())
  expect_match(code, "normal_lpdf")
})
