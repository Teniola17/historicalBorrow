# data-raw/simulate_example_data.R
# Generates `example_trial_data` saved in data/example_trial_data.rda
# Run once: source("data-raw/simulate_example_data.R")

set.seed(20240429)

n_hist <- 6L

# ── Binary outcome ────────────────────────────────────────────────────────────
true_p_hist <- 0.30
true_p_curr <- 0.32  # mild drift
tau_bin     <- 0.20  # between-study SD on logit scale

mu_true     <- qlogis(true_p_hist)
logit_p_i   <- rnorm(n_hist, mean = mu_true, sd = tau_bin)
p_i         <- plogis(logit_p_i)
n_bin_hist  <- c(100L, 120L, 90L, 110L, 95L, 105L)

hist_binary <- data.frame(
  study  = paste0("Historical_", seq_len(n_hist)),
  n      = n_bin_hist,
  r      = as.integer(rbinom(n_hist, size = n_bin_hist, prob = p_i)),
  true_p = round(p_i, 4)
)

curr_binary <- data.frame(
  study = "Current",
  n     = 80L,
  r     = as.integer(rbinom(1L, 80L, true_p_curr))
)

# ── Continuous outcome ────────────────────────────────────────────────────────
true_mean_hist <- 5.0
true_mean_curr <- 5.3
true_sd_cont   <- 1.1
tau_cont       <- 0.30

theta_cont_i  <- rnorm(n_hist, mean = true_mean_hist, sd = tau_cont)
n_cont_hist   <- c(50L, 60L, 55L, 45L, 58L, 52L)

hist_continuous <- data.frame(
  study = paste0("Historical_", seq_len(n_hist)),
  n     = n_cont_hist,
  mean  = round(rnorm(n_hist, theta_cont_i, true_sd_cont / sqrt(n_cont_hist)), 3),
  sd    = round(abs(rnorm(n_hist, true_sd_cont, 0.1)), 3)
)
hist_continuous$se <- round(hist_continuous$sd / sqrt(hist_continuous$n), 4)

curr_continuous <- data.frame(
  study = "Current",
  n     = 40L,
  mean  = round(rnorm(1L, true_mean_curr, true_sd_cont / sqrt(40)), 3),
  sd    = round(abs(rnorm(1L, true_sd_cont, 0.1)), 3)
)
curr_continuous$se <- round(curr_continuous$sd / sqrt(curr_continuous$n), 4)

# ── Count outcome ─────────────────────────────────────────────────────────────
true_rate_hist <- 0.060
true_rate_curr <- 0.065
tau_log        <- 0.25

log_lambda_i  <- rnorm(n_hist, log(true_rate_hist), tau_log)
person_time   <- c(200.0, 180.0, 220.0, 190.0, 210.0, 195.0)

hist_count <- data.frame(
  study      = paste0("Historical_", seq_len(n_hist)),
  n          = person_time,
  events     = as.integer(rpois(n_hist, person_time * exp(log_lambda_i))),
  true_rate  = round(exp(log_lambda_i), 5)
)

curr_count <- data.frame(
  study  = "Current",
  n      = 150.0,
  events = as.integer(rpois(1L, 150.0 * true_rate_curr))
)

# ── Bundle and save ───────────────────────────────────────────────────────────
example_trial_data <- list(
  binary = list(
    hist    = hist_binary,
    current = curr_binary
  ),
  continuous = list(
    hist    = hist_continuous,
    current = curr_continuous
  ),
  count = list(
    hist    = hist_count,
    current = curr_count
  )
)

usethis::use_data(example_trial_data, overwrite = TRUE)
message("Saved example_trial_data to data/example_trial_data.rda")
