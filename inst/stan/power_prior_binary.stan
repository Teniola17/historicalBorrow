// Power Prior — Binary Outcome
// Historical data raised to power a0 (fixed or estimated from data).
// When a0 = 1 this reduces to full pooling; a0 = 0 to no borrowing.

data {
  int<lower=1> J_hist;
  array[J_hist] int<lower=1> n_hist;
  array[J_hist] int<lower=0> r_hist;
  int<lower=0> n_curr;
  int<lower=0> r_curr;
  // Set to 0 to estimate a0; set to 1 and pass a fixed a0 in the model block
  int<lower=0, upper=1> a0_fixed_flag;
  real<lower=0, upper=1> a0_fixed_val;  // used only when a0_fixed_flag == 1
  real<lower=0> a0_alpha;
  real<lower=0> a0_beta;
}

parameters {
  real theta_curr;
  real<lower=0, upper=1> a0_raw;  // estimated only when a0_fixed_flag == 0
}

transformed parameters {
  real<lower=0, upper=1> a0;
  a0 = a0_fixed_flag == 1 ? a0_fixed_val : a0_raw;
}

model {
  theta_curr ~ normal(0, 2.5);
  if (a0_fixed_flag == 0)
    a0_raw ~ beta(a0_alpha, a0_beta);
  // Historical contribution scaled by a0
  for (j in 1:J_hist)
    target += a0 * binomial_logit_lpmf(r_hist[j] | n_hist[j], theta_curr);
  // Current trial
  r_curr ~ binomial_logit(n_curr, theta_curr);
}

generated quantities {
  real p_curr;
  p_curr = inv_logit(theta_curr);
  real log_lik;
  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);
}
