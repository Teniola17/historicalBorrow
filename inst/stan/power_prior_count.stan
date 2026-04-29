// Power Prior — Count (Poisson) Outcome

data {
  int<lower=1> J_hist;
  array[J_hist] int<lower=0> events_hist;
  vector[J_hist] log_n_hist;
  int<lower=0> events_curr;
  real log_n_curr;
  int<lower=0, upper=1> a0_fixed_flag;
  real<lower=0, upper=1> a0_fixed_val;
  real<lower=0> a0_alpha;
  real<lower=0> a0_beta;
}

parameters {
  real log_lambda_curr;
  real<lower=0, upper=1> a0_raw;
}

transformed parameters {
  real<lower=0, upper=1> a0;
  a0 = a0_fixed_flag == 1 ? a0_fixed_val : a0_raw;
}

model {
  log_lambda_curr ~ normal(0, 2.5);
  if (a0_fixed_flag == 0)
    a0_raw ~ beta(a0_alpha, a0_beta);
  for (j in 1:J_hist)
    target += a0 * poisson_log_lpmf(events_hist[j] | log_n_hist[j] + log_lambda_curr);
  events_curr ~ poisson_log(log_n_curr + log_lambda_curr);
}

generated quantities {
  real lambda_curr;
  lambda_curr = exp(log_lambda_curr);
  int events_rep;
  events_rep = poisson_log_rng(log_n_curr + log_lambda_curr);
  real log_lik;
  log_lik = poisson_log_lpmf(events_curr | log_n_curr + log_lambda_curr);
}
