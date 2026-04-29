// Power Prior — Continuous Outcome

data {
  int<lower=1> J_hist;
  vector[J_hist] y_hist;
  vector<lower=0>[J_hist] se_hist;
  int<lower=1> n_curr;
  real y_curr;
  real<lower=0> se_curr;
  int<lower=0, upper=1> a0_fixed_flag;
  real<lower=0, upper=1> a0_fixed_val;
  real<lower=0> a0_alpha;
  real<lower=0> a0_beta;
}

parameters {
  real theta_curr;
  real<lower=0, upper=1> a0_raw;
}

transformed parameters {
  real<lower=0, upper=1> a0;
  a0 = a0_fixed_flag == 1 ? a0_fixed_val : a0_raw;
}

model {
  theta_curr ~ normal(0, 100);
  if (a0_fixed_flag == 0)
    a0_raw ~ beta(a0_alpha, a0_beta);
  for (j in 1:J_hist)
    target += a0 * normal_lpdf(y_hist[j] | theta_curr, se_hist[j]);
  y_curr ~ normal(theta_curr, se_curr);
}

generated quantities {
  real y_rep;
  y_rep = normal_rng(theta_curr, se_curr);
  real log_lik;
  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);
}
