// Commensurate Prior — Continuous Outcome

data {
  int<lower=1> n_curr;
  real y_curr;
  real<lower=0> se_curr;
}

parameters {
  real theta_hist;
  real theta_curr;
  real<lower=0> tau_comm;
}

model {
  theta_hist ~ normal(0, 100);
  tau_comm   ~ cauchy(0, 1);
  theta_curr ~ normal(theta_hist, tau_comm);
  y_curr ~ normal(theta_curr, se_curr);
}

generated quantities {
  real y_rep;
  y_rep = normal_rng(theta_curr, se_curr);
  real log_lik;
  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);
}
