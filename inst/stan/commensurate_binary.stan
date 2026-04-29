// Commensurate Prior — Binary Outcome
// tau_comm controls the degree of borrowing:
//   small tau_comm => theta_curr close to theta_hist => strong borrowing
//   large tau_comm => weak borrowing, posterior driven by current data
// Cauchy(0,1) hyperprior on tau_comm allows the data to choose.

data {
  int<lower=0> n_curr;
  int<lower=0> r_curr;
}

parameters {
  real theta_hist;
  real theta_curr;
  real<lower=0> tau_comm;
}

model {
  theta_hist ~ normal(0, 2.5);
  tau_comm   ~ cauchy(0, 1);
  theta_curr ~ normal(theta_hist, tau_comm);
  r_curr ~ binomial_logit(n_curr, theta_curr);
}

generated quantities {
  real p_curr;
  p_curr = inv_logit(theta_curr);
  real log_lik;
  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);
}
