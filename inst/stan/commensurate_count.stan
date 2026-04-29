// Commensurate Prior — Count (Poisson) Outcome

data {
  int<lower=0> events_curr;
  real log_n_curr;
}

parameters {
  real log_lambda_hist;
  real log_lambda_curr;
  real<lower=0> tau_comm;
}

model {
  log_lambda_hist ~ normal(0, 2.5);
  tau_comm        ~ cauchy(0, 1);
  log_lambda_curr ~ normal(log_lambda_hist, tau_comm);
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
