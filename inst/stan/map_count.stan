// MAP Prior — Count (Poisson) Outcome
// Log-normal hierarchical model with person-time offset (log_n).
// Non-centered parameterisation on the log-rate scale.

data {
  int<lower=2> J;
  array[J] int<lower=0> events;
  vector[J] log_n;
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}

transformed parameters {
  vector[J] log_lambda;
  log_lambda = mu + tau * eta;
}

model {
  mu  ~ normal(0, 2.5);
  tau ~ normal(0, 0.5);
  eta ~ std_normal();
  for (j in 1:J)
    events[j] ~ poisson_log(log_n[j] + log_lambda[j]);
}

generated quantities {
  real mu_pred;
  mu_pred = normal_rng(mu, tau);
  real lambda_pred;
  lambda_pred = exp(mu_pred);
  vector[J] log_lik;
  for (j in 1:J)
    log_lik[j] = poisson_log_lpmf(events[j] | log_n[j] + log_lambda[j]);
}
