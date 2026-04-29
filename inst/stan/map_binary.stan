// MAP Prior — Binary Outcome
// Binomial-logit hierarchical model with non-centered parameterisation.
// mu_pred is the predictive logit(p) for a new (exchangeable) study.

data {
  int<lower=2> J;
  array[J] int<lower=1> n;
  array[J] int<lower=0> r;
}

parameters {
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}

transformed parameters {
  vector[J] theta;
  theta = mu + tau * eta;
}

model {
  mu  ~ normal(0, 2.5);
  tau ~ normal(0, 0.5);
  eta ~ std_normal();
  for (j in 1:J)
    r[j] ~ binomial_logit(n[j], theta[j]);
}

generated quantities {
  real mu_pred;
  mu_pred = normal_rng(mu, tau);
  real p_pred;
  p_pred = inv_logit(mu_pred);
  vector[J] log_lik;
  for (j in 1:J)
    log_lik[j] = binomial_logit_lpmf(r[j] | n[j], theta[j]);
}
