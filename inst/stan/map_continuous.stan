// MAP Prior — Continuous Outcome
// Normal-normal hierarchical model with non-centered parameterisation.
// Sufficient statistics: study-level means y[j] with standard errors se[j].

data {
  int<lower=2> J;
  array[J] int<lower=1> n;
  vector[J] y;
  vector<lower=0>[J] se;
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
  mu  ~ normal(0, 100);
  tau ~ normal(0, 10);
  eta ~ std_normal();
  y   ~ normal(theta, se);
}

generated quantities {
  real mu_pred;
  mu_pred = normal_rng(mu, tau);
  vector[J] log_lik;
  for (j in 1:J)
    log_lik[j] = normal_lpdf(y[j] | theta[j], se[j]);
}
