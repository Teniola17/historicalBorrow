// Borrowing Model — Continuous Outcome
// MAP/RMAP prior as K-component mixture of normals on the natural scale.

data {
  int<lower=1> K;
  simplex[K] mix_weights;
  vector[K] mix_means;
  vector<lower=0>[K] mix_sds;
  int<lower=1> n_curr;
  real y_curr;
  real<lower=0> se_curr;
}

parameters {
  real theta_curr;
}

model {
  if (K == 1) {
    theta_curr ~ normal(mix_means[1], mix_sds[1]);
  } else {
    vector[K] lp;
    for (k in 1:K)
      lp[k] = log(mix_weights[k]) +
               normal_lpdf(theta_curr | mix_means[k], mix_sds[k]);
    target += log_sum_exp(lp);
  }
  y_curr ~ normal(theta_curr, se_curr);
}

generated quantities {
  real y_rep;
  y_rep = normal_rng(theta_curr, se_curr);
  real log_lik;
  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);
}
