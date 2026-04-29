// Borrowing Model — Binary Outcome
// MAP/RMAP prior encoded as a K-component mixture of normals on the logit scale.
// log_sum_exp ensures numerical stability for K > 1.

data {
  int<lower=1> K;
  simplex[K] mix_weights;
  vector[K] mix_means;
  vector<lower=0>[K] mix_sds;
  int<lower=0> n_curr;
  int<lower=0> r_curr;
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
  r_curr ~ binomial_logit(n_curr, theta_curr);
}

generated quantities {
  real p_curr;
  p_curr = inv_logit(theta_curr);
  int r_rep;
  r_rep = binomial_rng(n_curr, p_curr);
  real log_lik;
  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);
}
