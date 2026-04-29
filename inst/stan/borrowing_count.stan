// Borrowing Model — Count (Poisson) Outcome
// MAP/RMAP prior as K-component mixture of normals on the log-rate scale.

data {
  int<lower=1> K;
  simplex[K] mix_weights;
  vector[K] mix_means;
  vector<lower=0>[K] mix_sds;
  int<lower=0> events_curr;
  real log_n_curr;
}

parameters {
  real log_lambda_curr;
}

model {
  if (K == 1) {
    log_lambda_curr ~ normal(mix_means[1], mix_sds[1]);
  } else {
    vector[K] lp;
    for (k in 1:K)
      lp[k] = log(mix_weights[k]) +
               normal_lpdf(log_lambda_curr | mix_means[k], mix_sds[k]);
    target += log_sum_exp(lp);
  }
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
