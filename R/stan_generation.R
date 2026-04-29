#' Generate Stan code for a given model configuration
#'
#' This is an internal function. Users do not call it directly; it is invoked
#' by [fit_map()] and [fit_borrowing_model()].
#'
#' @param outcome_type One of `"binary"`, `"continuous"`, `"count"`.
#' @param model_structure One of `"map_standard"`, `"borrowing_map"`,
#'   `"borrowing_rmap"`, `"borrowing_power"`, `"borrowing_commensurate"`.
#' @param prior_spec Named list of hyperprior values. If a value is missing, a
#'   sensible default is used.
#'
#' @return A character string of valid Stan code.
#' @noRd
.generate_stan_code <- function(outcome_type,
                                model_structure,
                                prior_spec = list()) {

  outcome_type    <- match.arg(outcome_type,    c("binary", "continuous", "count"))
  model_structure <- match.arg(model_structure, c(
    "map_standard",
    "borrowing_map",
    "borrowing_rmap",
    "borrowing_power",
    "borrowing_commensurate"
  ))

  switch(model_structure,
    map_standard          = .stan_map(outcome_type, prior_spec),
    borrowing_map         = ,
    borrowing_rmap        = .stan_borrowing_mixture(outcome_type, prior_spec),
    borrowing_power       = .stan_power_prior(outcome_type, prior_spec),
    borrowing_commensurate = .stan_commensurate(outcome_type, prior_spec)
  )
}

# ── Defaults ──────────────────────────────────────────────────────────────────

.default <- function(prior_spec, key, default) {
  val <- prior_spec[[key]]
  if (is.null(val)) default else val
}

# ── MAP hierarchical model ─────────────────────────────────────────────────────

.stan_map <- function(outcome_type, prior_spec) {
  mu_mean   <- .default(prior_spec, "mu_mean",   0)
  mu_sd     <- .default(prior_spec, "mu_sd",     2.5)
  tau_scale <- .default(prior_spec, "tau_scale", 0.5)

  data_block    <- .map_data_block(outcome_type)
  likelihood    <- .map_likelihood(outcome_type)
  log_lik_block <- .map_log_lik(outcome_type)
  pred_param    <- .map_pred_param(outcome_type)

  glue::glue(
'// MAP hierarchical model — {outcome_type} outcome
// Non-centered parameterisation to avoid Neal\'s funnel
data \\{
{data_block}
\\}

parameters \\{
  real mu;
  real<lower=0> tau;
  vector[J] eta;
\\}

transformed parameters \\{
  vector[J] theta;
  theta = mu + tau * eta;
\\}

model \\{
  mu  ~ normal({mu_mean}, {mu_sd});
  tau ~ normal(0, {tau_scale});
  eta ~ std_normal();
{likelihood}
\\}

generated quantities \\{
  real mu_pred;
  mu_pred = normal_rng(mu, tau);
  {pred_param}
  vector[J] log_lik;
{log_lik_block}
\\}
'
  )
}

.map_data_block <- function(outcome_type) {
  switch(outcome_type,
    binary = "  int<lower=2> J;\n  array[J] int<lower=1> n;\n  array[J] int<lower=0> r;",
    continuous = "  int<lower=2> J;\n  array[J] int<lower=1> n;\n  vector[J] y;\n  vector<lower=0>[J] se;",
    count = "  int<lower=2> J;\n  array[J] int<lower=0> events;\n  vector[J] log_n;"
  )
}

.map_likelihood <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  for (j in 1:J)\n    r[j] ~ binomial_logit(n[j], theta[j]);",
    continuous = "  y ~ normal(theta, se);",
    count      = "  for (j in 1:J)\n    events[j] ~ poisson_log(log_n[j] + theta[j]);"
  )
}

.map_log_lik <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  for (j in 1:J)\n    log_lik[j] = binomial_logit_lpmf(r[j] | n[j], theta[j]);",
    continuous = "  for (j in 1:J)\n    log_lik[j] = normal_lpdf(y[j] | theta[j], se[j]);",
    count      = "  for (j in 1:J)\n    log_lik[j] = poisson_log_lpmf(events[j] | log_n[j] + theta[j]);"
  )
}

.map_pred_param <- function(outcome_type) {
  switch(outcome_type,
    binary     = "real p_pred;\n  p_pred = inv_logit(mu_pred);",
    continuous = "",  # mu_pred is already on natural scale
    count      = "real lambda_pred;\n  lambda_pred = exp(mu_pred);"
  )
}

# ── Borrowing model with mixture prior ────────────────────────────────────────

.stan_borrowing_mixture <- function(outcome_type, prior_spec) {
  data_block <- .borrow_data_block(outcome_type)
  likelihood <- .borrow_likelihood(outcome_type)
  gq_block   <- .borrow_gq(outcome_type)
  param_name <- .borrow_param_name(outcome_type)

  glue::glue(
'// Borrowing model — {outcome_type} outcome
// MAP/RMAP prior encoded as mixture of normals (K components passed as data)
// Numerical stability via log_sum_exp
data \\{
  // MAP/RMAP mixture prior
  int<lower=1> K;
  simplex[K] mix_weights;
  vector[K] mix_means;
  vector<lower=0>[K] mix_sds;
  // Current trial data
{data_block}
\\}

parameters \\{
  real {param_name};
\\}

model \\{
  // Mixture prior
  if (K == 1) \\{
    {param_name} ~ normal(mix_means[1], mix_sds[1]);
  \\} else \\{
    vector[K] lp;
    for (k in 1:K)
      lp[k] = log(mix_weights[k]) +
               normal_lpdf({param_name} | mix_means[k], mix_sds[k]);
    target += log_sum_exp(lp);
  \\}
  // Likelihood
{likelihood}
\\}

generated quantities \\{
{gq_block}
\\}
'
  )
}

.borrow_param_name <- function(outcome_type) {
  switch(outcome_type,
    binary     = "theta_curr",
    continuous = "theta_curr",
    count      = "log_lambda_curr"
  )
}

.borrow_data_block <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  int<lower=0> n_curr;\n  int<lower=0> r_curr;",
    continuous = "  int<lower=1> n_curr;\n  real y_curr;\n  real<lower=0> se_curr;",
    count      = "  int<lower=0> events_curr;\n  real log_n_curr;"
  )
}

.borrow_likelihood <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  r_curr ~ binomial_logit(n_curr, theta_curr);",
    continuous = "  y_curr ~ normal(theta_curr, se_curr);",
    count      = "  events_curr ~ poisson_log(log_n_curr + log_lambda_curr);"
  )
}

.borrow_gq <- function(outcome_type) {
  switch(outcome_type,
    binary = paste0(
      "  real p_curr;\n",
      "  p_curr = inv_logit(theta_curr);\n",
      "  int r_rep;\n",
      "  r_rep = binomial_rng(n_curr, p_curr);\n",
      "  real log_lik;\n",
      "  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);"
    ),
    continuous = paste0(
      "  real y_rep;\n",
      "  y_rep = normal_rng(theta_curr, se_curr);\n",
      "  real log_lik;\n",
      "  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);"
    ),
    count = paste0(
      "  real lambda_curr;\n",
      "  lambda_curr = exp(log_lambda_curr);\n",
      "  int events_rep;\n",
      "  events_rep = poisson_log_rng(log_n_curr + log_lambda_curr);\n",
      "  real log_lik;\n",
      "  log_lik = poisson_log_lpmf(events_curr | log_n_curr + log_lambda_curr);"
    )
  )
}

# ── Power prior model ─────────────────────────────────────────────────────────

.stan_power_prior <- function(outcome_type, prior_spec) {
  a0_fixed  <- .default(prior_spec, "a0_fixed",  -1)    # -1 = estimate a0
  a0_alpha  <- .default(prior_spec, "a0_alpha",   1)
  a0_beta   <- .default(prior_spec, "a0_beta",    1)
  mu_mean   <- .default(prior_spec, "mu_mean",    0)
  mu_sd     <- .default(prior_spec, "mu_sd",      2.5)

  random_a0 <- (a0_fixed < 0)

  data_block <- .power_data_block(outcome_type)
  hist_ll    <- .power_hist_log_lik(outcome_type)
  curr_ll    <- .power_curr_likelihood(outcome_type)
  gq_block   <- .power_gq(outcome_type)
  param_name <- .borrow_param_name(outcome_type)

  a0_param_block <- if (random_a0) {
    "  real<lower=0, upper=1> a0;"
  } else {
    ""
  }

  a0_prior <- if (random_a0) {
    glue::glue("  a0 ~ beta({a0_alpha}, {a0_beta});")
  } else {
    glue::glue("  // a0 fixed at {a0_fixed}")
  }

  a0_value <- if (random_a0) "a0" else as.character(a0_fixed)

  glue::glue(
'// Power prior model — {outcome_type} outcome
// Historical likelihood raised to power a0 (0 <= a0 <= 1)
data \\{
{data_block}
\\}

parameters \\{
  real {param_name};
  {a0_param_block}
\\}

model \\{
  // Weakly informative prior on treatment parameter
  {param_name} ~ normal({mu_mean}, {mu_sd});
{a0_prior}
  // Historical data contribution scaled by a0
  target += {a0_value} * ({hist_ll});
  // Current data likelihood
{curr_ll}
\\}

generated quantities \\{
{gq_block}
\\}
'
  )
}

.power_data_block <- function(outcome_type) {
  switch(outcome_type,
    binary = paste0(
      "  // Historical\n",
      "  int<lower=1> J_hist;\n",
      "  array[J_hist] int<lower=1> n_hist;\n",
      "  array[J_hist] int<lower=0> r_hist;\n",
      "  // Current\n",
      "  int<lower=0> n_curr;\n",
      "  int<lower=0> r_curr;"
    ),
    continuous = paste0(
      "  // Historical\n",
      "  int<lower=1> J_hist;\n",
      "  vector[J_hist] y_hist;\n",
      "  vector<lower=0>[J_hist] se_hist;\n",
      "  // Current\n",
      "  int<lower=1> n_curr;\n",
      "  real y_curr;\n",
      "  real<lower=0> se_curr;"
    ),
    count = paste0(
      "  // Historical\n",
      "  int<lower=1> J_hist;\n",
      "  array[J_hist] int<lower=0> events_hist;\n",
      "  vector[J_hist] log_n_hist;\n",
      "  // Current\n",
      "  int<lower=0> events_curr;\n",
      "  real log_n_curr;"
    )
  )
}

.power_hist_log_lik <- function(outcome_type) {
  switch(outcome_type,
    binary     = "sum(binomial_logit_lpmf(r_hist | n_hist, rep_vector(theta_curr, J_hist)))",
    continuous = "normal_lpdf(y_hist | theta_curr, se_hist)",
    count      = "sum(poisson_log_lpmf(events_hist | log_n_hist + log_lambda_curr))"
  )
}

.power_curr_likelihood <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  r_curr ~ binomial_logit(n_curr, theta_curr);",
    continuous = "  y_curr ~ normal(theta_curr, se_curr);",
    count      = "  events_curr ~ poisson_log(log_n_curr + log_lambda_curr);"
  )
}

.power_gq <- function(outcome_type) {
  switch(outcome_type,
    binary = paste0(
      "  real p_curr;\n",
      "  p_curr = inv_logit(theta_curr);\n",
      "  real log_lik;\n",
      "  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);"
    ),
    continuous = paste0(
      "  real log_lik;\n",
      "  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);"
    ),
    count = paste0(
      "  real lambda_curr;\n",
      "  lambda_curr = exp(log_lambda_curr);\n",
      "  real log_lik;\n",
      "  log_lik = poisson_log_lpmf(events_curr | log_n_curr + log_lambda_curr);"
    )
  )
}

# ── Commensurate prior model ──────────────────────────────────────────────────

.stan_commensurate <- function(outcome_type, prior_spec) {
  mu_hist  <- .default(prior_spec, "mu_hist",    0)
  sig_hist <- .default(prior_spec, "sig_hist",   2.5)
  mu_mean  <- .default(prior_spec, "mu_mean",    0)
  mu_sd    <- .default(prior_spec, "mu_sd",      2.5)

  data_block <- .comm_data_block(outcome_type)
  curr_ll    <- .comm_curr_likelihood(outcome_type)
  gq_block   <- .comm_gq(outcome_type)

  param_hist <- paste0(.borrow_param_name(outcome_type), "_hist")
  param_curr <- .borrow_param_name(outcome_type)

  glue::glue(
'// Commensurate prior model — {outcome_type} outcome
// tau_comm controls borrowing: small tau_comm => strong borrowing
data \\{
{data_block}
\\}

parameters \\{
  real {param_hist};
  real {param_curr};
  real<lower=0> tau_comm;
\\}

model \\{
  // Prior on historical parameter
  {param_hist} ~ normal({mu_hist}, {sig_hist});
  // Commensurate prior: Cauchy hyperprior on tau_comm
  tau_comm ~ cauchy(0, 1);
  // Current parameter borrows from historical via commensurability
  {param_curr} ~ normal({param_hist}, tau_comm);
  // Likelihood
{curr_ll}
\\}

generated quantities \\{
{gq_block}
  real log_lik_tau_comm;
  log_lik_tau_comm = cauchy_lpdf(tau_comm | 0, 1);
\\}
'
  )
}

.comm_data_block <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  int<lower=0> n_curr;\n  int<lower=0> r_curr;",
    continuous = "  int<lower=1> n_curr;\n  real y_curr;\n  real<lower=0> se_curr;",
    count      = "  int<lower=0> events_curr;\n  real log_n_curr;"
  )
}

.comm_curr_likelihood <- function(outcome_type) {
  switch(outcome_type,
    binary     = "  r_curr ~ binomial_logit(n_curr, theta_curr);",
    continuous = "  y_curr ~ normal(theta_curr, se_curr);",
    count      = "  events_curr ~ poisson_log(log_n_curr + log_lambda_curr);"
  )
}

.comm_gq <- function(outcome_type) {
  switch(outcome_type,
    binary = paste0(
      "  real p_curr;\n",
      "  p_curr = inv_logit(theta_curr);\n",
      "  real log_lik;\n",
      "  log_lik = binomial_logit_lpmf(r_curr | n_curr, theta_curr);"
    ),
    continuous = paste0(
      "  real log_lik;\n",
      "  log_lik = normal_lpdf(y_curr | theta_curr, se_curr);"
    ),
    count = paste0(
      "  real lambda_curr;\n",
      "  lambda_curr = exp(log_lambda_curr);\n",
      "  real log_lik;\n",
      "  log_lik = poisson_log_lpmf(events_curr | log_n_curr + log_lambda_curr);"
    )
  )
}
