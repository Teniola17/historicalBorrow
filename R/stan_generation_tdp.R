#' Generate Stan code for Time-Dependent Prior (TDP) models
#'
#' Adapts the user-provided TDP Ornstein-Uhlenbeck model to different outcome
#' types. The core temporal structure remains the same; only the likelihood
#' changes (normal for continuous, binomial_logit for binary, poisson_log for count).
#'
#' @param outcome_type One of "binary", "continuous", "count".
#' @param hyperpriors Named list of hyperprior values.
#'
#' @return A character string of valid Stan code.
#' @noRd
.generate_stan_code_tdp <- function(outcome_type, hyperpriors = list()) {

  outcome_type <- match.arg(outcome_type, c("binary", "continuous", "count"))

  mu_mean     <- hyperpriors$mu_mean     %||% 0
  mu_sd       <- hyperpriors$mu_sd       %||% 2.5
  mu1_mean    <- hyperpriors$mu1_mean    %||% 0
  mu1_sd      <- hyperpriors$mu1_sd      %||% 2.5
  sigma1_mean <- hyperpriors$sigma1_mean %||% 0
  sigma1_sd   <- hyperpriors$sigma1_sd   %||% 1
  tau_mean    <- hyperpriors$tau_mean    %||% 0
  tau_sd      <- hyperpriors$tau_sd      %||% 0.5
  rho_alpha   <- hyperpriors$rho_alpha   %||% 2
  rho_beta    <- hyperpriors$rho_beta    %||% 2

  data_block  <- .tdp_data_block(outcome_type)
  likelihood  <- .tdp_likelihood(outcome_type)
  param_decl  <- .tdp_param_declarations(outcome_type)

  glue::glue(
'// Time-Dependent Prior (TDP) — {outcome_type} outcome
// Ornstein-Uhlenbeck continuous-time process with temporal decay

data \\{
{data_block}
  vector[H] time;
  real<lower=time[H]> t_star;
\\}

transformed data \\{
  vector[H - 1] delta;
  real<lower=0> delta_star;
  for (h in 2:H) \\{
    delta[h - 1] = time[h] - time[h - 1];
  \\}
  delta_star = t_star - time[H];
\\}

parameters \\{
  vector[H] theta;
  real mu;
  real mu1;
  real<lower=0> sigma1;
  real<lower=0> tau;
  real<lower=0, upper=1> rho;
\\}

model \\{
  // Hyperpriors
  mu ~ normal({mu_mean}, {mu_sd});
  mu1 ~ normal({mu1_mean}, {mu1_sd});
  sigma1 ~ normal({sigma1_mean}, {sigma1_sd});
  tau ~ normal({tau_mean}, {tau_sd});
  rho ~ beta({rho_alpha}, {rho_beta});

  // Evolution: theta[1] is the initial parameter value
  theta[1] ~ normal(mu1, sigma1);

  // Temporal evolution for h >= 2 via O-U process
  for (h in 2:H) \\{
    real mean_h;
    real var_h;
    mean_h = mu + pow(rho, delta[h - 1]) * (theta[h - 1] - mu);
    var_h  = square(tau) * (1 - pow(rho, 2 * delta[h - 1])) / (1 - square(rho));
    theta[h] ~ normal(mean_h, sqrt(var_h));
  \\}

  // Likelihood
{likelihood}
\\}

generated quantities \\{
  // Predictive parameter at current time
  real mean_star;
  real var_star;
  {param_decl}

  mean_star  = mu + pow(rho, delta_star) * (theta[H] - mu);
  var_star   = square(tau) * (1 - pow(rho, 2 * delta_star)) / (1 - square(rho));
  theta_star = normal_rng(mean_star, sqrt(var_star));
\\}
'
  )
}

.tdp_data_block <- function(outcome_type) {
  switch(outcome_type,
    binary = "  int<lower=2> H;\n  vector[H] ybar;\n  vector<lower=0>[H] se;",
    continuous = "  int<lower=2> H;\n  vector[H] ybar;\n  vector<lower=0>[H] se;",
    count = "  int<lower=2> H;\n  array[H] int<lower=0> events;\n  vector[H] log_n;"
  )
}

.tdp_likelihood <- function(outcome_type) {
  switch(outcome_type,
    binary = "  // theta[h] is logit(p) for study h\n  ybar ~ normal(theta, se);",
    continuous = "  ybar ~ normal(theta, se);",
    count = "  // theta[h] is log-rate for study h\n  for (h in 1:H)\n    events[h] ~ poisson_log(log_n[h] + theta[h]);"
  )
}

.tdp_param_declarations <- function(outcome_type) {
  switch(outcome_type,
    binary = "real p_star;\n  p_star = inv_logit(theta_star);",
    continuous = "",
    count = "real lambda_star;\n  lambda_star = exp(theta_star);"
  )
}
