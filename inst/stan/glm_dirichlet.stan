// Dirichlet-multinomial / Dirichlet noise model:
// counts ~ DirichletMultinomial(mu * concentration)
// proportions ~ Dirichlet(mu * concentration)
// No mean-dispersion association - scalar concentration is independent of mean

functions{

  #include common_functions.stan

  real partial_sum_dirichlet_multinomial_lpmf(
    array[] int idx_y,
    int start,
    int end,
    int is_proportion,
    array[,] int y,
    array[,] real y_proportion,
    matrix X,
    matrix Xa,
    matrix beta,
    vector alpha,
    array[] int ncol_X_random_eff,
    matrix X_random_effect,
    matrix X_random_effect_2,
    matrix random_effect,
    matrix random_effect_2,
    int M
  ){
    int N = end - start + 1;
    real lp = 0;

    // mu
    matrix[M, N] mu = (X[idx_y,] * beta)';
    if(ncol_X_random_eff[1] > 0)
      mu = mu + (X_random_effect[idx_y,] * random_effect)';
    if(ncol_X_random_eff[2] > 0)
      mu = mu + (X_random_effect_2[idx_y,] * random_effect_2)';

    for(n in 1:N) mu[,n] = softmax(mu[,n]);

    // Scalar concentration per sample
    vector[N] concentration = exp(Xa[idx_y,] * alpha);

    // Dirichlet(-multinomial): alpha_vec = mu .* concentration (> 0)
    for(n in 1:N) {
      vector[M] alpha_vec = mu[,n] .* concentration[n];
      if(is_proportion)
        lp += dirichlet_lpdf(to_vector(y_proportion[idx_y[n]]) | alpha_vec);
      else
        lp += dirichlet_multinomial_lupmf(y[idx_y[n]] | alpha_vec);
    }
    return lp;
  }

}
data{
  int<lower=0, upper=1> is_proportion;
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> C;
  int<lower=1> A;
  int<lower=1> A_intercept_columns;
  int<lower=1> B_intercept_columns;
  int<lower=1> Ar;
  array[N] int exposure;  // Unused for Dirichlet(-multinomial), kept for data compatibility
  array[N * !is_proportion, M] int<lower=0> y;
  array[N * is_proportion, M] real<lower=0, upper=1> y_proportion;
  matrix[N, C] X;
  matrix[Ar, A] XA;
  matrix[N, A] Xa;

  int<lower=1, upper=N*M> TNS;
  array[TNS] int<lower=1, upper=N*M> truncation_not_idx;
  int TNIM;
  array[TNIM,2] int<lower=1, upper=N*M> truncation_not_idx_minimal;

  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;   // Unused, kept for compatibility
  array[2] real prior_prec_sd;
  array[2] real prior_mean_intercept;
  array[2] real prior_mean_coefficients;

  int<lower=0, upper=1> exclude_priors;
  int<lower=0, upper=1> bimodal_mean_variability_association;  // Unused for Dirichlet-multinomial
  int<lower=0, upper=1> use_data;
  int<lower=1> grainsize;
  int<lower=0, upper=1> intercept_in_design;

  int is_random_effect;
  array[2] int ncol_X_random_eff;
  matrix[N, ncol_X_random_eff[1]] X_random_effect;
  matrix[N, ncol_X_random_eff[2]] X_random_effect_2;

  array[2] int n_groups;
  array[2] int how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;

  int<lower=0, upper=1> enable_loo;
}
transformed data{
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);

  array[N] int array_N;
  for(n in 1:N) array_N[n] = n;
}
parameters{
  array[C] sum_to_zero_vector[M] beta_raw;
  vector[A] alpha;  // Scalar concentration (log scale)

  array[ncol_X_random_eff[1] * (is_random_effect>0)] sum_to_zero_vector[M] random_effect_raw;
  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0)] sum_to_zero_vector[M] random_effect_raw_2;

  array[2 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[2 * (is_random_effect>0)] real random_effect_sigma_sigma;

  array[M * (is_random_effect>0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma_raw;
  array[M * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;

  array[M * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_raw_2;
  array[M * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;

  array[is_random_effect>0] real zero_random_effect;
}
transformed parameters{
  matrix[C,M] beta;
  for(c in 1:C) beta[c,] = to_row_vector(beta_raw[c]);

  array[M * (ncol_X_random_eff[1]> 0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma;
  if(ncol_X_random_eff[1]> 0) for(m in 1:M) random_effect_sigma[m] = random_effect_sigma_mu[1] + random_effect_sigma_sigma[1] * random_effect_sigma_raw[m];
  if(ncol_X_random_eff[1]> 0) for(m in 1:M) random_effect_sigma[m] = exp(random_effect_sigma[m]/3.0);

  array[M * (ncol_X_random_eff[2]> 0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_2;
  if(ncol_X_random_eff[2]> 0) for(m in 1:M) random_effect_sigma_2[m] = random_effect_sigma_mu[2] + random_effect_sigma_sigma[2] * random_effect_sigma_raw_2[m];
  if(ncol_X_random_eff[2]> 0) for(m in 1:M) random_effect_sigma_2[m] = exp(random_effect_sigma_2[m]/3.0);

  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M] random_effect;
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M] random_effect_2;

  if(ncol_X_random_eff[1]> 0){
    array[ncol_X_random_eff[1]] vector[M] random_effect_raw_vec;
    for(i in 1:ncol_X_random_eff[1]) random_effect_raw_vec[i] = to_vector(random_effect_raw[i]);
    random_effect = get_random_effect_matrix(
      M, n_groups[1], how_many_factors_in_random_design[1],
      is_random_effect, ncol_X_random_eff[1],
      group_factor_indexes_for_covariance,
      random_effect_raw_vec, random_effect_sigma, sigma_correlation_factor
    );
  }

  if(ncol_X_random_eff[2]>0){
    array[ncol_X_random_eff[2]] vector[M] random_effect_raw_2_vec;
    for(i in 1:ncol_X_random_eff[2]) random_effect_raw_2_vec[i] = to_vector(random_effect_raw_2[i]);
    random_effect_2 = get_random_effect_matrix(
      M, n_groups[2], how_many_factors_in_random_design[2],
      is_random_effect, ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2,
      random_effect_raw_2_vec, random_effect_sigma_2, sigma_correlation_factor_2
    );
  }
}
model{
  if(use_data == 1){
    target += reduce_sum(
      partial_sum_dirichlet_multinomial_lpmf,
      array_N,
      grainsize,
      is_proportion,
      y,
      y_proportion,
      X, Xa, beta, alpha,
      ncol_X_random_eff,
      X_random_effect, X_random_effect_2,
      random_effect, random_effect_2,
      M
    );
  }

  alpha ~ normal(prior_prec_intercept[1], prior_prec_sd[1]);

  for(c in 1:B_intercept_columns) beta_raw[c] ~ normal(prior_mean_intercept[1], prior_mean_intercept[2] * inv(sqrt(1 - inv(M))));
  if(C>B_intercept_columns) for(c in (B_intercept_columns+1):C) beta_raw[c] ~ normal(prior_mean_coefficients[1], prior_mean_coefficients[2] * inv(sqrt(1 - inv(M))));

  if(is_random_effect>0){
    for(m in 1:M) random_effect_raw[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for(m in 1:M) random_effect_sigma_raw[m] ~ std_normal();
    for(m in 1:M) sigma_correlation_factor[m] ~ lkj_corr_cholesky(2);
    random_effect_sigma_mu ~ std_normal();
    random_effect_sigma_sigma ~ std_normal();
    zero_random_effect ~ std_normal();
  }
  if(ncol_X_random_eff[2]>0){
    for(m in 1:M) random_effect_raw_2[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for(m in 1:M) random_effect_sigma_raw_2[m] ~ std_normal();
    for(m in 1:M) sigma_correlation_factor_2[m] ~ lkj_corr_cholesky(2);
  }
}
generated quantities {
  matrix[A, M] alpha_normalised = rep_matrix(alpha, M);  // No mean-dispersion adjustment

  // One log_lik per sample (full composition vector)
  vector[N] log_lik = rep_vector(0, N);

  if(enable_loo == 1){
    matrix[M, N] mu = (X * beta)';
    if(ncol_X_random_eff[1]> 0) mu = mu + (X_random_effect * random_effect)';
    if(ncol_X_random_eff[2]>0) mu = mu + (X_random_effect_2 * random_effect_2)';
    for(n in 1:N) mu[,n] = softmax(mu[,n]);

    vector[N] concentration = exp(Xa * alpha);

    for(n in 1:N){
      vector[M] alpha_vec = mu[,n] .* concentration[n];
      if(is_proportion)
        log_lik[n] = dirichlet_lpdf(to_vector(y_proportion[n]) | alpha_vec);
      else
        log_lik[n] = dirichlet_multinomial_lpmf(y[n] | alpha_vec);
    }
  }
}
