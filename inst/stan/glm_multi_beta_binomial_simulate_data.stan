data{
  int N_simulated; // Number of subjects
  int M_simulated; // Number of categories
  int C_simulated;
  int A_simulated;
  array[N_simulated] int exposure_simulated;
  matrix[N_simulated, C_simulated] X_simulated;
  matrix[A_simulated, A_simulated] XA_simulated;
  matrix[C_simulated,M_simulated] beta_simulated;
  real<lower=0> variability_multiplier;

	int M;
	int C;
	int A;
	array[2] int ncol_X_random_eff;
	int is_random_effect;
  array[2] int how_many_factors_in_random_design;


}

parameters{

  matrix[C, M-1] beta_raw_raw; // matrix with C rows and number of cells (-1) columns
  matrix[A, M] alpha; // Variability
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept // matrix with N_groupings rows and number of cells (-1) columns
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_effect_raw;
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M-1] random_effect_raw_2;
  
  // sd of random intercept
  array[is_random_effect>0] real random_effect_sigma_mu;
  array[is_random_effect>0] real random_effect_sigma_sigma;

	// Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[1]]  random_effect_sigma_raw;
	array[M-1 * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;

	// Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]]  random_effect_sigma_raw_2;
	array[M-1 * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;

  // If I have just one group
  array[is_random_effect>0] real zero_random_effect;
  
}

generated quantities{
  array[N_simulated, M_simulated] int counts_uncorrected;
  matrix[N_simulated, M_simulated] counts;
  matrix[A_simulated,M_simulated] alpha_simulated;
  matrix[M_simulated,N_simulated] mu = (X_simulated * beta_simulated)';
  matrix[M_simulated,N_simulated] precision;

  matrix[A_simulated,M_simulated] beta_intercept_slope;
  // Vector of the generated exposure_simulated
  array[N_simulated] real generated_exposure;

  // matrix[A_simulated,M_simulated] alpha_intercept_slope;

  // All this because if A_simulated ==1 we have ocnversion problems
  // This works only with two discrete groups
  if(A_simulated == 1) beta_intercept_slope = to_matrix(beta_simulated[A_simulated,], A_simulated, M_simulated, 0);
  else beta_intercept_slope = (XA_simulated * beta_simulated[1:A_simulated,]);
  // if(A_simulated == 1)  alpha_intercept_slope = alpha;
  // else alpha_intercept_slope = (XA_simulated * alpha);

  // PRECISION REGRESSION
  for(a in 1:A_simulated) for(m in 1:M_simulated) alpha_simulated[a,m] = normal_rng( beta_intercept_slope[a,m] * prec_coeff[2] + prec_coeff[1], prec_sd);

  precision = (X_simulated[,1:A_simulated] * alpha_simulated)';

  // Precision adjustment
  precision = precision - log(variability_multiplier);

  for(i in 1:N_simulated) mu[,i] = softmax(mu[,i]);
  for(i in 1:cols(mu)) {
    counts_uncorrected[i,] = beta_binomial_rng(
      exposure_simulated[i],
      mu[,i] .* exp(precision[,i]),
      (1.0 - mu[,i]) .* exp(precision[,i])
      );
  }

  // Calculate the generated exposure_simulated
  for(n in 1:N_simulated) generated_exposure[n] = sum(counts_uncorrected[n]);
  for(n in 1:N_simulated) counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure_simulated[n];

}
