functions{
  #include common_functions.stan
}

data{
  int N_simulated; // Number of subjects
  int M_simulated; // Number of categories
  int C_simulated;
  int A_simulated;
  array[N_simulated] int exposure_simulated;
  matrix[N_simulated, C_simulated] X_simulated;
  matrix[N_simulated, A_simulated] Xa_simulated; // Variability design matrix for simulated data
  matrix[A_simulated, A_simulated] XA_simulated; // Unique variability design (for compatibility with old code)

  real<lower=0> variability_multiplier;
  
  // Optional: provided coefficients to override posterior beta_raw
  int<lower=0, upper=1> user_provided_beta; // 1 if beta_simulated is provided, 0 to use posterior beta_raw
  // NOTE: Using vector[M_simulated] instead of sum_to_zero_vector[M_simulated] to avoid floating-point precision issues
  // Floating-point precision can cause sum_to_zero_vector to fail the strict sum-to-zero constraint when values are passed from R
  // The sum-to-zero constraint is a modeling constraint (compositional data must sum to 1, so log-ratios sum to 0)
  // but does not need to be enforced at the Stan type level for generate_quantities
  // Small precision errors are acceptable since Stan doesn't enforce strict constraints with vector types
  array[C_simulated] vector[M_simulated] beta_simulated_provided; // Provided coefficients
  
  // Optional: provided slope and intercept for mean-dispersion association
  int<lower=0, upper=1> user_provided_prec_coeff_slope; // 1 if prec_coeff_slope is provided, 0 to use posterior prec_coeff[2]
  real prec_coeff_slope_provided; // Provided slope for mean-dispersion association (alpha = beta * slope + intercept)
  int<lower=0, upper=1> user_provided_prec_coeff_intercept; // 1 if prec_coeff_intercept is provided, 0 to use posterior prec_coeff[1]
  real prec_coeff_intercept_provided; // Provided intercept for mean-dispersion association (alpha = beta * slope + intercept)

  int M;
	int C;
	int A;
	int A_intercept_columns; // How many intercept columns in variability design
	int intercept_in_design; // Whether intercept is in design
	int bimodal_mean_variability_association; // Whether to use bimodal association
	
	// Indices for subsetting beta and alpha (matching replicate_data approach)
	int length_X_which;
	int length_XA_which;
	array[length_X_which] int X_which;
	array[length_XA_which] int XA_which;
	array[2] int ncol_X_random_eff;
	int is_random_effect;
  array[2] int how_many_factors_in_random_design;
  
  // Random effects design matrices (if random effects exist, otherwise empty matrices)
  // Stan requires at least 1 column, so we use max(1, ncol) - the conditional logic will handle when ncol is 0
  matrix[N_simulated, max(1, ncol_X_random_eff[1])] X_random_effect_simulated;
  matrix[N_simulated, max(1, ncol_X_random_eff[2])] X_random_effect_2_simulated;
  
  // Covariance setup for random effects (if random effects exist)
  array[2] int n_groups;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
}

transformed data{
  // Issue 5: Ensure we have valid dimensions for random effect matrices
  // These are used in generated quantities to subset matrices correctly
  // ncol_re1 and ncol_re2 represent the actual matrix dimensions (always >= 1)
  // ncol_X_random_eff[1] and ncol_X_random_eff[2] represent the number of random effect parameters (can be 0)
  int ncol_re1 = max(1, ncol_X_random_eff[1]);
  int ncol_re2 = max(1, ncol_X_random_eff[2]);
}

parameters{
  // NOTE: Using vector[M] instead of sum_to_zero_vector[M] to avoid floating-point precision issues
  // Floating-point precision can cause sum_to_zero_vector to fail the strict sum-to-zero constraint
  // when posterior draws are reloaded from disk. The sum-to-zero constraint is a modeling constraint
  // (compositional data must sum to 1, so log-ratios sum to 0) but does not need to be enforced
  // at the Stan type level for generate_quantities. Small precision errors are acceptable.
  array[C] vector[M] beta_raw; // Each row is a vector of length M
  matrix[A, M] alpha; // Variability - kept in parameters for generate_quantities compatibility
  // Note: alpha from posterior is NOT used - alpha_simulated is computed from beta in generated quantities
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept - using vector instead of sum_to_zero_vector to avoid precision issues
  array[ncol_X_random_eff[1] * (is_random_effect>0)] vector[M] random_effect_raw;
  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0)] vector[M] random_effect_raw_2;
  
  // sd of random intercept
  array[2 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[2 * (is_random_effect>0)] real random_effect_sigma_sigma;

	// Covariance
  array[M * (is_random_effect>0)] vector[how_many_factors_in_random_design[1]]  random_effect_sigma_raw;
	array[M * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;

	// Covariance
  array[M * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]]  random_effect_sigma_raw_2;
	array[M * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;

  // If I have just one group
  array[is_random_effect>0] real zero_random_effect;
  
}

transformed parameters{
  // Convert vector to matrix
  // If user_provided_beta is 1, use provided coefficients
  // Otherwise, use beta_raw from posterior draws
  // Note: We use vector[M] instead of sum_to_zero_vector[M] to avoid floating-point precision issues
  // Small precision errors in sum-to-zero constraint are acceptable since Stan doesn't enforce strict constraints
  matrix[C,M] beta;
  if(user_provided_beta == 1) {
    // Use provided coefficients - map them to the full beta matrix using X_which
    // First, create beta matrix from beta_raw (needed for rows not in X_which)
    for(c in 1:C) {
      beta[c,] = to_row_vector(beta_raw[c]);
    }
    // Then override the rows specified by X_which with provided coefficients
    for(i in 1:length_X_which) {
      int c = X_which[i];
      // Provided coefficients are vector[M_simulated]
      // Map them to full M dimensions (no padding needed - X_which ensures correct mapping)
      for(m in 1:M_simulated) {
        beta[c, m] = beta_simulated_provided[i][m];
      }
      // If M_simulated < M, pad remaining with zeros
      if(M_simulated < M) {
        for(m in (M_simulated + 1):M) {
          beta[c, m] = 0.0;
        }
      }
    }
  } else {
    // Use beta_raw from posterior draws
    for(c in 1:C) {
      beta[c,] = to_row_vector(beta_raw[c]);
    }
  }
  
  // Non centered parameterisation SD of random effects (matching main model)
  array[M * (ncol_X_random_eff[1]> 0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma;
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M)) random_effect_sigma[m] = random_effect_sigma_mu[1] + random_effect_sigma_sigma[1] * random_effect_sigma_raw[m];
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M)) random_effect_sigma[m] = exp(random_effect_sigma[m]/3.0);
  
  // Non centered parameterisation SD of random effects 2
  array[M * (ncol_X_random_eff[2]> 0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_2;
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M)) random_effect_sigma_2[m] = random_effect_sigma_mu[2] + random_effect_sigma_sigma[2] * random_effect_sigma_raw_2[m];
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M)) random_effect_sigma_2[m] = exp(random_effect_sigma_2[m]/3.0);
    
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M] random_effect;
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M] random_effect_2;
  
  // random intercept
  if(ncol_X_random_eff[1]> 0){
    
    // Convert vector array
    array[ncol_X_random_eff[1]] vector[M] random_effect_raw_vec;
    for(i in 1:ncol_X_random_eff[1]) {
      random_effect_raw_vec[i] = random_effect_raw[i];
    }
    
    // Covariate setup
    random_effect =
    get_random_effect_matrix(
      M,
      n_groups[1],
      how_many_factors_in_random_design[1],
      is_random_effect,
      ncol_X_random_eff[1],
      group_factor_indexes_for_covariance,
      random_effect_raw_vec,
      random_effect_sigma,
      sigma_correlation_factor
      );
      
  }
  
  // random intercept
  if(ncol_X_random_eff[2]>0 ){

    // Convert vector array
    array[ncol_X_random_eff[2]] vector[M] random_effect_raw_2_vec;
    for(i in 1:ncol_X_random_eff[2]) {
      random_effect_raw_2_vec[i] = random_effect_raw_2[i];
    }

    // Covariate setup
    random_effect_2 =
    get_random_effect_matrix(
      M,
      n_groups[2],
      how_many_factors_in_random_design[2],
      is_random_effect,
      ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2,
      random_effect_raw_2_vec,
      random_effect_sigma_2,
      sigma_correlation_factor_2
      );

  }
}

generated quantities{
  array[N_simulated, M_simulated] int counts_uncorrected;
  matrix[N_simulated, M_simulated] counts;
  matrix[M_simulated,N_simulated] mu;
  matrix[M_simulated,N_simulated] precision;

  // Vector of the generated exposure_simulated
  array[N_simulated] real generated_exposure;

  // Use the transformed beta from transformed parameters
  // Subset beta using X_which (matching replicate_data approach) - no padding needed
  matrix[length_X_which, M_simulated] my_beta = beta[X_which, 1:M_simulated];
  
  // Compute alpha from beta using the association relationship (matching root model)
  // First compute full alpha_simulated from beta, then subset using XA_which
  matrix[A_simulated, M_simulated] alpha_simulated;
  
  // Use user-provided slope if available, otherwise use posterior prec_coeff[2]
  real prec_coeff_slope = user_provided_prec_coeff_slope ? prec_coeff_slope_provided : prec_coeff[2];
  // Use user-provided intercept if available, otherwise use posterior prec_coeff[1]
  real prec_coeff_intercept = user_provided_prec_coeff_intercept ? prec_coeff_intercept_provided : prec_coeff[1];
  
  // Build alpha from beta following the root model's association
  // Map A_simulated columns to corresponding beta columns via XA_which -> X_which
  if(A_simulated == 1) {
    // Simple case: single intercept column (variability ~ 1)
    // Find intercept column in X_which
    int beta_col = intercept_in_design && length_X_which > 0 ? X_which[1] : 1;
    for(m in 1:M_simulated) {
      alpha_simulated[1, m] = beta[beta_col, m] * prec_coeff_slope + prec_coeff_intercept;
    }
  } else {
    // Multiple columns: handle intercept and non-intercept columns separately
    int A_intercept_columns_sim = min(A_intercept_columns, A_simulated);
    
    // Intercept columns: alpha = beta * prec_coeff_slope + prec_coeff_intercept
    for(a in 1:A_intercept_columns_sim) {
      int alpha_col_idx = XA_which[a];
      // Find corresponding beta column - intercept columns map to first columns in X_which
      int beta_col = intercept_in_design && length_X_which > 0 ? X_which[min(a, length_X_which)] : 1;
      for(m in 1:M_simulated) {
        alpha_simulated[a, m] = beta[beta_col, m] * prec_coeff_slope + prec_coeff_intercept;
      }
    }
    
    // Non-intercept columns: alpha = beta * prec_coeff_slope
    if(A_simulated > A_intercept_columns_sim) {
      for(a in (A_intercept_columns_sim + 1):A_simulated) {
        int alpha_col_idx = XA_which[a];
        // Find corresponding beta column - map via XA_which to X_which
        int beta_col_idx = a;
        int beta_col = beta_col_idx <= length_X_which ? X_which[beta_col_idx] : X_which[length_X_which];
        for(m in 1:M_simulated) {
          alpha_simulated[a, m] = beta[beta_col, m] * prec_coeff_slope;
        }
      }
    }
  }
  
  // Subset alpha_simulated using XA_which to get my_alpha (matching replicate_data approach)
  matrix[length_XA_which, M_simulated] my_alpha = alpha_simulated[XA_which,];

  // Calculate mu using the subsetted beta (matching replicate_data approach)
  mu = (X_simulated * my_beta)';
  
  // Add random effects if present
  // Issue 5: Clarify dimension handling - ncol_X_random_eff[1] is the actual number of random effect parameters
  // The matrix X_random_effect_simulated has ncol_re1 columns (which is max(1, ncol_X_random_eff[1]))
  // We only use the first ncol_X_random_eff[1] columns if ncol_X_random_eff[1] > 0
  if(ncol_X_random_eff[1]> 0 && is_random_effect > 0) {
    // Extract the relevant subset of random effects (ncol_X_random_eff[1] rows, M_simulated columns)
    matrix[ncol_X_random_eff[1], M_simulated] random_effect_subset = random_effect[1:ncol_X_random_eff[1], 1:M_simulated];
    // Use only the first ncol_X_random_eff[1] columns of the design matrix
    mu = mu + (X_random_effect_simulated[,1:ncol_X_random_eff[1]] * random_effect_subset)';
  }
  
  if(ncol_X_random_eff[2]>0 && is_random_effect > 0) {
    // Extract the relevant subset of random effects (ncol_X_random_eff[2] rows, M_simulated columns)
    matrix[ncol_X_random_eff[2], M_simulated] random_effect_2_subset = random_effect_2[1:ncol_X_random_eff[2], 1:M_simulated];
    // Use only the first ncol_X_random_eff[2] columns of the design matrix
    mu = mu + (X_random_effect_2_simulated[,1:ncol_X_random_eff[2]] * random_effect_2_subset)';
  }

  // Calculate precision using the actual alpha from posterior
  // Issue 2: Use Xa_simulated (variability design matrix) with my_alpha (matching replicate_data approach)
  precision = (Xa_simulated * my_alpha)';

  // Precision adjustment (variability multiplier)
  precision = precision - log(variability_multiplier);

  // Convert to proportions
  for(i in 1:N_simulated) mu[,i] = softmax(mu[,i]);
  
  // Generate counts
  for(i in 1:N_simulated) {
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
