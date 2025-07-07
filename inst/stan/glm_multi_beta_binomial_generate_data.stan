functions{
    #include common_functions.stan
    
}
data {

	int N;
	int N_original;
	int M;
	int C;
	int A;
	array[N] int exposure;

	// Which column of design, coefficient matrices should be used to generate the data
	int length_X_which;
	int length_XA_which;
	array[length_X_which] int X_which;
	array[length_XA_which] int XA_which;
	matrix[N, length_X_which] X;
  matrix[N, length_XA_which] Xa; // The variability design
  
	matrix[N_original, C] X_original;
	int is_truncated;
	real<lower=1> truncation_ajustment;

	// Random intercept
	
  int is_random_effect;
  
	array[2] int length_X_random_effect_which;
	array[length_X_random_effect_which[1]] int X_random_effect_which;
	array[2] int ncol_X_random_eff;
	array[2] int  ncol_X_random_eff_new;
	matrix[N, ncol_X_random_eff_new[1]] X_random_effect;
	matrix[N, ncol_X_random_eff_new[2]] X_random_effect_2;
	array[length_X_random_effect_which[2]] int X_random_effect_which_2;
        
  // Should I create intercept for generate quantities
  int<lower=0, upper=1> create_intercept;
  int<lower=0> A_intercept_columns;
  
  // Covariance setup
  array[2] int n_groups;
  array[2] int how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;

  array[2] int<lower=0, upper = 1> unknown_grouping;
  
  // Dimensions for unseen random effects
  array[2] int ncol_X_random_eff_unseen;
  
  // Matrix for unseen random effects
  matrix[N, ncol_X_random_eff_unseen[1]] X_random_effect_unseen;
  matrix[N, ncol_X_random_eff_unseen[2]] X_random_effect_2_unseen;
}
transformed data{
  // If needed recreate the intercept
  matrix[N,1] X_intercept;
  
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);
  
  X_intercept = to_matrix(rep_vector(1, N));
}

parameters {

  array[C] vector[M] beta_raw; // Each row is a vector of length M
  matrix[A, M] alpha; // Variability
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept // array of sum_to_zero_vector for each random effect
  array[ncol_X_random_eff[1] * (is_random_effect>0)] sum_to_zero_vector[M] random_effect_raw;
  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0)] sum_to_zero_vector[M] random_effect_raw_2;
  
  // sd of random intercept
  array[2 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[2 * (is_random_effect>0)] real random_effect_sigma_sigma;

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

  array[N, M] int counts_uncorrected;

  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  array[N] real generated_exposure;

  matrix[C,M] beta;
  
  // Convert vectors to matrix and apply sum-to-zero constraint manually
  for(c in 1:C) {
    vector[M] temp_beta = beta_raw[c];
    // NOTE: Due to floating point precision, we must explicitly normalize to sum to zero
    // instead of declating the sum_to_zero variabe
    temp_beta = normalize_sum_to_zero(temp_beta);
    beta[c,] = to_row_vector(temp_beta);
  }
  
  // Subset for mean and deviation
  matrix[length_X_which,M] my_beta = beta[X_which,];
  matrix[length_XA_which,M] my_alpha = alpha[XA_which,];
  
  matrix[M,N] mu;
  matrix[M,N] precision;
  
  // If needed recreate the intercept
  if(create_intercept == 1){

    // Create mean and deviation
    mu = (
      append_col(
        to_matrix(rep_vector(1, N)), // Intercept
        X // Rest
        ) *
        append_row(
          average_by_col(beta[1:A_intercept_columns,]), // Intercept
          my_beta // Rest
        )
    )';

    precision = (
      append_col(
        to_matrix(rep_vector(1, N)), // Intercept
        Xa // Rest
        ) *
        append_row(
          average_by_col(alpha[1:A_intercept_columns,]), // Intercept
          my_alpha // Rest
          )
    )' /
    (is_truncated ? truncation_ajustment : 1);

  }
  else {
    // Create mean and deviation
    mu = (X * my_beta)';
    precision = (Xa * my_alpha)' / (is_truncated ? truncation_ajustment : 1);

  }

  // Non centered parameterisation SD of random effects
  array[M-1 * (ncol_X_random_eff[1]> 0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma;
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M-1)) random_effect_sigma[m] = random_effect_sigma_mu[1] + random_effect_sigma_sigma[1] * random_effect_sigma_raw[m];
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M-1)) random_effect_sigma[m] = exp(random_effect_sigma[m]/3.0);
  
  // Non centered parameterisation SD of random effects 2
  array[M-1 * (ncol_X_random_eff[2]> 0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_2;
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M-1)) random_effect_sigma_2[m] = random_effect_sigma_mu[2] + random_effect_sigma_sigma[2] * random_effect_sigma_raw_2[m];
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M-1)) random_effect_sigma_2[m] = exp(random_effect_sigma_2[m]/3.0);
    
  // Random intercept
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M] random_effect; 
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M] random_effect_2; 
  
  // For first random effect
  if(length_X_random_effect_which[1]>0) {

    // Convert sum_to_zero_vector array to vector array for function call
    array[ncol_X_random_eff[1]] vector[M] random_effect_raw_vec;
    for(i in 1:ncol_X_random_eff[1]) {
      random_effect_raw_vec[i] = to_vector(random_effect_raw[i]);
    }

    // Generate random effects matrix - either from fitted effects or random draws
    
    // Get transformed random effects
    random_effect = get_random_effect_matrix(
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
    
    // Apply random effects
    mu = mu + (X_random_effect * random_effect[X_random_effect_which,])';
    
    // Add random effects for unseen groups if they exist
    if(ncol_X_random_eff_unseen[1] > 0) {
      matrix[ncol_X_random_eff_unseen[1], M] unseen_random_effect = 
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[1] * M), ncol_X_random_eff_unseen[1], M);
      
      // Apply sum-to-zero constraint to unseen random effects
      for(i in 1:ncol_X_random_eff_unseen[1]) {
        unseen_random_effect[i,] = to_row_vector(normalize_sum_to_zero(to_vector(unseen_random_effect[i,])));
      }
      
      mu = mu + (X_random_effect_unseen * unseen_random_effect)';
    }
  }
  
  // For second random effect
  if(length_X_random_effect_which[2]>0) {

    // Convert sum_to_zero_vector array to vector array for function call
    array[ncol_X_random_eff[2]] vector[M] random_effect_raw_2_vec;
    for(i in 1:ncol_X_random_eff[2]) {
      random_effect_raw_2_vec[i] = to_vector(random_effect_raw_2[i]);
    }

    // Generate random effects matrix - either from fitted effects or random draws
    
    // Get transformed random effects
    random_effect_2 = get_random_effect_matrix(
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
    
    // Apply random effects
    mu = mu + (X_random_effect_2 * random_effect_2[X_random_effect_which_2,])';
    
    // Add random effects for unseen groups if they exist
    if(ncol_X_random_eff_unseen[2] > 0) {
      matrix[ncol_X_random_eff_unseen[2], M] unseen_random_effect_2 = 
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[2] * M), ncol_X_random_eff_unseen[2], M);
      
      // Apply sum-to-zero constraint to unseen random effects
      for(i in 1:ncol_X_random_eff_unseen[2]) {
        unseen_random_effect_2[i,] = to_row_vector(normalize_sum_to_zero(to_vector(unseen_random_effect_2[i,])));
      }
      
      mu = mu + (X_random_effect_2_unseen * unseen_random_effect_2)';
    }
  }

  // Calculate proportions
  for(i in 1:N) mu[,i] = softmax(mu[,i]);

	// Generate
	for(i in 1:N) {
    	counts_uncorrected[i,] = beta_binomial_rng(
    	  exposure[i],
    	  mu[,i] .* exp(precision[,i]),
    	  (1.0 - mu[,i]) .* exp(precision[,i])
    	 );
	}

	// Calculate the generated exposure
  for(n in 1:N) generated_exposure[n] = max( sum(counts_uncorrected[n]), 1); // avoid dividing by zero
  for(n in 1:N) counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure[n];
}


