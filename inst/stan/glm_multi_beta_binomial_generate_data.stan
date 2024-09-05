functions{

 vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }
  
   row_vector sum_to_zero_QR(row_vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    row_vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
  
    row_vector average_by_col(matrix beta){
    return
    rep_row_vector(1.0, rows(beta)) * beta / rows(beta);
  }
  
matrix get_random_effect_matrix(
		int M, 
		int n_groups, 
		int how_many_factors_in_random_design, 
		int is_random_effect,
		int ncol_X_random_eff,
		array[,] int group_factor_indexes_for_covariance,
	
		matrix random_effect_raw,
		
		array[] vector random_effect_sigma_raw,
		array[] real random_effect_sigma_mu,
		array[] real random_effect_sigma_sigma,
		matrix sigma_correlation_factor
	){
		
		matrix[ncol_X_random_eff * (is_random_effect>0), M-1] random_effect; 
		
		
		// PIVOT WIDER
		// increase of one dimension array[cell_type] matrix[group, factor]
		array[M-1] matrix[ how_many_factors_in_random_design, n_groups] matrix_of_random_effects_raw;
		
		for(w in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
			
			// If I don't have the factor for one group 
			if(group_factor_indexes_for_covariance[j,i] == 0)
			matrix_of_random_effects_raw[w, j,i] = 0;
			else 
			matrix_of_random_effects_raw[w, j,i] = random_effect_raw[group_factor_indexes_for_covariance[j,i], w];
		}
		
		// Design L
		array[M-1] matrix[how_many_factors_in_random_design, how_many_factors_in_random_design] L;
		array[M-1] matrix[how_many_factors_in_random_design, n_groups] matrix_of_random_effects;
		
		// Non centered parameterisation
		array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design] random_effect_sigma;
		for(w in 1:(M-1)) random_effect_sigma[w] = random_effect_sigma_mu[1] + random_effect_sigma_sigma[1] * random_effect_sigma_raw[w];
		for(w in 1:(M-1)) random_effect_sigma[w] = exp(random_effect_sigma[w]/3.0);
		
		
		for(w in 1:(M-1)) L[w] = diag_pre_multiply(random_effect_sigma[w], sigma_correlation_factor) ;
		for(w in 1:(M-1)) matrix_of_random_effects[w] = L[w] * matrix_of_random_effects_raw[w];
		
		// Pivot longer
		for(w in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
			
			// If I don't have the factor for one group 
			if(group_factor_indexes_for_covariance[j,i] > 0)
			random_effect[group_factor_indexes_for_covariance[j,i], w] = matrix_of_random_effects[w,j,i];
		}
		
		return(random_effect);
}

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


}
transformed data{
  matrix[C, C] R_ast_inverse;
  vector[2*M] Q_r;
  
  // thin and scale the QR decomposition
  R_ast_inverse = inverse(qr_thin_R(X_original) / sqrt(N_original - 1));
  Q_r = Q_sum_to_zero_QR(M);
  
  // If needed recreate the intercept
  matrix[N,1] X_intercept;
  
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);
  
  X_intercept = to_matrix(rep_vector(1, N));


  // If I get crazy diagonal matrix omit it
  if(is_random_effect>0) { 
    if(max(R_ast_inverse)>1000 )
      print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
  
}

parameters {

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
	cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;

	// Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]]  random_effect_sigma_raw_2;
	cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;

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
  matrix[C,M] beta_raw;
  for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);
  beta = R_ast_inverse * beta_raw; // coefficients on x
  
  // Subset for mean and deviation
  matrix[length_X_which,M] my_beta = beta[X_which,];
  matrix[length_XA_which,M] my_alpha = alpha[XA_which,];
  
  matrix[M,N] mu;
  matrix[M,N] precision;
  

//   // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
//   if(N_grouping==0) beta_random_effect[1] = rep_row_vector(0.0, M);
// 
//   // Random effect
//   else{
//       random_effect = get_random_effect_matrix(
// 			M, 
// 			how_many_groups, 
// 			how_many_factors_in_random_design, 
// 			N_random_effects,
// 			N_grouping,
// 			group_factor_indexes_for_covariance,
// 		
// 			random_effect_raw,
// 			
// 			random_effect_sigma_raw,
// 			random_effect_sigma_mu,
// 			random_effect_sigma_sigma,
// 			sigma_correlation_factor
// 		);
//       beta_random_effect[,1:(M-1)] = random_effect;
//       for(n in 1:N_grouping) beta_random_effect[n, M] = -sum(random_effect[n,]);
//   }
  
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

// Random intercept
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_effect; 
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M-1] random_effect_2; 
  
  
  if(length_X_random_effect_which[1]>0){
    
    
    random_effect = 
  	get_random_effect_matrix(
			M, 
			n_groups[1], 
			how_many_factors_in_random_design[1], 
			is_random_effect,
			ncol_X_random_eff[1],
			group_factor_indexes_for_covariance,
		
			random_effect_raw,
			
			random_effect_sigma_raw,
			random_effect_sigma_mu,
			random_effect_sigma_sigma,
			sigma_correlation_factor
		);
    
        mu = mu + append_row((X_random_effect * random_effect[X_random_effect_which,])', rep_row_vector(0, N));

  }
  
    // Random intercept 2
  if(length_X_random_effect_which[2]>0){
    
     // Covariate setup 
  random_effect_2 = 
  	get_random_effect_matrix(
			M, 
			n_groups[2], 
			how_many_factors_in_random_design[2], 
			is_random_effect,
			ncol_X_random_eff[2],
			group_factor_indexes_for_covariance_2,
		
			random_effect_raw_2,
			
			random_effect_sigma_raw_2,
			random_effect_sigma_mu,
			random_effect_sigma_sigma,
			sigma_correlation_factor_2
		);

    // Update with summing mu_random_effect
    mu = mu + append_row((X_random_effect_2 * random_effect_2[X_random_effect_which_2,])', rep_row_vector(0, N));


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


