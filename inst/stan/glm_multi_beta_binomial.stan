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

  vector sum_to_zero_QR_vector(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

  array[] int rep_each(array[] int x, int K) {
    int N = size(x);
    array[N * K] int y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }


  real partial_sum_lpmf(array[] int slice_y,

                        int start,
                        int end,
                        array[] int exposure_array,
                        vector mu_array,
                        vector precision_array
                        ) {

return beta_binomial_lupmf(
    slice_y |
    exposure_array[start:end],
    (mu_array[start:end] .* precision_array[start:end]),
    (1.0 - mu_array[start:end]) .* precision_array[start:end]
  ) ;

}



  row_vector average_by_col(matrix X) {
    int rows_X = rows(X);
    int cols_X = cols(X);
    row_vector[cols_X] means;
    
    for (j in 1:cols_X) {
      means[j] = mean(X[, j]);
    }
    
    return means;


  }

  real abundance_variability_regression(row_vector variability, row_vector abundance, array[] real prec_coeff, real prec_sd, int bimodal_mean_variability_association, real mix_p){

    real lp = 0;
    // If mean-variability association is bimodal such as for single-cell RNA use mixed model
    if(bimodal_mean_variability_association == 1){
      for(m in 1:cols(variability))
      lp += log_mix(mix_p,
      normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
      normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
      );

      // If no bimodal
    } else {
      lp =  normal_lpdf(variability | abundance * prec_coeff[2] + prec_coeff[1], prec_sd );
    }

    return(lp);
  }
  
	matrix get_random_effect_matrix(
		int M, 
		int n_groups, 
		int how_many_factors_in_random_design, 
		int is_random_effect,
		int ncol_X_random_eff,
		array[,] int group_factor_indexes_for_covariance,
	
		matrix random_intercept_raw,
		
		array[] vector random_intercept_sigma_raw,
		array[] real random_intercept_sigma_mu,
		array[] real random_intercept_sigma_sigma,
		matrix sigma_correlation_factor
	){
		
		matrix[ncol_X_random_eff * (is_random_effect>0), M-1] random_intercept; 
		
		
		// PIVOT WIDER
		// increase of one dimension array[cell_type] matrix[group, factor]
		array[M-1] matrix[ how_many_factors_in_random_design, n_groups] matrix_of_random_effects_raw;
		
		for(w in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
			
			// If I don't have the factor for one group 
			if(group_factor_indexes_for_covariance[j,i] == 0)
			matrix_of_random_effects_raw[w, j,i] = 0;
			else 
			matrix_of_random_effects_raw[w, j,i] = random_intercept_raw[group_factor_indexes_for_covariance[j,i], w];
		}
		
		// Design L
		array[M-1] matrix[how_many_factors_in_random_design, how_many_factors_in_random_design] L;
		array[M-1] matrix[how_many_factors_in_random_design, n_groups] matrix_of_random_effects;
		
		// Non centered parameterisation
		array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design] random_intercept_sigma;
		for(w in 1:(M-1)) random_intercept_sigma[w] = random_intercept_sigma_mu[1] + random_intercept_sigma_sigma[1] * random_intercept_sigma_raw[w];
		for(w in 1:(M-1)) random_intercept_sigma[w] = exp(random_intercept_sigma[w]/3.0);
		
		
		for(w in 1:(M-1)) L[w] = diag_pre_multiply(random_intercept_sigma[w], sigma_correlation_factor) ;
		for(w in 1:(M-1)) matrix_of_random_effects[w] = L[w] * matrix_of_random_effects_raw[w];
		
		// Pivot longer
		for(w in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
			
			// If I don't have the factor for one group 
			if(group_factor_indexes_for_covariance[j,i] > 0)
			random_intercept[group_factor_indexes_for_covariance[j,i], w] = matrix_of_random_effects[w,j,i];
		}
		
		return(random_intercept);
}


}
data{
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> C;
  int<lower=1> A; // How many column in variability design\
  int<lower=1> A_intercept_columns; // How many intercept column in varibility design
  int<lower=1> B_intercept_columns; // How many intercept column in varibility design
  int<lower=1> Ar; // Rows of unique variability design
  array[N] int exposure;
  array[N,M] int y;
  matrix[N, C] X;
  matrix[Ar, A] XA; // The unique variability design
  matrix[N, A] Xa; // The variability design

  // Truncation
  int is_truncated;
  array[N,M] int truncation_up;
  array[N,M] int truncation_down;
  int<lower=1, upper=N*M> TNS; // truncation_not_size
  array[TNS] int<lower=1, upper=N*M> truncation_not_idx;
  int<lower=0, upper=1> is_vb;

  // Prior info
  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;
  array[2] real prior_prec_sd;
  array[2] real prior_mean_intercept;
  array[2] real prior_mean_coefficients;

  // Exclude priors for testing purposes
  int<lower=0, upper=1> exclude_priors;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;


  // Parallel chain
  int<lower=1> grainsize;

  // Does the design icludes intercept
  int <lower=0, upper=1> intercept_in_design;

  // Random intercept
  
  // Is the parameters in random effect matrix, minus ther sub to zero parameters, for example if I have four groups, this will be 3
  int is_random_effect;

  // Is the parameters in random effect matrix
  array[2] int ncol_X_random_eff;
  matrix[N, ncol_X_random_eff[1]] X_random_intercept;
  matrix[N, ncol_X_random_eff[2]] X_random_intercept_2;

  // Covariance setup
  array[2] int n_groups;
  array[2] int how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;

  // LOO
  int<lower=0, upper=1> enable_loo;


}
transformed data{
  vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;
  array[N*M] int y_array;
  array[N*M] int truncation_down_array;
  array[N*M] int exposure_array;
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);


  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));
  
  // If I get crazy diagonal matrix omit it
  if(is_random_effect>0) { 
    if(max(R_ast_inverse)>1000 )
      print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
  
  // Data vectorised
  y_array =  to_array_1d(y);
  exposure_array = rep_each(exposure, M);
}
parameters{
  matrix[C, M-1] beta_raw_raw; // matrix with C rows and number of cells (-1) columns
  matrix[A, M] alpha; // Variability
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept // matrix with ncol_X_random_effs rows and number of cells (-1) columns
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_intercept_raw;
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M-1] random_intercept_raw_2;

  
  // sd of random intercept
  array[is_random_effect>0] real random_intercept_sigma_mu;
  array[is_random_effect>0] real random_intercept_sigma_sigma;

	// Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[1]]  random_intercept_sigma_raw;
	cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;

	// Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]]  random_intercept_sigma_raw_2;
	cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;


  // If I have just one group
  array[is_random_effect>0] real zero_random_intercept;
  
  
}
transformed parameters{

  // Initialisation
  matrix[C,M] beta_raw;
  matrix[M, N] precision = (Xa * alpha)';
  matrix[C,M] beta;

  // locations distribution
  matrix[M, N] mu;

  // vectorisation
  vector[N*M] mu_array;
  vector[N*M] precision_array;

  for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);
  beta = R_ast_inverse * beta_raw; // coefficients on x

  // Calculate locations distribution
  mu = (Q_ast * beta_raw)';

  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_intercept; 
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M-1] random_intercept_2; 

  // random intercept
  if(ncol_X_random_eff[1]> 0){
  	
  // Covariate setup 
  random_intercept = 
  	get_random_effect_matrix(
			M, 
			n_groups[1], 
			how_many_factors_in_random_design[1], 
			is_random_effect,
			ncol_X_random_eff[1],
			group_factor_indexes_for_covariance,
		
			random_intercept_raw,
			
			random_intercept_sigma_raw,
			random_intercept_sigma_mu,
			random_intercept_sigma_sigma,
			sigma_correlation_factor
		);

    // Update with summing mu_random_intercept
    mu = mu + append_row((X_random_intercept * random_intercept)', rep_row_vector(0, N));
  }


  // random intercept
  if(ncol_X_random_eff[2]>0 ){
  	
  // Covariate setup 
  random_intercept_2 = 
  	get_random_effect_matrix(
			M, 
			n_groups[2], 
			how_many_factors_in_random_design[2], 
			is_random_effect,
			ncol_X_random_eff[2],
			group_factor_indexes_for_covariance_2,
		
			random_intercept_raw_2,
			
			random_intercept_sigma_raw_2,
			random_intercept_sigma_mu,
			random_intercept_sigma_sigma,
			sigma_correlation_factor_2
		);

    // Update with summing mu_random_intercept
    mu = mu + append_row((X_random_intercept_2 * random_intercept_2)', rep_row_vector(0, N));
  }


  // Calculate proportions
  for(n in 1:N)  mu[,n] = softmax(mu[,n]);

  // Convert the matrix m to a column vector in column-major order.
  mu_array = to_vector(mu);
  precision_array = to_vector(exp(precision));
  

}
model{

  // Fit main distribution
  if(use_data == 1){
    
    target +=  reduce_sum(
      partial_sum_lupmf,
      y_array[truncation_not_idx],
      grainsize,
      exposure_array[truncation_not_idx],
      mu_array[truncation_not_idx],
      precision_array[truncation_not_idx]
    );
    

  }



  // Priors
  if(exclude_priors == 0){
    
    // If interceopt in design or I have complex variability design
    // This would include the models 
    // composition ~ 1 + ...; composition ~ 0 + ...; 
    // variability ~ 1
    if(A == 1){
      target += abundance_variability_regression(
        alpha[1],
        beta[1], // average_by_col(beta[1:B_intercept_columns,]),
        prec_coeff,
        prec_sd,
        bimodal_mean_variability_association,
        mix_p
        );
    }
    else {
      // Loop across the intercept columns in case of a intercept-less design (covariate are intercepts)
      for(a in 1:A_intercept_columns)
        target += abundance_variability_regression(
          alpha[a],
          beta[a],
          prec_coeff,
          prec_sd,
          bimodal_mean_variability_association,
          mix_p
          );
          
        // Variability effect if the formula is more complex
        if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) alpha[a] ~ normal(beta[a] * prec_coeff[2], 2 );
    }

  }
  
  // If I don't have priors for overdispersion
  else{
     // Priors variability
     if(intercept_in_design || A > 1){
       for(a in 1:A_intercept_columns) alpha[a]  ~ normal( prec_coeff[1], prec_sd );
        if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) to_vector(alpha[a]) ~ normal ( 0, 2 );
     }
     // if ~ 0 + covariuate
     else {
       alpha[1]  ~ normal( prec_coeff[1], prec_sd );
     }
  }

  // // Priors abundance
  for(c in 1:B_intercept_columns) beta_raw_raw[c] ~ normal ( prior_mean_intercept[1], prior_mean_intercept[2] );
  if(C>B_intercept_columns) for(c in (B_intercept_columns+1):C) to_vector(beta_raw_raw[c]) ~ normal ( prior_mean_coefficients[1], prior_mean_coefficients[2]);

  // Hyper priors
  mix_p ~ beta(1,5);
  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1],prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);
  prec_sd ~ std_normal();
  prec_coeff ~ std_normal();
  to_vector(beta_raw_raw) ~ std_normal();
  
  // Random intercept
  if(is_random_effect>0){
    for(m in 1:(M-1)) random_intercept_raw[,m] ~ std_normal();
    for(m in 1:(M-1)) random_intercept_sigma_raw[m] ~ std_normal();
    sigma_correlation_factor ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
    
    for(m in 1:(M-1)) random_intercept_raw_2[,m] ~ std_normal();
    for(m in 1:(M-1)) random_intercept_sigma_raw_2[m] ~ std_normal();
    if(is_random_effect>1) sigma_correlation_factor_2 ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
    
    random_intercept_sigma_mu ~ std_normal();
    random_intercept_sigma_sigma ~ std_normal();
    
    // If I have just one group
  	zero_random_intercept ~ std_normal();
  }
}
generated quantities {
  matrix[A, M] alpha_normalised = alpha;
  // Rondom effect
  matrix[ncol_X_random_eff_WINDOWS_BUG_FIX, M] beta_random_intercept;
  matrix[ncol_X_random_eff_WINDOWS_BUG_FIX_2, M] beta_random_intercept_2;

  // LOO
  vector[TNS] log_lik = rep_vector(0, TNS);

  // These instructions regress out the effect of mean proportion to the overdispersion
  // This adjustment provide A overdispersion value that can be tested for a hypotheses for example differences between two conditions
  if(intercept_in_design){
    if(A > 1) for(a in 2:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }
  else{
    for(a in 1:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }

  
  // Random effect
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  if(ncol_X_random_eff[1]==0) beta_random_intercept[1] = rep_row_vector(0.0, M);
  else{
     beta_random_intercept[,1:(M-1)] = random_intercept;
    for(n in 1:ncol_X_random_eff[1]) beta_random_intercept[n, M] = -sum(random_intercept[n,]);
  }
  
    // Random effect 2
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  if(ncol_X_random_eff[2]==0) beta_random_intercept_2[1] = rep_row_vector(0.0, M);
  else{
     beta_random_intercept_2[,1:(M-1)] = random_intercept_2;
    for(n in 1:ncol_X_random_eff[2]) beta_random_intercept_2[n, M] = -sum(random_intercept_2[n,]);
  }
  

  // LOO
  if(enable_loo==1)
    for (n in 1:TNS) {
      log_lik[n] = beta_binomial_lpmf(
        y_array[truncation_not_idx[n]] |
        exposure_array[truncation_not_idx[n]],
        (mu_array[truncation_not_idx[n]] .* precision_array[truncation_not_idx[n]]),
        ((1.0 - mu_array[truncation_not_idx[n]]) .* precision_array[truncation_not_idx[n]])
        ) ;
    }


}
