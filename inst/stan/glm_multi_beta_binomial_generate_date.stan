functions{
    row_vector average_by_col(matrix beta){
    return
    rep_row_vector(1.0, rows(beta)) * beta / rows(beta);
  }
  }
data {
	int N;
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

	int is_truncated;
	real<lower=1> truncation_ajustment;

	// Random intercept
	array[2] int length_X_random_intercept_which;
	array[length_X_random_intercept_which[1]] int X_random_intercept_which;
	array[2] int ncol_X_random_eff;
	array[2] int  ncol_X_random_eff_new;
	matrix[N, ncol_X_random_eff_new[1]] X_random_intercept;
	matrix[N, ncol_X_random_eff_new[2]] X_random_intercept_2;
	array[length_X_random_intercept_which[2]] int X_random_intercept_which_2;

  // Should I create intercept for generate quantities
  int<lower=0, upper=1> create_intercept;
  int<lower=0> A_intercept_columns;
}
transformed data{
  // If needed recreate the intercept
  matrix[N,1] X_intercept = to_matrix(rep_vector(1, N));

  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);
}
parameters {

	matrix[C,M] beta;
	matrix[A,M] alpha;

	// Random intercept
	matrix[ncol_X_random_eff_WINDOWS_BUG_FIX, M] beta_random_intercept;
	matrix[ncol_X_random_eff_WINDOWS_BUG_FIX_2, M] beta_random_intercept_2;

}
transformed parameters{
  // // If needed recreate the intercept
  // matrix[1,M] beta_intercept;
  // matrix[1,M] alpha_intercept;

}
generated quantities{

  array[N, M] int counts_uncorrected;


  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  array[N] real generated_exposure;

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

  // Random intercept
  if(length_X_random_intercept_which[1]>0){
      matrix[M, N] mu_random_intercept = (X_random_intercept * beta_random_intercept[X_random_intercept_which,])';
      mu = mu + mu_random_intercept;
  }
  
    // Random intercept 2
  if(length_X_random_intercept_which[2]>0){
      matrix[M, N] mu_random_intercept_2 = (X_random_intercept_2 * beta_random_intercept_2[X_random_intercept_which_2,])';
      mu = mu + mu_random_intercept_2;
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


