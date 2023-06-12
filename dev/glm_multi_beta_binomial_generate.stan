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
	int exposure[N];

	// Which column of design, coefficient matrices should be used to generate the data
	int length_X_which;
	int length_XA_which;
	int X_which[length_X_which];
	int XA_which[length_XA_which];
	matrix[N, length_X_which] X;
  matrix[N, length_XA_which] Xa; // The variability design

	int is_truncated;
	real<lower=1> truncation_ajustment;

 // Random intercept
 int length_X_random_intercept_which;
	int X_random_intercept_which[length_X_random_intercept_which];
  int N_random_intercepts;
  int N_minus_sum;
  int paring_cov_random_intercept[N_random_intercepts, 2];
  int N_grouping;
  matrix[N, N_grouping] X_random_intercept;
  int idx_group_random_intercepts[N_grouping, 2];

  // Should I create intercept for generate quantities
  int<lower=0, upper=1> create_intercept;
  int<lower=0> A_intercept_columns;
}
transformed data{
  // If needed recreate the intercept
  matrix[N,1] X_intercept = to_matrix(rep_vector(1, N));

  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int N_grouping_WINDOWS_BUG_FIX = max(N_grouping, 1);
  
  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;

  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));
  
  // If I get crazy diagonal matrix omit it
  if(max(R_ast_inverse)>1000 || N_random_intercepts>0){
   // print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
}

parameters {

  matrix[C, M-1] beta_raw_raw; // matrix with C rows and number of cells (-1) columns
  matrix[A, M] alpha; // Variability
  // To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept // matrix with N_groupings rows and number of cells (-1) columns
  matrix[N_random_intercepts * (N_random_intercepts>0), M-1] random_intercept_raw;
  
  // sd of random intercept
  real random_intercept_sigma_mu[N_random_intercepts>0];
  real random_intercept_sigma_sigma[N_random_intercepts>0];
  row_vector[(M-1) * (N_random_intercepts>0)] random_intercept_sigma_raw;
  
  // If I have just one group
  real zero_random_intercept[N_random_intercepts>0];
  
  // Initialisation
  matrix[C,M] beta_raw;

  // Random effects
  matrix[N_minus_sum, M-1] random_intercept_minus_sum;
  row_vector[M-1] random_intercept_sigma;
  matrix[N_grouping, M-1] beta_random_intercept_raw;

// Random intercept
	matrix[N_grouping_WINDOWS_BUG_FIX, M-1] beta_random_intercept;
	


}
transformed parameters{
  // // If needed recreate the intercept
  // matrix[1,M] beta_intercept;
  // matrix[1,M] alpha_intercept;

}
generated quantities{

  int counts_uncorrected[N, M];

  matrix[C,M] beta = R_ast_inverse * beta_raw; // coefficients on x


  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  real generated_exposure[N];

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
        Q_ast // Rest
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
    mu = (Q_ast * my_beta)';
    precision = (Xa * my_alpha)' / (is_truncated ? truncation_ajustment : 1);

  }

  // Random intercept
  if(length_X_random_intercept_which>0){
      matrix[M, N] mu_random_intercept = append_row((X_random_intercept[,X_random_intercept_which] * beta_random_intercept[X_random_intercept_which,])', rep_row_vector(0, N));
      mu = mu + mu_random_intercept;
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
  for(n in 1:N) generated_exposure[n] = sum(counts_uncorrected[n]);
  for(n in 1:N) counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure[n];
}


