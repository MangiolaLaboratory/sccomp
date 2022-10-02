data {
	int N;
	int M;
	int C;
	int A;
	int exposure[N];
	matrix[N, C] X;
  matrix[N, A] Xa; // The variability design

	int is_truncated;
	real<lower=1> truncation_ajustment;

	// Which column of design, coefficient matrices should be used to generate the data
	int length_X_which;
	int length_XA_which;
	int X_which[length_X_which];
	int XA_which[length_XA_which];

	// Random intercept
	int length_X_random_intercept_which;
	int X_random_intercept_which[length_X_random_intercept_which];
	int N_grouping;
	matrix[N, N_grouping] X_random_intercept;



}
parameters {

	matrix[C,M] beta;
	matrix[A,M] alpha;

	// Random intercept
	matrix[N_grouping, M-1] beta_random_intercept;

}
generated quantities{

  int counts_uncorrected[N, M];

  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  real generated_exposure[N];

  matrix[M,N] mu = (X[,X_which] * beta[X_which,])';

  // Random intercept
  if(N_grouping>1){
      matrix[M, N] mu_random_intercept = append_row((X_random_intercept[,X_random_intercept_which] * beta_random_intercept[X_random_intercept_which,])', rep_row_vector(0, N));
      mu = mu + mu_random_intercept;
  }

  matrix[M,N] precision = (Xa[,XA_which] * alpha[XA_which,])' / (is_truncated ? truncation_ajustment : 1);



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


