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

}
parameters {

	matrix[C,M] beta;
	matrix[A,M] alpha;

}
generated quantities{

  int counts_uncorrected[N, M];

  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  real generated_exposure[N];

  matrix[M,N] mu = (X * beta)';

  matrix[M,N] precision = (Xa * alpha)' / (is_truncated ? truncation_ajustment : 1);
  //matrix[M,N] precision = (X[,1:A] * alpha)'  / (is_truncated ? truncation_ajustment : 1);

	for(i in 1:N) mu[,i] = softmax(mu[,i]);
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


