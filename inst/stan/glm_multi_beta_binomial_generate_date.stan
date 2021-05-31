data {
	int N;
	int M;
	int C;
	int A;
	int exposure[N];
	matrix[N, C] X;

	int is_truncated;
	real<lower=1> truncation_ajustment;

}
parameters {

	matrix[C,M] beta;
	matrix[A,M] alpha;

}
generated quantities{

  int counts[N, M];
  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision = (X[,1:A] * alpha)'  / (is_truncated ? truncation_ajustment : 1);

	for(i in 1:N) mu[,i] = softmax(mu[,i]);
	for(i in 1:cols(mu)) {
    	counts[i,] = beta_binomial_rng(
    	  exposure[i],
    	  mu[,i] .* exp(precision[,i]),
    	  (1.0 - mu[,i]) .* exp(precision[,i])
    	 );
	}

}
