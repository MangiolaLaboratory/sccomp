functions{

  int[,] beta_binomial_regression_rng(int[] tot, matrix X, matrix beta, vector phi){

		int y[cols(beta),rows(X)];
    matrix[cols(beta),rows(X)] mu = (X * beta)';
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {

     	y[,i] = beta_binomial_rng( tot[i], (mu[,i] .* phi), ((1.0 - mu[,i]) .* phi) );

		}
		return (y);
	}

}
data {
	int N;
	int M;
	int C;
	int tot[N];
	matrix[N, C] X;
}
parameters {

	matrix[C,M] beta;
  vector[M] precision;

}
generated quantities{

  int counts[N, M];

	counts = beta_binomial_regression_rng(tot, X, beta, exp(precision)) ;

}
