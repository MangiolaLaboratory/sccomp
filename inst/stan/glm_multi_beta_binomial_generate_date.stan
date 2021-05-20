data {
	int N;
	int M;
	int C;
	int A;
	int exposure[N];
	matrix[N, C] X;
}
parameters {

	matrix[C,M] beta;
	matrix[A,M] alpha;

}
generated quantities{

  int counts[N, M];
  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision = (X[,1:A] * alpha)';

	for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
	for(i in 1:cols(mu)) {
    	counts[i,] = beta_binomial_rng( exposure[i], (mu[,i] .* precision[,i]), ((1.0 - mu[,i]) .* precision[,i]) );
	}

}
