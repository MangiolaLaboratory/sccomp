data {
	int<lower=0> N;
	int<lower=0> M;
	int<lower=0> C;
	array[N] int exposure;
	matrix[N,C] X;
}
parameters {

	matrix[C, M] beta;


}
generated quantities{

  array[N, M] int counts;

  matrix[N, M] alpha = (X * beta);

	for(n in 1:N)
	counts[n] = multinomial_rng(dirichlet_rng(to_vector(alpha[n])), exposure[n]) ;

}
