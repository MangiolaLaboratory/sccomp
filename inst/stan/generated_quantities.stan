data {
	int<lower=0> N;
	int<lower=0> M;
	int exposure[N];
}
parameters {

	matrix[N, M] alpha;

}
generated quantities{

  int counts[N, M];

	for(n in 1:N)
	counts[n] = multinomial_rng(dirichlet_rng(to_vector(alpha[n])), exposure[n]) ;

}
