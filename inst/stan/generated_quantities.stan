functions{

	int[] dirichlet_multinomial_rng(vector alpha, int exposure) {
	    return multinomial_rng(dirichlet_rng(alpha), exposure);
	}

}
data {
	int<lower=0> N;
	int<lower=0> M;
	int exposure[N];
}
parameters {

	matrix[M, N] alpha;

}
generated quantities{

  int counts[N, M];

	for(n in 1:N)
	counts[n] = dirichlet_multinomial_rng(alpha[,n], exposure[n]) ;

}
