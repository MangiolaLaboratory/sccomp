functions{

	int[] dirichlet_multinomial_rng(vector alpha, int exposure) {
	    return multinomial_rng(dirichlet_rng(alpha), exposure);
	}

matrix vector_array_to_matrix(vector[] x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
}



}
data{
	int N;
	int M;
	int C;
	int A;
	vector[M] y[N];
	matrix[N,C] X;

	// To exclude
	int<lower=0> how_namy_to_include;
	int to_include[how_namy_to_include, 2]; // Table with column N and M

  // RNG
	int I; // iterations
	vector[A] precision[I];
	int exposure[N];

}
parameters{
		matrix[C,M] beta;
	vector<lower=0>[A] sigma;

	// // To exclude
	// vector<lower=0, upper=1>[how_namy_to_exclude] excluded;

}
transformed parameters{
		matrix[ N,M] alpha = X * beta; // for generated quantities. It is cell types in the rows and samples as columns
}
model{

  for(i in 1:how_namy_to_include)
    y[to_include[i,1], to_include[i,2]] ~ normal(
      alpha[ to_include[i,1], to_include[i,2]],
      X[to_include[i,1], 1:A] * sigma
    );


	 sigma ~ normal(0,1);
	 for(i in 1:C) beta[i] ~ normal(0, 5);
}
generated quantities{
  vector[M] y_rng[N];
  vector[M] y_simplex[N];
  int counts[N, M];

  // Random precision
  int my_n;
  real unif = uniform_rng(0,1);
  real idx = unif * I;

  // Get proportions
  for(n in 1:N) for(m in 1:M) {
      y_rng[n, m] = normal_rng( alpha[n,m],  X[n, 1:A] * sigma  );
      y_simplex[n] = softmax(y_rng[n] );
    }

  // Random precision
  for (n in 1:I) if (n > idx) { my_n = n; break; }


	for(n in 1:N)
	  counts[n] = dirichlet_multinomial_rng((X[n,1:A] * precision[my_n] * 100) * y_simplex[n], exposure[n]) ;



}
