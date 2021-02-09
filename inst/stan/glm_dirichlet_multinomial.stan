functions{

  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    	 real alpha_plus = sum(alpha);

    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

matrix vector_array_to_matrix(vector[] x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
}

vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
}
data{
	int N;
	int M;
	int y[N,M];
	int X[N,2];

}
transformed data{
	int C = 2;
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	vector[M-1] beta_raw[C];
	real<lower=0> precision[2];

	// To exclude

}
transformed parameters{
		vector[M] beta[C];
		real precision_diff = precision[1] - precision[2];
		matrix[M, N] alpha; // for generated quantities. It is cell types in the rows and samples as columns
	  for(c in 1:C)	beta[c] =  sum_to_zero_QR(beta_raw[c], Q_r);

	  // For generated quantities
	  for(n in 1:N) alpha[,n] = (precision[X[n,2]+1] * 100) * softmax( beta[1] + beta[2] * X[n,2] );

}
model{

	 //for(n in 1:N) y[n] ~ dirichlet_multinomial( precision * softmax( vector_array_to_matrix(beta) )) );

	 for(n in 1:N) y[n] ~ dirichlet_multinomial( alpha[,n] );

	 precision ~ normal(0,1);
	 for(i in 1:C) beta_raw[i] ~ normal(0, x_raw_sigma);
}
