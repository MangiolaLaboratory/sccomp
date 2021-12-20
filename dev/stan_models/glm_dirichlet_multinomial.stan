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

  row_vector sum_to_zero_QR(row_vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    row_vector [N] x;
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
	int C;
	int A;
	int y[N,M];
	matrix[N,C] X;
}
transformed data{

	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	matrix[C, M-1] beta_raw;
	vector<upper=10>[A] precision;

	// To exclude

}
transformed parameters{
  real buffer;
  real plateau = 0.5;
	matrix[C,M] beta;
	vector[N] precision_per_sample = X[,1:A] * precision;

	//real precision_diff = precision[1] - precision[2];
	matrix[ N, M] alpha; // for generated quantities. It is cell types in the rows and samples as columns
  for(c in 1:C)	beta[c] =  sum_to_zero_QR(beta_raw[c], Q_r);

  // For generated quantities
  alpha = (X * beta); // N x M

	for(n in 1:N){
	  alpha[n] = to_row_vector(softmax(to_vector(alpha[n])));
		buffer = 1.0/min(alpha[n]) * plateau;
	  alpha[n] =  alpha[n] *  exp(precision_per_sample[n]) ;
	}


}
model{

	 //for(n in 1:N) y[n] ~ dirichlet_multinomial( precision * softmax( vector_array_to_matrix(beta) )) );

	 for(n in 1:N) y[n] ~ dirichlet_multinomial( to_vector(alpha[n] ));

	 precision ~ normal(0,5);
	 for(i in 1:C) beta_raw[i] ~ normal(0, x_raw_sigma );
}
