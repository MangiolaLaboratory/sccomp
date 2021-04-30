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

	vector[] beta_regression_rng( matrix X, matrix beta, vector phi){

		vector[cols(beta)] p[rows(X)];


    matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';



		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {

     	p[i] = to_vector(beta_rng((mu[,i] .* phi), ((1.0 - mu[,i]) .* phi)));

		}

		return (p);
	}


}
data{
	int N;
	int M;
	int C;
	matrix[N, 2] X;
	vector[M] beta[C];
	vector<lower=0>[M] precision;

}
transformed data{

	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
generated quantities{

  vector[M] y[N];

  y = beta_regression_rng(  X, vector_array_to_matrix(beta),  precision);
}
