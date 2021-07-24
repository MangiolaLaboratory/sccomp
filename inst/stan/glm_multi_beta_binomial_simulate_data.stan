functions{
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
	int N; // Number of subjects
	int M; // Number of categories
	int C;
	int A;
	int exposure[N];
	matrix[N, C] X;
	matrix[A, A] XA;

  matrix[C,M] beta;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;

}
transformed data{
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	matrix[A, M] alpha;


}
transformed parameters{

		matrix[A,M] beta_intercept_slope;
  	matrix[A,M] alpha_intercept_slope;

    // All this because if A ==1 we have ocnversion problems
    if(A == 1) beta_intercept_slope = to_matrix(beta[A,], A, M, 0);
    else beta_intercept_slope = (XA * beta[1:A,])';
		if(A == 1)  alpha_intercept_slope = alpha;
		else alpha_intercept_slope = (XA * alpha)';

}
model{

  // PRECISION REGRESSION
  // to_vector(alpha_intercept_slope) ~ student_t( 8, to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd);
  to_vector(alpha_intercept_slope) ~ normal( to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd);


}
generated quantities{

  int counts[N, M];
  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision = (X[,1:A] * alpha)';

	for(i in 1:N) mu[,i] = softmax(mu[,i]);
	for(i in 1:cols(mu)) {
    	counts[i,] = beta_binomial_rng(
    	  exposure[i],
    	  mu[,i] .* exp(precision[,i]),
    	  (1.0 - mu[,i]) .* exp(precision[,i])
    	 );
	}

}
