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
	int N;
	int M;
	int C;
	int exposure[N];
	int y[N,M];
	matrix[N, C] X;

	// Truncation
	int is_truncated;
	int truncation_up[N,M];
	int truncation_down[N,M];


}
transformed data{
  int A = 2;
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	matrix[C, M-1] beta_raw;
	matrix[A, M] alpha;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;
}
transformed parameters{
		matrix[C,M] beta;
		matrix[A,M] beta_intercept_slope;
		matrix[A,M] alpha_intercept_slope;
    matrix[num_elements(y[1]), num_elements(y[,1])] precision = (X[,1:A] * alpha)';

	  for(c in 1:C)	beta[c,] =  sum_to_zero_QR(beta_raw[c,], Q_r);

    beta_intercept_slope[1] = beta[1];
    beta_intercept_slope[2] = beta[1] + beta[2];
    alpha_intercept_slope[1] = alpha[1];
    alpha_intercept_slope[2] = alpha[1] + alpha[2];

}
model{

  // Calculate MU
  matrix[num_elements(y[1]), num_elements(y[,1])] mu = (X * beta)';

  for(i in 1:cols(mu)) { mu[,i] = softmax(mu[,i]); }

  // NON TRUNCATION
  if(is_truncated == 0){
    for(i in 1:cols(mu))
      y[i,] ~ beta_binomial( exposure[i], (mu[,i] .* exp(precision[,i])), ((1.0 - mu[,i]) .* exp(precision[,i])) );
  }

  // YES TRUNCATION
  else{
    for(i in 1:cols(mu)) { // SAMPLE
      for(j in 1:rows(mu)){ // CATEGORY
        if(truncation_down[i, j] >=0)
          y[i,j] ~ beta_binomial( exposure[i], (mu[j,i] .* exp(precision[j,i])), ((1.0 - mu[j,i]) .* exp(precision[j,i])) );
      }
    }
  }

  for(i in 1:C) to_vector(beta_raw[i]) ~ student_t (8, 0, x_raw_sigma );

  // PRECISION REGRESSION
  to_vector(alpha_intercept_slope) ~ normal( to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd);
  prec_sd ~ normal(0,2);
  prec_coeff ~ normal(0,5);



}

