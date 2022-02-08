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


}data {
	int N;
	int M;
	int C;
	int A;
	int exposure[N];
	matrix[N, C] X;

	int is_truncated;
	real<lower=1> truncation_ajustment;

}
transformed data{
  matrix[N, C] Q_ast;
  matrix[C, C] R_ast_inverse;
  vector[2*M] Q_r = Q_sum_to_zero_QR(M);

  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));
}
parameters {

	matrix[C, M-1] beta_raw_raw;
	matrix[A, M] alpha;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;

  real<lower=0, upper=1> mix_p;

}
transformed parameters{
  	matrix[C,M] beta_raw;
    matrix[C,M] beta;

    for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);
    beta = R_ast_inverse * beta_raw; // coefficients on x
}
generated quantities{

  int counts_uncorrected[N, M];

  // Matrix for correcting for exposure
  matrix[N, M] counts;

  // Vector of the generated exposure
  real generated_exposure[N];

  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision = (X[,1:A] * alpha)'  / (is_truncated ? truncation_ajustment : 1);

	for(i in 1:N) mu[,i] = softmax(mu[,i]);
	for(i in 1:N) {
    	counts_uncorrected[i,] = beta_binomial_rng(
    	  exposure[i],
    	  mu[,i] .* exp(precision[,i]),
    	  (1.0 - mu[,i]) .* exp(precision[,i])
    	 );
	}

	// Calculate the generated exposure
  for(n in 1:N) generated_exposure[n] = sum(counts_uncorrected[n]);
  for(n in 1:N) counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure[n];
}


