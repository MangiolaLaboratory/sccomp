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

  real beta_binomial_regression_lpmf(int[,] p, int[] exposure, matrix X, matrix beta, vector phi){

		real lp = 0;
    matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {

     	lp += beta_binomial_lpmf(p[i] | exposure[i], (mu[,i] .* phi), ((1.0 - mu[,i]) .* phi) );

		}
		return (lp);
	}

	 real beta_binomial_truncated_lpmf(int p, int exposure, real alpha, real beta, int lower, int upper){

		real lp = 0;

    lp += beta_binomial_lpmf(p | exposure, alpha, beta );
    if (p < lower || p > upper)
      lp += negative_infinity();
    else
      lp += -log_sum_exp(
        beta_binomial_lpmf(lower | exposure, alpha, beta),
        log_diff_exp(beta_binomial_lcdf(upper | exposure, alpha, beta),
        beta_binomial_lcdf(lower | exposure, alpha, beta))
      );

		return (lp);
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
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
}
parameters{
	matrix[C, M-1] beta_raw;
	vector[M] precision;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;
}
transformed parameters{
		matrix[C,M] beta;
	  for(c in 1:C)	beta[c,] =  sum_to_zero_QR(beta_raw[c,], Q_r);


}
model{

   // y ~ beta_binomial_regression(exposure, X, beta, exp(precision) );

//     matrix[num_elements(y[1]), num_elements(y[,1])] mu = (X * beta)';
// 		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
// 		for(i in 1:cols(mu)) {
//
//      	target += beta_binomial_lpmf(y[i] | exposure[i], (mu[,i] .* exp(precision)), ((1.0 - mu[,i]) .* exp(precision)) );
//
// 		}

  matrix[num_elements(y[1]), num_elements(y[,1])] mu = (X * beta)';

	for(i in 1:cols(mu)) {

    mu[,i] = softmax(mu[,i]);

    for(j in 1:rows(mu)){

      // If I have outlier exclusion
      if(is_truncated == 1 && truncation_down[i, j] != -99) {
        if(truncation_down[i, j] > 0) // If it is < 0 it is an outlier
        target += beta_binomial_truncated_lpmf(
          y[i,j] |
          exposure[i],
          (mu[j,i] .* exp(precision[j])),
          ((1.0 - mu[j,i]) .* exp(precision[j])) ,
          truncation_down[i, j],
          truncation_up[i, j]
        );

      }

      // If I don't have any outlier elimination
      else
       	target += beta_binomial_lpmf(y[i,j] | exposure[i], (mu[j,i] .* exp(precision[j])), ((1.0 - mu[j,i]) .* exp(precision[j])) );


    }

	}

    precision ~ normal( beta[1] * prec_coeff[2] + prec_coeff[1], prec_sd);
    prec_sd ~ normal(0,2);
    prec_coeff ~ normal(0,5);

    for(i in 1:C) to_vector(beta_raw[i]) ~ normal(0, x_raw_sigma );


}
