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

}
data{
	int N;
	int M;
	int C;
	int exposure[N];
	int y[N,M];
	matrix[N, C] X;
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
	vector[M] precision_log;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;
}
transformed parameters{
		matrix[C,M] beta;
		vector[M] precision ;
	  for(c in 1:C)	beta[c,] =  sum_to_zero_QR(beta_raw[c,], Q_r);

  // non centrred parametrisation of regression
   precision = exp(precision_log * prec_sd + to_vector(beta[1] * prec_coeff[2] + prec_coeff[1]) );
}
model{

    matrix[num_elements(y[1]), num_elements(y[,1])] mu = (X * beta)';

	  y ~ beta_binomial_regression(exposure, X, beta, precision) ; // truncation );

	 precision_log ~ std_normal();
   prec_sd ~ normal(0,2);
   prec_coeff ~ normal(0,5);

	 for(c in 1:C) to_vector(beta_raw[c]) ~ normal(0, x_raw_sigma * 5);

    for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {

     //	y[i] ~ beta_binomial( exposure[i], (mu[,i] .* precision), ((1.0 - mu[,i]) .* precision) );

		}

// 	for(i in 1:cols(mu)) {
//
//     mu[,i] = softmax(mu[,i]);
//
//     for(j in 1:rows(mu)){
//
//       // If I have outlier exclusion
//       if(is_truncated == 1) {
//         if(truncation_down[i, j] > 0) // If it is < 0 it is an outlier
//         y[i,j] ~ beta_binomial(
//           exposure[i],
//           (mu[j,i] .* precision[j]),
//           ((1.0 - mu[j,i]) .* precision[j])
//         ) T[truncation_down[i, j],truncation_up[i, j]];
//       }
//
//       // If I don't have any outlier elimination
//       else
//         y[i,j] ~ beta_binomial(
//           exposure[i],
//           (mu[j,i] .* precision[j]),
//           ((1.0 - mu[j,i]) .* precision[j])
//         ) ;
//
//     }
//
// 	}
}
