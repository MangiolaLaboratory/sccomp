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
	int A;
	array[N] int exposure;
	array[N,M] int y;
	matrix[N, C] X;
	matrix[A, A] XA;

	// Truncation
	int is_truncated;
	array[N,M] int truncation_up;
	array[N,M] int truncation_down;

	int<lower=0, upper=1> is_vb;

  // Prior info
  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;
  array[2] real prior_prec_sd;

  // Exclude priors for testing purposes
  int<lower=0, upper=1> exclude_priors;

}
transformed data{
	vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));

  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast = qr_thin_R(X) / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
}
parameters{
	matrix[C, M] beta_raw;
	matrix[A, M] alpha;

	// To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;

  real<lower=0, upper=1> mix_p;
}
transformed parameters{
		matrix[A,M] beta_intercept_slope;
		matrix[A,M] alpha_intercept_slope;
    matrix[M, N] precision = (X[,1:A] * alpha)';

    // All this because if A ==1 we have ocnversion problems
    // This works only with two discrete groups
    if(A == 1) beta_intercept_slope = to_matrix(beta_raw[A,], A, M, 0);
    else beta_intercept_slope = (XA * beta_raw[1:A,]);
		if(A == 1)  alpha_intercept_slope = alpha;
		else alpha_intercept_slope = (XA * alpha);

}
model{

  // Calculate MU
  matrix[M, N] mu = (Q_ast * beta_raw)';

  for(n in 1:N) { mu[,n] = inv_logit(mu[,n]); }

  // NON TRUNCATION
  if(is_truncated == 0){
    for(i in 1:N)
      target += beta_binomial_lpmf(y[i,] | exposure[i], (mu[,i] .* exp(precision[,i])), ((1.0 - mu[,i]) .* exp(precision[,i])) );
  }

  // YES TRUNCATION
  else{
    for(i in 1:cols(mu)) { // SAMPLE
      for(j in 1:rows(mu)){ // CATEGORY
        if(truncation_down[i, j] >=0)
          target += beta_binomial_lpmf(y[i,j] | exposure[i], (mu[j,i] .* exp(precision[j,i])), ((1.0 - mu[j,i]) .* exp(precision[j,i])) );
      }
    }
  }

  // Priors
  for(i in 1:C) to_vector(beta_raw[i]) ~ normal ( 0, 5 );

  if(exclude_priors == 0){


    // PRECISION REGRESSION
    // to_vector(alpha_intercept_slope) ~ student_t( 8, to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd);
    // if(is_vb==0){
    for (a in 1:A) for(m in 1:M)
      target += log_mix(mix_p,
                      normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
                      normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
                    );
  } else {
    for (a in 1:A) alpha_intercept_slope[a]  ~ normal( 5, 5 );
  }

  // Hyper priors
  mix_p ~ beta(1,5);

  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1], prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);

}

generated quantities {
  matrix[C,M] beta;
  beta = R_ast_inverse * beta_raw; // coefficients on x
}
