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

real partial_sum_lpmf(int[,] y,
                        int start, int end,
                        matrix Q_ast,
                        matrix beta_raw,
                        matrix precision,

                        int[] exposure,
                        int N,
                        int M,
                        int is_truncated
                        ) {


    return neg_binomial_2_log_lupmf(slice_Y | mu[start:end], 1.0 ./ exp(shape[start:end]) );

  real lp = 0;
  vector[N] exp_precision;

  // Calculate MU
  matrix[M, N] mu = (Q_ast * beta_raw)';
  for(n in 1:N) { mu[,n] = softmax(mu[,n]); }



  // NON TRUNCATION
  if(is_truncated == 0){
     exp_precision = exp(precision[,i]);

    for(i in 1:N)
      lp += beta_binomial_lupmf(y[i,] | exposure[i], (mu[,i] .* exp_precision), ((1.0 - mu[,i]) .* exp_precision) );
  }

  // YES TRUNCATION
  else{
    for(i in 1:cols(mu)) { // SAMPLE
      for(j in 1:rows(mu)){ // CATEGORY
        if(truncation_down[i, j] >=0)
         exp_precision = exp(precision[j,i]);

          lp += beta_binomial_lupmf(y[i,j] | exposure[i], (mu[j,i] .* exp_precision), ((1.0 - mu[j,i]) .* exp_precision) );
      }
    }
  }

  return(lp);

}

}
data{
	int<lower=1> N;
	int<lower=1> M;
	int<lower=1> C;
	int<lower=1> A;
	int<lower=1> Ar; // Rows of unique variability design
	int exposure[N];
	int y[N,M];
	matrix[N, C] X;
	matrix[Ar, A] XA; // The unique variability design
	matrix[N, A] Xa; // The variability design

	// Truncation
	int is_truncated;
	int truncation_up[N,M];
	int truncation_down[N,M];

	int<lower=0, upper=1> is_vb;

  // Prior info
  real prior_prec_intercept[2] ;
  real prior_prec_slope[2] ;
  real prior_prec_sd[2] ;

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
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));

  // If I get crazy diagonal matrix omit it
  if(max(R_ast_inverse)>1000){
    print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
}
parameters{
	matrix[C, M-1] beta_raw_raw;
	matrix[A, M] alpha;

	// To exclude
  real prec_coeff[2];
  real<lower=0> prec_sd;

  real<lower=0, upper=1> mix_p;
}
transformed parameters{
		matrix[C,M] beta_raw;
		matrix[A,M] beta_intercept_slope;
		matrix[A,M] alpha_intercept_slope;
    matrix[M, N] precision = (Xa * alpha)';
    matrix[C,M] beta;

	  for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);


		// Beta
	  beta = R_ast_inverse * beta_raw; // coefficients on x

    // All this because if A ==1 we have ocnversion problems
    // This works only with two discrete groups
    if(A == 1) beta_intercept_slope = to_matrix(beta[A,], A, M, 0);
    else beta_intercept_slope = (XA * beta[1:A,]);
		if(A == 1)  alpha_intercept_slope = alpha;
		else alpha_intercept_slope = (XA * alpha);

}
model{

   target += reduce_sum(
    partial_sum_lupmf,
    y,
    grainsize,
    Q_ast,
    beta_raw,
    precision,

    exposure,
    N,
    M,
    is_truncated
  );

  // Priors

  if(exclude_priors == 0){

    for(i in 1:C) to_vector(beta_raw_raw[i]) ~ normal ( 0, x_raw_sigma );

    // PRECISION REGRESSION
    // to_vector(alpha_intercept_slope) ~ student_t( 8, to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd);
    // if(is_vb==0){
    for (a in 1:A) for(m in 1:M)
      target += log_mix(mix_p,
                      normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
                      normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
                    );
  } else {
    for(i in 1:C) to_vector(beta_raw_raw[i]) ~ normal ( 0, 5 );
    for (a in 1:A) alpha[a]  ~ normal( 5, 5 );
  }

  // Hyper priors
  mix_p ~ beta(1,5);

  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1], prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);

}

generated quantities {
  matrix[A, M] alpha_normalised = alpha;


  if(A > 1) for(a in 2:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );

}

