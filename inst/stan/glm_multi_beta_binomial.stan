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

  int[] rep_each(int[] x, int K) {
    int N = size(x);
    int y[N * K];
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }

  matrix average_by_col(matrix beta){
    return to_matrix(
      rep_row_vector(1.0, rows(beta)) * beta / rows(beta),
      1, cols(beta), 0
    );


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
  int<lower=1, upper=N*M> TNS; // truncation_not_size
  int<lower=1, upper=N*M> truncation_not_idx[TNS];

  int<lower=0, upper=1> is_vb;

  // Prior info
  real prior_prec_intercept[2] ;
  real prior_prec_slope[2] ;
  real prior_prec_sd[2] ;

  // Exclude priors for testing purposes
  int<lower=0, upper=1> exclude_priors;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;

  // Does the design icludes intercept
  int <lower=0, upper=1> intercept_in_design;

}
transformed data{
  vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));

  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;

  int y_array[N*M];
  int truncation_down_array[N*M];
  int exposure_array[N*M];

  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));

  // If I get crazy diagonal matrix omit it
  if(max(R_ast_inverse)>1000){
    print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }

  // Data vectorised
  y_array =  to_array_1d(y);
  truncation_down_array = to_array_1d(truncation_down);
  exposure_array = rep_each(exposure, M);
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

  // All this because if A ==1 we have conversion problems
  // This works only with two discrete groups
  if(A == 1) {

    // If I have an intercept in the design matrix
    if(intercept_in_design) beta_intercept_slope = to_matrix(beta[1,], 1, M, 0);
    else beta_intercept_slope = average_by_col(beta);
  }
  else beta_intercept_slope = (XA * beta[1:A,]);
  if(A == 1)  alpha_intercept_slope = alpha;
  else alpha_intercept_slope = (XA * alpha);

}
model{

  // Calculate MU
  matrix[M, N] mu = (Q_ast * beta_raw)';
  vector[N*M] mu_array;
  vector[N*M] precision_array;

  for(n in 1:N) { mu[,n] = softmax(mu[,n]); }

  // Convert the matrix m to a column vector in column-major order.
  mu_array = to_vector(mu);
  precision_array = to_vector(exp(precision));


  if(use_data == 1){
    target += beta_binomial_lpmf(
      y_array[truncation_not_idx] |
      exposure_array[truncation_not_idx],
      (mu_array[truncation_not_idx] .* precision_array[truncation_not_idx]),
      ((1.0 - mu_array[truncation_not_idx]) .* precision_array[truncation_not_idx])
    ) ;
  }

  // Priors
  if(exclude_priors == 0){

    for(i in 1:C) to_vector(beta_raw_raw[i]) ~ normal ( 0, x_raw_sigma );


    // If mean-variability association is bimodal such as for single-cell RNA use mixed model
    if(bimodal_mean_variability_association == 1){
      for (a in 1:A)
      for(m in 1:M)
        target += log_mix(mix_p,
                        normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
                        normal_lpdf(alpha_intercept_slope[a,m] | beta_intercept_slope[a,m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
                      );

    // If no bimodal
    } else {
      to_vector(alpha_intercept_slope) ~ normal(to_vector(beta_intercept_slope) * prec_coeff[2] + prec_coeff[1], prec_sd );
    }

  // If no priors
  } else {
    for(i in 1:C) to_vector(beta_raw_raw[i]) ~ normal ( 0, 2 );
    for (a in 1:A) alpha[a]  ~ normal( 5, 2 );
  }

  // Hyper priors
  mix_p ~ beta(1,5);
  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1], prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);

}

generated quantities {
  matrix[A, M] alpha_normalised = alpha;

  if(intercept_in_design){
    if(A > 1) for(a in 2:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }
  else{
    for(a in 1:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }

}

