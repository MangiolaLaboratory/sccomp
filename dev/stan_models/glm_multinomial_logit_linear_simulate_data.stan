data{
  int N; // Number of subjects
  int M; // Number of categories
  int C;
  int A;
  int exposure[N];
  matrix[N, C] X;
  matrix[A, A] XA;

  matrix[C,M] beta;

  real variability_multiplier;

}

parameters{

  real prec_coeff[2];
  real<lower=0> prec_sd;
}

generated quantities{
  int counts[N, M];
  matrix[M, N] normal_draws;

  matrix[A,M] alpha;
  matrix[M,N] mu = (X * beta)';
  matrix[M,N] precision;

  matrix[A,M] beta_intercept_slope;


  // All this because if A ==1 we have ocnversion problems
  // This works only with two discrete groups
  if(A == 1) beta_intercept_slope = to_matrix(beta[A,], A, M, 0);
  else beta_intercept_slope = (XA * beta[1:A,]);

  // PRECISION REGRESSION
  for(a in 1:A) for(m in 1:M) alpha[a,m] = normal_rng( beta_intercept_slope[a,m] * -prec_coeff[2] + prec_coeff[1], prec_sd);

  precision = (X[,1:A] * alpha)';

   for(i in 1:cols(mu)) {
      normal_draws[,i] = to_vector( normal_rng(mu[,i],  precision[,i] / variability_multiplier )  ); // sd decreased because different representation from beta binomial
    }

  for(i in 1:cols(normal_draws)) {
    counts[i,] = multinomial_rng(  softmax(normal_draws[,i]), exposure[i]);
  }



}
