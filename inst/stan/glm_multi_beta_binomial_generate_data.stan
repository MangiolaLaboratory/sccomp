functions{
    #include common_functions.stan

}
data {
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int N;
  int N_original;
  int M;
  int C;
  int A;
  array[N] int exposure;

  // Which columns of fixed / variability design enter the linear predictors
  int length_X_which;
  int length_XA_which;
  array[length_X_which] int X_which;
  array[length_XA_which] int XA_which;
  matrix[N, length_X_which] X;
  matrix[N, length_XA_which] Xa;

  matrix[N_original, C] X_original;
  int is_truncated;
  real<lower=1> truncation_ajustment;

  int is_random_effect;

  // Four uniform random-effect slots (same layout as glm_multi_beta_binomial.stan)
  array[4] int ncol_X_random_eff;
  array[4] int ncol_X_random_eff_new;

  matrix[N, ncol_X_random_eff_new[1]] X_random_effect_1;
  matrix[N, ncol_X_random_eff_new[2]] X_random_effect_2;
  matrix[N, ncol_X_random_eff_new[3]] X_random_effect_3;
  matrix[N, ncol_X_random_eff_new[4]] X_random_effect_4;

  array[4] int n_groups;
  array[4] int how_many_factors_in_random_design;

  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance_1;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
  array[how_many_factors_in_random_design[3], n_groups[3]] int group_factor_indexes_for_covariance_3;
  array[how_many_factors_in_random_design[4], n_groups[4]] int group_factor_indexes_for_covariance_4;

  // Per-slot column counts for X_random_effect_which_* arrays (from R: ncol per slot)
  array[4] int length_X_random_effect_which;
  array[length_X_random_effect_which[1]] int X_random_effect_which_1;
  array[length_X_random_effect_which[2]] int X_random_effect_which_2;
  array[length_X_random_effect_which[3]] int X_random_effect_which_3;
  array[length_X_random_effect_which[4]] int X_random_effect_which_4;

  int<lower=0, upper=1> create_intercept;
  int<lower=0> A_intercept_columns;

  array[4] int<lower=0, upper=1> unknown_grouping;

  array[4] int ncol_X_random_eff_unseen;
  matrix[N, ncol_X_random_eff_unseen[1]] X_random_effect_1_unseen;
  matrix[N, ncol_X_random_eff_unseen[2]] X_random_effect_2_unseen;
  matrix[N, ncol_X_random_eff_unseen[3]] X_random_effect_3_unseen;
  matrix[N, ncol_X_random_eff_unseen[4]] X_random_effect_4_unseen;
}
transformed data {
  matrix[N,1] X_intercept;
  X_intercept = to_matrix(rep_vector(1, N));

  int ncol_X_random_eff_safe_1 = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_safe_2 = max(ncol_X_random_eff[2], 1);
  int ncol_X_random_eff_safe_3 = max(ncol_X_random_eff[3], 1);
  int ncol_X_random_eff_safe_4 = max(ncol_X_random_eff[4], 1);
}
parameters {
  // Unconstrained vectors so posterior CSVs (rounded sig_figs) still validate;
  // we apply sum-to-zero in transformed parameters (same idea as the old GQ model).
  array[C] vector[M] beta_raw;
  matrix[A, M] alpha;
  array[A] ordered[1 + bimodal_mean_variability_association] prec_intercept;
  array[A] real prec_slope_1;
  array[A * bimodal_mean_variability_association] real prec_slope_2;
  array[A] real log_prec_sd;
  real<lower=0, upper=1> mix_p;

  array[ncol_X_random_eff[1] * (ncol_X_random_eff[1]>0)] vector[M] random_effect_raw_1;
  array[M * (ncol_X_random_eff[1]>0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma_raw_1;
  array[M * (ncol_X_random_eff[1]>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (ncol_X_random_eff[1]>0)] sigma_correlation_factor_1;

  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0)] vector[M] random_effect_raw_2;
  array[M * (ncol_X_random_eff[2]>0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_raw_2;
  array[M * (ncol_X_random_eff[2]>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (ncol_X_random_eff[2]>0)] sigma_correlation_factor_2;

  array[ncol_X_random_eff[3] * (ncol_X_random_eff[3]>0)] vector[M] random_effect_raw_3;
  array[M * (ncol_X_random_eff[3]>0)] vector[how_many_factors_in_random_design[3]] random_effect_sigma_raw_3;
  array[M * (ncol_X_random_eff[3]>0)] cholesky_factor_corr[how_many_factors_in_random_design[3] * (ncol_X_random_eff[3]>0)] sigma_correlation_factor_3;

  array[ncol_X_random_eff[4] * (ncol_X_random_eff[4]>0)] vector[M] random_effect_raw_4;
  array[M * (ncol_X_random_eff[4]>0)] vector[how_many_factors_in_random_design[4]] random_effect_sigma_raw_4;
  array[M * (ncol_X_random_eff[4]>0)] cholesky_factor_corr[how_many_factors_in_random_design[4] * (ncol_X_random_eff[4]>0)] sigma_correlation_factor_4;

  array[4 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[4 * (is_random_effect>0)] real random_effect_sigma_sigma;
  array[is_random_effect>0] real zero_random_effect;
}
transformed parameters {
  array[A] real prec_intercept_1;
  array[A * bimodal_mean_variability_association] real prec_intercept_2;
  array[A] real<lower=0> prec_sd;
  for (a in 1:A) {
    prec_intercept_1[a] = prec_intercept[a][1];
    prec_sd[a] = exp(log_prec_sd[a]);
    if (bimodal_mean_variability_association == 1)
      prec_intercept_2[a] = prec_intercept[a][2];
  }

  matrix[C,M] beta;
  for (c in 1:C) {
    beta[c,] = to_row_vector(normalize_sum_to_zero(beta_raw[c]));
  }

  real mix_p_scalar = bimodal_mean_variability_association == 1 ? mix_p : 0.5;

  matrix[ncol_X_random_eff_safe_1 * (is_random_effect>0), M] random_effect_1;
  matrix[ncol_X_random_eff_safe_2 * (is_random_effect>0), M] random_effect_2;
  matrix[ncol_X_random_eff_safe_3 * (is_random_effect>0), M] random_effect_3;
  matrix[ncol_X_random_eff_safe_4 * (is_random_effect>0), M] random_effect_4;

  if (ncol_X_random_eff[1] > 0) {
    array[ncol_X_random_eff[1]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[1]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_1[i]);
    random_effect_1 = build_re_block(
      M, n_groups[1], how_many_factors_in_random_design[1], ncol_X_random_eff[1],
      group_factor_indexes_for_covariance_1, raw_vec,
      random_effect_sigma_mu[1], random_effect_sigma_sigma[1],
      random_effect_sigma_raw_1, sigma_correlation_factor_1
    );
  }
  if (ncol_X_random_eff[2] > 0) {
    array[ncol_X_random_eff[2]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[2]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_2[i]);
    random_effect_2 = build_re_block(
      M, n_groups[2], how_many_factors_in_random_design[2], ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2, raw_vec,
      random_effect_sigma_mu[2], random_effect_sigma_sigma[2],
      random_effect_sigma_raw_2, sigma_correlation_factor_2
    );
  }
  if (ncol_X_random_eff[3] > 0) {
    array[ncol_X_random_eff[3]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[3]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_3[i]);
    random_effect_3 = build_re_block(
      M, n_groups[3], how_many_factors_in_random_design[3], ncol_X_random_eff[3],
      group_factor_indexes_for_covariance_3, raw_vec,
      random_effect_sigma_mu[3], random_effect_sigma_sigma[3],
      random_effect_sigma_raw_3, sigma_correlation_factor_3
    );
  }
  if (ncol_X_random_eff[4] > 0) {
    array[ncol_X_random_eff[4]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[4]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_4[i]);
    random_effect_4 = build_re_block(
      M, n_groups[4], how_many_factors_in_random_design[4], ncol_X_random_eff[4],
      group_factor_indexes_for_covariance_4, raw_vec,
      random_effect_sigma_mu[4], random_effect_sigma_sigma[4],
      random_effect_sigma_raw_4, sigma_correlation_factor_4
    );
  }
}
model {
}
generated quantities {
  array[N, M] int counts_uncorrected;
  matrix[N, M] counts;
  array[N] real generated_exposure;

  matrix[length_X_which, M] my_beta = beta[X_which, ];
  matrix[length_XA_which, M] my_alpha = alpha[XA_which, ];

  matrix[M, N] mu;
  matrix[M, N] precision;

  if (create_intercept == 1) {
    mu = (
      append_col(
        to_matrix(rep_vector(1, N)),
        X
        ) *
        append_row(
          average_by_col(beta[1:A_intercept_columns, ]),
          my_beta
        )
    )';

    precision = (
      append_col(
        to_matrix(rep_vector(1, N)),
        Xa
        ) *
        append_row(
          average_by_col(alpha[1:A_intercept_columns, ]),
          my_alpha
          )
    )' /
    (is_truncated ? truncation_ajustment : 1);
  } else {
    mu = (X * my_beta)';
    precision = (Xa * my_alpha)' / (is_truncated ? truncation_ajustment : 1);
  }

  if (ncol_X_random_eff[1] > 0) {
    mu = mu + (X_random_effect_1 * random_effect_1[X_random_effect_which_1, ])';
    if (ncol_X_random_eff_unseen[1] > 0) {
      matrix[ncol_X_random_eff_unseen[1], M] unseen =
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[1] * M),
                 ncol_X_random_eff_unseen[1], M);
      for (i in 1:ncol_X_random_eff_unseen[1])
        unseen[i, ] = to_row_vector(normalize_sum_to_zero(to_vector(unseen[i, ])));
      mu = mu + (X_random_effect_1_unseen * unseen)';
    }
  }
  if (ncol_X_random_eff[2] > 0) {
    mu = mu + (X_random_effect_2 * random_effect_2[X_random_effect_which_2, ])';
    if (ncol_X_random_eff_unseen[2] > 0) {
      matrix[ncol_X_random_eff_unseen[2], M] unseen =
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[2] * M),
                 ncol_X_random_eff_unseen[2], M);
      for (i in 1:ncol_X_random_eff_unseen[2])
        unseen[i, ] = to_row_vector(normalize_sum_to_zero(to_vector(unseen[i, ])));
      mu = mu + (X_random_effect_2_unseen * unseen)';
    }
  }
  if (ncol_X_random_eff[3] > 0) {
    mu = mu + (X_random_effect_3 * random_effect_3[X_random_effect_which_3, ])';
    if (ncol_X_random_eff_unseen[3] > 0) {
      matrix[ncol_X_random_eff_unseen[3], M] unseen =
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[3] * M),
                 ncol_X_random_eff_unseen[3], M);
      for (i in 1:ncol_X_random_eff_unseen[3])
        unseen[i, ] = to_row_vector(normalize_sum_to_zero(to_vector(unseen[i, ])));
      mu = mu + (X_random_effect_3_unseen * unseen)';
    }
  }
  if (ncol_X_random_eff[4] > 0) {
    mu = mu + (X_random_effect_4 * random_effect_4[X_random_effect_which_4, ])';
    if (ncol_X_random_eff_unseen[4] > 0) {
      matrix[ncol_X_random_eff_unseen[4], M] unseen =
        to_matrix(rep_vector(std_normal_rng(), ncol_X_random_eff_unseen[4] * M),
                 ncol_X_random_eff_unseen[4], M);
      for (i in 1:ncol_X_random_eff_unseen[4])
        unseen[i, ] = to_row_vector(normalize_sum_to_zero(to_vector(unseen[i, ])));
      mu = mu + (X_random_effect_4_unseen * unseen)';
    }
  }

  matrix[M, N] mu_unconstrained = mu;

  for (i in 1:N) {
    mu[, i] = softmax(mu[, i]);
  }

  for (i in 1:N) {
    counts_uncorrected[i, ] = beta_binomial_rng(
      exposure[i],
      mu[, i] .* exp(precision[, i]),
      (1.0 - mu[, i]) .* exp(precision[, i])
    );
  }

  for (n in 1:N)
    generated_exposure[n] = fmax(sum(to_vector(counts_uncorrected[n])), 1);
  for (n in 1:N)
    counts[n] = to_row_vector(counts_uncorrected[n]) / generated_exposure[n] * exposure[n];
}
