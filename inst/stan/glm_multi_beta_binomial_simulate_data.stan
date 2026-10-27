functions{
    #include common_functions.stan

    // Linear predictor contributions for random-effect design blocks; mu is M × N (types × samples).
    matrix add_seen_random_effect_contribution(matrix mu_mn,
                                               matrix x_n_times_k,
                                               matrix re_all_rows_times_m,
                                               array[] int re_row_index) {
      return mu_mn + (x_n_times_k * re_all_rows_times_m[re_row_index, ])';
    }

    matrix add_unseen_random_effect_contribution_rng(matrix mu_mn,
                                                     matrix x_unseen_n_times_u,
                                                     int ncol_u,
                                                     int M_param) {
      matrix[ncol_u, M_param] unseen =
          to_matrix(rep_vector(std_normal_rng(), ncol_u * M_param),
                    ncol_u, M_param);
      for (i in 1:ncol_u)
        unseen[i, ] = to_row_vector(normalize_sum_to_zero(to_vector(unseen[i, ])));
      return mu_mn + (x_unseen_n_times_u * unseen)';
    }

}
data {
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> exclude_mean_variability_association;
  int N;
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

  // Composition column each variability column is associated with
  array[A] int<lower=1, upper=C> variability_to_composition_map;

  // Optional overrides. Rows of beta_provided follow X_which.
  int<lower=0, upper=1> user_provided_beta;
  matrix[length_X_which, M] beta_provided;
  array[A] int<lower=0, upper=1> user_provided_slope;
  array[A] real slope_provided;
  array[A] int<lower=0, upper=1> user_provided_intercept;
  array[A] real intercept_provided;

  real<lower=0> variability_multiplier;
  int<lower=0, upper=1> variability_scatter;

  int is_random_effect;

  // Six uniform random-effect slots (same layout as glm_multi_beta_binomial.stan)
  array[6] int ncol_X_random_eff;
  array[6] int ncol_X_random_eff_new;

  matrix[N, ncol_X_random_eff_new[1]] X_random_effect_1;
  matrix[N, ncol_X_random_eff_new[2]] X_random_effect_2;
  matrix[N, ncol_X_random_eff_new[3]] X_random_effect_3;
  matrix[N, ncol_X_random_eff_new[4]] X_random_effect_4;
  matrix[N, ncol_X_random_eff_new[5]] X_random_effect_5;
  matrix[N, ncol_X_random_eff_new[6]] X_random_effect_6;

  array[6] int n_groups;
  array[6] int how_many_factors_in_random_design;

  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance_1;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
  array[how_many_factors_in_random_design[3], n_groups[3]] int group_factor_indexes_for_covariance_3;
  array[how_many_factors_in_random_design[4], n_groups[4]] int group_factor_indexes_for_covariance_4;
  array[how_many_factors_in_random_design[5], n_groups[5]] int group_factor_indexes_for_covariance_5;
  array[how_many_factors_in_random_design[6], n_groups[6]] int group_factor_indexes_for_covariance_6;

  array[6] int length_X_random_effect_which;
  array[length_X_random_effect_which[1]] int X_random_effect_which_1;
  array[length_X_random_effect_which[2]] int X_random_effect_which_2;
  array[length_X_random_effect_which[3]] int X_random_effect_which_3;
  array[length_X_random_effect_which[4]] int X_random_effect_which_4;
  array[length_X_random_effect_which[5]] int X_random_effect_which_5;
  array[length_X_random_effect_which[6]] int X_random_effect_which_6;

  array[6] int ncol_X_random_eff_unseen;
  matrix[N, ncol_X_random_eff_unseen[1]] X_random_effect_1_unseen;
  matrix[N, ncol_X_random_eff_unseen[2]] X_random_effect_2_unseen;
  matrix[N, ncol_X_random_eff_unseen[3]] X_random_effect_3_unseen;
  matrix[N, ncol_X_random_eff_unseen[4]] X_random_effect_4_unseen;
  matrix[N, ncol_X_random_eff_unseen[5]] X_random_effect_5_unseen;
  matrix[N, ncol_X_random_eff_unseen[6]] X_random_effect_6_unseen;
}
transformed data {
  int ncol_X_random_eff_safe_1 = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_safe_2 = max(ncol_X_random_eff[2], 1);
  int ncol_X_random_eff_safe_3 = max(ncol_X_random_eff[3], 1);
  int ncol_X_random_eff_safe_4 = max(ncol_X_random_eff[4], 1);
  int ncol_X_random_eff_safe_5 = max(ncol_X_random_eff[5], 1);
  int ncol_X_random_eff_safe_6 = max(ncol_X_random_eff[6], 1);
}
parameters {
  // Must match the parameter block of glm_multi_beta_binomial.stan, with
  // unconstrained vectors so posterior CSVs (rounded sig_figs) still validate.
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

  array[ncol_X_random_eff[5] * (ncol_X_random_eff[5]>0)] vector[M] random_effect_raw_5;
  array[M * (ncol_X_random_eff[5]>0)] vector[how_many_factors_in_random_design[5]] random_effect_sigma_raw_5;
  array[M * (ncol_X_random_eff[5]>0)] cholesky_factor_corr[how_many_factors_in_random_design[5] * (ncol_X_random_eff[5]>0)] sigma_correlation_factor_5;

  array[ncol_X_random_eff[6] * (ncol_X_random_eff[6]>0)] vector[M] random_effect_raw_6;
  array[M * (ncol_X_random_eff[6]>0)] vector[how_many_factors_in_random_design[6]] random_effect_sigma_raw_6;
  array[M * (ncol_X_random_eff[6]>0)] cholesky_factor_corr[how_many_factors_in_random_design[6] * (ncol_X_random_eff[6]>0)] sigma_correlation_factor_6;

  array[6 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[6 * (is_random_effect>0)] real random_effect_sigma_sigma;
  array[is_random_effect>0] real zero_random_effect;
}
transformed parameters {
  matrix[C,M] beta;
  for (c in 1:C) {
    beta[c,] = to_row_vector(normalize_sum_to_zero(beta_raw[c]));
  }

  matrix[ncol_X_random_eff_safe_1 * (is_random_effect>0), M] random_effect_1;
  matrix[ncol_X_random_eff_safe_2 * (is_random_effect>0), M] random_effect_2;
  matrix[ncol_X_random_eff_safe_3 * (is_random_effect>0), M] random_effect_3;
  matrix[ncol_X_random_eff_safe_4 * (is_random_effect>0), M] random_effect_4;
  matrix[ncol_X_random_eff_safe_5 * (is_random_effect>0), M] random_effect_5;
  matrix[ncol_X_random_eff_safe_6 * (is_random_effect>0), M] random_effect_6;

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
  if (ncol_X_random_eff[5] > 0) {
    array[ncol_X_random_eff[5]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[5]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_5[i]);
    random_effect_5 = build_re_block(
      M, n_groups[5], how_many_factors_in_random_design[5], ncol_X_random_eff[5],
      group_factor_indexes_for_covariance_5, raw_vec,
      random_effect_sigma_mu[5], random_effect_sigma_sigma[5],
      random_effect_sigma_raw_5, sigma_correlation_factor_5
    );
  }
  if (ncol_X_random_eff[6] > 0) {
    array[ncol_X_random_eff[6]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[6]) raw_vec[i] = normalize_sum_to_zero(random_effect_raw_6[i]);
    random_effect_6 = build_re_block(
      M, n_groups[6], how_many_factors_in_random_design[6], ncol_X_random_eff[6],
      group_factor_indexes_for_covariance_6, raw_vec,
      random_effect_sigma_mu[6], random_effect_sigma_sigma[6],
      random_effect_sigma_raw_6, sigma_correlation_factor_6
    );
  }
}
model {
}
generated quantities {
  array[N, M] int counts_uncorrected;
  matrix[N, M] counts;
  array[N] real generated_exposure;

  matrix[C, M] beta_simulated = beta;
  if (user_provided_beta == 1)
    for (i in 1:length_X_which)
      beta_simulated[X_which[i], ] = beta_provided[i, ];

  // Variability follows the abundance it is associated with, through the
  // first (main) component of the mean-variability association. The posterior
  // alpha is not used, which is what lets user-provided coefficients move it.
  // With scatter, each cell group deviates from the trend by the residual SD
  // of the association, as fitted cell groups do.
  matrix[A, M] alpha_simulated;
  for (a in 1:A) {
    real slope =
      user_provided_slope[a] == 1 ? slope_provided[a] :
      (exclude_mean_variability_association == 1 ? 0 : prec_slope_1[a]);
    real intercept =
      user_provided_intercept[a] == 1 ? intercept_provided[a] : prec_intercept[a][1];
    alpha_simulated[a, ] = beta_simulated[variability_to_composition_map[a], ] * slope + intercept;
    if (variability_scatter == 1)
      for (m in 1:M)
        alpha_simulated[a, m] = normal_rng(alpha_simulated[a, m], exp(log_prec_sd[a]));
  }

  matrix[M, N] mu = (X * beta_simulated[X_which, ])';
  matrix[M, N] precision = (Xa * alpha_simulated[XA_which, ])' - log(variability_multiplier);

  if (ncol_X_random_eff[1] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_1, random_effect_1,
                                             X_random_effect_which_1);
    if (ncol_X_random_eff_unseen[1] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_1_unseen,
                                                    ncol_X_random_eff_unseen[1], M);
  }
  if (ncol_X_random_eff[2] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_2, random_effect_2,
                                             X_random_effect_which_2);
    if (ncol_X_random_eff_unseen[2] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_2_unseen,
                                                    ncol_X_random_eff_unseen[2], M);
  }
  if (ncol_X_random_eff[3] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_3, random_effect_3,
                                             X_random_effect_which_3);
    if (ncol_X_random_eff_unseen[3] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_3_unseen,
                                                    ncol_X_random_eff_unseen[3], M);
  }
  if (ncol_X_random_eff[4] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_4, random_effect_4,
                                             X_random_effect_which_4);
    if (ncol_X_random_eff_unseen[4] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_4_unseen,
                                                    ncol_X_random_eff_unseen[4], M);
  }
  if (ncol_X_random_eff[5] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_5, random_effect_5,
                                             X_random_effect_which_5);
    if (ncol_X_random_eff_unseen[5] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_5_unseen,
                                                    ncol_X_random_eff_unseen[5], M);
  }
  if (ncol_X_random_eff[6] > 0) {
    mu = add_seen_random_effect_contribution(mu, X_random_effect_6, random_effect_6,
                                             X_random_effect_which_6);
    if (ncol_X_random_eff_unseen[6] > 0)
      mu = add_unseen_random_effect_contribution_rng(mu, X_random_effect_6_unseen,
                                                    ncol_X_random_eff_unseen[6], M);
  }

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
