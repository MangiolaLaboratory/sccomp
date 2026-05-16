

functions{

  #include common_functions.stan

  array[] int rep_each(array[] int x, int K) {
    int N = size(x);
    array[N * K] int y;
    int pos = 1;
    for (n in 1:N) {
      for (k in 1:K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }


    real abundance_variability_regression(
      row_vector variability,
      row_vector abundance,
      real prec_intercept_1,
      real prec_slope_1,
      real prec_slope_2,
      real prec_intercept_2,
      real prec_sd,
      int bimodal_mean_variability_association,
      real mix_p,
      // 0: variability prior depends on abundance via the regression slope.
      // 1: drop the abundance term; the prior reduces to an intercept-only
      //    Normal (or two-component mixture in the bimodal case), which is
      //    what callers want when the mean-variability association is to be
      //    disabled without abandoning the rest of the hierarchical structure.
      //    Same name and polarity as the user-facing flag, so no inversion
      //    happens at the call site.
      int exclude_mean_variability_association
    ){

      // Zeroing the slopes when the association is excluded preserves the rest
      // of the lpdf expression (including the bimodal mixture) unchanged, so we
      // do not duplicate the likelihood for the two cases.
      real eff_slope_1 = exclude_mean_variability_association == 1 ? 0 : prec_slope_1;
      real eff_slope_2 = exclude_mean_variability_association == 1 ? 0 : prec_slope_2;

      real lp = 0;
      // If mean-variability association is bimodal such as for single-cell RNA use mixed model
      if(bimodal_mean_variability_association == 1){
        for(m in 1:cols(variability))
        lp += log_mix(mix_p,
        normal_lpdf(variability[m] |
                    abundance[m] * eff_slope_1 + prec_intercept_1,
                    prec_sd),
        normal_lpdf(variability[m] |
                    abundance[m] * eff_slope_2 + prec_intercept_2,
                    prec_sd)
        );

        // If no bimodal
      } else {
        lp = normal_lpdf(variability |
                         abundance * eff_slope_1 + prec_intercept_1,
                         prec_sd);
      }

      return(lp);
    }

      real partial_sum_2_lpmf(
        // Parallel
        array[] int idx_y,
        int start,
        int end,

        // General
        int is_proportion,
        array[,] int y,
        array[,] real y_proportion,
        array[] int exposure,  // Sliced

        // Precision
        matrix Xa,                   // Sliced
        matrix alpha,

        // Fixed effects
        matrix X,                   // Sliced
        matrix beta,
        int M,

      // Random effects (up to 4 uniform slots)
      array[] int ncol_X_random_eff,
      matrix X_random_effect_1,   // Sliced
      matrix X_random_effect_2,   // Sliced
      matrix X_random_effect_3,   // Sliced
      matrix X_random_effect_4,   // Sliced
      matrix random_effect_1,
      matrix random_effect_2,
      matrix random_effect_3,
      matrix random_effect_4,

      // truncation
      array[,] int truncation_not_idx_minimal

      ){

        int N = end-start+1;

        // mu = fixed effects + each non-empty random-effect slot
        matrix[M, N] mu = (X[idx_y,] * beta)';
        if(ncol_X_random_eff[1]>0) mu = mu + (X_random_effect_1[idx_y,] * random_effect_1)';
        if(ncol_X_random_eff[2]>0) mu = mu + (X_random_effect_2[idx_y,] * random_effect_2)';
        if(ncol_X_random_eff[3]>0) mu = mu + (X_random_effect_3[idx_y,] * random_effect_3)';
        if(ncol_X_random_eff[4]>0) mu = mu + (X_random_effect_4[idx_y,] * random_effect_4)';

          for(n in 1:N)  mu[,n] = softmax(mu[,n]);

          // Precision
          matrix[M, N] precision = (Xa[idx_y,] * alpha)';

          // vectorisation
          vector[N*M] mu_array = to_vector(mu);
          vector[N*M] precision_array = to_vector(exp(precision));
          int W = count_filtered_indices(truncation_not_idx_minimal, idx_y);

          // truncation
          if(W == 0){

            // If input is proportions
            if(is_proportion)
              return beta_lupdf(
                to_array_1d(y_proportion[idx_y,]) |
                (mu_array .* precision_array),
                (1.0 - mu_array) .* precision_array
                ) ;

              // If input is counts
              else
                return beta_binomial_lupmf(
                  to_array_1d(y[idx_y,]) |
                  rep_each(exposure[idx_y], M),
                  (mu_array .* precision_array),
                  (1.0 - mu_array) .* precision_array
                ) ;

          }
          else{

            // If truncation is null for my chunk
            // Get non missing, invert the missing, this will be a big array
            array[N * M - W] int non_missing_indices =
            get_non_missing_indices(
              N,
              M,
              filter_missing_indices(truncation_not_idx_minimal, idx_y)
            );

            // If input is proportions
             if(is_proportion)
              return beta_lupdf(
              to_array_1d(y_proportion[idx_y,])[non_missing_indices] |
              (mu_array[non_missing_indices] .* precision_array[non_missing_indices]),
              (1.0 - mu_array[non_missing_indices]) .* precision_array[non_missing_indices]
              ) ;

             // If input is counts
             else
             return beta_binomial_lupmf(
              to_array_1d(y[idx_y,])[non_missing_indices] |
              rep_each(exposure[idx_y], M)[non_missing_indices],
              (mu_array[non_missing_indices] .* precision_array[non_missing_indices]),
              (1.0 - mu_array[non_missing_indices]) .* precision_array[non_missing_indices]
              ) ;

          }


        }

        /**
        * Counts the number of rows in missing_indices where the first column matches any value in idx_y.
        *
        * @param missing_indices  A two-dimensional integer array of size [TNIM, 2], containing {row, col} pairs of missing elements.
        * @param idx_y            An integer array containing the row indices to filter on.
        * @return                 An integer representing the number of rows in missing_indices where the first column matches any value in idx_y.
        *
        * @details
        * This function iterates over missing_indices and counts how many times the first column (row index) matches any value in idx_y.
        */
        int count_filtered_indices(array[,] int missing_indices, array[] int idx_y) {
          int num_missing = dims(missing_indices)[1];  // Number of rows in missing_indices
          int num_idx_y = num_elements(idx_y);         // Number of elements in idx_y
          int num_filtered = 0;                        // Initialize the count of filtered rows

          for (i in 1:num_missing) {
            for (j in 1:num_idx_y) {
              if (missing_indices[i, 1] == idx_y[j]) {
                num_filtered += 1;
                break;  // Exit the inner loop once a match is found
              }
            }
          }

          return num_filtered;
        }

        /**
        * Filters the missing_indices matrix to include only rows where the first column matches values in idx_y.
        *
        * @param missing_indices  A two-dimensional integer array of size [TNIM, 2], containing {row, col} pairs of missing elements.
        * @param idx_y            An integer array containing the row indices to filter on.
        * @return                 A two-dimensional integer array containing only the rows from missing_indices where the first column matches any value in idx_y.
        *
        * @details
        * This function uses the count_filtered_indices function to determine the size of the output array.
        * It then iterates over missing_indices to collect the matching rows.
        */
array[,] int filter_missing_indices(array[,] int missing_indices, array[] int idx_y) {
  int num_missing = dims(missing_indices)[1];  // Number of rows in missing_indices
  int num_idx_y = num_elements(idx_y);         // Number of elements in idx_y

  // Use the count_filtered_indices function to get the number of filtered rows
  int num_filtered = count_filtered_indices(missing_indices, idx_y);

  // Allocate the output array with the determined size
  array[num_filtered, 2] int missing_indices_filtered;

  // Fill the output array with matching rows
  int count = 0;
  for (i in 1:num_missing) {
    for (j in 1:num_idx_y) {
      if (missing_indices[i, 1] == idx_y[j]) {
        count += 1;

        // Adjust the row index in the filtered array
        missing_indices_filtered[count, 1] = j;  // Set to the relative position in idx_y
        missing_indices_filtered[count, 2] = missing_indices[i, 2];  // Keep the column index

        break;  // Exit the inner loop once a match is found
      }
    }
  }

  return missing_indices_filtered;
}
        /**
        * Compute indices of non-missing elements in a matrix when flattened in column-major order.
        * (Existing function; included here for completeness)
        */
        array[] int get_non_missing_indices(int n_rows, int n_cols, array[,] int missing_indices) {
          // Total number of elements in the matrix
          int N = n_rows * n_cols;

          // Number of missing elements
          int num_missing = dims(missing_indices)[1];  // Assuming missing_indices is [num_missing, 2]

          // Initialize a matrix to track missing data (0 = not missing, 1 = missing)
           array[n_rows, n_cols] int is_missing = rep_array(0, n_rows, n_cols);

          //
          // Mark the missing positions in the is_missing matrix
          for (i in 1:num_missing) {
            int row = missing_indices[i, 1];  // Row index of missing element
            int col = missing_indices[i, 2];  // Column index of missing element
            is_missing[row, col] = 1;         // Mark as missing
          }

          // Preallocate an array to hold the indices of non-missing elements
          int max_non_missing = N - num_missing;       // Maximum possible non-missing elements
          array[max_non_missing] int non_missing_indices;    // Array to store indices
          int count = 0;                               // Counter for non-missing elements

          // **Iterate over the matrix in row-major order**
          for (row in 1:n_rows) {
            for (col in 1:n_cols) {
              if (is_missing[row, col] == 0) {         // If the element is not missing
              count += 1;
              // **Use row-major indexing**
              non_missing_indices[count] = (row - 1) * n_cols + col;
              }
            }
          }

          // Return the array of non-missing indices (trimmed to the actual count)
          return non_missing_indices[1:count];

        }
}
data{
  int<lower=0, upper=1> is_proportion;
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> C;
  int<lower=1> A; // How many column in variability design
  int<lower=1> A_intercept_columns; // How many intercept column in varibility design
  int<lower=1> B_intercept_columns; // How many intercept column in varibility design
  int<lower=1> Ar; // Rows of unique variability design
  array[N] int exposure;
  array[N * !is_proportion,M] int<lower=0> y;
  array[N * is_proportion,M] real<lower=0, upper=1> y_proportion;
  matrix[N, C] X;
  matrix[Ar, A] XA; // The unique variability design
  matrix[N, A] Xa; // The variability design
  array[A] int<lower=1, upper=C> variability_to_composition_map;

  // Truncation
  int is_truncated;
  array[N,M] int truncation_up;
  array[N,M] int truncation_down;
  int<lower=1, upper=N*M> TNS; // truncation_not_size
  array[TNS] int<lower=1, upper=N*M> truncation_not_idx;

  int TNIM; // truncation_not_size
  array[TNIM,2] int<lower=1, upper=N*M> truncation_not_idx_minimal;

  int<lower=0, upper=1> is_vb;

  // Prior info
  array[2] real prior_prec_intercept;
  array[2] real prior_prec_slope;
  array[2] real prior_prec_sd;
  array[2] real prior_mean_intercept;
  array[2] real prior_mean_coefficients;

  // Exclude priors for testing purposes
  int<lower=0, upper=1> exclude_mean_variability_association;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;

  // Parallel chain
  int<lower=1> grainsize;

  // Does the design icludes intercept
  int <lower=0, upper=1> intercept_in_design;

  // ----------------------------------------------------------------------
  // Random effect blocks: up to 4 uniform "slots", one block per slot.
  // Each slot has its own n_factors (K) so that no padding is needed across
  // slots. A slot with ncol_X_random_eff[k] == 0 is unused (its arrays have
  // length 0 and no parameters get sampled).
  //
  //   slot 1 -> *_1     slot 2 -> *_2     slot 3 -> *_3     slot 4 -> *_4
  //
  // To go past 4 slots, paste another "slot 4" block in this file and bump
  // the array length below from 4 to 5. There is no other architectural
  // limit.
  // ----------------------------------------------------------------------
  int is_random_effect;

  array[4] int ncol_X_random_eff;
  matrix[N, ncol_X_random_eff[1]] X_random_effect_1;
  matrix[N, ncol_X_random_eff[2]] X_random_effect_2;
  matrix[N, ncol_X_random_eff[3]] X_random_effect_3;
  matrix[N, ncol_X_random_eff[4]] X_random_effect_4;

  // Covariance setup (per slot)
  array[4] int n_groups;
  array[4] int how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance_1;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
  array[how_many_factors_in_random_design[3], n_groups[3]] int group_factor_indexes_for_covariance_3;
  array[how_many_factors_in_random_design[4], n_groups[4]] int group_factor_indexes_for_covariance_4;

  // LOO
  int<lower=0, upper=1> enable_loo;


}
transformed data{
  // Floors of 1 so the per-slot random_effect matrices can be declared
  // without a 0-row size on Windows (Stan generates degenerate code there).
  // The matrices remain logically empty when the slot is unused.
  int ncol_X_random_eff_safe_1 = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_safe_2 = max(ncol_X_random_eff[2], 1);
  int ncol_X_random_eff_safe_3 = max(ncol_X_random_eff[3], 1);
  int ncol_X_random_eff_safe_4 = max(ncol_X_random_eff[4], 1);

  // For parallelisation
  array[N] int array_N;
  for(n in 1:N) array_N[n] = n;

  // Data vectorised
  // y_array =  to_array_1d(y);
  // exposure_array = rep_each(exposure, M);
}
parameters{
  // Use the new sum_to_zero_vector type instead of QR decomposition
  array[C] sum_to_zero_vector[M] beta_raw; // Each row is a sum_to_zero_vector of length M
  matrix[A, M] alpha; // Variability
  // Mean-variability intercepts: length 1 if unimodal, length 2 (strictly increasing) if bimodal.
  array[A] ordered[1 + bimodal_mean_variability_association] prec_intercept;
  // Mean-variability slopes
  array[A] real prec_slope_1; // s1, always present
  array[A * bimodal_mean_variability_association] real prec_slope_2; // s2, only for bimodal
  // Log-scale residual SD avoids a hard boundary at 0 and reduces funnel neck pathologies.
  array[A] real log_prec_sd;
  real<lower=0, upper=1> mix_p;

  // ----------------------------------------------------------------------
  // Random effect parameters - 4 uniform slots.
  //
  // For each slot k:
  //   * random_effect_raw_k       : sum_to_zero_vector[M] per design column
  //   * random_effect_sigma_raw_k : per-category raw SD vector (length n_factors[k])
  //   * sigma_correlation_factor_k: per-category Cholesky of correlation matrix
  //                                  (n_factors[k] x n_factors[k]; 1x1 = no LKJ work)
  //
  // Hyperprior scalars sigma_mu / sigma_sigma are shared in length-4 arrays.
  // ----------------------------------------------------------------------

  // Slot 1
  array[ncol_X_random_eff[1] * (ncol_X_random_eff[1]>0)] sum_to_zero_vector[M] random_effect_raw_1;
  array[M * (ncol_X_random_eff[1]>0)] vector[how_many_factors_in_random_design[1]]  random_effect_sigma_raw_1;
  array[M * (ncol_X_random_eff[1]>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (ncol_X_random_eff[1]>0)] sigma_correlation_factor_1;

  // Slot 2
  array[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0)] sum_to_zero_vector[M] random_effect_raw_2;
  array[M * (ncol_X_random_eff[2]>0)] vector[how_many_factors_in_random_design[2]]  random_effect_sigma_raw_2;
  array[M * (ncol_X_random_eff[2]>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (ncol_X_random_eff[2]>0)] sigma_correlation_factor_2;

  // Slot 3
  array[ncol_X_random_eff[3] * (ncol_X_random_eff[3]>0)] sum_to_zero_vector[M] random_effect_raw_3;
  array[M * (ncol_X_random_eff[3]>0)] vector[how_many_factors_in_random_design[3]]  random_effect_sigma_raw_3;
  array[M * (ncol_X_random_eff[3]>0)] cholesky_factor_corr[how_many_factors_in_random_design[3] * (ncol_X_random_eff[3]>0)] sigma_correlation_factor_3;

  // Slot 4
  array[ncol_X_random_eff[4] * (ncol_X_random_eff[4]>0)] sum_to_zero_vector[M] random_effect_raw_4;
  array[M * (ncol_X_random_eff[4]>0)] vector[how_many_factors_in_random_design[4]]  random_effect_sigma_raw_4;
  array[M * (ncol_X_random_eff[4]>0)] cholesky_factor_corr[how_many_factors_in_random_design[4] * (ncol_X_random_eff[4]>0)] sigma_correlation_factor_4;

  // Shared hyperprior scalars (one mu, one sigma per slot)
  array[4 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[4 * (is_random_effect>0)] real random_effect_sigma_sigma;

  // For models with a single group (kept from the original design)
  array[is_random_effect>0] real zero_random_effect;


}
transformed parameters{

  array[A] real prec_intercept_1;
  array[A * bimodal_mean_variability_association] real prec_intercept_2;
  array[A] real<lower=0> prec_sd;  // residual scale per effect for mean-variability association
  for (a in 1:A) {
    prec_intercept_1[a] = prec_intercept[a][1];
    prec_sd[a] = exp(log_prec_sd[a]);
    if (bimodal_mean_variability_association == 1)
      prec_intercept_2[a] = prec_intercept[a][2];
  }

  // Initialisation
  matrix[C,M] beta;
  matrix[M, N] precision = (Xa * alpha)';

  // Convert sum_to_zero_vector to regular matrix
  for(c in 1:C) {
    beta[c,] = to_row_vector(beta_raw[c]);
  }

  real mix_p_scalar = bimodal_mean_variability_association == 1 ? mix_p : 0.5;

  // ----------------------------------------------------------------------
  // Random effect contributions: one matrix per slot.
  //
  // Per slot, build_re_block() does the per-category SD non-centered build
  // and the LKJ-Cholesky x SD covariance construction. Slots with
  // ncol_X_random_eff[k] == 0 keep an unused 1xM placeholder (Windows guard)
  // and never get touched.
  // ----------------------------------------------------------------------
  matrix[ncol_X_random_eff_safe_1 * (is_random_effect>0), M] random_effect_1;
  matrix[ncol_X_random_eff_safe_2 * (is_random_effect>0), M] random_effect_2;
  matrix[ncol_X_random_eff_safe_3 * (is_random_effect>0), M] random_effect_3;
  matrix[ncol_X_random_eff_safe_4 * (is_random_effect>0), M] random_effect_4;

  if (ncol_X_random_eff[1] > 0) {
    array[ncol_X_random_eff[1]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[1]) raw_vec[i] = to_vector(random_effect_raw_1[i]);
    random_effect_1 = build_re_block(
      M, n_groups[1], how_many_factors_in_random_design[1], ncol_X_random_eff[1],
      group_factor_indexes_for_covariance_1, raw_vec,
      random_effect_sigma_mu[1], random_effect_sigma_sigma[1],
      random_effect_sigma_raw_1, sigma_correlation_factor_1
    );
  }

  if (ncol_X_random_eff[2] > 0) {
    array[ncol_X_random_eff[2]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[2]) raw_vec[i] = to_vector(random_effect_raw_2[i]);
    random_effect_2 = build_re_block(
      M, n_groups[2], how_many_factors_in_random_design[2], ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2, raw_vec,
      random_effect_sigma_mu[2], random_effect_sigma_sigma[2],
      random_effect_sigma_raw_2, sigma_correlation_factor_2
    );
  }

  if (ncol_X_random_eff[3] > 0) {
    array[ncol_X_random_eff[3]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[3]) raw_vec[i] = to_vector(random_effect_raw_3[i]);
    random_effect_3 = build_re_block(
      M, n_groups[3], how_many_factors_in_random_design[3], ncol_X_random_eff[3],
      group_factor_indexes_for_covariance_3, raw_vec,
      random_effect_sigma_mu[3], random_effect_sigma_sigma[3],
      random_effect_sigma_raw_3, sigma_correlation_factor_3
    );
  }

  if (ncol_X_random_eff[4] > 0) {
    array[ncol_X_random_eff[4]] vector[M] raw_vec;
    for (i in 1:ncol_X_random_eff[4]) raw_vec[i] = to_vector(random_effect_raw_4[i]);
    random_effect_4 = build_re_block(
      M, n_groups[4], how_many_factors_in_random_design[4], ncol_X_random_eff[4],
      group_factor_indexes_for_covariance_4, raw_vec,
      random_effect_sigma_mu[4], random_effect_sigma_sigma[4],
      random_effect_sigma_raw_4, sigma_correlation_factor_4
    );
  }
}
model{
  // Fit main distribution
  if(use_data == 1){

   target += reduce_sum(
      partial_sum_2_lupmf,
      array_N,
      grainsize,

      // General
      is_proportion,
      y,
      y_proportion,
      exposure,

      // Precision
      Xa,
      alpha,

      // Fixed effects
      X,
      beta,
      M,

      // Random effects (4 uniform slots)
      ncol_X_random_eff,
      X_random_effect_1,
      X_random_effect_2,
      X_random_effect_3,
      X_random_effect_4,
      random_effect_1,
      random_effect_2,
      random_effect_3,
      random_effect_4,

      //truncation
      truncation_not_idx_minimal

      );


  }

  // Variability prior via the mean-variability regression.
  // When `exclude_mean_variability_association = 1` the regression drops the
  // abundance term but keeps the hierarchical hyperpriors on
  // `prec_intercept_*`, `prec_slope_*` and `prec_sd` identical across both
  // modes — only the slope contribution to alpha's prior gets zeroed. This
  // avoids a parallel ad-hoc prior block.
  for(a in 1:A){
    target += abundance_variability_regression(
      alpha[a],
      beta[variability_to_composition_map[a]],
      prec_intercept_1[a],
      prec_slope_1[a],
      bimodal_mean_variability_association == 1 ? prec_slope_2[a] : 0,
      bimodal_mean_variability_association == 1 ? prec_intercept_2[a] : 0,
      prec_sd[a],
      bimodal_mean_variability_association,
      mix_p_scalar,
      exclude_mean_variability_association
    );
  }

  // Hyper priors: i1/s1 (and prec_sd) shared; bimodal adds i2/s2 and mix_p shape
  if(bimodal_mean_variability_association == 1)
    mix_p ~ beta(1, 5);
  else
    mix_p ~ beta(1, 1);

  for(a in 1:A){
    // If design has intercept, first column gets intercept-centred prior, others are centred at 0.
    if(intercept_in_design == 1 && a == 1){
      prec_intercept[a][1] ~ student_t(3, prior_prec_intercept[1], prior_prec_intercept[2]);
      if(bimodal_mean_variability_association == 1)
        prec_intercept[a][2] ~ student_t(3, prior_prec_intercept[1], prior_prec_intercept[2]);
    } else {
      prec_intercept[a][1] ~ student_t(3, 0, prior_prec_intercept[2]);
      if(bimodal_mean_variability_association == 1)
        prec_intercept[a][2] ~ student_t(3, 0, prior_prec_intercept[2]);
    }
    prec_slope_1[a] ~ student_t(3, prior_prec_slope[1], prior_prec_slope[2]); // s1
    if(bimodal_mean_variability_association == 1){
      prec_slope_2[a] ~ student_t(3, prior_prec_slope[1], prior_prec_slope[2]); // s2
    }
  }
  for(a in 1:A) log_prec_sd[a] ~ normal(prior_prec_sd[1], prior_prec_sd[2]);

  // // Priors abundance - use correct scale for sum_to_zero_vector
  for(c in 1:B_intercept_columns) beta_raw[c] ~ normal ( prior_mean_intercept[1], prior_mean_intercept[2] * inv(sqrt(1 - inv(M))) );
  if(C>B_intercept_columns) for(c in (B_intercept_columns+1):C) beta_raw[c] ~ normal ( prior_mean_coefficients[1], prior_mean_coefficients[2] * inv(sqrt(1 - inv(M))) );
  
  // ----------------------------------------------------------------------
  // Random effect priors. Per-slot block, same pattern; the `if` keeps an
  // unused slot's empty-length arrays out of any prior call. The shared
  // hyperprior arrays are sampled once, outside the per-slot loop.
  // ----------------------------------------------------------------------
  if (is_random_effect > 0) {
    random_effect_sigma_mu ~ std_normal();
    random_effect_sigma_sigma ~ std_normal();
    zero_random_effect ~ std_normal();
  }

  if (ncol_X_random_eff[1] > 0) {
    for (m in 1:M) random_effect_raw_1[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for (m in 1:M) random_effect_sigma_raw_1[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_1[m] ~ lkj_corr_cholesky(2);
  }

  if (ncol_X_random_eff[2] > 0) {
    for (m in 1:M) random_effect_raw_2[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for (m in 1:M) random_effect_sigma_raw_2[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_2[m] ~ lkj_corr_cholesky(2);
  }

  if (ncol_X_random_eff[3] > 0) {
    for (m in 1:M) random_effect_raw_3[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for (m in 1:M) random_effect_sigma_raw_3[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_3[m] ~ lkj_corr_cholesky(2);
  }

  if (ncol_X_random_eff[4] > 0) {
    for (m in 1:M) random_effect_raw_4[,m] ~ normal(0, inv(sqrt(1 - inv(M))));
    for (m in 1:M) random_effect_sigma_raw_4[m] ~ std_normal();
    for (m in 1:M) sigma_correlation_factor_4[m] ~ lkj_corr_cholesky(2);
  }
  }
generated quantities {
  // LOO
  vector[TNS] log_lik = rep_vector(0, TNS);


  // LOO
  if(enable_loo==1){

    matrix[M, N] mu;
    vector[N*M] mu_array;
    vector[N*M] precision_array;

    mu = (X * beta)';

    // Each non-empty random-effect slot contributes additively
    if(ncol_X_random_eff[1]>0) mu = mu + (X_random_effect_1 * random_effect_1)';
    if(ncol_X_random_eff[2]>0) mu = mu + (X_random_effect_2 * random_effect_2)';
    if(ncol_X_random_eff[3]>0) mu = mu + (X_random_effect_3 * random_effect_3)';
    if(ncol_X_random_eff[4]>0) mu = mu + (X_random_effect_4 * random_effect_4)';


    // Calculate proportions
    for(n in 1:N)  mu[,n] = softmax(mu[,n]);

    // Convert the matrix m to a column vector in column-major order.
    mu_array = to_vector(mu);
    precision_array = to_vector(exp(precision));

    if(is_proportion)
          for (n in 1:TNS) {
      log_lik[n] = beta_lpdf(
        to_array_1d(y_proportion)[truncation_not_idx[n]] |
        (mu_array[truncation_not_idx[n]] .* precision_array[truncation_not_idx[n]]),
        ((1.0 - mu_array[truncation_not_idx[n]]) .* precision_array[truncation_not_idx[n]])
        ) ;
      }
    else
      for (n in 1:TNS) {
       log_lik[n] = beta_binomial_lpmf(
        to_array_1d(y)[truncation_not_idx[n]] |
        rep_each(exposure, M)[truncation_not_idx[n]],
        (mu_array[truncation_not_idx[n]] .* precision_array[truncation_not_idx[n]]),
        ((1.0 - mu_array[truncation_not_idx[n]]) .* precision_array[truncation_not_idx[n]])
        ) ;
    }

  }
}

