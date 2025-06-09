

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
  
    
    real abundance_variability_regression(row_vector variability, row_vector abundance, array[] real prec_coeff, real prec_sd, int bimodal_mean_variability_association, real mix_p){
      
      real lp = 0;
      // If mean-variability association is bimodal such as for single-cell RNA use mixed model
      if(bimodal_mean_variability_association == 1){
        for(m in 1:cols(variability))
        lp += log_mix(mix_p,
        normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + prec_coeff[1], prec_sd ),
        normal_lpdf(variability[m] | abundance[m] * prec_coeff[2] + 1, prec_sd)  // -0.73074903 is what we observe in single-cell dataset Therefore it is safe to fix it for this mixture model as it just want to capture few possible outlier in the association
        );
        
        // If no bimodal
      } else {
        lp =  normal_lpdf(variability | abundance * prec_coeff[2] + prec_coeff[1], prec_sd );
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
        
        // Random effects
        array[] int ncol_X_random_eff,
        matrix X_random_effect,   // Sliced
        matrix X_random_effect_2,  // Sliced
        matrix random_effect,
        matrix random_effect_2,
        
        // truncation
        array[,] int truncation_not_idx_minimal
        
        ){
          
          int N = end-start+1;
          
          // mu
          matrix[M, N] mu = (X[idx_y,] * beta)';
          if(ncol_X_random_eff[1]> 0)
          mu = mu + append_row((X_random_effect[idx_y,] * random_effect)', rep_row_vector(0, N));
          
          if(ncol_X_random_eff[2]>0 )
          mu = mu + append_row((X_random_effect_2[idx_y,] * random_effect_2)', rep_row_vector(0, N));
          
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
  int<lower=1> A; // How many column in variability design\
  int<lower=1> A_intercept_columns; // How many intercept column in varibility design
  int<lower=1> B_intercept_columns; // How many intercept column in varibility design
  int<lower=1> Ar; // Rows of unique variability design
  array[N] int exposure;
  array[N * !is_proportion,M] int<lower=0> y;
  array[N * is_proportion,M] real<lower=0, upper=1> y_proportion;
  matrix[N, C] X;
  matrix[Ar, A] XA; // The unique variability design
  matrix[N, A] Xa; // The variability design
  
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
  int<lower=0, upper=1> exclude_priors;
  int<lower=0, upper=1> bimodal_mean_variability_association;
  int<lower=0, upper=1> use_data;
  
  // Parallel chain
  int<lower=1> grainsize;
  
  // Does the design icludes intercept
  int <lower=0, upper=1> intercept_in_design;
  
  // Random intercept
  
  // Is the parameters in random effect matrix, minus ther sub to zero parameters, for example if I have four groups, this will be 3
  int is_random_effect;
  
  // Is the parameters in random effect matrix
  array[2] int ncol_X_random_eff;
  matrix[N, ncol_X_random_eff[1]] X_random_effect;
  matrix[N, ncol_X_random_eff[2]] X_random_effect_2;
  
  // Covariance setup
  array[2] int n_groups;
  array[2] int how_many_factors_in_random_design;
  array[how_many_factors_in_random_design[1], n_groups[1]] int group_factor_indexes_for_covariance;
  array[how_many_factors_in_random_design[2], n_groups[2]] int group_factor_indexes_for_covariance_2;
  
  // LOO
  int<lower=0, upper=1> enable_loo;
  
  
}
transformed data{
  vector[2*M] Q_r = Q_sum_to_zero_QR(M);
  real x_raw_sigma = inv_sqrt(1 - inv(M));
  matrix[N, C] Q_ast;
  matrix[C, C] R_ast;
  matrix[C, C] R_ast_inverse;
  // array[N*M] real y_array;
  array[N*M] int truncation_down_array;
  // array[N*M] int exposure_array;
  // EXCEPTION MADE FOR WINDOWS GENERATE QUANTITIES IF RANDOM EFFECT DO NOT EXIST
  int ncol_X_random_eff_WINDOWS_BUG_FIX = max(ncol_X_random_eff[1], 1);
  int ncol_X_random_eff_WINDOWS_BUG_FIX_2 = max(ncol_X_random_eff[2], 1);
  
  // thin and scale the QR decomposition
  Q_ast = qr_thin_Q(X) * sqrt(N - 1);
  R_ast_inverse = inverse(qr_thin_R(X) / sqrt(N - 1));
  
  // If I get crazy diagonal matrix omit it
  if(is_random_effect>0) { 
    if(max(R_ast_inverse)>1000 )
    print("sccomp says: The QR deconposition resulted in extreme values, probably for the correlation structure of your design matrix. Omitting QR decomposition.");
    Q_ast = X;
    R_ast_inverse = diag_matrix(rep_vector(1.0, C));
  }
  
  // For parallelisation
  array[N] int array_N;
  for(n in 1:N) array_N[n] = n;
    
  // Data vectorised
  // y_array =  to_array_1d(y);
  // exposure_array = rep_each(exposure, M);
}
parameters{
  matrix[C, M-1] beta_raw_raw; // matrix with C rows and number of cells (-1) columns
  matrix[A, M] alpha; // Variability
  // To exclude
  array[2] real prec_coeff;
  real<lower=0> prec_sd;
  real<lower=0, upper=1> mix_p;
  
  // Random intercept // matrix with ncol_X_random_effs rows and number of cells (-1) columns
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_effect_raw;
  matrix[ncol_X_random_eff[2] * (ncol_X_random_eff[2]>0), M-1] random_effect_raw_2;
  
  // sd of random intercept
  array[2 * (is_random_effect>0)] real random_effect_sigma_mu;
  array[2 * (is_random_effect>0)] real random_effect_sigma_sigma;
  
  // Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[1]]  random_effect_sigma_raw;
  array[M-1 * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[1] * (is_random_effect>0)] sigma_correlation_factor;
  
  // Covariance
  array[M-1 * (is_random_effect>0)] vector[how_many_factors_in_random_design[2]]  random_effect_sigma_raw_2;
  array[M-1 * (is_random_effect>0)] cholesky_factor_corr[how_many_factors_in_random_design[2] * (is_random_effect>0)] sigma_correlation_factor_2;
  
  // If I have just one group
  array[is_random_effect>0] real zero_random_effect;
  
  
}
transformed parameters{
  
  // Initialisation
  matrix[C,M] beta_raw;
  matrix[M, N] precision = (Xa * alpha)';
  matrix[C,M] beta;
  
  // // // locations distribution
  // matrix[M, N] mu;
  // // vectorisation
  // vector[N*M] mu_array;
  // vector[N*M] precision_array;
  
  for(c in 1:C)	beta_raw[c,] =  sum_to_zero_QR(beta_raw_raw[c,], Q_r);
  beta = R_ast_inverse * beta_raw; // coefficients on x
  
  // // JUST FOR SANITY CHECK
  //  mu = (Q_ast * beta_raw)';
  
  // Non centered parameterisation SD of random effects
  array[M-1 * (ncol_X_random_eff[1]> 0)] vector[how_many_factors_in_random_design[1]] random_effect_sigma;
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M-1)) random_effect_sigma[m] = random_effect_sigma_mu[1] + random_effect_sigma_sigma[1] * random_effect_sigma_raw[m];
  if(ncol_X_random_eff[1]> 0) for(m in 1:(M-1)) random_effect_sigma[m] = exp(random_effect_sigma[m]/3.0);
  
  // Non centered parameterisation SD of random effects 2
  array[M-1 * (ncol_X_random_eff[2]> 0)] vector[how_many_factors_in_random_design[2]] random_effect_sigma_2;
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M-1)) random_effect_sigma_2[m] = random_effect_sigma_mu[2] + random_effect_sigma_sigma[2] * random_effect_sigma_raw_2[m];
  if(ncol_X_random_eff[2]> 0) for(m in 1:(M-1)) random_effect_sigma_2[m] = exp(random_effect_sigma_2[m]/3.0);
    
  
  matrix[ncol_X_random_eff[1] * (is_random_effect>0), M-1] random_effect;
  matrix[ncol_X_random_eff[2] * (is_random_effect>0), M-1] random_effect_2;
  
  // random intercept
  if(ncol_X_random_eff[1]> 0){
    
    // Covariate setup
    random_effect =
    get_random_effect_matrix(
      M,
      n_groups[1],
      how_many_factors_in_random_design[1],
      is_random_effect,
      ncol_X_random_eff[1],
      group_factor_indexes_for_covariance,
      
      random_effect_raw,
      random_effect_sigma,
      sigma_correlation_factor
      );
      
  }
  
  // random intercept
  if(ncol_X_random_eff[2]>0 ){

    // Covariate setup
    random_effect_2 =
    get_random_effect_matrix(
      M,
      n_groups[2],
      how_many_factors_in_random_design[2],
      is_random_effect,
      ncol_X_random_eff[2],
      group_factor_indexes_for_covariance_2,

      random_effect_raw_2,
      random_effect_sigma_2,
      sigma_correlation_factor_2
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
      Q_ast,                   
      beta_raw, 
      M, 
      
      // Random effects
      ncol_X_random_eff,
      X_random_effect,   
      X_random_effect_2, 
      random_effect,
      random_effect_2,
      
      //truncation
      truncation_not_idx_minimal
      
      );
      
      // print("2---", reduce_sum(
        //   partial_sum_lupmf,
        //   y_array[truncation_not_idx],
        //   grainsize,
        //   exposure_array[truncation_not_idx],
        //   mu_array[truncation_not_idx],
        //   precision_array[truncation_not_idx]
        //   ));
        
        
  }
  
  // Priors
  if(exclude_priors == 0){
    
    // If interceopt in design or I have complex variability design
    // This would include the models 
    // composition ~ 1 + ...; composition ~ 0 + ...; 
    // variability ~ 1
    if(A == 1){
      target += abundance_variability_regression(
        alpha[1],
        beta[1], // average_by_col(beta[1:B_intercept_columns,]),
        prec_coeff,
        prec_sd,
        bimodal_mean_variability_association,
        mix_p
        );
    }
    else {
      // Loop across the intercept columns in case of a intercept-less design (covariate are intercepts)
      for(a in 1:A_intercept_columns)
      target += abundance_variability_regression(
        alpha[a],
        beta[a],
        prec_coeff,
        prec_sd,
        bimodal_mean_variability_association,
        mix_p
        );
        
        // Variability effect if the formula is more complex
        if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) alpha[a] ~ normal(beta[a] * prec_coeff[2], 2 );
    }
    
  }
  
  // If I don't have priors for overdispersion
  else{
    // Priors variability
    if(intercept_in_design || A > 1){
      for(a in 1:A_intercept_columns) alpha[a]  ~ normal( prec_coeff[1], prec_sd );
      if(A>A_intercept_columns) for(a in (A_intercept_columns+1):A) to_vector(alpha[a]) ~ normal ( 0, 2 );
    }
    // if ~ 0 + covariuate
    else {
      alpha[1]  ~ normal( prec_coeff[1], prec_sd );
    }
  }
  
  // // Priors abundance
  for(c in 1:B_intercept_columns) beta_raw_raw[c] ~ normal ( prior_mean_intercept[1], prior_mean_intercept[2] );
  if(C>B_intercept_columns) for(c in (B_intercept_columns+1):C) to_vector(beta_raw_raw[c]) ~ normal ( prior_mean_coefficients[1], prior_mean_coefficients[2]);
  
  // Hyper priors
  mix_p ~ beta(1,5);
  prec_coeff[1] ~ normal(prior_prec_intercept[1], prior_prec_intercept[2]);
  prec_coeff[2] ~ normal(prior_prec_slope[1],prior_prec_slope[2]);
  prec_sd ~ gamma(prior_prec_sd[1],prior_prec_sd[2]);
  prec_coeff ~ std_normal();
  to_vector(beta_raw_raw) ~ std_normal();
  
  // Random intercept
  if(is_random_effect>0){
    for(m in 1:(M-1)) random_effect_raw[,m] ~ std_normal(); 
    for(m in 1:(M-1)) random_effect_sigma_raw[m] ~ std_normal();
    for(m in 1:(M-1)) sigma_correlation_factor[m] ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
    
    random_effect_sigma_mu ~ std_normal();
    random_effect_sigma_sigma ~ std_normal();
    
    // If I have just one group
    zero_random_effect ~ std_normal();
  }
  if(ncol_X_random_eff[2]>0){
    for(m in 1:(M-1)) random_effect_raw_2[,m] ~ std_normal();
    for(m in 1:(M-1)) random_effect_sigma_raw_2[m] ~ std_normal();
    for(m in 1:(M-1)) sigma_correlation_factor_2[m] ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
    
  }
}
generated quantities {
  matrix[A, M] alpha_normalised = alpha;
  
  // // Rondom effect
  // matrix[ncol_X_random_eff_WINDOWS_BUG_FIX, M] beta_random_effect;
  // matrix[ncol_X_random_eff_WINDOWS_BUG_FIX_2, M] beta_random_effect_2;
  
  // LOO
  vector[TNS] log_lik = rep_vector(0, TNS);
  
  // These instructions regress out the effect of mean proportion to the overdispersion
  // This adjustment provide A overdispersion value that can be tested for a hypotheses for example differences between two conditions
  if(intercept_in_design){
    if(A > 1) for(a in 2:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }
  else{
    for(a in 1:A) alpha_normalised[a] = alpha[a] - (beta[a] * prec_coeff[2] );
  }
  
  // LOO
  if(enable_loo==1){

    matrix[M, N] mu;
    vector[N*M] mu_array;
    vector[N*M] precision_array;

    mu = (Q_ast * beta_raw)';

    // random intercept
    if(ncol_X_random_eff[1]> 0)
    mu = mu + append_row((X_random_effect * random_effect)', rep_row_vector(0, N));
    if(ncol_X_random_eff[2]>0 )
    mu = mu + append_row((X_random_effect_2 * random_effect_2)', rep_row_vector(0, N));


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

