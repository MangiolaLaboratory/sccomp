    array[] matrix reshape_to_3d_matrix(
      int M,
      int n_groups,
      int how_many_factors_in_random_design,
      matrix input_matrix,
      array[,] int group_factor_indexes_for_covariance
    ) {
      // PIVOT WIDER
      // increase of one dimension array[cell_type] matrix[group, factor]
      array[M-1] matrix[how_many_factors_in_random_design, n_groups] matrix_of_random_effects_raw;
      
      for(m in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
        
        // If I don't have the factor for one group 
        if(group_factor_indexes_for_covariance[j,i] == 0)
          matrix_of_random_effects_raw[m, j,i] = 0;
        else 
          matrix_of_random_effects_raw[m, j,i] = input_matrix[group_factor_indexes_for_covariance[j,i], m];
      }
      
      return matrix_of_random_effects_raw;
    }
    
    matrix reshape_to_2d_matrix(
      int M,
      int n_groups,
      int how_many_factors_in_random_design,
      array[] matrix matrix_of_random_effects,
      array[,] int group_factor_indexes_for_covariance,
      int ncol_X_random_eff
    ) {
      matrix[ncol_X_random_eff , M-1] random_effect;
      
      // Pivot longer
      for(m in 1:(M-1)) for(i in 1:n_groups) for(j in 1:how_many_factors_in_random_design)  {
        
        // If I don't have the factor for one group 
        if(group_factor_indexes_for_covariance[j,i] > 0)
          random_effect[group_factor_indexes_for_covariance[j,i], m] = matrix_of_random_effects[m,j,i];
      }
      
      return random_effect;
    }
    
    matrix get_random_effect_matrix(
      int M,                           // Number of categories/outcomes
      int n_groups,                    // Number of groups in the random effects design
      int how_many_factors_in_random_design,  // Number of factors in the random effects design
      int is_random_effect,            // Flag indicating if random effects are used (0/1)
      int ncol_X_random_eff,           // Number of columns in the random effects design matrix
      array[,] int group_factor_indexes_for_covariance,  // 2D array mapping factors to groups for covariance structure
      
      matrix random_effect_raw,      // Raw random effects matrix before transformation
      array[] vector random_effect_sigma,  // Standard deviations for each random effect
      array[] matrix sigma_correlation_factor  // Correlation matrices for random effects
      ){
        
        // PIVOT WIDER, as my columns should be covariates, not groups
        array[M-1] matrix[how_many_factors_in_random_design, n_groups] matrix_of_random_effects_raw = 
          reshape_to_3d_matrix(
            M, 
            n_groups, 
            how_many_factors_in_random_design, 
            random_effect_raw, 
            group_factor_indexes_for_covariance
          );
        
        // Design L
        array[M-1] matrix[how_many_factors_in_random_design, how_many_factors_in_random_design] L;
        array[M-1] matrix[how_many_factors_in_random_design, n_groups] matrix_of_random_effects;
        
        for(m in 1:(M-1)) L[m] = diag_pre_multiply(random_effect_sigma[m], sigma_correlation_factor[m]) ;
        for(m in 1:(M-1)) matrix_of_random_effects[m] = L[m] * matrix_of_random_effects_raw[m];
        
        // PIVOT LONGER 
        return reshape_to_2d_matrix(
          M, 
          n_groups, 
          how_many_factors_in_random_design, 
          matrix_of_random_effects, 
          group_factor_indexes_for_covariance,
          ncol_X_random_eff
        );
      }
      
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
  
      row_vector average_by_col(matrix X) {
      int rows_X = rows(X);
      int cols_X = cols(X);
      row_vector[cols_X] means;
      
      for (j in 1:cols_X) {
        means[j] = mean(X[, j]);
      }
      
      return means;
      
      
    }