#' @importFrom purrr when
#' @importFrom stats model.matrix
#' @importFrom tidyr expand_grid
#' @importFrom stringr str_detect
#' @importFrom stringr str_remove_all
#' @importFrom purrr reduce
#' @importFrom purrr map_int
#' @importFrom stats as.formula
#'
#' @keywords internal
#' @noRd
#'
data_spread_to_model_input =
  function(
    .data_spread, formula, .sample, .cell_type, .count,
    truncation_ajustment = 1, approximate_posterior_inference,
    formula_variability = ~ 1,
    contrasts = NULL,
    bimodal_mean_variability_association = FALSE,
    use_data = TRUE,
    random_effect_elements,
    accept_NA_as_average_effect = FALSE){
    
    # Define the variables as NULL to avoid CRAN NOTES
    exposure <- NULL
    design <- NULL
    mat <- NULL
    factor___numeric <- NULL
    mean_idx <- NULL
    design_matrix <- NULL
    minus_sum <- NULL
    group___numeric <- NULL
    idx <- NULL
    group___label <- NULL
    parameter <- NULL
    group <- NULL
    design_matrix_col <- NULL
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .count = enquo(.count)
    .grouping_for_random_effect =
      random_effect_elements |>
      pull(grouping) |>
      unique() 
    
    if (length(.grouping_for_random_effect)==0 ) .grouping_for_random_effect = "random_effect"
    
    
    X  =
      
      .data_spread |>
      get_design_matrix(
        # Drop random intercept
        formula |>
          as.character() |>
          str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
          paste(collapse="") |>
          as.formula(),
        !!.sample,
        accept_NA_as_average_effect = accept_NA_as_average_effect
      )
    
    Xa  =
      .data_spread |>
      get_design_matrix(
        # Drop random intercept
        formula_variability |>
          as.character() |>
          str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
          paste(collapse="") |>
          as.formula() ,
        !!.sample,
        accept_NA_as_average_effect = accept_NA_as_average_effect
      )
    
    XA = Xa %>%
      as_tibble() %>%
      distinct()
    
    A = ncol(XA);
    Ar = nrow(XA);
    
    factor_names = parse_formula(formula)
    factor_names_variability = parse_formula(formula_variability)
    cell_cluster_names = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% colnames()
    
    # Random intercept
    if(nrow(random_effect_elements)>0 ) {
      
      
      #check_random_effect_design(.data_spread, any_of(factor_names), random_effect_elements, formula, X)
      random_effect_grouping = formula |>
        formula_to_random_effect_formulae() |>
        mutate(design = map2(
          formula, grouping,
          ~ {
            get_random_effect_design3(.data_spread, .x, .y, !!.sample )
          }))
      
      # Actual parameters, excluding for the sum to one parameters
      is_random_effect = 1
      
      random_effect_grouping =
        random_effect_grouping |>
        mutate(design_matrix = map(
          design,
          ~ ..1 |>
            select(!!.sample, group___label, value) |>
            pivot_wider(names_from = group___label, values_from = value) |>
            mutate(across(everything(), ~ .x |> replace_na(0)))
        )) 
      
      
      X_random_effect = 
        random_effect_grouping |> 
        pull(design_matrix) |> 
        _[[1]]  |>  
        as_matrix(rownames = quo_name(.sample))
      
      # Separate NA group column into X_random_effect_unseen
      X_random_effect_unseen = X_random_effect[, colnames(X_random_effect) |> str_detect("___NA$"), drop = FALSE]
      X_random_effect = X_random_effect[, !colnames(X_random_effect) |> str_detect("___NA$"), drop = FALSE]
      
      # For now that stan does not have tuples, I just allow max two levels
      if(random_effect_grouping |> nrow() > 2) stop("sccomp says: at the moment sccomp allow max two groupings")
      # This will be modularised with the new stan
      if(random_effect_grouping |> nrow() > 1){
        X_random_effect_2 =   
          random_effect_grouping |> 
          pull(design_matrix) |> 
          _[[2]] |>  
          as_matrix(rownames = quo_name(.sample))
        
        # Separate NA group column into X_random_effect_2_unseen  
        X_random_effect_2_unseen = X_random_effect_2[, colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
        X_random_effect_2 = X_random_effect_2[, !colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
      }
      
      else X_random_effect_2 =  X_random_effect[,0,drop=FALSE]
      
      n_random_eff = random_effect_grouping |> nrow()
      
      ncol_X_random_eff =
        random_effect_grouping |>
        mutate(n = map_int(design, ~.x |> distinct(group___numeric) |> nrow())) |>
        pull(n) 
      
      if(ncol_X_random_eff |> length() < 2) ncol_X_random_eff[2] = 0
      
      # TEMPORARY
      group_factor_indexes_for_covariance = 
        X_random_effect |> 
        colnames() |> 
        enframe(value = "parameter", name = "order")  |> 
        separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
        complete(factor, group, fill = list(order=0)) |> 
        select(-parameter) |> 
        pivot_wider(names_from = group, values_from = order)  |> 
        as_matrix(rownames = "factor")
      
      n_groups = group_factor_indexes_for_covariance |> ncol()
      
      # This will be modularised with the new stan
      if(random_effect_grouping |> nrow() > 1)
        group_factor_indexes_for_covariance_2 = 
        X_random_effect_2 |> 
        colnames() |> 
        enframe(value = "parameter", name = "order")  |> 
        separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
        complete(factor, group, fill = list(order=0)) |> 
        select(-parameter) |> 
        pivot_wider(names_from = group, values_from = order)  |> 
        as_matrix(rownames = "factor")
      else group_factor_indexes_for_covariance_2 = matrix()[0,0, drop=FALSE]
      
      n_groups = n_groups |> c(group_factor_indexes_for_covariance_2 |> ncol())
      
      how_many_factors_in_random_design = list(group_factor_indexes_for_covariance, group_factor_indexes_for_covariance_2) |> map_int(nrow)
      
      
    } else {
      X_random_effect = matrix(rep(1, nrow(.data_spread)))[,0, drop=FALSE]
      X_random_effect_2 = matrix(rep(1, nrow(.data_spread)))[,0, drop=FALSE] # This will be modularised with the new stan
      is_random_effect = 0
      ncol_X_random_eff = c(0,0)
      n_random_eff = 0
      n_groups = c(0,0)
      how_many_factors_in_random_design = c(0,0)
      group_factor_indexes_for_covariance = matrix()[0,0, drop=FALSE]
      group_factor_indexes_for_covariance_2 = matrix()[0,0, drop=FALSE] # This will be modularised with the new stan
    }
    
    
    y = .data_spread %>% select(-any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% as_matrix(rownames = quo_name(.sample))
    
    # If proportion ix 0 issue
    is_proportion = y |> as.numeric() |> max()  |> between(0,1) |> all()
    if(is_proportion){
      y_proportion = y
      y = y[0,,drop = FALSE]
    }
    else{
      y = y
      y_proportion = y[0,,drop = FALSE]
    }
    
    data_for_model =
      list(
        N = .data_spread %>% nrow(),
        M = .data_spread %>% select(-!!.sample, -any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% ncol(),
        exposure = .data_spread$exposure,
        is_proportion = is_proportion,
        y = y,
        y_proportion = y_proportion,
        X = X,
        XA = XA,
        Xa = Xa,
        C = ncol(X),
        A = A,
        Ar = Ar,
        truncation_ajustment = truncation_ajustment,
        is_vb = as.integer(approximate_posterior_inference),
        bimodal_mean_variability_association = bimodal_mean_variability_association,
        use_data = use_data,
        
        # Random intercept
        is_random_effect = is_random_effect,
        ncol_X_random_eff = ncol_X_random_eff,
        n_random_eff = n_random_eff,
        n_groups  = n_groups,
        X_random_effect = X_random_effect,
        X_random_effect_2 = X_random_effect_2,
        group_factor_indexes_for_covariance = group_factor_indexes_for_covariance,
        group_factor_indexes_for_covariance_2 = group_factor_indexes_for_covariance_2,
        how_many_factors_in_random_design = how_many_factors_in_random_design,
        
        # For parallel chains
        grainsize = 1,
        
        ## LOO
        enable_loo = FALSE
      )
    
    # Add censoring
    data_for_model$is_truncated = 0
    data_for_model$truncation_up = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_down = matrix(rep(-1, data_for_model$M * data_for_model$N), ncol = data_for_model$M)
    data_for_model$truncation_not_idx = seq_len(data_for_model$M*data_for_model$N)
    data_for_model$TNS = length(data_for_model$truncation_not_idx)
    data_for_model$truncation_not_idx_minimal = matrix(c(1,1), nrow = 1)[0,,drop=FALSE]
    data_for_model$TNIM = 0
    
    # Add parameter factor dictionary
    data_for_model$factor_parameter_dictionary = tibble()
    
    if(.data_spread  |> select(any_of(parse_formula(formula))) |> lapply(class) %in% c("factor", "character") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |> bind_rows(
        # For discrete
        .data_spread  |>
          select(any_of(parse_formula(formula)))  |>
          distinct()  |>
          
          # Drop numerical
          select_if(function(x) !is.numeric(x)) |>
          pivot_longer(everything(), names_to =  "factor", values_to = "parameter") %>%
          unite("design_matrix_col", c(`factor`, parameter), sep="", remove = FALSE)  |>
          select(-parameter) |>
          filter(design_matrix_col %in% colnames(data_for_model$X)) %>%
          distinct()
        
      )
    
    # For continuous
    if(.data_spread  |> select(all_of(parse_formula(formula))) |> lapply(class) |> equals("numeric") |> any())
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |>
      bind_rows(
        tibble(
          design_matrix_col =  .data_spread  |>
            select(all_of(parse_formula(formula)))  |>
            distinct()  |>
            
            # Drop numerical
            select_if(function(x) is.numeric(x)) |>
            names()
        ) |>
          mutate(`factor` = design_matrix_col)
      )
    
    # If constrasts is set it is a bit more complicated
    if(! is.null(contrasts))
      data_for_model$factor_parameter_dictionary =
      data_for_model$factor_parameter_dictionary |>
      distinct() |>
      expand_grid(parameter=contrasts) |>
      filter(str_detect(parameter, design_matrix_col )) |>
      select(-design_matrix_col) |>
      rename(design_matrix_col = parameter) |>
      distinct()
    
    data_for_model$intercept_in_design = X[,1] |> unique() |> identical(1)
    
    
    if (data_for_model$intercept_in_design | length(factor_names_variability) == 0) {
      data_for_model$A_intercept_columns = 1
    } else {
      data_for_model$A_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    
    if (data_for_model$intercept_in_design ) {
      data_for_model$B_intercept_columns = 1
    } else {
      data_for_model$B_intercept_columns = 
        .data_spread |> 
        select(any_of(factor_names[1])) |> 
        distinct() |> 
        nrow()
    }
    
    # Default all grouping known. This is used for data generation to estimate unknown groupings.
    data_for_model$unknown_grouping = c(FALSE, FALSE)
    
    
    # Return
    data_for_model
  }

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_simulation_to_model_input =
  function(.data, formula, .sample, .cell_type, .exposure, .coefficients, truncation_ajustment = 1, approximate_posterior_inference ){
    
    # Define the variables as NULL to avoid CRAN NOTES
    sd <- NULL
    . <- NULL
    
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .exposure = enquo(.exposure)
    .coefficients = enquo(.coefficients)
    
    factor_names = parse_formula(formula)
    
    sample_data =
      .data %>%
      select(!!.sample, any_of(factor_names)) %>%
      distinct() %>%
      arrange(!!.sample)
    X =
      sample_data %>%
      model.matrix(formula, data=.) %>%
      apply(2, function(x) {
        
        if(sd(x)==0 ) x
        else x |> scale(scale=FALSE)
        
      } ) %>%
      {
        .x = (.)
        rownames(.x) = sample_data %>% pull(!!.sample)
        .x
      }
    
    if(factor_names == "1") XA = X[,1, drop=FALSE]
    else XA = X[,c(1,2), drop=FALSE]
    
    XA = XA |> 
      as_tibble()  |> 
      distinct()
    
    cell_cluster_names =
      .data %>%
      distinct(!!.cell_type) %>%
      arrange(!!.cell_type) %>%
      pull(!!.cell_type)
    
    coefficients =
      .data %>%
      select(!!.cell_type, !!.coefficients) %>%
      unnest(!!.coefficients) %>%
      distinct() %>%
      arrange(!!.cell_type) %>%
      as_matrix(rownames = quo_name(.cell_type)) %>%
      t()
    
    list(
      N = .data %>% distinct(!!.sample) %>% nrow(),
      M = .data %>% distinct(!!.cell_type) %>% nrow(),
      exposure = .data %>% distinct(!!.sample, !!.exposure) %>% arrange(!!.sample) %>% pull(!!.exposure),
      X = X,
      XA = XA,
      C = ncol(X),
      A =  ncol(XA),
      beta = coefficients
    )
    
  }
