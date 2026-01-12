#' sccomp_simulate
#'
#' @description This function simulates data from a fitted model.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null sym
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#' @importFrom lifecycle deprecate_warn
#' @importFrom tidyr crossing expand_grid cross_join
#'
#' @param fit The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment
#' @param new_data A tibble including sample-specific columns (sample identifier and factor columns from formula). If coefficients are provided separately, cell_group column is optional. Otherwise, should include cell_group column.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment
#' @param coefficients A data frame/tibble with cell-type specific coefficients. Must contain a column matching the cell_group column name, and columns matching the design matrix column names (e.g., "(Intercept)", "typeB"). If NULL, posterior beta_raw will be used.
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#' @param .coefficients (Deprecated) The column names for coefficients in new_data, for example, c(b_0, b_1). Use the 'coefficients' parameter instead. This parameter is placed last for backward compatibility.
#' 
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item \strong{sample} - A character column representing the sample name.
#'   \item \strong{type} - A factor column representing the type of the sample.
#'   \item \strong{phenotype} - A factor column representing the phenotype in the data.
#'   \item \strong{count} - An integer column representing the original cell counts.
#'   \item \strong{cell_group} - A character column representing the cell group identifier.
#'   \item \strong{b_0} - A numeric column representing the first coefficient used for simulation.
#'   \item \strong{b_1} - A numeric column representing the second coefficient used for simulation.
#'   \item \strong{generated_proportions} - A numeric column representing the generated proportions from the simulation.
#'   \item \strong{generated_counts} - An integer column representing the generated cell counts from the simulation.
#'   \item \strong{replicate} - An integer column representing the replicate number for each draw from the posterior distribution.
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
#'
#' @examples
#'
#' # print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' # \donttest{
#' #   if (instantiate::stan_cmdstan_exists()) {
#' #     data("counts_obj")
#' #     library(dplyr)
#'
#' #     estimate = sccomp_estimate(
#' #       counts_obj,
#' #       ~ type, ~1, "sample", "cell_group", "count",
#' #       cores = 1
#' #     )
#'
#' #     # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
#' #     counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)
#'
#' #     # Simulate data
#' #     sccomp_simulate(estimate, ~type, ~1, counts_obj, sample, cell_group, c(b_0, b_1))
#' #   }
#' # }
sccomp_simulate <- function(fit,
                          formula_composition,
                          formula_variability = NULL,
                          new_data = NULL,
                          coefficients = NULL,
                          .sample = NULL,
                          .cell_group = NULL,
                          variability_multiplier = 5,
                          number_of_draws = 1,
                          mcmc_seed = sample_seed(),
                          cores = detectCores(),
                          sig_figs = 9,
                          cache_stan_model = sccomp_stan_models_cache_dir,
                          .coefficients = NULL) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_simulate", fit)
}

#' @export
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom SingleCellExperiment counts
#' @importFrom purrr map_dbl
#' @importFrom purrr pmap
#' @importFrom readr read_file
#' @importFrom tibble column_to_rownames
#'
sccomp_simulate.sccomp_tbl = function(fit,
                             formula_composition,
                             formula_variability = NULL,
                                     new_data = NULL,
                                     coefficients = NULL,
                             .sample = NULL,
                             .cell_group = NULL,
                             variability_multiplier = 5,
                             number_of_draws = 1,
                             mcmc_seed = sample_seed(),
                             cores = detectCores(),
                             sig_figs = 9,
                                     cache_stan_model = sccomp_stan_models_cache_dir,
                                     .coefficients = NULL) {
  
  
  model_data = attr(fit, "model_input")
  original_data = attr(fit, "model_input")
  
  # Get sample and cell_group from fit attributes if not provided
  if(is.null(.sample)) {
    .sample_attr = attr(fit, ".sample")
    if(!is.null(.sample_attr)) {
      .sample = .sample_attr
    } else {
      stop("sccomp says: .sample must be provided or available in fit attributes")
    }
  } else {
  .sample = enquo(.sample)
  }
  
  if(is.null(.cell_group)) {
    .cell_group_attr = attr(fit, ".cell_group")
    if(!is.null(.cell_group_attr)) {
      .cell_group = .cell_group_attr
    } else {
      stop("sccomp says: .cell_group must be provided or available in fit attributes")
    }
  } else {
  .cell_group = enquo(.cell_group)
  }
  
  # Handle coefficients parameter
  # If coefficients table is provided, use it; otherwise check for .coefficients column names (backward compatibility)
  .coefficients_quo = enquo(.coefficients)
  
  # Deprecate .coefficients parameter in favor of coefficients table
  if(!rlang::quo_is_null(.coefficients_quo)) {
    lifecycle::deprecate_warn(
      "2.2.0",
      "sccomp_simulate(.coefficients)",
      details = "sccomp says: .coefficients parameter is deprecated. Please use the 'coefficients' parameter with a data frame/tibble instead. See ?sccomp_simulate for details."
    )
  }
  
  # Use new_data if provided, otherwise use count_data from fit
  if(is.null(new_data)) {
    .data = attr(fit, "count_data")
  } else {
    .data = new_data
  }
  
  # Check column class
  # If coefficients table is provided separately, new_data might not have cell_group column
  # In that case, skip cell_group validation for new_data
  if(is.null(coefficients) || quo_name(.cell_group) %in% colnames(.data)) {
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)
  }
  # If coefficients are provided separately, we'll get cell_group info from coefficients table
  
  # Get formulas if not provided
  if(is.null(formula_variability)) {
    formula_variability = attr(fit, "formula_variability")
  }
  original_formula_composition = attr(fit, "formula_composition")
  
  # Validate dimensions before proceeding
  # Use prepare_replicate_data machinery to get X_which for proper beta subsetting
  # This avoids padding beta with zeros and matches the approach used in sccomp_predict
  original_X = original_data$X
  original_Xa = original_data$Xa
  
  # Prepare data using the same machinery as replicate_data
  prepared_data = prepare_replicate_data(
    X = original_X,
    Xa = original_Xa,
    N = nrow(.data |> distinct(!!.sample)),
    intercept_in_design = original_data$intercept_in_design,
    X_random_effect = original_data$X_random_effect,
    X_random_effect_2 = original_data$X_random_effect_2,
    .sample = !!.sample,
    .cell_group = !!.cell_group,
    .count = quo(count),  # Dummy - not used for simulation
    formula_composition = formula_composition,
    original_formula_composition = original_formula_composition,
    formula_variability = formula_variability,
    new_data = .data |>
      nest(data___ = -!!.sample) |>
      mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) |>
      unnest(data___) |>
      select(!!.sample, any_of(c(parse_formula(formula_composition), parse_formula(formula_variability)))) |>
      distinct(),
    original_count_data = attr(fit, "count_data") |>
      distinct(!!.sample) |>
      mutate(dummy = 1)
  )
  
  # Get X_which for beta subsetting
  X_which = prepared_data$X_which
  XA_which = prepared_data$XA_which
  
  # Create data_for_model with simulated dimensions
  # Use prepared_data$X and prepared_data$Xa (which match X_which) instead of creating new ones
  # If coefficients are provided separately, expand .data to include all cell_groups
  if(!is.null(coefficients) && !quo_name(.cell_group) %in% colnames(.data)) {
    cell_group_colname = quo_name(.cell_group)
    all_cell_groups = coefficients[[cell_group_colname]]
    # Create a tibble with all cell_groups and cross join with .data
    cell_group_tibble = tibble(!!.cell_group := all_cell_groups)
    .data_expanded = 
      .data |>
      cross_join(cell_group_tibble)
  } else {
    .data_expanded = .data
  }
  
  # If coefficients table is provided, don't pass coefficients to data_to_simulation_covariates
  # (we'll extract them separately later)
  # Use an empty quosure if coefficients table is provided
  if(!is.null(coefficients)) {
    coefficients_quo_for_data_prep = quo(NULL)
  } else {
    coefficients_quo_for_data_prep = .coefficients_quo
  }
  
  data_for_model =
    .data_expanded |>
    nest(data___ = -!!.sample) |>
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) |>
    unnest(data___) |>
    data_to_simulation_covariates(
      formula_composition,
      formula_variability,
      !!.sample, !!.cell_group, .exposure, !!coefficients_quo_for_data_prep
    )
  names(data_for_model)  = names(data_for_model) |> stringr::str_c("_simulated")
  
  # Update X_simulated and Xa_simulated to use the prepared matrices (which match X_which)
  data_for_model$X_simulated = prepared_data$X
  data_for_model$Xa_simulated = prepared_data$Xa
  data_for_model$C_simulated = ncol(prepared_data$X)
  data_for_model$A_simulated = ncol(prepared_data$Xa)
  
  # Extract coefficients based on prepared_data$X columns (which match X_which)
  # Coefficients can be provided as:
  # 1. A separate table (coefficients parameter) - cell-type specific
  # 2. Column names in new_data (.coefficients parameter) - backward compatibility
  # 3. NULL - use posterior beta_raw
  
  # Get design matrix column names from prepared_data$X (which matches X_which)
  X_simulated_colnames = colnames(prepared_data$X)
  
  # Get cell group names
  # If coefficients table is provided, get cell_group names from there
  # Otherwise, get from .data_expanded (which includes cell_group)
  if(!is.null(coefficients)) {
    cell_group_colname = quo_name(.cell_group)
    
    # Validate coefficients table has cell_group column
    if(!cell_group_colname %in% colnames(coefficients)) {
      stop("sccomp says: coefficients table must contain a column matching the cell_group column name (", cell_group_colname, ")")
    }
    
    cell_group_names = 
      coefficients |>
      pull(!!sym(cell_group_colname)) |>
      unique() |>
      sort()
  } else {
    cell_group_names = 
      .data_expanded |>
      distinct(!!.cell_group) |>
      arrange(!!.cell_group) |>
      pull(!!.cell_group)
  }
  
  # Initialize beta_simulated_aligned
  beta_simulated_aligned = NULL
  
  # Case 1: coefficients table is provided
  if(!is.null(coefficients)) {
    # Create beta matrix: rows = cell_types, cols = design columns (matching X_simulated_colnames)
    beta_simulated_aligned = matrix(0, nrow = length(cell_group_names), ncol = length(X_simulated_colnames))
    rownames(beta_simulated_aligned) = cell_group_names
    colnames(beta_simulated_aligned) = X_simulated_colnames
    
    # Extract coefficients from coefficients table
    # Use column name directly instead of quosure
    coefficients_data = 
      coefficients |>
      select(all_of(cell_group_colname), any_of(X_simulated_colnames)) |>
      arrange(!!sym(cell_group_colname))
    
    # Map coefficients to design matrix columns
    # Coefficient column names should match design matrix column names
    for(col_name in X_simulated_colnames) {
      if(col_name %in% colnames(coefficients_data)) {
        # Match by cell_group and extract coefficient values
        for(i in 1:length(cell_group_names)) {
          cell_name = cell_group_names[i]
          matching_row = coefficients_data[[cell_group_colname]] == cell_name
          if(any(matching_row)) {
            beta_simulated_aligned[i, col_name] = coefficients_data[matching_row, col_name][[1]]
          }
        }
      }
    }
    
    # Note: Normalization to sum-to-zero happens after transpose (see below)
    # beta_simulated_aligned has [cell_types, design_columns], we need to normalize columns
    # but we'll normalize rows after transpose to [design_columns, cell_types]
    data_for_model$beta_simulated = beta_simulated_aligned
    
  } else if(!rlang::quo_is_null(.coefficients_quo)) {
    # Case 2: .coefficients column names provided in new_data (backward compatibility)
    # Get coefficient column names
    coeff_names = quo_names(.coefficients_quo)
    
    # Create beta matrix: rows = cell_types, cols = design columns (matching X_simulated_colnames)
    beta_simulated_aligned = matrix(0, nrow = length(cell_group_names), ncol = length(X_simulated_colnames))
    rownames(beta_simulated_aligned) = cell_group_names
    colnames(beta_simulated_aligned) = X_simulated_colnames
    
    # Extract coefficients from data
    coefficients_data = 
      .data |>
      select(!!.cell_group, all_of(coeff_names)) |>
      distinct() |>
      arrange(!!.cell_group)
    
    # Map coefficients to design matrix columns
    # If coefficient names match design column names, use them directly
    # Otherwise, map by position (assuming same order)
    for(i in 1:length(X_simulated_colnames)) {
      col_name = X_simulated_colnames[i]
      # Try to find matching coefficient column
      if(col_name %in% coeff_names) {
        beta_simulated_aligned[, col_name] = coefficients_data[[col_name]]
      } else if(i <= length(coeff_names)) {
        # Map by position if no name match
        beta_simulated_aligned[, i] = coefficients_data[[coeff_names[i]]]
      }
    }
    
    data_for_model$beta_simulated = beta_simulated_aligned
  }
  # Case 3: coefficients is NULL and .coefficients is NULL - use posterior beta_raw (no action needed)
  
  # Issue 3: Validate dimensions
  if(data_for_model$C_simulated > original_data$C) {
    stop("sccomp says: C_simulated (", data_for_model$C_simulated, ") cannot be larger than C (", original_data$C, ") from the fitted model. The simulated design matrix has more columns than the original model.")
  }
  if(data_for_model$M_simulated > original_data$M) {
    stop("sccomp says: M_simulated (", data_for_model$M_simulated, ") cannot be larger than M (", original_data$M, ") from the fitted model. The simulated data has more cell groups than the original model.")
  }
  if(data_for_model$A_simulated > original_data$A) {
    stop("sccomp says: A_simulated (", data_for_model$A_simulated, ") cannot be larger than A (", original_data$A, ") from the fitted model. The simulated variability design has more columns than the original model.")
  }
  
  # Extract necessary data from original model - reuse pattern from replicate_data
  # These fields are always present in model_input from sccomp_estimate (set in data_spread_to_model_input)
  original_data_subset = original_data[(names(original_data) %in% c("C", "M", "A", "A_intercept_columns", "intercept_in_design", "bimodal_mean_variability_association", "ncol_X_random_eff", "is_random_effect", "how_many_factors_in_random_design", "n_groups", "group_factor_indexes_for_covariance", "group_factor_indexes_for_covariance_2"))]
  
  # Create proper random effect design matrices
  # Stan requires matrices to have at least 1 column, so we use max(1, ncol) to ensure valid dimensions
  ncol_re1 = max(1, original_data_subset$ncol_X_random_eff[1])
  ncol_re2 = max(1, original_data_subset$ncol_X_random_eff[2])
  
  # Initialize with empty matrices
  X_random_effect_simulated = matrix(0, nrow = nrow(data_for_model$X_simulated), ncol = ncol_re1)
  X_random_effect_2_simulated = matrix(0, nrow = nrow(data_for_model$X_simulated), ncol = ncol_re2)
  
  # Create proper random effect design matrices if they exist in the original model
  if(original_data_subset$is_random_effect > 0 && 
     !is.null(original_data$X_random_effect) && 
     original_data_subset$ncol_X_random_eff[1] > 0) {
    
    # Get random effect elements from formula
    random_effect_elements = parse_formula_random_effect(formula_composition)
    original_grouping_names = original_formula_composition |> formula_to_random_effect_formulae() |> pull(grouping)
    
    if(length(original_grouping_names) > 0 && 
       (random_effect_elements$grouping %in% original_grouping_names[1]) |> any()) {
      
      # Create random effect design for simulated data
      random_effect_grouping = 
        formula_composition |>
        formula_to_random_effect_formulae() |>
        mutate(design = map2(
          formula, grouping,
          ~  get_random_effect_design3(.data, .x, .y, !!.sample,
                                       accept_NA_as_average_effect = TRUE )
        ))
      
      if((random_effect_grouping$grouping %in% original_grouping_names[1]) |> any()) {
        X_random_effect_simulated =
          random_effect_grouping |>
          filter(grouping==original_grouping_names[1]) |> 
          mutate(design_matrix = map(
            design,
            ~ ..1 |>
              select(!!.sample, group___label, value) |>
              filter(group___label %in% colnames(original_data$X_random_effect)) |> 
              pivot_wider(names_from = group___label, values_from = value) |>
              mutate(across(everything(), ~ .x |> replace_na(0)))
          )) |>
          pull(design_matrix) |> 
          _[[1]] |> 
          column_to_rownames(quo_name(.sample))
        
        # Ensure matrix has correct dimensions (pad with zeros if needed)
        if(ncol(X_random_effect_simulated) < ncol_re1) {
          padding = matrix(0, nrow = nrow(X_random_effect_simulated), 
                          ncol = ncol_re1 - ncol(X_random_effect_simulated))
          X_random_effect_simulated = cbind(X_random_effect_simulated, padding)
        } else if(ncol(X_random_effect_simulated) > ncol_re1) {
          X_random_effect_simulated = X_random_effect_simulated[, 1:ncol_re1, drop = FALSE]
        }
      }
    }
  }
  
  # Create second random effect design matrix if it exists
  if(original_data_subset$is_random_effect > 0 && 
     !is.null(original_data$X_random_effect_2) && 
     original_data_subset$ncol_X_random_eff[2] > 0) {
    
    random_effect_elements = parse_formula_random_effect(formula_composition)
    original_grouping_names = original_formula_composition |> formula_to_random_effect_formulae() |> pull(grouping)
    
    if(length(original_grouping_names) > 1 && 
       (random_effect_elements$grouping %in% original_grouping_names[2]) |> any()) {
      
      random_effect_grouping = 
        formula_composition |>
        formula_to_random_effect_formulae() |>
        mutate(design = map2(
          formula, grouping,
          ~  get_random_effect_design3(.data, .x, .y, !!.sample,
                                       accept_NA_as_average_effect = TRUE )
        ))
      
      if((random_effect_grouping$grouping %in% original_grouping_names[2]) |> any()) {
        X_random_effect_2_simulated =
          random_effect_grouping |>
          filter(grouping==original_grouping_names[2]) |> 
          mutate(design_matrix = map(
            design,
            ~ ..1 |>
              select(!!.sample, group___label, value) |>
              filter(group___label %in% colnames(original_data$X_random_effect_2)) |> 
              pivot_wider(names_from = group___label, values_from = value) |>
              mutate(across(everything(), ~ .x |> replace_na(0)))
          )) |>
          pull(design_matrix) |> 
          _[[1]] |> 
          column_to_rownames(quo_name(.sample))
        
        # Ensure matrix has correct dimensions
        if(ncol(X_random_effect_2_simulated) < ncol_re2) {
          padding = matrix(0, nrow = nrow(X_random_effect_2_simulated), 
                          ncol = ncol_re2 - ncol(X_random_effect_2_simulated))
          X_random_effect_2_simulated = cbind(X_random_effect_2_simulated, padding)
        } else if(ncol(X_random_effect_2_simulated) > ncol_re2) {
          X_random_effect_2_simulated = X_random_effect_2_simulated[, 1:ncol_re2, drop = FALSE]
        }
      }
    }
  }
  
  mod_rng = load_model("glm_multi_beta_binomial_simulate_data", threads = cores, cache_dir = cache_stan_model)
  
  # Check if coefficients are provided and validate/normalize them
  user_provided_beta = 0
  # When coefficients are not provided, pass empty array with correct dimensions for Stan
  # Stan expects array[C_simulated] vector[M_simulated], so we create a list of length C_simulated
  # Each element is an empty vector of length M_simulated
  beta_simulated_provided = lapply(1:data_for_model$C_simulated, function(i) {
    rep(0.0, data_for_model$M_simulated)
  })
  
  if(!is.null(data_for_model$beta_simulated) && nrow(data_for_model$beta_simulated) > 0 && ncol(data_for_model$beta_simulated) > 0) {
    # Coefficients are provided - validate dimensions
    # Note: beta_simulated has been aligned with prepared_data$X columns above
    # beta_simulated has rows=cell_types, cols=coefficients (matching X_simulated_colnames)
    # Stan expects rows=design_columns (C), cols=cell_types (M), so we need to transpose
    beta_provided = t(data_for_model$beta_simulated)
    
    # Validate dimensions after transpose
    # beta_provided should match length_X_which (number of design columns in X_which)
    expected_rows = length(X_which)
    if(nrow(beta_provided) != expected_rows) {
      stop("sccomp says: beta_simulated must have ", expected_rows, " design columns (matching X_which), but has ", nrow(beta_provided), ". This may happen if the design matrix columns don't match between the formula and the original model.")
    }
    if(ncol(beta_provided) != data_for_model$M_simulated) {
      stop("sccomp says: beta_simulated must have ", data_for_model$M_simulated, " cell groups (after transpose), but has ", ncol(beta_provided))
    }
    
    # Validate that coefficients sum to zero for each design column
    # User-provided coefficients must sum to exactly zero (within floating-point tolerance)
    # This is required for compositional models where log-ratios must sum to zero
    max_abs_val = max(abs(beta_provided))
    tolerance = .Machine$double.eps * 100 * max(max_abs_val, 1)
    
    # Check row sums using tidyverse approach
    row_sums = beta_provided |>
      as_tibble(.name_repair = "minimal") |>
      rowwise() |>
      mutate(row_sum = sum(c_across(everything()))) |>
      pull(row_sum)
    
    # Find columns that don't sum to zero
    invalid_cols = which(abs(row_sums) > tolerance)
    if(length(invalid_cols) > 0) {
      col_name = X_simulated_colnames[invalid_cols[1]]
      row_sum = row_sums[invalid_cols[1]]
      stop(
        "sccomp says: User-provided coefficients for design column '", col_name, 
        "' do not sum to zero (sum = ", round(row_sum, 10), "). ",
        "Coefficients in compositional models must sum to zero. ",
        "Please normalize your coefficients (e.g., set the last element to -sum(others)) ",
        "or subtract the mean from each coefficient vector."
      )
    }
    
    # Convert to list format for Stan (array[length_X_which] vector[M_simulated])
    # NOTE: We use vector[M_simulated] instead of sum_to_zero_vector[M_simulated] in Stan
    # to avoid floating-point precision issues. Small precision errors in sum-to-zero constraint
    # are acceptable since Stan doesn't enforce strict constraints with vector types.
    # However, user-provided coefficients must sum to zero within reasonable tolerance (checked above).
    # Stan expects a list/array where each element is a vector of length M_simulated
    # In R, we pass this as a list of vectors
    # The length matches X_which, not C_simulated
    beta_simulated_provided = lapply(1:nrow(beta_provided), function(c) {
      as.numeric(beta_provided[c, ])
    })
    
    user_provided_beta = 1
  }
  
  # Issue 2: Xa_simulated is already in data_for_model, so we don't need to add it again
  # Combine all data for Stan model
  stan_data = data_for_model |> 
    c(original_data_subset) |> 
    c(list(
      variability_multiplier = variability_multiplier,
      X_random_effect_simulated = X_random_effect_simulated,
      X_random_effect_2_simulated = X_random_effect_2_simulated,
      user_provided_beta = user_provided_beta,
      beta_simulated_provided = beta_simulated_provided,
      length_X_which = length(X_which),
      length_XA_which = length(XA_which),
      X_which = X_which,
      XA_which = XA_which
    ))
  
  # Get posterior draws - reuse pattern from replicate_data
  number_of_draws_in_the_fit = attr(fit, "fit") |> get_output_samples()
  number_of_draws = min(number_of_draws, number_of_draws_in_the_fit)
  
  draws_matrix = attr(fit, "fit")$draws(format = "matrix")
  
  if(number_of_draws > nrow(draws_matrix)) {
    number_of_draws = nrow(draws_matrix)
  }
  
  # Sample draws if needed (reuse pattern from replicate_data)
  if(number_of_draws < nrow(draws_matrix)) {
    draws_matrix = draws_matrix[sample(seq_len(nrow(draws_matrix)), size = number_of_draws),, drop = FALSE]
  }
  
  # Note: We no longer normalize sum-to-zero parameters since we use vector[M] instead of sum_to_zero_vector[M]
  # Small precision errors in sum-to-zero constraint are acceptable since Stan doesn't enforce strict constraints
  
  # Generate quantities - reuse pattern from replicate_data
  fit = mod_rng |> sample_safe(
    generate_quantities_fx,
    draws_matrix,
    data = stan_data,
    seed = mcmc_seed,
    parallel_chains = attr(fit, "fit")$metadata()$threads_per_chain, 
    threads_per_chain = cores,
    sig_figs = sig_figs
  )
  
  # Parse generated quantities - reuse parse_generated_quantities from replicate_data
  # Get cell group names from the simulated data (same as used in data_to_simulation_covariates)
  # Use .data_expanded which includes cell_group (either from original data or expanded from coefficients)
  cell_group_names = 
    .data_expanded |>
    distinct(!!.cell_group) |>
    arrange(!!.cell_group) |>
    pull(!!.cell_group)
  
  sample_names = rownames(data_for_model$X_simulated)
  
  parsed_fit =
    fit |>
    parse_generated_quantities(number_of_draws = number_of_draws) |>
    
    # Get sample name - reuse pattern from sccomp_replicate
    nest(data = -N) |>
    arrange(N) |>
    mutate(!!.sample := sample_names) |>
    unnest(data) |>
    
    # get cell type name - reuse pattern from sccomp_replicate
    nest(data = -M) |>
    mutate(!!.cell_group := cell_group_names) |>
    unnest(data) |>
    
    select(-N, -M)
  
  # Join with original data - reuse pattern from sccomp_predict
  # Use .data_expanded which includes cell_group (either from original data or expanded from coefficients)
  .data_expanded |>
    left_join(
      parsed_fit,
      by = c(quo_name(.sample), quo_name(.cell_group))
    )
  
  
}

#' DEPRECATED: simulate_data
#'
#' @description This function is DEPRECATED. Please use \code{\link{sccomp_simulate}} instead.
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients (Deprecated) The column names for coefficients in new_data, for example, c(b_0, b_1). Use the 'coefficients' parameter in sccomp_simulate() instead.
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#'
#' @export
#'
simulate_data <- function(.data,
                          .estimate_object,
                          formula_composition,
                          formula_variability = NULL,
                          .sample = NULL,
                          .cell_group = NULL,
                          .coefficients = NULL,
                          variability_multiplier = 5,
                          number_of_draws = 1,
                          mcmc_seed = sample_seed(),
                          cores = detectCores(),
                          sig_figs = 9,
                          cache_stan_model = sccomp_stan_models_cache_dir) {
  
  lifecycle::deprecate_warn(
    "2.2.0",
    "sccomp::simulate_data()",
    details = "sccomp says: simulate_data is deprecated. Please use sccomp_simulate() instead."
  )
  
  sccomp_simulate(
    fit = .estimate_object,
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = .data,
    .sample = .sample,
    .cell_group = .cell_group,
    .coefficients = .coefficients,
    variability_multiplier = variability_multiplier,
    number_of_draws = number_of_draws,
    mcmc_seed = mcmc_seed,
    cores = cores,
    sig_figs = sig_figs,
    cache_stan_model = cache_stan_model
  )
}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_to_simulation_covariates =
  function(.data, formula, formula_variability = ~ 1, .sample, .cell_type, .exposure, .coefficients, truncation_ajustment = 1, approximate_posterior_inference ){
    
    # Define the variables as NULL to avoid CRAN NOTES
    sd <- NULL
    . <- NULL
    
    
    # Prepare column same enquo
    .sample = enquo(.sample)
    .cell_type = enquo(.cell_type)
    .exposure = enquo(.exposure)
    .coefficients = enquo(.coefficients)
    
    factor_names = parse_formula(formula)
    factor_names_variability = parse_formula(formula_variability)
    
    sample_data =
      .data %>%
      select(!!.sample, any_of(c(factor_names, factor_names_variability))) %>%
      distinct() %>%
      arrange(!!.sample)
    
    # Create composition design matrix
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
    
    # Create variability design matrix (Issue 2: Use formula_variability)
    Xa =
      sample_data %>%
      model.matrix(formula_variability, data=.) %>%
      apply(2, function(x) {
        
        if(sd(x)==0 ) x
        else x |> scale(scale=FALSE)
        
      } ) %>%
      {
        .x = (.)
        rownames(.x) = sample_data %>% pull(!!.sample)
        .x
      }
    
    # Unique variability design (for compatibility)
    XA = Xa %>%
      as_tibble() %>%
      distinct()
    
    cell_cluster_names =
      .data %>%
      distinct(!!.cell_type) %>%
      arrange(!!.cell_type) %>%
      pull(!!.cell_type)
    
    # Extract coefficients
    # .coefficients is a quosure pointing to column names like c(b_0, b_1)
    # Use quo_names to extract the actual column names from the quosure
    # If .coefficients is NULL quosure, skip coefficient extraction
    if(rlang::quo_is_null(.coefficients)) {
      coeff_names = character(0)
      coefficients = matrix(nrow = 0, ncol = 0)
    } else {
      coeff_names = quo_names(.coefficients)
      
      if(length(coeff_names) > 0) {
        # Pivot to long format: cell_type | coefficient_name | value
        coefficients_long =
      .data %>%
          select(!!.cell_type, all_of(coeff_names)) %>%
      distinct() %>%
      arrange(!!.cell_type) %>%
          pivot_longer(cols = all_of(coeff_names), names_to = "coefficient_name", values_to = "value")
        
        # Pivot to wide format: coefficient_name | (cell_type_1) | (cell_type_2) | ...
        coefficients_wide = coefficients_long %>%
          pivot_wider(names_from = quo_name(.cell_type), values_from = value) %>%
          column_to_rownames("coefficient_name") %>%
          as.matrix()
        
        # Transpose to get: rows = coefficients, columns = cell_types
        coefficients = t(coefficients_wide)
      } else {
        coefficients = matrix(nrow = 0, ncol = 0)
      }
    }
    
    list(
      N = .data %>% distinct(!!.sample) %>% nrow(),
      M = .data %>% distinct(!!.cell_type) %>% nrow(),
      exposure = .data %>% distinct(!!.sample, !!.exposure) %>% arrange(!!.sample) %>% pull(!!.exposure),
      X = X,
      Xa = Xa,  # Issue 2: Add Xa for variability design
      XA = XA,
      C = ncol(X),
      A =  ncol(XA),
      beta = coefficients
    )
    
  }
