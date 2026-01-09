#' simulate_data
#'
#' @description This function simulates data from a fitted model.
#'
#' @import dplyr
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @importFrom rlang quo_is_null
#' @importFrom SingleCellExperiment colData
#' @importFrom parallel detectCores
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column
#' @param .estimate_object The result of sccomp_estimate execution. This is used for sampling from real-data properties.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment
#' @param .sample A column name as symbol. The sample identifier
#' @param .cell_group A column name as symbol. The cell_group identifier
#' @param .coefficients The column names for coefficients, for example, c(b_0, b_1)
#' @param variability_multiplier A real scalar. This can be used for artificially increasing the variability of the simulation for benchmarking purposes.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' @param sig_figs Number of significant figures to use for Stan model output. Default is 9.
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
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
#' #     simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
#' #   }
#' # }
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
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("simulate_data", .data)
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
simulate_data.tbl = function(.data,
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
  
  
  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)
  .coefficients = enquo(.coefficients)
  
  #Check column class
  check_if_columns_right_class(.data, !!.sample, !!.cell_group)
  
  model_data = attr(.estimate_object, "model_input")
  original_data = attr(.estimate_object, "model_input")
  
  # Get formulas if not provided
  if(is.null(formula_variability)) {
    formula_variability = attr(.estimate_object, "formula_variability")
  }
  original_formula_composition = attr(.estimate_object, "formula_composition")
  
  # Validate dimensions before proceeding
  # Issue 3: Add dimension validation
  data_for_model =
    .data |>
    nest(data___ = -!!.sample) |>
    mutate(.exposure = sample(model_data$exposure, size = n(), replace = TRUE )) |>
    unnest(data___) |>
    data_simulation_to_model_input(
      formula_composition,
      formula_variability,  # Issue 2: Use formula_variability
      !!.sample, !!.cell_group, .exposure, !!.coefficients
    )
  names(data_for_model)  = names(data_for_model) |> stringr::str_c("_simulated")
  
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
  
  # Issue 2: Xa_simulated is already in data_for_model, so we don't need to add it again
  # Combine all data for Stan model
  stan_data = data_for_model |> 
    c(original_data_subset) |> 
    c(list(
      variability_multiplier = variability_multiplier,
      X_random_effect_simulated = X_random_effect_simulated,
      X_random_effect_2_simulated = X_random_effect_2_simulated
    ))
  
  # Get posterior draws - reuse pattern from replicate_data
  number_of_draws_in_the_fit = attr(.estimate_object, "fit") |> get_output_samples()
  number_of_draws = min(number_of_draws, number_of_draws_in_the_fit)
  
  draws_matrix = attr(.estimate_object, "fit")$draws(format = "matrix")
  
  if(number_of_draws > nrow(draws_matrix)) {
    number_of_draws = nrow(draws_matrix)
  }
  
  # Sample draws if needed (reuse pattern from replicate_data)
  if(number_of_draws < nrow(draws_matrix)) {
    draws_matrix = draws_matrix[sample(seq_len(nrow(draws_matrix)), size = number_of_draws),, drop = FALSE]
  }
  
  # Normalize all sum_to_zero_vector parameters to ensure they sum to exactly zero
  # This fixes floating-point precision issues when using generate_quantities
  # Uses the utility function from utilities.R
  draws_matrix = normalize_sum_to_zero_params(draws_matrix, "^beta_raw\\[")
  draws_matrix = normalize_sum_to_zero_params(draws_matrix, "^random_effect_raw\\[")
  draws_matrix = normalize_sum_to_zero_params(draws_matrix, "^random_effect_raw_2\\[")
  
  # Generate quantities - reuse pattern from replicate_data
  fit = mod_rng |> sample_safe(
    generate_quantities_fx,
    draws_matrix,
    data = stan_data,
    seed = mcmc_seed,
    parallel_chains = attr(.estimate_object, "fit")$metadata()$threads_per_chain, 
    threads_per_chain = cores,
    sig_figs = sig_figs
  )
  
  # Parse generated quantities - reuse parse_generated_quantities from replicate_data
  # Get cell group names from the simulated data (same as used in data_simulation_to_model_input)
  cell_group_names = 
    .data |>
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
  .data |>
    left_join(
      parsed_fit,
      by = c(quo_name(.sample), quo_name(.cell_group))
    )
  
  
}

#' @importFrom purrr when
#' @importFrom stats model.matrix
#'
#' @keywords internal
#' @noRd
#'
data_simulation_to_model_input =
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
    coeff_names = quo_names(.coefficients)
    
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
