#' sccomp_replicate
#'
#' @description This function replicates counts from a real-world dataset.
#'
#' @param fit The result of sccomp_estimate.
#' @param formula_composition A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param formula_variability A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.
#' @param number_of_draws An integer. How may copies of the data you want to draw from the model joint posterior distribution.
#' @param mcmc_seed An integer. Used for Markov-chain Monte Carlo reproducibility. By default a random number is sampled from 1 to 999999. This itself can be controlled by set.seed()
#' @param cache_stan_model A character string specifying the cache directory for compiled Stan models. 
#'                        The sccomp version will be automatically appended to ensure version isolation.
#'                        Default is `sccomp_stan_models_cache_dir` which points to `~/.sccomp_models`.
#'
#' @return A tibble `tbl` with cell_group-wise statistics
#'
#' @return A tibble (`tbl`), with the following columns:
#' \itemize{
#'   \item \strong{cell_group} - A character column representing the cell group being tested.
#'   \item \strong{sample} - A factor column representing the sample name from which data was generated.
#'   \item \strong{generated_proportions} - A numeric column representing the proportions generated from the model.
#'   \item \strong{generated_counts} - An integer column representing the counts generated from the model.
#'   \item \strong{replicate} - An integer column representing the replicate number, where each row corresponds to a different replicate of the data.
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @importFrom tibble deframe
#' @importFrom tibble column_to_rownames
#' @export
#'
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
#'     data("counts_obj")
#'
#'     sccomp_estimate(
#'       counts_obj,
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_replicate()
#'   }
#' }
sccomp_replicate <- function(fit,
                             formula_composition = NULL,
                             formula_variability = NULL,
                             number_of_draws = 1,
                             mcmc_seed = sample_seed(),
                             cache_stan_model = sccomp_stan_models_cache_dir) {
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_replicate", fit)
}

#' @export
#'
sccomp_replicate.sccomp_tbl = function(fit,
                                       formula_composition = NULL,
                                       formula_variability = NULL,
                                       number_of_draws = 1,
                                       mcmc_seed = sample_seed(),
                                       cache_stan_model = sccomp_stan_models_cache_dir){
  
  .sample = attr(fit, ".sample")
  .cell_group = attr(fit, ".cell_group")
  
  sample_names =
    fit |>
    attr("count_data") |>
    distinct(!!.sample) |> 
    pull(!!.sample)
  
  rng =
    replicate_data(
      fit,
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      number_of_draws = number_of_draws,
      mcmc_seed = mcmc_seed,
      cache_stan_model = cache_stan_model
    )
  
  model_input = attr(fit, "model_input")
  
  # mean generated
  rng |>
    
    # Parse
    parse_generated_quantities(number_of_draws = number_of_draws) |>
    
    # Get sample name
    nest(data = -N) |>
    arrange(N) |>
    mutate(!!.sample := sample_names) |>
    unnest(data) |>
    
    # get cell type name
    nest(data = -M) |>
    mutate(!!.cell_group := colnames(model_input$y)) |>
    unnest(data) |>
    
    select(-N, -M) |>
    select(!!.cell_group, !!.sample, everything())
  
}


#' Prepare data for replicate_data function
#' 
#' @description
#' Internal function that prepares the data for the replicate_data function. It handles:
#' - Data validation and preprocessing
#' - Design matrix creation
#' - Random effects handling
#' - Unseen group handling
#' 
#' @param X Original design matrix
#' @param Xa Original variability design matrix
#' @param N Original number of samples
#' @param intercept_in_design Whether intercept is in design
#' @param X_random_effect Original random effect design matrix
#' @param X_random_effect_2 Original second random effect design matrix
#' @param .sample Quosure for the sample identifier column
#' @param .cell_group Quosure for the cell group column
#' @param .count Quosure for the count column
#' @param formula_composition Formula for the composition model
#' @param formula_variability Formula for the variability model
#' @param new_data New data to generate predictions for. If NULL, uses the original data
#' @param original_count_data Original count data from the model
#' 
#' @return A list containing:
#' - model_input: The prepared model input data
#' - X_which: Indices for the composition design matrix
#' - XA_which: Indices for the variability design matrix
#' - X_random_effect_which: Indices for the first random effect design matrix
#' - X_random_effect_which_2: Indices for the second random effect design matrix
#' - create_intercept: Boolean indicating if intercept should be created
#' 
#' @noRd
prepare_replicate_data = function(X,
                                  Xa,
                                  N,
                                  intercept_in_design,
                                  X_random_effect,
                                  X_random_effect_2,
                                  .sample,
                                  .cell_group,
                                  .count,
                                  formula_composition,
                                  original_formula_composition,
                                  formula_variability,
                                  new_data = NULL,
                                  original_count_data) {
  
  
  .sample = enquo(.sample)
  .count = enquo(.count)
  
  # New data
  if(new_data |> is.null())
    new_data = original_count_data 
  
  # Check that original_count_data dont have duplicated sample names
  if(original_count_data |> count(!!.sample) |> pull(n) |> max() > 1)
    stop("sccomp says: your original_count_data has duplicated sample names. Please ensure that the sample names are unique.") 
  
  # Check that new_data dont have duplicated sample names
  if(new_data |> count(!!.sample) |> pull(n) |> max() > 1)
    stop("sccomp says: your new_data has duplicated sample names. Please ensure that the sample names are unique.") 
  
  
  # If seurat
  else if(new_data |> is("Seurat")) new_data = new_data[[]]
  
  # Check if the input new data is not suitable
  if(!parse_formula(formula_composition) %in% colnames(new_data) |> all())
    stop("sccomp says: your `new_data` might be malformed. It might have the covariate columns with multiple values for some element of the \"%s\" column. As a generic example, a sample identifier (\"Sample_123\") might be associated with multiple treatment values, or age values.")
  
  # Match factors with old data
  nrow_new_data = nrow(new_data)
  
  new_exposure = new_data |>
    nest(data = -!!.sample) |>
    mutate(exposure = map_dbl(
      data,
      ~{
        if (quo_name(.count) %in% colnames(.x)) .x |> pull(!!.count) |> sum()
        else 5000
      })) |>
    select(!!.sample, exposure) |>
    deframe() |>
    as.array()
  
  # Update data, merge with old data because
  # I need the same ordering of the design matrix
  old_data = original_count_data |>
    
    # Change sample names to make unique
    mutate(dummy = "OLD") |>
    tidyr::unite(!!.sample, c(!!.sample, dummy), sep="___")
  
  # Harmonise factors
  new_data = new_data |> as_tibble() |> harmonise_factor_levels(old_data)
  
  # check that for each column of the old data. The new data has values that are found in the old data, omit NA , ignore !!.sample column
  map(
    old_data %>%
      # Just apply to categorical
      select(where(~ is.factor(.) || is.character(.))) %>%
      names() %>%
      setdiff(quo_name(.sample)),
    ~ if (any(!new_data[[.x]][!is.na(new_data[[.x]])] %in% old_data[[.x]])) {
      stop(
        paste(
          "sccomp says: The values for column `",.x,"` in the new data must be among the factor levels in the original data. Please ensure that the factor levels are consistent between the two datasets."
        )
      )
    }
  )
  
  new_data =  old_data |> bind_rows( new_data ) 
  
  new_X =
    new_data |>
    get_design_matrix(
      # Drop random intercept
      formula_composition |>
        as.character() |>
        str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
        paste(collapse="") |>
        as.formula(),
      !!.sample, 
      accept_NA_as_average_effect = TRUE
    ) |>
    tail(nrow_new_data) %>%
    # Remove columns that are not in the original design matrix
    .[,colnames(.) %in% colnames(X), drop=FALSE]
  
  # Check that all effect combination were present when the model was fitted
  check_missing_parameters(
    new_X |> colnames(), 
    X |> colnames()
  ) 
  
  X_which =
    colnames(new_X) |>
    match(
      X %>%
        colnames()
    ) |>
    na.omit() |>
    as.array()
  
  # Variability
  new_Xa =
    new_data |>
    get_design_matrix(
      # Drop random intercept
      formula_variability |>
        as.character() |>
        str_remove_all("\\+ ?\\(.+\\|.+\\)") |>
        paste(collapse="") |>
        as.formula(),
      !!.sample, 
      accept_NA_as_average_effect = TRUE
    ) |>
    tail(nrow_new_data) %>%
    # Remove columns that are not in the original design matrix
    .[,colnames(.) %in% colnames(Xa), drop=FALSE]
  
  XA_which =
    colnames(new_Xa) |>
    match(
      Xa %>%
        colnames()
    ) |>
    na.omit() |>
    as.array()
  
  # If I want to replicate data with intercept and I don't have intercept in my fit
  create_intercept =
    !intercept_in_design &
    "(Intercept)" %in% colnames(new_X)
  if(create_intercept) warning("sccomp says: your estimated model is intercept free, while your desired replicated data do have an intercept term. The intercept estimate will be calculated averaging your first factor in your formula ~ 0 + <factor>. If you don't know the meaning of this warning, this is likely undesired, and please reconsider your formula for replicate_data()")
  
  # Original grouping
  original_grouping_names = original_formula_composition |> formula_to_random_effect_formulae() |> pull(grouping)
  
  # Random intercept
  random_effect_elements = parse_formula_random_effect(formula_composition)
  
  # Initialize unseen random effect variables
  new_X_random_effect_unseen = matrix(rep(0, nrow_new_data))[,0, drop=FALSE]
  new_X_random_effect_2_unseen = matrix(rep(0, nrow_new_data))[,0, drop=FALSE]
  
  # Set default random intercept
  X_random_effect_which = array()[0]
  new_X_random_effect = matrix(rep(0, nrow_new_data))[,0, drop=FALSE]
  
  # setup default unknown_grouping variable for generated quantities
  unknown_grouping = c(FALSE, FALSE)
  
  #check_random_effect_design(.data_spread, any_of(factor_names), random_effect_elements, formula, X)
  random_effect_grouping = 
    formula_composition |>
    formula_to_random_effect_formulae() |>
    mutate(design = map2(
      formula, grouping,
      ~  get_random_effect_design3(new_data, .x, .y, !!.sample,
                                   accept_NA_as_average_effect = TRUE )
    ))
  
  
  
  if((random_effect_grouping$grouping %in% original_grouping_names[1]) |> any() && !unknown_grouping[1]) {
    new_X_random_effect =
      random_effect_grouping |>
      filter(grouping==original_grouping_names[1]) |> 
      mutate(design_matrix = map(
        design,
        ~ ..1 |>
          select(!!.sample, group___label, value) |>
          
          # Some combinations might not have been present in a specific group so the parameter does not exist
          filter(group___label %in% colnames(X_random_effect)) |> 
          
          pivot_wider(names_from = group___label, values_from = value) |>
          mutate(across(everything(), ~ .x |> replace_na(0)))
      )) |>
      # Merge
      pull(design_matrix) |> 
      _[[1]] |> 
      column_to_rownames(quo_name(.sample))  |>
      tail(nrow_new_data) 
    
    # Separate NA group column into new_X_random_effect_unseen
    new_X_random_effect_unseen = new_X_random_effect[, colnames(new_X_random_effect) |> str_detect("___NA$"), drop = FALSE]
    new_X_random_effect = new_X_random_effect[, !colnames(new_X_random_effect) |> str_detect("___NA$"), drop = FALSE]
    
    # Check that all effect combination were present when the model was fitted
    check_missing_parameters(
      new_X_random_effect |> colnames(), 
      X_random_effect |> colnames()
    ) 
    
    X_random_effect_which =
      colnames(new_X_random_effect) |>
      match(
        X_random_effect %>%
          colnames()
      ) |>
      as.array()
  }
  
  # Set default X random intercept
  X_random_effect_which_2 = array()[0]
  new_X_random_effect_2 = matrix(rep(0, nrow_new_data))[,0, drop=FALSE]
  
  if((random_effect_grouping$grouping %in% original_grouping_names[2]) |> any()){
    new_X_random_effect_2 =
      random_effect_grouping |>
      filter(grouping==original_grouping_names[2]) |> 
      mutate(design_matrix = map(
        design,
        ~ ..1 |>
          select(!!.sample, group___label, value) |>
          
          # Some combinations might not have been present in a specific group so the parameter does not exist
          filter(group___label %in% colnames(X_random_effect_2)) |> 
          
          pivot_wider(names_from = group___label, values_from = value) |>
          mutate(across(everything(), ~ .x |> replace_na(0)))
      )) |>
      # Merge
      pull(design_matrix) |> 
      _[[1]] |> 
      column_to_rownames(quo_name(.sample))  |>
      tail(nrow_new_data)
    
    # Separate NA group column into new_X_random_effect_2_unseen
    new_X_random_effect_2_unseen = new_X_random_effect_2[, colnames(new_X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
    new_X_random_effect_2 = new_X_random_effect_2[, !colnames(new_X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
    
    # Check that all effect combination were present when the model was fitted
    check_missing_parameters(
      new_X_random_effect_2 |> colnames(), 
      X_random_effect_2 |> colnames()
    )
    
    X_random_effect_which_2 =
      colnames(new_X_random_effect_2) |>
      match(
        X_random_effect_2 %>%
          colnames()
      ) |>
      as.array()
  }
  
  # Prepare the list of inputs to the model
  list(
    X = new_X,
    Xa = new_Xa,
    N = nrow_new_data,  
    exposure = new_exposure,
    X_random_effect = new_X_random_effect,
    X_random_effect_2 = new_X_random_effect_2,
    X_random_effect_unseen = new_X_random_effect_unseen,
    X_random_effect_2_unseen = new_X_random_effect_2_unseen,
    ncol_X_random_eff_new = c(ncol(new_X_random_effect), ncol(new_X_random_effect_2)),  
    unknown_grouping = unknown_grouping,
    ncol_X_random_eff_unseen = c(ncol(new_X_random_effect_unseen), ncol(new_X_random_effect_2_unseen)),
    
    X_which = X_which,
    XA_which = XA_which,
    X_random_effect_which = X_random_effect_which,
    X_random_effect_which_2 = X_random_effect_which_2,
    create_intercept = create_intercept
  )
  
}

#' Generate replicated data from a fitted sccomp model
#' 
#' @description
#' This function generates replicated data from a fitted sccomp model. It can be used to:
#' - Generate predictions for new data
#' - Simulate data from the fitted model
#' - Generate posterior predictive samples
#' 
#' @param .data A sccomp model object
#' @param formula_composition Formula for the composition model. If NULL, uses the formula from .data
#' @param formula_variability Formula for the variability model. If NULL, uses the formula from .data
#' @param new_data New data to generate predictions for. If NULL, uses the original data
#' @param number_of_draws Number of posterior draws to use for generation
#' @param mcmc_seed Random seed for reproducibility
#' @param cores Number of CPU cores to use
#' 
#' @return A cmdstanr model object containing the generated quantities
#' 
#' @examples
#' # Generate predictions for new data
#' replicate_data(fitted_model, new_data = new_samples)
#' 
#' # Generate 1000 posterior predictive samples
#' replicate_data(fitted_model, number_of_draws = 1000)
#' 
#' @noRd
replicate_data = function(.data,
                          formula_composition = NULL,
                          formula_variability = NULL,
                          new_data = NULL,
                          number_of_draws = 1,
                          mcmc_seed = sample(1e5, 1),
                          cores = detectCores(),
                          cache_stan_model = sccomp_stan_models_cache_dir){
  
  # Extract required components from .data
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .count = attr(.data, ".count")
  
  # Get formulas if not provided
  if(is.null(formula_composition)) formula_composition = attr(.data, "formula_composition")
  if(is.null(formula_variability)) formula_variability = attr(.data, "formula_variability")
  
  # create model input 
  model_input = attr(.data, "model_input")
  
  # Prepare data
  prepared_data = prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!.sample,
    .cell_group = !!.cell_group,
    .count = !!.count,
    formula_composition = formula_composition,
    original_formula_composition = .data |> attr("formula_composition"),
    formula_variability = formula_variability,
    new_data = new_data,
    original_count_data =
      .data |>
      attr("count_data") |>
      .subset(!!.sample) 
  )
  
  # Original input 
  model_input$X_original = model_input$X
  model_input$N_original = model_input$N
  
  # New input
  model_input$X = prepared_data$X
  model_input$Xa = prepared_data$Xa 
  model_input$N = prepared_data$N
  model_input$exposure = prepared_data$exposure
  model_input$X_random_effect = prepared_data$X_random_effect
  model_input$X_random_effect_2 = prepared_data$X_random_effect_2
  model_input$X_random_effect_unseen = prepared_data$X_random_effect_unseen
  model_input$X_random_effect_2_unseen = prepared_data$X_random_effect_2_unseen
  model_input$ncol_X_random_eff_new = prepared_data$ncol_X_random_eff_new
  model_input$unknown_grouping = prepared_data$unknown_grouping
  model_input$ncol_X_random_eff_unseen = prepared_data$ncol_X_random_eff_unseen
  
  # Add subset of coefficients
  # Add subset of coefficients
  model_input$length_X_which = length(prepared_data$X_which)
  model_input$length_XA_which = length(prepared_data$XA_which)
  model_input$X_which = prepared_data$X_which
  model_input$XA_which = prepared_data$XA_which
  
  # Add random effect coefficients
  model_input$X_random_effect_which = prepared_data$X_random_effect_which
  model_input$X_random_effect_which_2 = prepared_data$X_random_effect_which_2
  model_input$length_X_random_effect_which = c(
    length(prepared_data$X_random_effect_which),
    length(prepared_data$X_random_effect_which_2)
  )
  
  # Should I create an intercept for generate quantities?
  model_input$create_intercept = prepared_data$create_intercept
  
  number_of_draws_in_the_fit = attr(.data, "fit") |>  get_output_samples()
  
  # To avoid error in case of a NULL posterior sample
  number_of_draws = min(number_of_draws, number_of_draws_in_the_fit)
  
  # Load model
  mod_rng = load_model("glm_multi_beta_binomial_generate_data", threads = cores, cache_dir = cache_stan_model)

  draws_matrix <- attr(.data, "fit")$draws(format = "matrix")
  
  if(number_of_draws > nrow(draws_matrix)) {
    number_of_draws <- nrow(draws_matrix)
  }
  
  draws_matrix <- draws_matrix[sample(seq_len(nrow(draws_matrix)), size=number_of_draws),, drop=FALSE]

  # Generate quantities
  mod_rng |> sample_safe(
    generate_quantities_fx,
    draws_matrix,
    data = model_input,
    seed = mcmc_seed, 
    threads_per_chain = 1
  )
}

#' @keywords internal
#' @noRd
#'
parse_generated_quantities = function(rng, number_of_draws = 1){
  
  # Define the variables as NULL to avoid CRAN NOTES
  .draw <- NULL
  N <- NULL
  .value <- NULL
  generated_counts <- NULL
  M <- NULL
  generated_proportions <- NULL
  
  draws_to_tibble_x_y(rng, "counts", "N", "M", number_of_draws) %>%
    with_groups(c(.draw, N), ~ .x %>% mutate(generated_proportions = .value/max(1, sum(.value)))) %>%
    filter(.draw<= number_of_draws) %>%
    rename(generated_counts = .value, replicate = .draw) %>%
    
    mutate(generated_counts = as.integer(generated_counts)) %>%
    select(M, N, generated_proportions, generated_counts, replicate)
  
}

#' harmonise_factor_levels
#'
#' @description
#' A helper function to make sure that factor levels in new data match the old data
#'
#' @param new_data A data frame containing potential factor variables
#' @param old_data A data frame containing the reference factor levels
#'
#' @return A data frame with the same dimensions as `new_data`
#' @noRd
#'
harmonise_factor_levels <- function(new_data, old_data) {
  if (!is.data.frame(new_data) || !is.data.frame(old_data)) {
    stop("Both new_data and old_data must be data frames")
  }
  
  # Ensure arguments are in the correct order
  if (!identical(names(formals(harmonise_factor_levels))[1:2], c("new_data", "old_data"))) {
    warning("Arguments appear to be in the wrong order. 'new_data' should be the data to be harmonized, 'old_data' is the reference.")
  }
  
  # Original implementation follows
  f_old <- sapply(old_data, is.factor)
  f_new <- sapply(new_data, is.factor)
  
  # Get the common factor column names
  common_f <- intersect(names(old_data)[f_old], names(new_data)[f_new])
  
  # If there are no common factor columns, return the new_data as is
  if (length(common_f) == 0) {
    return(new_data)
  }
  
  # For each common factor column
  for (col in common_f) {
    # Get all unique levels from both datasets
    all_levels <- unique(c(levels(old_data[[col]]), levels(new_data[[col]])))
    
    # Check if there are new levels in new_data that don't exist in old_data
    new_levels <- setdiff(levels(new_data[[col]]), levels(old_data[[col]]))
    if (length(new_levels) > 0) {
      message(paste("sccomp says: New levels found in column", col, ":", paste(new_levels, collapse = ", ")))
    }
    
    # Set levels of new_data to match old_data, adding any missing levels
    new_data[[col]] <- factor(new_data[[col]], levels = levels(old_data[[col]]))
  }
  
  # Also check for character columns in new_data that are factors in old_data
  char_new <- sapply(new_data, is.character)
  factor_cols_in_old <- names(old_data)[f_old]
  
  char_to_factor <- intersect(names(new_data)[char_new], factor_cols_in_old)
  
  for (col in char_to_factor) {
    # Convert character to factor with the same levels as in old_data
    new_data[[col]] <- factor(new_data[[col]], levels = levels(old_data[[col]]))
    
    # Check if there are values in new_data that don't exist in old_data's levels
    invalid_values <- setdiff(unique(as.character(new_data[[col]])), levels(old_data[[col]]))
    if (length(invalid_values) > 0 && !all(is.na(invalid_values))) {
      warning(paste("sccomp says: Values in column", col, "not found in reference levels:", paste(invalid_values[!is.na(invalid_values)], collapse = ", ")))
    }
  }
  
  return(new_data)
}

#' Get Output Samples from a Stan Fit Object
#'
#' This function retrieves the number of output samples from a Stan fit object, 
#' supporting different methods (MHC and Variational) based on the available data within the object.
#'
#' @param fit A `stanfit` object, which is the result of fitting a model via Stan.
#' @return The number of output samples used in the Stan model. 
#'         Returns from MHC if available, otherwise from Variational inference.
#' @examples
#' # Assuming 'fit' is a stanfit object obtained from running a Stan model
#' print("samples_count = get_output_samples(fit)")
#'
#' @noRd
#' 
get_output_samples = function(fit){
  
  # Check if the output_samples field is present in the metadata of the fit object
  # This is generally available when the model is fit using MHC (Markov chain Monte Carlo)
  if(!is.null(fit$metadata()$output_samples)) {
    # Return the output_samples from the metadata
    fit$metadata()$output_samples
  }
  
  # If the output_samples field is not present, check for iter_sampling
  # This occurs typically when the model is fit using Variational inference methods
  else if(!is.null(fit$metadata()$iter_sampling)) {
    # Return the iter_sampling from the metadata
    fit$metadata()$iter_sampling
  }
  else
    fit$metadata()$num_psis_draws
}