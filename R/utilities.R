# Define global variable
# sccomp_stan_models_cache_dir = file.path(path.expand("~"), ".sccomp_models", packageVersion("sccomp"))



# Define global variable (without version - version will be added in load_model when needed)
#' Default cache directory for Stan models
#' 
#' A global variable that defines the default cache directory for Stan models used by the sccomp package.
#' This directory is used to store compiled Stan models to avoid recompilation on subsequent runs.
#' 
#' @details
#' The cache directory is set to \code{~/.sccomp_models} by default. This location is used by
#' various sccomp functions to store and retrieve compiled Stan models, improving performance
#' by avoiding unnecessary recompilation of models that have already been compiled.
#' 
#' Users can override this default by specifying a different cache directory in function calls
#' that accept a \code{cache_stan_model} parameter.
#' 
#' @return A character string containing the path to the default cache directory.
#' 
#' @examples
#' # View the default cache directory
#' sccomp_stan_models_cache_dir
#' 
#' # Use a custom cache directory in a function call
#' # sccomp_estimate(data, cache_stan_model = "/path/to/custom/cache")
#' 
#' @seealso
#' \code{\link{sccomp_estimate}}, \code{\link{sccomp_replicate}}, \code{\link{clear_stan_model_cache}}
#' 
#' @export
sccomp_stan_models_cache_dir = file.path(path.expand("~"), ".sccomp_models")

#' Add attribute to abject
#'
#' @keywords internal
#' @noRd
#'
#'
#' @param var A tibble
#' @param attribute An object
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_attr = function(var, attribute, name) {
  attr(var, name) <- attribute
  var
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom magrittr extract2
#' @importFrom stats terms
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula <- function(fm) {
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  as.character(attr(terms(fm), "variables")) |>
    str_subset("\\|", negate = TRUE) %>%
    
    # Does not work the following
    # |>
    # extract2(-1)
    .[-1]
}


#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#' @importFrom stringr str_trim
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
formula_to_random_effect_formulae <- function(fm) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  random_effect_elements =
    as.character(attr(terms(fm), "variables")) |>
    
    # Select random intercept part
    str_subset("\\|")
  
  if(length(random_effect_elements) > 0){
    
    random_effect_elements |>
      
      # Divide grouping from factors
      str_split("\\|") |>
      
      # Set name
      map_dfr(~ .x |> set_names(c("formula", "grouping"))) |>
      
      # Create formula
      mutate(formula = map(formula, ~ formula(glue("~ {.x}")))) |>
      mutate(grouping = grouping |> str_trim())
    
  }
  
  else
    tibble(`formula` = list(), grouping = character())
  
}

#' Formula parser
#'
#' @param fm A formula
#'
#' @importFrom stringr str_subset
#' @importFrom stringr str_split
#' @importFrom stringr str_remove_all
#' @importFrom rlang set_names
#' @importFrom purrr map_dfr
#'
#' @importFrom magrittr extract2
#'
#' @return A character vector
#'
#' @keywords internal
#' @noRd
parse_formula_random_effect <- function(fm) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  formula <- NULL
  
  stopifnot("The formula must be of the kind \"~ factors\" " = attr(terms(fm), "response") == 0)
  
  random_effect_elements =
    as.character(attr(terms(fm), "variables")) |>
    
    # Select random intercept part
    str_subset("\\|")
  
  if(length(random_effect_elements) > 0){
    
    formula_to_random_effect_formulae(fm) |>
      
      # Divide factors
      mutate(factor = map(
        formula,
        ~
          # Attach intercept
          .x |>
          terms() |>
          attr("intercept") |>
          str_replace("^1$", "(Intercept)") |>
          str_subset("0", negate = TRUE) |>
          
          # Attach variables
          c(
            .x |>
              terms() |>
              attr("variables") |>
              as.character() |>
              str_split("\\+") |>
              as.character() %>%
              .[-1]
          )
      )) |>
      unnest(factor)
    
  }
  
  else
    tibble(factor = character(), grouping = character())
  
}

#' Get matrix from tibble
#'
#' @import dplyr
#'
#' @keywords internal
#' @noRd
#'
#' @import dplyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#' @return A tibble
# REMOVED: ifelse_pipe function, conditional_apply function, and as_matrix function




#' draws_to_tibble_x_y
#'
#' @importFrom tidyr pivot_longer
#' @importFrom rlang :=
#'
#' @param fit A fit object
#' @param par A character vector. The parameters to extract.
#' @param x A character. The first index.
#' @param y A character. The first index.
#'
#' @keywords internal
#' @noRd
draws_to_tibble_x_y = function(fit, par, x, y, number_of_draws = NULL) {
  
  # Define the variables as NULL to avoid CRAN NOTES
  dummy <- NULL
  .variable <- NULL
  .chain <- NULL
  .iteration <- NULL
  .draw <- NULL
  .value <- NULL

  
  fit$draws(variables = par, format = "draws_df") %>%
    mutate(.iteration = seq_len(n())) %>%
    
    pivot_longer(
      names_to = "parameter", # c( ".chain", ".variable", x, y),
      cols = contains(par),
      #names_sep = "\\.?|\\[|,|\\]|:",
      # names_ptypes = list(
      #   ".variable" = character()),
      values_to = ".value"
    ) %>%
    tidyr::extract(parameter, c(".chain", ".variable", x, y), "([1-9]+)?\\.?([a-zA-Z0-9_\\.]+)\\[([0-9]+),([0-9]+)") |> 
    
    # Warning message:
    # Expected 5 pieces. Additional pieces discarded
    suppressWarnings() %>%
    
    mutate(
      !!as.symbol(x) := as.integer(!!as.symbol(x)),
      !!as.symbol(y) := as.integer(!!as.symbol(y))
    ) %>%
    arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
    group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
    mutate(.draw = seq_len(n())) %>%
    ungroup() %>%
    select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value) %>%
    filter(.variable == par)
  
}


#' @importFrom tidyr separate
#' @importFrom purrr when
#' 
#' @param fit A fit object from a statistical model from 'cmdstanr' package.
#' @param par A character vector specifying the parameters to extract from the fit object.
#' @param x A character string specifying the first index variable name in the parameter names.
#' @param y A character string specifying the second index variable name in the parameter names (optional). If NULL, the function will attempt to parse both indices from the parameter names.
#' @param probs A numerical vector specifying the quantiles to extract. Default is c(0.025, 0.25, 0.50, 0.75, 0.975).
#' @param robust A logical value indicating whether to use robust statistics (median and MAD) instead of mean and standard deviation. Default is FALSE.
#' 
#' @return A tibble containing the summary statistics for the specified parameters, with separate columns for the parsed indices.
#' 
#' @details This function extracts summary statistics from a fitted model and parses parameter names 
#' that contain indices (e.g., "beta[1,2]") into separate columns. It supports both classical 
#' (mean/sd) and robust (median/mad) summary statistics.
#' 
#' @keywords internal
#' @noRd
summary_to_tibble = function(fit, par, x, y = NULL, probs = c(0.025, 0.25, 0.50, 0.75, 0.975), robust = FALSE) {
  
  # Avoid bug
  #if(fit@stan_args[[1]]$method %>% is.null) fit@stan_args[[1]]$method = "hmc"
  
  # Choose between robust and classical summary statistics
  if(!robust)
    # Use classical statistics: mean, standard deviation, and quantiles
    summary = 
      fit$summary(variables = par, "mean", "sd", ~ quantile(.x, probs = probs,  na.rm=TRUE), "rhat", "ess_bulk", "ess_tail") %>%
      rename(.variable = variable ) 
  else
    # Use robust statistics: median, median absolute deviation, and quantiles
    summary = 
      fit$summary(variables = par, "median", "mad", ~ quantile(.x, probs = probs,  na.rm=TRUE), "rhat", "ess_bulk", "ess_tail") %>%
      rename(.variable = variable ) 
  
  # Parse parameter names that contain indices (e.g., "beta[1,2]")
  if(is.null(y))
    # If y is not provided, extract both x and y indices from parameter names
    summary |> 
      tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop")
  else
    # If y is provided, still extract both indices (y parameter is ignored in this case)
    # This maintains consistency in the output format regardless of whether y is provided
    summary |> 
      tidyr::separate(col = .variable,  into = c(".variable", x, y), sep="\\[|,|\\]", convert = TRUE, extra="drop")

  
}


#' Get Random Intercept Design 3
#'
#' This function processes the formula composition elements in the data and creates design matrices
#' for random intercept models.
#'
#' @param .data_ A data frame containing the data.
#' @param .sample A quosure representing the sample variable.
#' @param formula_composition A data frame containing the formula composition elements.
#' 
#' @return A data frame with the processed design matrices for random intercept models.
#' 
#' @importFrom glue glue
#' @importFrom magrittr subtract
#' @importFrom purrr map2
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr pull
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_all
#' @importFrom dplyr mutate_if
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom tidyselect all_of
#' @importFrom readr type_convert
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr not
#' @importFrom tibble deframe
#' @importFrom tibble column_to_rownames
#' @noRd
get_random_effect_design3 = function(
  .data_, formula, grouping, .sample, 
  accept_NA_as_average_effect = FALSE 
){
  
  # Define the variables as NULL to avoid CRAN NOTES
  .sample = enquo(.sample)
  
  mydesign = .data_ |> get_design_matrix(formula, !!.sample, accept_NA_as_average_effect = accept_NA_as_average_effect)
  
  # Create a matrix of group assignments
  group_matrix = .data_ |> 
    select(all_of(grouping)) |> 
    pull(1) |> 
    rep(ncol(mydesign)) |> 
    matrix(ncol = ncol(mydesign))
  
  # Create a mask where design matrix has non-zero values
  mask = mydesign != 0
  
  # Apply mask to group matrix (setting non-masked values to "")
  group_matrix[!mask] = ""
  
  # Handle NAs in group matrix
  if(accept_NA_as_average_effect) {
    # For rows with NA in grouping, set group to "NA"
    group_matrix[is.na(group_matrix)] = "NA"
  }
  
  colnames(group_matrix) = colnames(mydesign)
  rownames(group_matrix) = rownames(mydesign)
  
  # Create the result matrix
  result = matrix(0, nrow = nrow(mydesign), ncol = ncol(mydesign))
  rownames(result) = rownames(mydesign)
  colnames(result) = colnames(mydesign)
  
  # For each column in the design matrix
  for(col in seq_len(ncol(mydesign))) {
    # Get the non-zero values and their corresponding groups
    non_zero = mydesign[, col] != 0
    if(any(non_zero)) {
      # Set the values in the result matrix
      result[non_zero, col] = mydesign[non_zero, col]
    }
  }
  
  # Convert to long format with the correct column names
  result_long = result |>
    as.data.frame() |>
    rownames_to_column(quo_name(.sample)) |>
    pivot_longer(-!!.sample, names_to = "factor", values_to = "value") |>
    filter(value != 0) |>
    mutate(
      group___label = paste0(factor, "___", group_matrix[cbind(!!.sample, factor)]),
      group___numeric = as.integer(factor(group___label)),
      factor___numeric = as.integer(factor(factor))
    ) |>
    select(!!.sample, group___label, group___numeric, value)
  
  result_long
}

#' @importFrom glue glue
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr pull
#' @importFrom dplyr where
#' @importFrom rlang enquo
#' @noRd
get_design_matrix = function(.data_spread, formula, .sample, accept_NA_as_average_effect = FALSE){
  
  .sample = enquo(.sample)
  
  .data_spread = .data_spread %>%
    
    select(!!.sample, parse_formula(formula)) |>
    mutate(across(where(is.numeric),  scale)) 

  # Check for NAs in the data
  has_na = any(is.na(.data_spread |> select(parse_formula(formula))))
  
  # If NAs are present and we don't accept them as average effects, throw an error
  if(has_na && !accept_NA_as_average_effect){
    stop("sccomp says: NA values are present in the design matrix factors. Set accept_NA_as_average_effect = TRUE to handle NAs as average effects across factor levels.")
  }

  # Check if we should handle NAs as average effects
  if(accept_NA_as_average_effect && has_na){
    return(get_design_matrix_with_na_handling(.data_spread, formula, !!.sample))
  }
  
  design_matrix =
    .data_spread |>
    model.matrix(formula, data=_)
  
  rownames(design_matrix) = .data_spread |> pull(!!.sample)
  
  design_matrix
}

#' @description
#' This function handles NA values in a special way for design matrices.
#' When NA values are present in categorical variables (factors or characters), they are treated as a separate
#' level "NA" initially, but then this level is removed from the design matrix. For rows with NA values,
#' the function sets equal weights (0.5) across all remaining factor levels. This approach allows for modeling
#' the average effect across all levels of a factor when the specific level is unknown or missing, rather than
#' excluding these observations or treating NA as its own distinct category. This is particularly useful when
#' you want to include observations with missing factor values but don't want to bias the model toward any
#' particular level of that factor.
#' 
#' @param .data_spread A data frame containing the data.
#' @param formula The formula to use for creating the design matrix.
#' @param .sample A quosure representing the sample variable.
#' 
#' @return A design matrix with rows named by the sample variable.
#' 
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr pull
#' @importFrom dplyr where
#' @importFrom rlang enquo
#' @noRd
get_design_matrix_with_na_handling = function(.data_spread, formula, .sample){
  
  .sample = enquo(.sample)
  
  # Convert NA to a factor level "NA" for categorical variables
  data_with_na <- handle_missing_values(.data_spread)
  
  # Create design matrix
  design_matrix = model.matrix(formula, data = data_with_na, na.action = NULL)
  
  # Here I get
 # design_matrix |> tail()
 #  (Intercept) group2__GROUP22 typehealthy typeNA group2__GROUP22:typehealthy group2__GROUP22:typeNA
 # 16           1               0           0      0                           0                      0
 # 17           1               1           0      0                           0                      0
 # 18           1               0           0      0                           0                      0
 # 19           1               0           1      0                           0                      0
 # 20           1               0           1      0                           0                      0
 # 21           1               0           0      1                           0                      0

# Split the design matrix into two parts:
# 1. The part that is not affected by NA values
# 2. The columns that are affected by NA values


na_rows = which(rowSums(design_matrix[, grep("NA$|NA:", colnames(design_matrix)), drop=FALSE]) > 0) |> unique()

# Loop over na_rows and calculate the fraction contribution for each row
for(row in na_rows) { 

  # Select the columns that are affected by NA values and also have a non-zero value in the design matrix
  na_cols = grep("NA$|NA:", colnames(design_matrix[row, , drop = FALSE]))
  na_cols = na_cols[design_matrix[row, na_cols] != 0]
  
  if(length(na_cols) == 0) next
  
  my_design_matrix = design_matrix[row,, drop = FALSE]
 
  fraction_contribution = calculate_na_fraction_contribution(my_design_matrix, na_cols, design_matrix, data_with_na)

  # Now for every row in fraction_contribution, I want to add the fraction contribution to the design matrix according to the factor_name
  for(fraction_contribution_row in 1:nrow(fraction_contribution)) {
    my_col = which(colnames(design_matrix) == (fraction_contribution[fraction_contribution_row,,drop=FALSE] |> select(all_of("factor_name")) |> pull(1)))
    design_matrix[row, my_col] = design_matrix[row, my_col] + (fraction_contribution[fraction_contribution_row,,drop=FALSE] |> select(all_of("fraction_contribution")) |> pull(1))
  }
}

rownames(design_matrix) = .data_spread |> pull(!!.sample)
design_matrix[,grep("NA$|NA:", colnames(design_matrix), invert = TRUE), drop=FALSE]

}

#' Calculate fraction contribution for NA values in design matrix
#'
#' This function calculates how to distribute the effect of NA values across all possible factor levels.
#' For each column in the design matrix that contains NA values (either directly or in interactions),
#' it:
#' 1. Splits the column name into its components
#' 2. For each component containing "NA", replaces it with all possible factor levels
#' 3. Creates all possible combinations of these levels
#' 4. Calculates the fraction of the effect that should be attributed to each combination
#' 5. Filters to only keep combinations that exist in the original design matrix
#'
#' @param design_with_na A matrix containing columns with NA values
#' @param design_matrix The original design matrix
#' @param data_with_na The data frame containing the original data with NA values handled
#'
#' @importFrom stringr str_remove
#' @importFrom stringr str_split
#' @importFrom stringr str_replace_all
#' @importFrom tidyr crossing
#' @importFrom tibble enframe
#' @importFrom dplyr rowwise
#' @importFrom dplyr c_across
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr matches
#' @importFrom dplyr all_of
#' @importFrom purrr map
#' @importFrom purrr reduce
#' 
#' @return A list of tibbles, each containing the factor combinations and their fraction contributions
#' @noRd
calculate_na_fraction_contribution = function(my_design_matrix, na_cols, design_matrix, data_with_na) {
  
  design_with_na = my_design_matrix[, na_cols, drop = FALSE]
  
  map(colnames(design_with_na), function(col) {
    
    list_of_contributions = 
      str_split(col, ":", simplify = TRUE) |>
      map(function(x) {

        factor_name = x[1] |> str_remove("NA$")
        
        # Now I want to count the possible values of the factor in the input data
        if(data_with_na[[factor_name]] |> is.numeric())
          factor_levels = factor_name
        else{
          factor_levels = levels(factor(data_with_na[[factor_name]]))
          factor_levels = factor_levels[factor_levels != "NA"]
          if(factor_levels |> length() == 0) factor_levels = ""
        }
       
        
        str_replace_all(x, "NA", factor_levels) |> 
          enframe(value = "factor_name") |> 
          select(-name) |> 
          mutate(fraction_contribution = 1/n()) |> 
          mutate(fraction_contribution = fraction_contribution * my_design_matrix[,x ])
      }) 
    
    # Parse contribution
    list_of_contributions |>
      
      purrr::reduce(crossing,  .name_repair = "unique")  |> 
      suppressMessages() |> 
      rowwise() |> 
      mutate(
        fraction_contribution = prod( c_across(matches("fraction_contribution")) ),
        factor_name = paste0( c_across(matches("factor_name")), collapse=":" )
      ) |> 
      
      filter(factor_name %in% colnames(design_matrix)) |> 
      select(factor_name, fraction_contribution)
  }) |>
    
    bind_rows()
}



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
    .data_spread, formula, .sample, .cell_group, .count,
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
    .cell_group = enquo(.cell_group)
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
        column_to_rownames(quo_name(.sample))
      
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
          column_to_rownames(quo_name(.sample))
        
        # Separate NA group column into X_random_effect_2_unseen  
        X_random_effect_2_unseen = X_random_effect_2[, colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
        X_random_effect_2 = X_random_effect_2[, !colnames(X_random_effect_2) |> str_detect("___NA$"), drop = FALSE]
      }
      
      else X_random_effect_2 =  X_random_effect[,0,drop=FALSE]
      
      n_random_eff = random_effect_grouping |> nrow()
      
      ncol_X_random_eff = c(ncol(X_random_effect), ncol(X_random_effect_2))
      
      # TEMPORARY
      group_factor_indexes_for_covariance = 
        X_random_effect |> 
        colnames() |> 
        enframe(value = "parameter", name = "order")  |> 
        separate(parameter, c("factor", "group"), "___", remove = FALSE) |> 
        complete(factor, group, fill = list(order=0)) |> 
        select(-parameter) |> 
        pivot_wider(names_from = group, values_from = order)  |> 
        column_to_rownames("factor") |> as.matrix()

      
      
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
        column_to_rownames("factor") |> as.matrix()
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
    
    
    y = .data_spread %>% select(-any_of(factor_names), -exposure, -!!.grouping_for_random_effect) %>% column_to_rownames(quo_name(.sample)) %>% as.matrix()
    
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



contrasts_to_parameter_list = function(contrasts, drop_back_quotes = TRUE){
  
  if(is.null(names(contrasts)))
    names(contrasts) <- contrasts
  
  contrast_list <-
    contrasts |>
    
    # Remove fractions used in multiplication before or after '*'
    str_remove_all_ignoring_if_inside_backquotes("([0-9]+/[0-9]+ ?\\* ?)|(\\* ?[0-9]+/[0-9]+)") |>
    
    # Remove decimals used in multiplication before or after '*'
    str_remove_all_ignoring_if_inside_backquotes("([-+]?[0-9]+\\.[0-9]+ ?\\* ?)|(\\* ?[-+]?[0-9]+\\.[0-9]+)") |>
    
    # Remove fractions used in divisions
    str_remove_all_ignoring_if_inside_backquotes("/ ?[0-9]+") |>
    
    # Remove standalone numerical constants not inside backquotes
    str_remove_all_ignoring_if_inside_backquotes("\\b[0-9]+\\b") |>
    
    # Split by "+", "-", "*", "/"
    str_split_ignoring_if_inside_backquotes("\\+|-|\\*|/") |>
    unlist() |>
    
    # Remove parentheses and spaces
    str_remove_all_ignoring_if_inside_backquotes("[\\(\\) ]") |>
    
    # Remove empty strings
    {\(x) x[x != ""]}()
  
  if(drop_back_quotes)
    contrast_list <-
    contrast_list |>
    str_remove_all("`")
  
  contrast_list |> unique()
}








#' chatGPT - Remove Specified Regex Pattern from Each String in a Vector
#'
#' This function takes a vector of strings and a regular expression pattern.
#' It removes occurrences of the pattern from each string, except where the pattern
#' is found inside backticks. The function returns a vector of cleaned strings.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern to be removed
#' from the text.
#'
#' @return A character vector with the regex pattern removed from each string.
#' Occurrences of the pattern inside backticks are not removed.
#'
#' @examples
#' texts <- c("A string with (some) parentheses and `a (parenthesis) inside` backticks",
#'            "Another string with (extra) parentheses")
#' cleaned_texts <- str_remove_all_ignoring_if_inside_backquotes(texts, "\\(")
#' print(cleaned_texts)
#' 
#' @noRd
str_remove_all_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex removal for a single string
  remove_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- ""
    skip <- 0
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      if (skip > 0) {
        skip <- skip - 1
        next
      }
      
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        result <- paste0(result, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          next
        } else {
          result <- paste0(result, char)
        }
      } else {
        result <- paste0(result, char)
      }
    }
    
    return(result)
  }
  
  # Apply the function to each element in the vector
  sapply(text_vector, remove_regex_chars, regex)
}


#' chatGPT - Split Each String in a Vector by a Specified Regex Pattern
#'
#' This function takes a vector of strings and a regular expression pattern. It splits
#' each string based on the pattern, except where the pattern is found inside backticks.
#' The function returns a list, with each element being a vector of the split segments
#' of the corresponding input string.
#'
#' @param text_vector A character vector with the strings to be processed.
#' @param regex A character string containing a regular expression pattern used for splitting
#' the text.
#'
#' @return A list of character vectors. Each list element corresponds to an input string
#' from `text_vector`, split according to `regex`, excluding occurrences inside backticks.
#'
#' @examples
#' texts <- c("A string with, some, commas, and `a, comma, inside` backticks",
#'            "Another string, with, commas")
#' split_texts <- split_regex_chars_from_vector(texts, ",")
#' print(split_texts)
#' 
#' @noRd
str_split_ignoring_if_inside_backquotes <- function(text_vector, regex) {
  # Nested function to handle regex split for a single string
  split_regex_chars <- function(text, regex) {
    inside_backticks <- FALSE
    result <- c()
    current_segment <- ""
    
    chars <- strsplit(text, "")[[1]]
    for (i in seq_along(chars)) {
      char <- chars[i]
      if (char == "`") {
        inside_backticks <- !inside_backticks
        current_segment <- paste0(current_segment, char)
      } else if (!inside_backticks) {
        # Check the remaining text against the regex
        remaining_text <- paste(chars[i:length(chars)], collapse = "")
        match <- regexpr(regex, remaining_text)
        
        if (attr(match, "match.length") > 0 && match[1] == 1) {
          # Add current segment to result and start a new segment
          result <- c(result, current_segment)
          current_segment <- ""
          # Skip the length of the matched text
          skip <- attr(match, "match.length") - 1
          i <- i + skip
        } else {
          current_segment <- paste0(current_segment, char)
        }
      } else {
        current_segment <- paste0(current_segment, char)
      }
    }
    
    # Add the last segment to the result
    result <- c(result, current_segment)
    return(result)
  }
  
  # Apply the function to each element in the vector
  lapply(text_vector, split_regex_chars, regex)
}

#' Convert array of quosure (e.g. c(col_a, col_b)) into character vector
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom rlang quo_name
#' @importFrom rlang quo_squash
#'
#' @param v A array of quosures (e.g. c(col_a, col_b))
#'
#' @return A character vector
quo_names <- function(v) {
  
  v = quo_name(quo_squash(v))
  gsub('^c\\(|`|\\)$', '', v) |>
    strsplit(', ') |>
    unlist()
}

#' Add class to abject
#'
#' @keywords internal
#' @noRd
#'
#' @param var A tibble
#' @param name A character name of the attribute
#'
#' @return A tibble with an additional attribute
add_class = function(var, name) {
  
  if(!name %in% class(var)) class(var) <- append(class(var),name, after = 0)
  
  var
}

drop_environment <- function(obj) {
  # Check if the object has an environment
  if (!is.null(environment(obj))) {
    environment(obj) <- new.env()
  }
  return(obj)
}


#' Print Tibble in Red
#'
#' This function captures the console output of printing a tibble,
#' colours it in red and returns the coloured text.
#'
#' @param tbl A data frame or tibble to be printed and coloured in red.
#'
#' @return A character string containing the coloured tibble output.
#'
#' @importFrom crayon red
#' @noRd
print_red_tibble <- function(tbl) {
  # Capture the console output of printing the tibble
  example_text <- capture.output(print(tbl))
  
  # Combine all lines into one block and colour it in red
  red(paste(example_text, collapse = "\n"))
}


# This is needed because I have a `sample` argument, that when it is not defined affects the sample() function
sample_seed = function(){
  sample(1e5, 1)
}

#' Handle missing values in a data frame
#' 
#' @param data A data frame to process
#' @return A data frame with missing values handled:
#'   - For categorical variables (factor/character): NA values are converted to a factor level "NA"
#'   - For numeric variables: NA values are replaced with the mean of the column
handle_missing_values <- function(data) {
  data %>%
    mutate(
      # 1) Handle factor/character: turn NA into a real "NA" level, last in the ordering
      across(
        where(~ is.factor(.x) || is.character(.x)),
        ~{
          if (any(is.na(.x))) {
            x_char <- as.character(.x)
            x_char[is.na(x_char)] <- "NA"
            f <- factor(x_char)
            non_na_levels <- setdiff(levels(f), "NA")
            new_levels <- c(non_na_levels, "NA")
            factor(f, levels = new_levels)
          } else {
            .x
          }
        }
      ),
      # 2) Handle numeric: replace NA by the (non-NA) mean of that column
      across(
        where(is.numeric),
        ~{
          if (any(is.na(.x))) {
            mean_val <- mean(.x, na.rm = TRUE) 
            ifelse(is.na(.x), mean_val, .x)  # Only replace NAs, keep original values
          } else {
            .x
          }
        }
      )
    )
}




