#' Calculate Proportional Fold Change from sccomp Estimated Effects
#'
#' @description 
#' This function calculates the proportional fold change between two conditions using the estimated effects from a sccomp model.
#' The fold changes are derived from the model's posterior predictions rather than raw counts, providing a more robust estimate
#' that accounts for the model's uncertainty and covariate effects.
#' 
#' Note! This statistic is descriptive and should not be used to define significance - use sccomp_test() for that. 
#' While fold changes in proportions are easier to interpret than changes in logit space, they are not linear 
#' (the same proportional change has different meaning for rare vs abundant cell types). In contrast, 
#' the logit scale used internally by sccomp provides linear effects that are more appropriate for statistical inference.
#'
#' @param .data A sccomp estimate object (of class 'sccomp_tbl') obtained from running sccomp_estimate().
#'              This object contains the fitted model and estimated effects.
#' @param formula_composition The formula specifying which model effects to use for calculating fold changes.
#'                          This should match or be a subset of the formula used in the original sccomp_estimate() call.
#' @param from Character string specifying the reference/control condition (e.g., "benign").
#' @param to Character string specifying the comparison condition (e.g., "cancer").
#'
#' @return A tibble with the following columns:
#' \itemize{
#'   \item cell_group - The cell group identifier
#'   \item proportion_fold_change - The estimated fold change in proportions between conditions.
#'                                 Positive values indicate increases, negative values indicate decreases.
#'   \item average_uncertainty - The average uncertainty in the fold change estimate, derived from the credible intervals
#'   \item statement - A text description of the fold change, including the direction and the estimated proportions
#' }
#'
#' @examples
#' 
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#'
#' \donttest{
#' if (instantiate::stan_cmdstan_exists()) {
#'   # Load example data
#'   data("counts_obj")
#'
#'   # First estimate the composition effects
#'   estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'   )
#'  
#'   # Calculate proportional fold changes from the estimated effects
#'   estimate |> 
#'   sccomp_proportional_fold_change(
#'     formula_composition = ~  type, 
#'     from = "benign", 
#'     to = "cancer"
#'   ) 
#' }
#' }
#' @export
sccomp_proportional_fold_change <- function(.data, formula_composition, from, to) {
  UseMethod("sccomp_proportional_fold_change", .data)
}

#' Convert character columns to factors with proper levels
#'
#' @description
#' Internal helper function that converts character columns in new_data to factors
#' with the same levels as in the original training data. This is critical for
#' predictions to work correctly, as the model expects factors with specific levels.
#' 
#' Without this conversion, character values are treated as the reference level,
#' leading to incorrect predictions where all conditions produce the same output.
#'
#' @param new_data A data frame with character columns to convert
#' @param factor_columns Character vector of column names that should be converted to factors
#' @param original_data The original training data containing the factor levels
#'
#' @return The new_data with character columns converted to factors with proper levels
#' 
#' @keywords internal
#' @noRd
convert_to_factors_with_levels <- function(new_data, factor_columns, original_data) {
  
  for (factor_col in factor_columns) {
    
    # Extract factor levels from the original data
    if (is.factor(original_data[[factor_col]])) {
      factor_levels <- levels(original_data[[factor_col]])
    } else {
      # For character columns, use unique values as levels
      factor_levels <- unique(original_data[[factor_col]])
    }
    
    # Store original values before conversion to detect invalids
    original_values <- new_data[[factor_col]]
    
    # Convert new_data column to factor with the same levels
    new_data[[factor_col]] <- factor(new_data[[factor_col]], levels = factor_levels)
    
    # Check if the conversion created NAs (indicating invalid factor levels)
    na_mask <- is.na(new_data[[factor_col]])
    if (any(na_mask)) {
      invalid_values <- original_values[na_mask]
      stop(sprintf(
        "sccomp says: Invalid level(s) '%s' provided for factor '%s'. Valid levels are: %s",
        paste(invalid_values, collapse = "', '"),
        factor_col,
        paste(factor_levels, collapse = ", ")
      ))
    }
  }
  
  return(new_data)
}

#' Match interaction term parts to their corresponding factors
#'
#' @description
#' Internal helper function that matches values from interaction specifications
#' (e.g., "control:baseline") to their correct factor columns, regardless of the
#' order in which they are provided. This allows users to specify interactions
#' as either "treatment:timepoint" or "timepoint:treatment".
#'
#' @param parts Character vector of factor level values (e.g., c("control", "baseline"))
#' @param factors Character vector of factor column names (e.g., c("treatment", "timepoint"))
#' @param data The original data frame containing the factor columns
#'
#' @return A named list where names are factor column names and values are the matched levels
#' 
#' @keywords internal
#' @noRd
match_interaction_parts_to_factors <- function(parts, factors, data) {
  
  # Initialize result list
  matched <- list()
  
  # For each factor, find which part matches
  for (factor_name in factors) {
    # Get valid levels for this factor
    if (is.factor(data[[factor_name]])) {
      valid_levels <- levels(data[[factor_name]])
    } else {
      valid_levels <- unique(data[[factor_name]])
    }
    
    # Find which part matches this factor's valid levels
    match_idx <- which(parts %in% valid_levels)
    
    if (length(match_idx) == 0) {
      stop(sprintf(
        "sccomp says: None of the provided values '%s' match valid levels for factor '%s'. Valid levels are: %s",
        paste(parts, collapse = "', '"),
        factor_name,
        paste(valid_levels, collapse = ", ")
      ))
    } else if (length(match_idx) > 1) {
      stop(sprintf(
        "sccomp says: Multiple values '%s' match valid levels for factor '%s'. Each factor should have exactly one matching value.",
        paste(parts[match_idx], collapse = "', '"),
        factor_name
      ))
    }
    
    # Store the matched value
    matched[[factor_name]] <- parts[match_idx]
    
    # Remove this part from consideration for other factors
    parts <- parts[-match_idx]
  }
  
  # Check if there are leftover parts that didn't match any factor
  if (length(parts) > 0) {
    stop(sprintf(
      "sccomp says: The following values could not be matched to any factor: '%s'. Available factors are: %s",
      paste(parts, collapse = "', '"),
      paste(factors, collapse = ", ")
    ))
  }
  
  return(matched)
}

#' @export
#' 
#' @importFrom glue glue
#' @importFrom stringr str_split
#' 
sccomp_proportional_fold_change.sccomp_tbl = function(.data, formula_composition, from, to){
  
  my_factors = parse_formula(formula_composition)
  
  # Get the sample column name from the original data
  .sample = attr(.data, ".sample")
  
  # Handle interaction categories by parsing the from/to strings
  if (length(my_factors) > 1) {
    # For interactions, parse the category strings
    from_parts <- str_split(from, ":")[[1]]
    to_parts <- str_split(to, ":")[[1]]
    
    # # Validate that the number of parts matches the number of factors
    # if (length(from_parts) != length(my_factors) || length(to_parts) != length(my_factors)) {
    #   stop(sprintf(
    #     "sccomp says: Interaction term requires %d factor levels separated by ':', but got %d in 'from' and %d in 'to'",
    #     length(my_factors),
    #     length(from_parts),
    #     length(to_parts)
    #   ))
    # }
    
    # Get the original count data to determine which value belongs to which factor
    original_data <- attr(.data, "count_data")
    
    # Match each part to the correct factor by checking which factor contains that value
    from_matched <- match_interaction_parts_to_factors(from_parts, my_factors, original_data)
    to_matched <- match_interaction_parts_to_factors(to_parts, my_factors, original_data)
    
    # Create new_data with individual factor columns
    # Sample column needs dummy identifiers, not the interaction string
    new_data <- tibble(
      !!quo_name(.sample) := c("to", "from")
    )
    
    # Add each factor column with matched values
    for (i in seq_along(my_factors)) {
      new_data <- new_data %>%
        mutate(!!my_factors[i] := c(to_matched[[my_factors[i]]], from_matched[[my_factors[i]]]))
    }
  } else {   
    # For single factor, use the original approach
    # Sample column needs dummy identifiers, not the factor levels
    new_data <- tibble(
      !!quo_name(.sample) := c("to", "from"), 
      !!my_factors := c(to, from)
    )
  }
  
  # Convert character columns to factors with proper levels matching the original data
  # This is critical for the prediction to work correctly
  new_data <- convert_to_factors_with_levels(
    new_data, 
    my_factors, 
    attr(.data, "count_data")
  )
  
  # Predict the composition for the specified conditions
  predictions <- .data |> 
    sccomp_predict(
      formula_composition = formula_composition, 
      new_data = new_data
    )

  # Sample IDs for filtering
  sample_to <- "to"
  sample_from <- "from"
  
  predictions |> 
    # Nest the predicted data by cell group
    nest(data = -!!.data |> attr(".cell_group")) |> 
    
    
    # Calculate the ratio of proportions between 'to' and 'from' conditions
    mutate(
      ratio_mean = map_dbl(
        data, 
        ~ {
          x = .x |> arrange(!!.sample != sample_to) |> pull(proportion_mean); 
          x[2]/x[1] })
    ) |> 
    mutate(
      proportion_from = map_dbl(data, ~.x |> filter(!!.sample == sample_from) |> pull(proportion_mean)),
      proportion_to = map_dbl(data, ~.x |> filter(!!.sample == sample_to) |> pull(proportion_mean))
    ) |> 
    
    # Calculate the proportional fold change
    mutate(proportion_fold_change = if_else(ratio_mean<1, (-1/ratio_mean) , ratio_mean)) |> 
    
    # Calculate the ratio of credible interval between 'to' and 'from' conditions
    mutate(
      ratio_upper = map_dbl(
        data,  
        ~ {
          x = .x |> arrange(!!.sample != sample_to) |> pull(proportion_upper); 
          x[2]/x[1] }),
      ratio_lower = map_dbl(
        data, 
        ~ {
          x = .x |> arrange(!!.sample != sample_to) |> pull(proportion_lower); 
          x[2]/x[1] })
    ) |> 
    
    # Calculate the proportional fold change
    mutate(
      difference_proportion_upper_fold_change = if_else(ratio_mean<1, (-1/ratio_upper) , ratio_upper) - proportion_fold_change,
      difference_proportion_lower_fold_change = if_else(ratio_mean<1, (-1/ratio_lower) , ratio_lower) - proportion_fold_change
    ) |> 
    mutate(average_uncertainty = (abs(difference_proportion_upper_fold_change) + abs(difference_proportion_lower_fold_change))/2) |> 
    
    # Print expression
    mutate(increase_decrease = if_else(proportion_fold_change>0, "increase", "decrease")) |> 
    mutate(statement = glue("{round(abs(proportion_fold_change),1)}-fold {increase_decrease} (from {round(proportion_from, 4)} to {round(proportion_to, 4)})")) |> 
    
    # Select and return the relevant columns
    select(
      !!.data |> attr(".cell_group"), 
      proportion_fold_change, 
      average_uncertainty, 
      statement
    ) 
  
  
}

#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).


