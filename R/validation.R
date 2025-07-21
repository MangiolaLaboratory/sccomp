check_columns_exist = function(.data, ...){

  tryCatch({
    .data |> select(...)   # Try selecting the column
    TRUE  # If successful, the column exists
  }, error = function(e) {
    # Print a warning with additional context
    warning(
      "sccomp says: please check typos in your formula_composition, formula_variability (if applicable), \n",
      "and your .sample, .cell_group, .count (if applicable) arguments \n",
      e$message
    )
   
    .data |> select(...)
  })
  
  

}

check_if_count_integer = function(.data, .count){
  .count = enquo(.count)

  # Check if the counts column is an integer
  if (!.data %>% pull(!!.count) %>% is("integer"))
    stop(
      sprintf(
        "The column %s must be of class integer. You can do as mutate(`%s` = as.integer(`%s`))",
        quo_name(.count),
        quo_name(.count),
        quo_name(.count)
      )
    )
}

#' Check if NA
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr enquo
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | factors column
#' @param columns Columns to check
#'
#' @keywords internal
#' @noRd
check_if_any_NA = function(.data, ...){


  if(

    .data |>
    drop_na(...) |>
    nrow() < ( .data |> nrow() )

  )
    stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
}

check_if_within_posterior = function(.data, my_df, .do_check, .count){

  # Define the variables as NULL to avoid CRAN NOTES
  .lower <- NULL
  .upper <- NULL
  ppc <- NULL
  
  writeLines(sprintf("executing %s", "check_if_within_posterior"))

  .data %>%
    left_join(my_df, by = c("S", "G")) %>%
    filter((!!.do_check)) %>% # Filter only DE genes
    rowwise() %>%
    mutate(`ppc` = !!.count %>% between(`.lower`, `.upper`)) %>%
    mutate(`is higher than mean` = (!`ppc`) &
             (!!.count > mean)) %>%
    ungroup
}

check_if_columns_right_class = function(.data, .sample, .cell_group){

  .sample = enquo(.sample)
  .cell_group = enquo(.cell_group)

  # Check that sample and cell_type is chr of factor
  if(
    !(
      .data %>% pull(!!.sample) %>% is("character") |
      .data %>% pull(!!.sample) %>% is("factor")
    ) |
    !(
      .data %>% pull(!!.cell_group) %>% is("character") |
      .data %>% pull(!!.cell_group) %>% is("factor")
    )
  ) stop(sprintf("sccomp says: the columns %s and %s must be of class character or factor", quo_name(.sample), quo_name(.cell_group)))


}

#' Check if a Sample Column is a Unique Identifier
#'
#' This function checks if the `.sample` column in a wide dataset is truly
#' a unique identifier. If not, it throws an error containing the problematic
#' rows in red text.
#'
#' @param data_wide A data frame or tibble in wide format.
#' @param .sample   An unquoted column name indicating the sample column to check.
#'
#' @return Returns the original `data_wide` if `.sample` is unique. Otherwise,
#'   throws an error showing the problematic rows in red.
#'
#' @importFrom rlang enquo
#' @importFrom rlang quo_name
#' @importFrom dplyr count
#' @importFrom dplyr add_count
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom glue glue
#' @noRd
check_if_sample_is_a_unique_identifier <- function(data_wide, .sample) {
  .sample <- enquo(.sample)
  
  if (
    data_wide |>
    count(!!.sample) |>
    pull(n) |>
    max() > 1
  ) {
    stop(
      paste(
        glue("sccomp says: .sample column `{quo_name(.sample)}` should be a unique identifier, with a unique combination of factors. For example Sample_A cannot have both treated and untreated conditions in your input"),
        data_wide |>
          add_count(!!.sample, name = "n___") |>
          filter(n___ > 1) |>
          select(-n___) |>
          print_red_tibble(),
        sep = "\n\n"
      )
    )
  } else {
    return(data_wide)
  }
}

check_missing_parameters <- function(effects, model_effects) {
  # Find missing parameters
  
  missing_parameters <- 
    effects |> 
    setdiff(model_effects)
  
  # If there are any missing parameters, stop and show an error message
  if (length(missing_parameters) > 0) {
    stop(
      "sccomp says: Some of the parameters present in the data provided were not present when the model was fitted. For example:\n",
      paste(missing_parameters[1:min(3, length(missing_parameters))], collapse = "\n"),
      if (length(missing_parameters) > 3) "\n..."
    )
  }
  
}

#' Check Sample Consistency of Factors
#'
#' This function checks for each sample in the provided data frame if the number of unique
#' covariate values from a specified formula matches the number of samples. It is useful for
#' verifying data consistency before statistical analysis. The function stops and throws an
#' error if inconsistencies are found.
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr distinct
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr n_distinct
#' @importFrom dplyr all_of
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map_lgl
#'
#' @param .data A data frame containing the samples and covariates.
#' @param my_formula A formula object specifying the covariates to check (e.g., ~Timepoint*Response_binary+Patient_ID).
#' @param sample A character string specifying the sample column name.
#' @param cell_group A character string specifying the cell group column name.
#'
#' @details The function selects the sample and covariates based on `my_formula`, pivots
#' the data longer so each row represents a unique sample-covariate combination, nests
#' the data by covariate name, and checks if the number of unique sample-covariate
#' pairs matches the number of samples for each covariate.
#'
#' @return This function does not return a value; it stops with an error message if any
#' inconsistencies are found.
#'
#' @noRd
#' @keywords internal
check_sample_consistency_of_factors = function(.data, my_formula, sample, cell_group){
  
  # Convert character parameters to quos like the main methods
  .sample = sym(sample)
  .cell_group = sym(cell_group)
  
  # Check that I have one set of covariates per sample
  first_cell_group = .data |> pull(!!.cell_group) |> _[[1]]
  
  # If the formula is intercept only -> ~ 1 this test does not apply
  if(my_formula |> parse_formula() |> length() == 0)
    return(TRUE)
  
  # Get the formula variables
  formula_vars = parse_formula(my_formula)
  
  # Check for inconsistencies by finding samples with multiple values for the same factor
  # Convert all factor variables to character to avoid type conflicts in pivot_longer
  data_for_check = 
    .data |> 
    filter(!!.cell_group == first_cell_group) |> 
    select(!!.sample, all_of(formula_vars)) |>
    mutate(across(all_of(formula_vars), as.character))
  
  inconsistent_samples = 
    data_for_check |>
    pivot_longer(-!!.sample, names_to = "factor_name", values_to = "factor_value") |>
    group_by(!!.sample, factor_name) |>
    summarise(
      n_unique_values = n_distinct(factor_value),
      unique_values = paste(sort(unique(factor_value)), collapse = " vs "),
      .groups = "drop"
    ) |>
    filter(n_unique_values > 1) |>
    arrange(!!.sample, factor_name)
  
  if(inconsistent_samples |> nrow() > 0) {
    
    # Create a more informative error message
    error_msg = "sccomp says: your factors are mismatched across samples. Each sample should have consistent factor values across all cell groups.\n\n"
    error_msg = paste0(error_msg, "Problematic samples and their conflicting values:\n")
    
    # Group by sample for better readability
    for(sample_name in unique(inconsistent_samples |> pull(!!.sample))) {
      sample_issues = inconsistent_samples |> filter(!!.sample == sample_name)
      error_msg = paste0(error_msg, "\nSample: ", sample_name, "\n")
      for(i in seq_len(nrow(sample_issues))) {
        row = sample_issues[i, ]
        error_msg = paste0(error_msg, "  - ", row$factor_name, ": ", row$unique_values, "\n")
      }
    }
    
    error_msg = paste0(error_msg, "\nTo fix this issue, ensure each sample has consistent factor values across all cell groups.")
    
    stop(error_msg)
  }
  
}

#' chatGPT - Check for Valid Column Names in Tidyverse Context
#'
#' This function checks if each given column name in a vector contains only valid characters 
#' (letters, numbers, periods, and underscores) and does not start with a digit 
#' or an underscore, which are the conditions for a valid column name in `tidyverse`.
#'
#' @param column_names A character vector representing the column names to be checked.
#'
#' @return A logical vector: `TRUE` for each column name that contains only valid characters 
#' and does not start with a digit or an underscore; `FALSE` otherwise.
#'
#' @examples
#' contains_only_valid_chars_for_column(c("valid_column", "invalid column", "valid123", 
#' "123startWithNumber", "_startWithUnderscore"))
#'
#' @noRd
contains_only_valid_chars_for_column <- function(column_names) {
  # Function to check a single column name
  check_validity <- function(column_name) {
    # Regex pattern for valid characters (letters, numbers, periods, underscores)
    valid_char_pattern <- "[A-Za-z0-9._]"
    
    # Check if all characters in the string match the valid pattern
    all_chars_valid <- stringr::str_detect(column_name, paste0("^", valid_char_pattern, "+$"))
    
    # Check for leading digits or underscores
    starts_with_digit_or_underscore <- stringr::str_detect(column_name, "^[0-9_]")
    
    return(all_chars_valid && !starts_with_digit_or_underscore)
  }
  
  # Apply the check to each element of the vector
  sapply(column_names, check_validity)
}

#' function that check is there are NAs in the count column
#'
#' @param .data A tibble containing the data
#' @param .count Column containing count data
#'
#' @return A tibble with NAs in the count column
#' @keywords internal
#' @noRd
check_if_NAs_in_count = function(.data, .count){
  .count = enquo(.count)
  if(.data |> pull(!!.count) |> is.na() |> any())
    stop("sccomp says: the input data frame has NAs in the count column")
}