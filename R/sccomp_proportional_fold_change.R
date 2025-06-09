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

#' @export
#' 
#' @importFrom glue glue
#' 
sccomp_proportional_fold_change.sccomp_tbl = function(.data, formula_composition, from, to){
  
  my_factor = parse_formula(formula_composition)
  
  # Predict the composition for the specified conditions
  .data |> 
    sccomp_predict(
      formula_composition = formula_composition, 
      new_data = 
        tibble(sample=as.character(c(to, from)), factor = c(to, from)) |> 
        rename(!!my_factor := `factor`)
    ) |> 
    
    # Nest the predicted data by cell group
    nest(data = -!!.data |> attr(".cell_group")) |> 
    
    
    # Calculate the ratio of proportions between 'to' and 'from' conditions
    mutate(
      ratio_mean = map_dbl(
        data, 
        ~ {
          x = .x |> arrange(sample != !!from) |> pull(proportion_mean); 
          x[2]/x[1] })
    ) |> 
    mutate(
      proportion_from = map_dbl(data, ~.x |> filter(sample==from) |> pull(proportion_mean)),
      proportion_to = map_dbl(data, ~.x |> filter(sample!=from) |> pull(proportion_mean))
    ) |> 
    
    # Calculate the proportional fold change
    mutate(proportion_fold_change = if_else(ratio_mean<1, (-1/ratio_mean) , ratio_mean)) |> 
    
    # Calculate the ratio of credible interval between 'to' and 'from' conditions
    mutate(
      ratio_upper = map_dbl(
        data,  
        ~ {
          x = .x |> arrange(sample != !!from) |> pull(proportion_upper); 
          x[2]/x[1] }),
      ratio_lower = map_dbl(
        data, 
        ~ {
          x = .x |> arrange(sample != !!from) |> pull(proportion_lower); 
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


