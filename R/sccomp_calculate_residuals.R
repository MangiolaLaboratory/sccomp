#' Calculate Residuals Between Observed and Predicted Proportions
#'
#' @description
#' `sccomp_calculate_residuals` computes the residuals between observed cell group proportions and the predicted proportions from a fitted `sccomp` model. This function is useful for assessing model fit and identifying cell groups or samples where the model may not adequately capture the observed data. The residuals are calculated as the difference between the observed proportions and the predicted mean proportions from the model.
#'
#' @param .data A tibble of class `sccomp_tbl`, which is the result of `sccomp_estimate()`. This tibble contains the fitted model and associated data necessary for calculating residuals.
#'
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item \strong{sample} - A character column representing the sample identifiers.
#'   \item \strong{cell_group} - A character column representing the cell group identifiers.
#'   \item \strong{residuals} - A numeric column representing the residuals, calculated as the difference between observed and predicted proportions.
#'   \item \strong{exposure} - A numeric column representing the total counts (sum of counts across cell groups) for each sample.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the predicted mean proportions for each cell group and sample using `sccomp_predict()`.
#'   \item Calculates the observed proportions from the original count data.
#'   \item Computes residuals by subtracting the predicted proportions from the observed proportions.
#'   \item Returns a tibble containing the sample, cell group, residuals, and exposure (total counts per sample).
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr distinct
#' @importFrom dplyr with_groups
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @importFrom magrittr %$%
#' @importFrom sccomp sccomp_predict
#'
#' @export
#'
#' @examples
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
#' # Load example data
#' data("counts_obj")
#'
#' # Fit the sccomp model
#' estimates <- sccomp_estimate(
#'   counts_obj,
#'   formula_composition = ~ type,
#'   formula_variability = ~1,
#'   sample = "sample",
#'   cell_group = "cell_group",
#'   abundance = "count",
#'   approximate_posterior_inference = "all",
#'   cores = 1
#' )
#'
#' # Calculate residuals
#' residuals <- sccomp_calculate_residuals(estimates)
#'
#' # View the residuals
#' print(residuals)
#' }}
#'
sccomp_calculate_residuals <- function(.data) {
  UseMethod("sccomp_calculate_residuals", .data)
}

#' @export
sccomp_calculate_residuals.sccomp_tbl = function(.data){
  
  
  
  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .grouping_for_random_intercept = attr(.data, ".grouping_for_random_intercept")
  .count = attr(.data, ".count")
  
  # Residuals
  residual = 
    .data |>
    sccomp_predict( 
      number_of_draws = 
        .data |>
        attr("fit") |>
        dim() |>
        _[1] |> 
        min(500) 
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) 
  
  if(.data |> attr("model_input") %$% is_proportion)
    y = .data |> attr("model_input") %$% y_proportion
  else
    y = .data |> attr("model_input") %$% y
  
  y = 
    y |> 
    as_tibble(rownames = quo_name(.sample)) |>
    pivot_longer(-!!.sample, names_to = quo_name(.cell_group), values_to = quo_name(.count)) |>
    with_groups(!!.sample,  ~ .x |> mutate(observed_proportion := !!.count / sum(!!.count ))) |>
    with_groups(!!.sample,  ~ .x |>  mutate(exposure := sum(!!.count))  )
  
  
  residual |>
    left_join(
      y,
      by = c(quo_name(.sample), quo_name(.cell_group))
    ) |>
    mutate(residuals = observed_proportion - proportion_mean) |>
    select(!!.sample, !!.cell_group, residuals, exposure)
  
}