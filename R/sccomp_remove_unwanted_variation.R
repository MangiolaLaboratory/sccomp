

#' sccomp_remove_unwanted_variation
#'
#' @description This function uses the model to remove unwanted variation from a dataset using the estimates of the model. For example, if you fit your data with the formula `~ factor_1 + factor_2` and use the formula `~ factor_1` to remove unwanted variation, the `factor_2` effect will be factored out.
#'
#' @param .data A tibble. The result of `sccomp_estimate`.
#' @param formula_composition_keep A formula. The formula describing the model for differential abundance, for example `~type`. In this case, only the effect of the `type` factor will be preserved, while all other factors will be factored out.
#' @param formula_composition DEPRECATED. Use `formula_composition_keep` instead.
#' @param formula_variability DEPRECATED. Use `formula_variability_keep` instead.
#' @param cores Integer, the number of cores to be used for parallel calculations.
#' 
#' @return A tibble (`tbl`) with the following columns:
#' \itemize{
#'   \item \strong{sample} - A character column representing the sample name for which data was adjusted.
#'   \item \strong{cell_group} - A character column representing the cell group being tested.
#'   \item \strong{adjusted_proportion} - A numeric column representing the adjusted proportion after removing unwanted variation.
#'   \item \strong{adjusted_counts} - A numeric column representing the adjusted counts after removing unwanted variation.
#'   \item \strong{logit_residuals} - A numeric column representing the logit residuals calculated after adjustment.
#' }
#'
#' @export
#'
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimates = sccomp_estimate(
#'       counts_obj,
#'       ~ type, ~1, sample, cell_group, count,
#'       cores = 1
#'     ) |>
#'     sccomp_remove_unwanted_variation()
#'   }
#' }
#'
sccomp_remove_unwanted_variation <- function(.data,
                                             formula_composition_keep = NULL,
                                             formula_composition = NULL,
                                             formula_variability = NULL,
                                             cores = detectCores()) {
  
  # Check for deprecated arguments
  if (!is.null(formula_composition)) {
    warning("The argument 'formula_composition' is deprecated. Please use 'formula_composition_keep' instead.", call. = FALSE)
    formula_composition_keep <- formula_composition
  }
  
  if (!is.null(formula_variability)) {
    warning("The argument 'formula_variability' is deprecated as not used.", call. = FALSE)
  }
  
  # Run the function
  check_and_install_cmdstanr()
  
  UseMethod("sccomp_remove_unwanted_variation", .data)
}

#' @importFrom readr write_file
#' @export
#'
sccomp_remove_unwanted_variation.sccomp_tbl = function(.data,
                                                       formula_composition_keep = NULL,
                                                       formula_composition = NULL,
                                                       formula_variability = NULL,
                                                       cores = detectCores()){
  
  
  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .grouping_for_random_effect = attr(.data, ".grouping_for_random_effect")
  .count = attr(.data, ".count")
  
  # Residuals
  message("sccomp says: calculating residuals")
  residuals = .data |> sccomp_calculate_residuals()
  
  
  message("sccomp says: regressing out unwanted factors")
  
  
  # Generate quantities
  .data |>
    sccomp_predict(
      formula_composition = formula_composition_keep,
      number_of_draws = attr(.data, "fit") |>  get_output_samples() |> min(500)
      
    ) |>
    distinct(!!.sample, !!.cell_group, proportion_mean) |>
    # mutate(proportion_mean =
    #          proportion_mean |>
    #          # compress_zero_one() |>
    #          boot::logit()
    # ) |>
    left_join(residuals,  by = c(quo_name(.sample), quo_name(.cell_group))) |>
    mutate(adjusted_proportion = proportion_mean + residuals) |>
    mutate(adjusted_proportion = adjusted_proportion |> pmax(0)) |> 
    # mutate(adjusted_proportion = adjusted_proportion |> boot::inv.logit()) |>
    # with_groups(!!.sample,  ~ .x |> mutate(adjusted_proportion := adjusted_proportion / sum(adjusted_proportion ))) |>
    
    # Recostituite counts
    mutate(adjusted_counts = adjusted_proportion * exposure) |>
    
    select(!!.sample, !!.cell_group, adjusted_proportion, adjusted_counts, residuals)
  
  
  
}