#' Remove unwanted effects from a sccomp_tbl object
#'
#' @description This method removes unwanted effects from a dataset using the model estimates. For example, if you fit your data with the formula `~ factor_1 + factor_2` and use the formula `~ factor_1` to remove unwanted variation, the `factor_2` effect will be factored out.
#'
#' @param .data A `sccomp_tbl` object. The result of `sccomp_estimate`.
#' @param formula_composition_keep A formula. The formula describing the model for differential abundance, for example `~type`. In this case, only the effect of the `type` factor will be preserved, while all other factors will be factored out.
#' @param formula_composition DEPRECATED. Use `formula_composition_keep` instead.
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
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_remove_unwanted_effects()
#'   }
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
sccomp_remove_unwanted_effects <- function(.data,
                                           formula_composition_keep = NULL,
                                           formula_composition = NULL,
                                           cores = detectCores()) {
  # Check for deprecated arguments
  if (!is.null(formula_composition)) {
    warning("The argument 'formula_composition' is deprecated. Please use 'formula_composition_keep' instead.", call. = FALSE)
    formula_composition_keep <- formula_composition
  }

  check_and_install_cmdstanr()
  UseMethod("sccomp_remove_unwanted_effects", .data)
}

#' @importFrom readr write_file
#' @export
#'
sccomp_remove_unwanted_effects.sccomp_tbl = function(.data,
                                                       formula_composition_keep = NULL,
                                                       formula_composition = NULL,
                                                       cores = detectCores()){
  
  
  model_input = attr(.data, "model_input")
  .sample = attr(.data, ".sample")
  .cell_group = attr(.data, ".cell_group")
  .grouping_for_random_effect = attr(.data, ".grouping_for_random_effect")
  .count = attr(.data, ".count")
  
  # Calculate residuals from the model
  message("sccomp says: calculating residuals")
  residuals = .data |> sccomp_calculate_residuals()

  # Predict proportions using only the factors we want to keep
  message("sccomp says: regressing out unwanted factors")
  .data |>

    # Generate predictions using the specified formula
    sccomp_predict(
      formula_composition = formula_composition_keep,
      # Use up to 500 samples from the fit for prediction
      number_of_draws = attr(.data, "fit") |>  get_output_samples() |> min(500)
    ) |>
    # Keep only unique sample-cell_group combinations with their mean proportions
    distinct(!!.sample, !!.cell_group, proportion_mean) |>

    # Join with residuals to combine predicted proportions with residual effects
    left_join(residuals,  by = c(quo_name(.sample), quo_name(.cell_group))) |>

    # Calculate adjusted proportions by adding residuals to predicted means
    mutate(adjusted_proportion = proportion_mean + residuals) |>

    # Ensure proportions are non-negative
    mutate(adjusted_proportion = adjusted_proportion |> pmax(0)) |> 

    # Calculate adjusted counts based on adjusted proportions and exposure
    mutate(adjusted_counts = adjusted_proportion * exposure) |>

    # Select final columns for output
    select(!!.sample, !!.cell_group, adjusted_proportion, adjusted_counts, residuals)
}