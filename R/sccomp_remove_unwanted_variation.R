#' DEPRECATED: Remove Unwanted Variation from sccomp Estimates
#'
#' @description This function is DEPRECATED. Please use \code{\link{sccomp_remove_unwanted_effects}} instead.
#' This function uses the model to remove unwanted variation from a dataset using the estimates of the model. 
#' For example, if you fit your data with the formula `~ factor_1 + factor_2` and use the formula `~ factor_1` 
#' to remove unwanted variation, the `factor_2` effect will be factored out.
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
#'       ~ type, ~1, "sample", "cell_group", "count",
#'       cores = 1
#'     ) |>
#'     sccomp_remove_unwanted_variation()
#'   }
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
sccomp_remove_unwanted_variation <- function(.data,
                                             formula_composition_keep = NULL,
                                             formula_composition = NULL,
                                             formula_variability = NULL,
                                             cores = detectCores()) {
  lifecycle::deprecate_warn(
    "1.99.20",
    "sccomp::sccomp_remove_unwanted_variation()",
    details = "sccomp says: sccomp_remove_unwanted_variation is deprecated. Please use sccomp_remove_unwanted_effects() instead."
  )
  
  sccomp_remove_unwanted_effects(
    .data = .data,
    formula_composition_keep = formula_composition_keep,
    formula_composition = formula_composition,
    cores = cores
  )
}

