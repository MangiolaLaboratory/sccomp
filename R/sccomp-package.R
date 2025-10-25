#' @keywords internal
"_PACKAGE"

#' sccomp: Differential Composition and Variability Analysis for Single-Cell Data
#'
#' @description
#' The sccomp package provides comprehensive tools for differential composition 
#' and variability analysis in single-cell RNA sequencing, CyTOF, and microbiome data. 
#' It implements robust Bayesian modeling with outlier detection, random effects, 
#' and advanced statistical methods for cell type proportion analysis.
#'
#' @details
#' The main functions are:
#'
#' \itemize{
#' \item \code{\link{sccomp_estimate}} - Perform differential composition and variability analysis
#' \item \code{\link{sccomp_test}} - Test for differential composition with significance testing
#' \item \code{\link{sccomp_remove_outliers}} - Identify and remove outlier samples
#' \item \code{\link{sccomp_predict}} - Predict cell type proportions for new samples
#' \item \code{\link{sccomp_remove_unwanted_variation}} - Remove unwanted variation from data
#' \item \code{\link{sccomp_proportional_fold_change}} - Calculate proportional fold changes
#' \item Plotting functions: \code{\link{plot.sccomp_tbl}}, \code{\link{sccomp_boxplot}}, \code{\link{plot_1D_intervals}}, \code{\link{plot_2D_intervals}}
#' }
#'
#' For detailed information on usage, see the package vignettes:
#'
#' \code{browseVignettes("sccomp")}
#'
#' All software-related questions should be posted to the GitHub Issues page:
#'
#' \url{https://github.com/MangiolaLaboratory/sccomp/issues}
#'
#' The code can be viewed at the GitHub repository:
#'
#' \url{https://github.com/MangiolaLaboratory/sccomp}
#'
#' @references
#' Mangiola, S., Roth-Schulze, A.J., Trussart, M., Zozaya-Valdés, E., Ma, M., 
#' Gao, Z., Rubin, A.F., Speed, T.P., Shim, H., & Papenfuss, A.T. (2023).
#' sccomp: Robust differential composition and variability analysis for single-cell data.
#' Proceedings of the National Academy of Sciences, 120(33), e2203828120.
#' \doi{10.1073/pnas.2203828120}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/MangiolaLaboratory/sccomp}
#'   \item Report bugs at \url{https://github.com/MangiolaLaboratory/sccomp/issues}
#' }
#'
#' @docType package
#' @name sccomp-package
#' @aliases sccomp
NULL
