#' The 'sccomp' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name sccomp-package
#' @aliases sccomp
#' @useDynLib sccomp, .registration = TRUE
#' @import methods
#' @import RcppParallel
#' @import rstantools
#' @import SeuratObject
#' @import lifecycle
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
NULL
if (FALSE) {
  lifecycle::deprecate_soft()
  RcppParallel::CxxFlags()
  rstantools::loo_R2()
  SeuratObject::subset()
}
