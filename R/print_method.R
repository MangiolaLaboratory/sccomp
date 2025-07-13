#' Print method for sccomp objects
#'
#' @description
#' Print method for sccomp objects.
#' 
#' The print method for sccomp objects provides a summary of the model specifications,
#' data dimensions, and convergence diagnostics.
#' 
#' The output is formatted to be easy to read and understand.
#' 
#' @importFrom magrittr "%$%"
#' @param x A sccomp object
#' @param ... Additional arguments passed to print
#'
#' @return The printed object
#' @export
#'
#' @examples
#' 
#' # Note: Before running the example, ensure that the 'cmdstanr' package is installed:
#' # install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))
#'
#' \donttest{
#' if (instantiate::stan_cmdstan_exists()) {
#'   # Create a sccomp object
#'   data("counts_obj") 
#'   estimate <- sccomp_estimate(
#'     counts_obj,
#'     ~ type,
#'     ~1,
#'     "sample",
#'     "cell_group",
#'     "count",
#'     cores = 1
#'   )
#'   print(estimate)
#' }
#' }
#' 
print.sccomp_tbl <- function(x, ...) {
  
  # Get model information from attributes
  formula_composition <- attr(x, "formula_composition")
  formula_variability <- attr(x, "formula_variability")
  inference_method <- attr(x, "inference_method")
  noise_model <- attr(x, "noise_model")
  
  # Get data dimensions
  n_samples <- x |> attr("model_input") %$% N
  n_cell_groups <- x |> attr("model_input") %$% M 
  
  # Print header with model information
  cat("sccomp model\n")
  cat("============\n\n")
  
  cat("Model specifications:\n")
  cat("  Family:", noise_model, "\n")
  cat("  Composition formula:", deparse(formula_composition), "\n")
  cat("  Variability formula:", deparse(formula_variability), "\n")
  cat("  Inference method:", inference_method, "\n\n")
  
  cat("Data: Samples:", n_samples, "  Cell groups:", n_cell_groups, "\n\n")
  
  cat("Column prefixes: c_ -> composition parameters  v_ -> variability parameters\n\n")
  
  cat("Convergence diagnostics:\n")
  cat("  For each parameter, n_eff is the effective sample size and R_k_hat is the potential\n")
  cat("  scale reduction factor on split chains (at convergence, R_k_hat = 1).\n\n")
  
  # Print the tibble using pillar
  print(tibble::as_tibble(x), ...)
  
  # Return the object invisibly
  invisible(x)
} 