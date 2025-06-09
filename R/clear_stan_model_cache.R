#' Clear Stan Model Cache
#'
#' @description
#' This function removes the Stan model cache directory and all its contents. The cache is used by cmdstanr to store compiled Stan models.
#' Clearing the cache can be useful when:
#' \itemize{
#'   \item You're experiencing issues with cached models
#'   \item You want to force recompilation of models
#'   \item You're switching between different versions of Stan
#' }
#'
#' @details
#' The function uses `unlink()` with `recursive = TRUE` to remove the directory and all its contents.
#' If the cache directory doesn't exist, the function will simply print a message and continue.
#' This is a safe operation as cmdstanr will recreate the cache directory when needed.
#'
#' @param cache_dir A character string representing the path of the cache directory to delete.
#'   Defaults to `sccomp_stan_models_cache_dir`, which is the default cache location used by sccomp.
#'
#' @return NULL invisibly. The function is called for its side effect of removing the cache directory.
#'
#' @examples
#' # Clear the default sccomp cache directory
#' clear_stan_model_cache()
#'
#' # Clear a specific cache directory
#' clear_stan_model_cache("path/to/custom/cache")
#'
#' # Example of when you might want to clear the cache
#' \dontrun{
#' # If you're experiencing issues with model compilation
#' clear_stan_model_cache()
#' # Then try running your model again
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link[cmdstanr]{cmdstan_path}} for information about cmdstanr's cache management
#'   \item \code{\link[base]{unlink}} for details about the underlying file removal function
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
clear_stan_model_cache <- function(cache_dir = sccomp_stan_models_cache_dir) {
  
  # Check if the directory exists
  if (dir.exists(cache_dir)) {
    # Attempt to delete the directory and its contents
    unlink(cache_dir, recursive = TRUE)
    message("Cache deleted: ", cache_dir)
  } else {
    message("Cache does not exist: ", cache_dir)
  }
}