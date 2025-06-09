#' Clear Stan Model Cache
#'
#' This function removes the Stan model cache directory and its contents. This is useful when:
#' \itemize{
#'   \item You need to force recompilation of Stan models
#'   \item You're experiencing issues with cached models
#'   \item You want to ensure you're using the latest model versions
#' }
#'
#' @param cache_dir A character string representing the path of the cache directory to delete. 
#'   Defaults to `sccomp_stan_models_cache_dir`. The cache directory stores compiled Stan models
#'   to avoid recompilation on subsequent runs.
#' 
#' @return NULL, invisibly. The function prints messages to indicate the status of the operation.
#' 
#' @examples
#' # Basic usage with default cache directory
#' clear_stan_model_cache()
#' 
#' # Clear a specific cache directory
#' clear_stan_model_cache("path/to/custom/cache")
#' 
#' # Example of when you might need to clear the cache:
#' # 1. After updating the Stan model code
#' # 2. When switching between different versions of cmdstanr
#' # 3. If you're experiencing compilation issues
#' 
#' @export
clear_stan_model_cache <- function(cache_dir = sccomp_stan_models_cache_dir) {
  
  # Check if the directory exists
  if (dir.exists(cache_dir)) {
    # Attempt to delete the directory and its contents
    unlink(cache_dir, recursive = TRUE)
    message("Cache successfully deleted: ", cache_dir)
  } else {
    message("Cache directory not found: ", cache_dir)
  }
  
  # Return invisibly
  invisible(NULL)
}

