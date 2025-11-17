#' Clear Stan Draw Files
#'
#' @description
#' This function removes Stan MCMC draw CSV files from the output directory. These files can accumulate and consume significant disk space (often GBs).
#' Draw files are created during model fitting but are typically only needed during the analysis session.
#' Once results are extracted into R objects, the CSV files can be safely deleted.
#'
#' @details
#' The function can delete files based on:
#' \itemize{
#'   \item All files in the directory (default)
#'   \item Files older than a specified number of days
#'   \item Files larger than a specified size in MB
#'   \item Files matching a specific pattern
#' }
#'
#' @param output_directory A character string specifying the directory containing draw files.
#'   Defaults to "sccomp_draws_files", which is the default output directory used by sccomp.
#' @param older_than_days Numeric. If specified, only delete files older than this many days. Default is NULL (delete all).
#' @param larger_than_mb Numeric. If specified, only delete files larger than this size in MB. Default is NULL (no size filter).
#' @param pattern Character string. Regular expression pattern to match files. Default is "\\.csv$" (all CSV files).
#' @param dry_run Logical. If TRUE, only report what would be deleted without actually deleting. Default is FALSE.
#'
#' @return A list invisibly containing:
#'   \itemize{
#'     \item `files_deleted`: Number of files deleted
#'     \item `space_freed_mb`: Approximate disk space freed in MB
#'     \item `files`: Vector of file paths. When `dry_run = TRUE`, these are files that would be deleted; when `dry_run = FALSE`, these are files that were successfully deleted.
#'   }
#'
#' @examples
#' \dontrun{
#' # Clear all draw files in default directory
#' clear_draw_files()
#'
#' # See what would be deleted without actually deleting
#' clear_draw_files(dry_run = TRUE)
#'
#' # Delete only files older than 7 days
#' clear_draw_files(older_than_days = 7)
#'
#' # Delete only large files (>10 MB)
#' clear_draw_files(larger_than_mb = 10)
#'
#' # Delete from custom directory
#' clear_draw_files(output_directory = "my_custom_draws")
#'
#' # Combine filters: files older than 30 days AND larger than 50 MB
#' clear_draw_files(older_than_days = 30, larger_than_mb = 50)
#' }
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{clear_stan_model_cache}} for clearing compiled model cache
#'   \item \code{\link{sccomp_estimate}} which creates these draw files
#' }
#'
#' @references
#' S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z. Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust differential composition and variability analysis for single-cell data, Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120, https://doi.org/10.1073/pnas.2203828120 (2023).
#'
#' @export
clear_draw_files <- function(output_directory = "sccomp_draws_files",
                              older_than_days = NULL,
                              larger_than_mb = NULL,
                              pattern = "\\.csv$",
                              dry_run = FALSE) {
  
  # Check if directory exists
  if (!dir.exists(output_directory)) {
    message(sprintf("Draw file directory '%s' does not exist. Nothing to clean.", output_directory))
    return(invisible(list(files_deleted = 0, space_freed_mb = 0, files = character(0))))
  }
  
  # Get all files matching pattern
  all_files <- list.files(output_directory, pattern = pattern, full.names = TRUE)
  
  if (length(all_files) == 0) {
    message(sprintf("No files matching pattern '%s' found in '%s'.", pattern, output_directory))
    return(invisible(list(files_deleted = 0, space_freed_mb = 0, files = character(0))))
  }
  
  # Filter by age if specified
  if (!is.null(older_than_days)) {
    file_info <- file.info(all_files)
    age_days <- as.numeric(difftime(Sys.time(), file_info$mtime, units = "days"))
    all_files <- all_files[age_days > older_than_days]
    
    if (length(all_files) == 0) {
<<<<<<< HEAD
      message(sprintf("No files older than %.1f days found.", older_than_days))
=======
      message(sprintf("No files older than %d days found.", older_than_days))
>>>>>>> 3af1137 (Update version to 2.1.22 and add automatic cleanup for Stan draw files)
      return(invisible(list(files_deleted = 0, space_freed_mb = 0, files = character(0))))
    }
  }
  
  # Filter by size if specified
  if (!is.null(larger_than_mb)) {
    file_sizes_mb <- file.size(all_files) / (1024^2)
    all_files <- all_files[file_sizes_mb > larger_than_mb]
    
    if (length(all_files) == 0) {
<<<<<<< HEAD
      message(sprintf("No files larger than %.1f MB found.", larger_than_mb))
=======
      message(sprintf("No files larger than %d MB found.", larger_than_mb))
>>>>>>> 3af1137 (Update version to 2.1.22 and add automatic cleanup for Stan draw files)
      return(invisible(list(files_deleted = 0, space_freed_mb = 0, files = character(0))))
    }
  }
  
  # Calculate space to be freed
  space_to_free_mb <- sum(file.size(all_files)) / (1024^2)
  
  # Dry run - just report what would be deleted
  if (dry_run) {
    message(sprintf("DRY RUN: Would delete %d files (%.2f MB) from '%s'", 
                    length(all_files), space_to_free_mb, output_directory))
    if (length(all_files) <= 10) {
      message("Files that would be deleted:")
      message(paste("  -", basename(all_files), collapse = "\n"))
    } else {
      message(sprintf("First 10 files that would be deleted:"))
      message(paste("  -", basename(all_files[1:10]), collapse = "\n"))
      message(sprintf("  ... and %d more files", length(all_files) - 10))
    }
    return(invisible(list(files_deleted = length(all_files), 
                          space_freed_mb = space_to_free_mb, 
                          files = all_files)))
  }
  
  # Actually delete the files
  success <- file.remove(all_files)
  files_deleted <- sum(success)
  
  if (files_deleted > 0) {
    message(sprintf("Deleted %d draw files (%.2f MB freed) from '%s'", 
                    files_deleted, space_to_free_mb, output_directory))
  } else {
    warning("Failed to delete any files.")
  }
  
  invisible(list(
    files_deleted = files_deleted,
    space_freed_mb = space_to_free_mb,
    files = all_files[success]
  ))
}

