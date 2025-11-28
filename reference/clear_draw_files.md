# Clear Stan Draw Files

This function removes Stan MCMC draw CSV files from the output
directory. These files can accumulate and consume significant disk space
(often GBs). Draw files are created during model fitting but are
typically only needed during the analysis session. Once results are
extracted into R objects, the CSV files can be safely deleted.

## Usage

``` r
clear_draw_files(
  output_directory = "sccomp_draws_files",
  older_than_days = NULL,
  larger_than_mb = NULL,
  pattern = "\\.csv$",
  dry_run = FALSE
)
```

## Arguments

- output_directory:

  A character string specifying the directory containing draw files.
  Defaults to "sccomp_draws_files", which is the default output
  directory used by sccomp.

- older_than_days:

  Numeric. If specified, only delete files older than this many days.
  Default is NULL (delete all).

- larger_than_mb:

  Numeric. If specified, only delete files larger than this size in MB.
  Default is NULL (no size filter).

- pattern:

  Character string. Regular expression pattern to match files. Default
  is "\\csv\$" (all CSV files).

- dry_run:

  Logical. If TRUE, only report what would be deleted without actually
  deleting. Default is FALSE.

## Value

A list invisibly containing:

- `files_deleted`: Number of files deleted

- `space_freed_mb`: Approximate disk space freed in MB

- `files`: Vector of deleted file paths (if dry_run = FALSE)

## Details

The function can delete files based on:

- All files in the directory (default)

- Files older than a specified number of days

- Files larger than a specified size in MB

- Files matching a specific pattern

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## See also

- [`clear_stan_model_cache`](https://mangiolalaboratory.github.io/sccomp/reference/clear_stan_model_cache.md)
  for clearing compiled model cache

- [`sccomp_estimate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
  which creates these draw files

## Examples

``` r
if (FALSE) { # \dontrun{
# Clear all draw files in default directory
clear_draw_files()

# See what would be deleted without actually deleting
clear_draw_files(dry_run = TRUE)

# Delete only files older than 7 days
clear_draw_files(older_than_days = 7)

# Delete only large files (>10 MB)
clear_draw_files(larger_than_mb = 10)

# Delete from custom directory
clear_draw_files(output_directory = "my_custom_draws")

# Combine filters: files older than 30 days AND larger than 50 MB
clear_draw_files(older_than_days = 30, larger_than_mb = 50)
} # }
```
