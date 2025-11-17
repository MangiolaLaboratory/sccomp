# Clear Stan Model Cache

This function removes the Stan model cache directory and all its
contents. The cache is used by cmdstanr to store compiled Stan models.
Clearing the cache can be useful when:

- You're experiencing issues with cached models

- You want to force recompilation of models

- You're switching between different versions of Stan

## Usage

``` r
clear_stan_model_cache(cache_dir = sccomp_stan_models_cache_dir)
```

## Arguments

- cache_dir:

  A character string representing the path of the cache directory to
  delete. Defaults to `sccomp_stan_models_cache_dir`, which is the
  default cache location used by sccomp.

## Value

NULL invisibly. The function is called for its side effect of removing
the cache directory.

## Details

The function uses [`unlink()`](https://rdrr.io/r/base/unlink.html) with
`recursive = TRUE` to remove the directory and all its contents. If the
cache directory doesn't exist, the function will simply print a message
and continue. This is a safe operation as cmdstanr will recreate the
cache directory when needed.

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## See also

- [`cmdstan_path`](https://mc-stan.org/cmdstanr/reference/set_cmdstan_path.html)
  for information about cmdstanr's cache management

- [`unlink`](https://rdrr.io/r/base/unlink.html) for details about the
  underlying file removal function

## Examples

``` r
# Clear the default sccomp cache directory
clear_stan_model_cache()
#> Cache does not exist: /home/runner/.sccomp_models

# Clear a specific cache directory
clear_stan_model_cache("path/to/custom/cache")
#> Cache does not exist: path/to/custom/cache

# Example of when you might want to clear the cache
if (FALSE) { # \dontrun{
# If you're experiencing issues with model compilation
clear_stan_model_cache()
# Then try running your model again
} # }
```
