# Default cache directory for Stan models

A global variable that defines the default cache directory for Stan
models used by the sccomp package. This directory is used to store
compiled Stan models to avoid recompilation on subsequent runs.

## Usage

``` r
sccomp_stan_models_cache_dir
```

## Format

An object of class `character` of length 1.

## Value

A character string containing the path to the default cache directory.

## Details

The cache directory is set to `~/.sccomp_models` by default. This
location is used by various sccomp functions to store and retrieve
compiled Stan models, improving performance by avoiding unnecessary
recompilation of models that have already been compiled.

Users can override this default by specifying a different cache
directory in function calls that accept a `cache_stan_model` parameter.

## See also

[`sccomp_estimate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md),
[`sccomp_replicate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_replicate.md),
[`clear_stan_model_cache`](https://mangiolalaboratory.github.io/sccomp/reference/clear_stan_model_cache.md)

## Examples

``` r
# View the default cache directory
sccomp_stan_models_cache_dir
#> [1] "/home/runner/.sccomp_models"

# Use a custom cache directory in a function call
# sccomp_estimate(data, cache_stan_model = "/path/to/custom/cache")
```
