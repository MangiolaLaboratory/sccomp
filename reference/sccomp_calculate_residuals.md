# Calculate Residuals Between Observed and Predicted Proportions

`sccomp_calculate_residuals` computes the residuals between observed
cell group proportions and the predicted proportions from a fitted
`sccomp` model. This function is useful for assessing model fit and
identifying cell groups or samples where the model may not adequately
capture the observed data. The residuals are calculated as the
difference between the observed proportions and the predicted mean
proportions from the model.

## Usage

``` r
sccomp_calculate_residuals(.data)
```

## Arguments

- .data:

  A tibble of class `sccomp_tbl`, which is the result of
  [`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md).
  This tibble contains the fitted model and associated data necessary
  for calculating residuals.

## Value

A tibble (`tbl`) with the following columns:

- **sample** - A character column representing the sample identifiers.

- **cell_group** - A character column representing the cell group
  identifiers.

- **residuals** - A numeric column representing the residuals,
  calculated as the difference between observed and predicted
  proportions.

- **exposure** - A numeric column representing the total counts (sum of
  counts across cell groups) for each sample.

## Details

The function performs the following steps:

1.  Extracts the predicted mean proportions for each cell group and
    sample using
    [`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md).

2.  Calculates the observed proportions from the original count data.

3.  Computes residuals by subtracting the predicted proportions from the
    observed proportions.

4.  Returns a tibble containing the sample, cell group, residuals, and
    exposure (total counts per sample).

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## Examples

``` r
# \donttest{
  if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
# Load example data
data("counts_obj")

# Fit the sccomp model
estimates <- sccomp_estimate(
  counts_obj,
  formula_composition = ~ type,
  formula_variability = ~1,
  sample = "sample",
  cell_group = "cell_group",
  abundance = "count",
  approximate_posterior_inference = "all",
  cores = 1
)

# Calculate residuals
residuals <- sccomp_calculate_residuals(estimates)

# View the residuals
print(residuals)
}# }
#> Warning: The `approximate_posterior_inference` argument of `sccomp_estimate()` is
#> deprecated as of sccomp 1.7.7.
#> ℹ The argument approximate_posterior_inference is now deprecated. Please use
#>   inference_method. By default, inference_method value is inferred from
#>   approximate_posterior_inference.
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Init values were only set for a subset of parameters. 
#> Missing init values for the following parameters:
#> random_effect_raw, random_effect_raw_2, random_effect_sigma_mu, random_effect_sigma_sigma, random_effect_sigma_raw, sigma_correlation_factor, random_effect_sigma_raw_2, sigma_correlation_factor_2, zero_random_effect
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> ------------------------------------------------------------ 
#> EXPERIMENTAL ALGORITHM: 
#>   This procedure has not been thoroughly tested and may be unstable 
#>   or buggy. The interface is subject to change. 
#> ------------------------------------------------------------ 
#> Gradient evaluation took 0.000343 seconds 
#> 1000 transitions using 10 leapfrog steps per transition would take 3.43 seconds. 
#> Adjust your expectations accordingly! 
#> Begin eta adaptation. 
#> Iteration:   1 / 250 [  0%]  (Adaptation) 
#> Iteration:  50 / 250 [ 20%]  (Adaptation) 
#> Iteration: 100 / 250 [ 40%]  (Adaptation) 
#> Iteration: 150 / 250 [ 60%]  (Adaptation) 
#> Iteration: 200 / 250 [ 80%]  (Adaptation) 
#> Success! Found best value [eta = 1] earlier than expected. 
#> Begin stochastic gradient ascent. 
#>   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes  
#>    100        -4472.862             1.000            1.000 
#>    200        -3754.237             0.596            1.000 
#>    300        -3713.096             0.401            0.191 
#>    400        -3719.384             0.301            0.191 
#>    500        -3709.501             0.241            0.011 
#>    600        -3709.158             0.201            0.011 
#>    700        -3707.014             0.173            0.003   MEDIAN ELBO CONVERGED 
#> Drawing a sample of size 4000 from the approximate posterior...  
#> COMPLETED. 
#> Finished in  2.6 seconds.
#> Warning: Unknown or uninitialised column: `c_R_k_hat`.
#> Warning: no non-missing arguments to max; returning -Inf
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 4
#>    sample cell_group residuals exposure
#>    <chr>  <chr>          <dbl>    <int>
#>  1 10x_6K B1           0.0107      5030
#>  2 10x_6K B2          -0.0272      5030
#>  3 10x_6K B3          -0.00488     5030
#>  4 10x_6K BM           0.00204     5030
#>  5 10x_6K CD4 1        0.00210     5030
#>  6 10x_6K CD4 2       -0.0162      5030
#>  7 10x_6K CD4 3        0.0613      5030
#>  8 10x_6K CD4 4       -0.00105     5030
#>  9 10x_6K CD4 5       -0.00139     5030
#> 10 10x_6K CD8 1       -0.0283      5030
#> # ℹ 710 more rows
```
