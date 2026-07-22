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

- **residuals_unconstrained** - A numeric column representing the
  residuals on the unconstrained scale, calculated as the difference
  between the inverse softmax transform of observed proportions and the
  predicted unconstrained predictors.

## Details

The function performs the following steps:

1.  Extracts the predicted mean proportions and unconstrained predictors
    for each cell group and sample using
    [`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md).

2.  Calculates the observed proportions from the original count data.

3.  Computes residuals on the proportion scale by subtracting the
    predicted proportions from the observed proportions.

4.  Computes residuals on the unconstrained (log-ratio) scale by: (1)
    applying inverse softmax (log-ratio transform with sum-to-zero
    normalization) to observed proportions, (2) subtracting the
    predicted unconstrained predictors.

5.  Returns a tibble containing the sample, cell group, residuals
    (proportion scale), residuals_unconstrained (log-ratio scale), and
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
#> prec_slope_2, random_effect_raw_1, random_effect_sigma_raw_1, sigma_correlation_factor_1, random_effect_raw_2, random_effect_sigma_raw_2, sigma_correlation_factor_2, random_effect_raw_3, random_effect_sigma_raw_3, sigma_correlation_factor_3, random_effect_raw_4, random_effect_sigma_raw_4, sigma_correlation_factor_4, random_effect_sigma_mu, random_effect_sigma_sigma, zero_random_effect
#> 
#> To disable this message use options(cmdstanr_warn_inits = FALSE).
#> ------------------------------------------------------------ 
#> EXPERIMENTAL ALGORITHM: 
#>   This procedure has not been thoroughly tested and may be unstable 
#>   or buggy. The interface is subject to change. 
#> ------------------------------------------------------------ 
#> Gradient evaluation took 0.000349 seconds 
#> 1000 transitions using 10 leapfrog steps per transition would take 3.49 seconds. 
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
#>    100        -4668.733             1.000            1.000 
#>    200        -3777.042             0.618            1.000 
#>    300        -3704.740             0.419            0.236 
#>    400        -3698.475             0.314            0.236 
#>    500        -3694.261             0.252            0.020 
#>    600        -3694.743             0.210            0.020 
#>    700        -3699.562             0.180            0.002   MEDIAN ELBO CONVERGED 
#> Drawing a sample of size 4000 from the approximate posterior...  
#> COMPLETED. 
#> Finished in  2.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.558 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 5
#>    sample cell_group residuals exposure residuals_unconstrained
#>    <chr>  <chr>          <dbl>    <int>                   <dbl>
#>  1 10x_6K B1          0.0158       5030                 0.286  
#>  2 10x_6K B2         -0.0289       5030                -1.23   
#>  3 10x_6K B3         -0.00521      5030                -0.466  
#>  4 10x_6K BM          0.00159      5030                 0.248  
#>  5 10x_6K CD4 1      -0.000806     5030                 0.00530
#>  6 10x_6K CD4 2      -0.0128       5030                -0.272  
#>  7 10x_6K CD4 3       0.0588       5030                 0.586  
#>  8 10x_6K CD4 4      -0.00107      5030                -0.991  
#>  9 10x_6K CD4 5       0.00235      5030                 0.120  
#> 10 10x_6K CD8 1      -0.0204       5030                -0.179  
#> # ℹ 710 more rows
```
