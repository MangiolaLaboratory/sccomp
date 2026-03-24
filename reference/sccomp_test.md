# sccomp_test

This function test contrasts from a sccomp result.

## Usage

``` r
sccomp_test(
  .data,
  contrasts = NULL,
  percent_false_positive = 5,
  test_composition_above_logit_fold_change = 0.1,
  pass_fit = TRUE
)
```

## Arguments

- .data:

  A tibble. The result of sccomp_estimate.

- contrasts:

  A vector of character strings. For example if your formula is
  `~ 0 + treatment` and the factor treatment has values `yes` and `no`,
  your contrast could be "constrasts = c(treatmentyes - treatmentno)".

- percent_false_positive:

  A real between 0 and 100 non included. This used to identify outliers
  with a specific false positive rate.

- test_composition_above_logit_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the logit uncontrained
  scale.

- pass_fit:

  A boolean. Whether to pass the Stan fit as attribute in the output.
  Because the Stan fit can be very large, setting this to FALSE can be
  used to lower the memory imprint to save the output.

## Value

A tibble (`tbl`), with the following columns:

- cell_group - The cell groups being tested.

- parameter - The parameter being estimated from the design matrix
  described by the input formula_composition and formula_variability.

- factor - The covariate factor in the formula, if applicable (e.g., not
  present for Intercept or contrasts).

- c_lower - Lower (2.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_effect - Mean of the posterior distribution for a composition (c)
  parameter.

- c_upper - Upper (97.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_pH0 - Probability of the c_effect being smaller or bigger than the
  `test_composition_above_logit_fold_change` argument.

- c_FDR - False-discovery rate of the null hypothesis (no difference)
  for a composition (c).

- c_n_eff - Effective sample size - the number of independent draws in
  the sample, the higher the better.

- c_R_k_hat - R statistic, a measure of chain equilibrium, should be
  within 0.05 of 1.0.

- v_lower - Lower (2.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_effect - Mean of the posterior distribution for a variability (v)
  parameter.

- v_upper - Upper (97.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_pH0 - Probability of the null hypothesis (no difference) for a
  variability (v).

- v_FDR - False-discovery rate of the null hypothesis (no difference)
  for a variability (v).

- v_n_eff - Effective sample size for a variability (v) parameter.

- v_R_k_hat - R statistic for a variability (v) parameter.

- count_data - Nested input count data.

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimates = sccomp_estimate(
      counts_obj,
      ~ 0 + type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_test("typecancer - typebenign")
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: typebenign, typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481406.184653 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              41      -4.788e+05      8.490e-03   2.359e-01    1.000e+00  1.000e+00      2160 -3.710e+03 -3.712e+03                   
#> Path [1] :Best Iter: [34] ELBO (-3710.284388) evaluations: (2160) 
#> Path [2] :Initial log joint density = -481780.242506 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.098e-02   1.733e-01    1.000e+00  1.000e+00      2924 -3.711e+03 -3.713e+03                   
#> Path [2] :Best Iter: [42] ELBO (-3710.778947) evaluations: (2924) 
#> Path [3] :Initial log joint density = -481204.020824 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.142e-03   1.482e-01    9.023e-01  9.023e-01      2495 -3.709e+03 -3.715e+03                   
#> Path [3] :Best Iter: [36] ELBO (-3709.468091) evaluations: (2495) 
#> Path [4] :Initial log joint density = -481598.336702 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.391e-02   2.225e-01    9.623e-01  9.623e-01      2324 -3.711e+03 -3.713e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3711.352642) evaluations: (2324) 
#> Path [5] :Initial log joint density = -481329.092713 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.036e-02   2.238e-01    1.000e+00  1.000e+00      2560 -3.711e+03 -3.712e+03                   
#> Path [5] :Best Iter: [43] ELBO (-3710.881619) evaluations: (2560) 
#> Path [6] :Initial log joint density = -481718.265303 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.749e-02   2.071e-01    1.000e+00  1.000e+00      2516 -3.712e+03 -3.718e+03                   
#> Path [6] :Best Iter: [36] ELBO (-3711.543660) evaluations: (2516) 
#> Path [7] :Initial log joint density = -481423.066296 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.267e-02   1.936e-01    1.000e+00  1.000e+00      2754 -3.711e+03 -3.712e+03                   
#> Path [7] :Best Iter: [49] ELBO (-3710.527526) evaluations: (2754) 
#> Path [8] :Initial log joint density = -481783.822809 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.352e-02   2.033e-01    1.000e+00  1.000e+00      2773 -3.710e+03 -3.715e+03                   
#> Path [8] :Best Iter: [37] ELBO (-3710.357598) evaluations: (2773) 
#> Path [9] :Initial log joint density = -481595.510861 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.947e-02   2.473e-01    1.000e+00  1.000e+00      3150 -3.709e+03 -3.708e+03                   
#> Path [9] :Best Iter: [56] ELBO (-3708.088165) evaluations: (3150) 
#> Path [10] :Initial log joint density = -487902.679669 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      8.252e-03   2.121e-01    1.000e+00  1.000e+00      2791 -3.711e+03 -3.711e+03                   
#> Path [10] :Best Iter: [50] ELBO (-3711.054782) evaluations: (2791) 
#> Path [11] :Initial log joint density = -483381.449596 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      9.687e-03   1.338e-01    1.000e+00  1.000e+00      2560 -3.710e+03 -3.713e+03                   
#> Path [11] :Best Iter: [34] ELBO (-3710.369466) evaluations: (2560) 
#> Path [12] :Initial log joint density = -482567.496230 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      6.014e-03   2.116e-01    9.577e-01  9.577e-01      2723 -3.710e+03 -3.714e+03                   
#> Path [12] :Best Iter: [48] ELBO (-3710.171260) evaluations: (2723) 
#> Path [13] :Initial log joint density = -481232.180624 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      5.538e-03   1.214e-01    1.000e+00  1.000e+00      2378 -3.710e+03 -3.711e+03                   
#> Path [13] :Best Iter: [42] ELBO (-3709.867885) evaluations: (2378) 
#> Path [14] :Initial log joint density = -481517.923381 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.025e-02   2.392e-01    7.809e-01  7.809e-01      2589 -3.708e+03 -3.715e+03                   
#> Path [14] :Best Iter: [48] ELBO (-3708.121378) evaluations: (2589) 
#> Path [15] :Initial log joint density = -481180.385233 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      6.812e-03   1.675e-01    1.000e+00  1.000e+00      2295 -3.711e+03 -3.721e+03                   
#> Path [15] :Best Iter: [34] ELBO (-3710.509398) evaluations: (2295) 
#> Path [16] :Initial log joint density = -484042.303680 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.585e-03   1.921e-01    1.000e+00  1.000e+00      2789 -3.712e+03 -3.719e+03                   
#> Path [16] :Best Iter: [41] ELBO (-3711.579992) evaluations: (2789) 
#> Path [17] :Initial log joint density = -481393.426309 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      2.395e-02   3.235e-01    1.000e+00  1.000e+00      2800 -3.711e+03 -3.714e+03                   
#> Path [17] :Best Iter: [33] ELBO (-3710.545915) evaluations: (2800) 
#> Path [18] :Initial log joint density = -481671.923063 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      2.972e-03   2.162e-01    7.766e-01  7.766e-01      2635 -3.709e+03 -3.716e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3708.595833) evaluations: (2635) 
#> Path [19] :Initial log joint density = -481524.213331 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      1.619e-02   2.753e-01    9.195e-01  9.195e-01      2316 -3.711e+03 -3.713e+03                   
#> Path [19] :Best Iter: [29] ELBO (-3710.673958) evaluations: (2316) 
#> Path [20] :Initial log joint density = -481686.301390 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.461e-02   2.214e-01    1.000e+00  1.000e+00      2550 -3.711e+03 -3.714e+03                   
#> Path [20] :Best Iter: [33] ELBO (-3711.281616) evaluations: (2550) 
#> Path [21] :Initial log joint density = -481986.323442 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      9.567e-03   1.661e-01    9.824e-01  9.824e-01      2808 -3.711e+03 -3.713e+03                   
#> Path [21] :Best Iter: [43] ELBO (-3710.672959) evaluations: (2808) 
#> Path [22] :Initial log joint density = -481628.797837 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.504e-03   1.768e-01    9.631e-01  9.631e-01      2719 -3.712e+03 -3.721e+03                   
#> Path [22] :Best Iter: [34] ELBO (-3712.362526) evaluations: (2719) 
#> Path [23] :Initial log joint density = -483892.353935 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.137e-03   2.138e-01    1.000e+00  1.000e+00      2538 -3.708e+03 -3.712e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3707.734542) evaluations: (2538) 
#> Path [24] :Initial log joint density = -481637.876455 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.765e-02   2.464e-01    1.000e+00  1.000e+00      2993 -3.709e+03 -3.712e+03                   
#> Path [24] :Best Iter: [42] ELBO (-3709.447405) evaluations: (2993) 
#> Path [25] :Initial log joint density = -481563.913742 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      6.308e-03   1.770e-01    1.000e+00  1.000e+00      2576 -3.710e+03 -3.711e+03                   
#> Path [25] :Best Iter: [37] ELBO (-3709.919996) evaluations: (2576) 
#> Path [26] :Initial log joint density = -481418.204449 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.562e-02   1.878e-01    1.000e+00  1.000e+00      2587 -3.710e+03 -3.713e+03                   
#> Path [26] :Best Iter: [42] ELBO (-3709.941523) evaluations: (2587) 
#> Path [27] :Initial log joint density = -481295.297919 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      1.267e-02   3.348e-01    9.913e-01  9.913e-01      2295 -3.710e+03 -3.715e+03                   
#> Path [27] :Best Iter: [42] ELBO (-3710.175715) evaluations: (2295) 
#> Path [28] :Initial log joint density = -481329.992260 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      7.471e-03   1.656e-01    1.000e+00  1.000e+00      2379 -3.712e+03 -3.719e+03                   
#> Path [28] :Best Iter: [30] ELBO (-3711.548794) evaluations: (2379) 
#> Path [29] :Initial log joint density = -481396.404204 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.306e-02   2.222e-01    1.000e+00  1.000e+00      2725 -3.711e+03 -3.716e+03                   
#> Path [29] :Best Iter: [37] ELBO (-3711.123785) evaluations: (2725) 
#> Path [30] :Initial log joint density = -481576.896286 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.649e-02   2.039e-01    1.000e+00  1.000e+00      2812 -3.713e+03 -3.716e+03                   
#> Path [30] :Best Iter: [46] ELBO (-3712.898166) evaluations: (2812) 
#> Path [31] :Initial log joint density = -481518.794679 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.140e-02   1.945e-01    1.000e+00  1.000e+00      2961 -3.709e+03 -3.707e+03                   
#> Path [31] :Best Iter: [53] ELBO (-3706.681963) evaluations: (2961) 
#> Path [32] :Initial log joint density = -481248.005383 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.122e-02   1.838e-01    1.000e+00  1.000e+00      2469 -3.710e+03 -3.711e+03                   
#> Path [32] :Best Iter: [46] ELBO (-3710.023418) evaluations: (2469) 
#> Path [33] :Initial log joint density = -482382.903016 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.772e-02   2.831e-01    9.719e-01  9.719e-01      2298 -3.709e+03 -3.721e+03                   
#> Path [33] :Best Iter: [36] ELBO (-3708.789932) evaluations: (2298) 
#> Path [34] :Initial log joint density = -481318.395503 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.952e-03   1.459e-01    8.683e-01  8.683e-01      2408 -3.710e+03 -3.711e+03                   
#> Path [34] :Best Iter: [39] ELBO (-3709.500862) evaluations: (2408) 
#> Path [35] :Initial log joint density = -481598.740341 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.602e-02   2.767e-01    8.236e-01  8.236e-01      2554 -3.711e+03 -3.727e+03                   
#> Path [35] :Best Iter: [38] ELBO (-3710.712903) evaluations: (2554) 
#> Path [36] :Initial log joint density = -481546.640820 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.809e-03   1.186e-01    9.937e-01  9.937e-01      3166 -3.712e+03 -3.711e+03                   
#> Path [36] :Best Iter: [54] ELBO (-3711.309546) evaluations: (3166) 
#> Path [37] :Initial log joint density = -481308.912001 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      6.567e-03   1.358e-01    1.000e+00  1.000e+00      2144 -3.710e+03 -3.713e+03                   
#> Path [37] :Best Iter: [41] ELBO (-3710.041954) evaluations: (2144) 
#> Path [38] :Initial log joint density = -481468.810688 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.590e-02   3.147e-01    1.000e+00  1.000e+00      2407 -3.710e+03 -3.713e+03                   
#> Path [38] :Best Iter: [40] ELBO (-3710.031663) evaluations: (2407) 
#> Path [39] :Initial log joint density = -481432.051051 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      3.882e-03   2.741e-01    5.586e-01  5.586e-01      2771 -3.709e+03 -3.714e+03                   
#> Path [39] :Best Iter: [43] ELBO (-3709.360986) evaluations: (2771) 
#> Path [40] :Initial log joint density = -481444.943240 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.526e-03   1.523e-01    1.000e+00  1.000e+00      2408 -3.710e+03 -3.722e+03                   
#> Path [40] :Best Iter: [32] ELBO (-3709.765093) evaluations: (2408) 
#> Path [41] :Initial log joint density = -481435.622888 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              43      -4.788e+05      1.754e-02   2.331e-01    8.535e-01  8.535e-01      2190 -3.711e+03 -3.714e+03                   
#> Path [41] :Best Iter: [40] ELBO (-3710.521580) evaluations: (2190) 
#> Path [42] :Initial log joint density = -481816.847685 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      7.381e-03   2.524e-01    1.000e+00  1.000e+00      2405 -3.711e+03 -3.712e+03                   
#> Path [42] :Best Iter: [38] ELBO (-3710.934678) evaluations: (2405) 
#> Path [43] :Initial log joint density = -481330.230157 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.246e-02   2.371e-01    1.000e+00  1.000e+00      2296 -3.709e+03 -3.714e+03                   
#> Path [43] :Best Iter: [37] ELBO (-3709.258887) evaluations: (2296) 
#> Path [44] :Initial log joint density = -481315.512809 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      8.059e-03   1.993e-01    8.550e-01  8.550e-01      2444 -3.709e+03 -3.714e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3709.243692) evaluations: (2444) 
#> Path [45] :Initial log joint density = -481364.627046 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.284e-02   2.035e-01    1.000e+00  1.000e+00      2905 -3.709e+03 -3.716e+03                   
#> Path [45] :Best Iter: [51] ELBO (-3709.282303) evaluations: (2905) 
#> Path [46] :Initial log joint density = -481842.815747 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.216e-02   1.627e-01    1.000e+00  1.000e+00      3051 -3.711e+03 -3.710e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3710.418957) evaluations: (3051) 
#> Path [47] :Initial log joint density = -481468.680146 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.210e-02   3.238e-01    1.000e+00  1.000e+00      2691 -3.712e+03 -3.717e+03                   
#> Path [47] :Best Iter: [48] ELBO (-3711.784298) evaluations: (2691) 
#> Path [48] :Initial log joint density = -482155.125325 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              42      -4.788e+05      1.159e-02   2.768e-01    4.581e-01  1.000e+00      2154 -3.714e+03 -3.727e+03                   
#> Path [48] :Best Iter: [40] ELBO (-3713.859141) evaluations: (2154) 
#> Path [49] :Initial log joint density = -481968.308474 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.932e-02   3.092e-01    1.000e+00  1.000e+00      2922 -3.709e+03 -3.713e+03                   
#> Path [49] :Best Iter: [51] ELBO (-3709.073018) evaluations: (2922) 
#> Path [50] :Initial log joint density = -481925.514947 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.081e-02   2.127e-01    9.688e-01  9.688e-01      2368 -3.711e+03 -3.722e+03                   
#> Path [50] :Best Iter: [34] ELBO (-3710.524516) evaluations: (2368) 
#> Finished in  11.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
