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
#> Path [1] :Initial log joint density = -481220.950047 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.198e-02   1.678e-01    1.000e+00  1.000e+00      2447 -3.699e+03 -3.702e+03                   
#> Path [1] :Best Iter: [41] ELBO (-3699.209178) evaluations: (2447) 
#> Path [2] :Initial log joint density = -483113.145737 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.255e-02   3.050e-01    7.817e-01  7.817e-01      2688 -3.697e+03 -3.701e+03                   
#> Path [2] :Best Iter: [48] ELBO (-3697.154719) evaluations: (2688) 
#> Path [3] :Initial log joint density = -481788.226291 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      7.833e-03   1.594e-01    1.000e+00  1.000e+00      4075 -3.694e+03 -3.695e+03                   
#> Path [3] :Best Iter: [60] ELBO (-3694.464571) evaluations: (4075) 
#> Path [4] :Initial log joint density = -485739.237042 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      2.410e-02   2.388e-01    1.000e+00  1.000e+00      3027 -3.699e+03 -3.700e+03                   
#> Path [4] :Best Iter: [49] ELBO (-3699.292779) evaluations: (3027) 
#> Path [5] :Initial log joint density = -481469.067837 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      9.511e-03   1.789e-01    1.000e+00  1.000e+00      2559 -3.701e+03 -3.702e+03                   
#> Path [5] :Best Iter: [30] ELBO (-3701.061534) evaluations: (2559) 
#> Path [6] :Initial log joint density = -481510.727719 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      5.969e-03   1.221e-01    8.807e-01  8.807e-01      4133 -3.694e+03 -3.697e+03                   
#> Path [6] :Best Iter: [61] ELBO (-3694.315253) evaluations: (4133) 
#> Path [7] :Initial log joint density = -481373.714188 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      5.368e-03   1.551e-01    1.000e+00  1.000e+00      2416 -3.700e+03 -3.707e+03                   
#> Path [7] :Best Iter: [33] ELBO (-3699.563701) evaluations: (2416) 
#> Path [8] :Initial log joint density = -481368.675290 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      8.060e-03   1.362e-01    1.000e+00  1.000e+00      3824 -3.696e+03 -3.694e+03                   
#> Path [8] :Best Iter: [63] ELBO (-3694.381211) evaluations: (3824) 
#> Path [9] :Initial log joint density = -481402.983298 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.620e-02   2.261e-01    1.000e+00  1.000e+00      2531 -3.700e+03 -3.700e+03                   
#> Path [9] :Best Iter: [34] ELBO (-3699.974092) evaluations: (2531) 
#> Path [10] :Initial log joint density = -481539.991953 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      2.184e-02   2.334e-01    8.868e-01  8.868e-01      2573 -3.697e+03 -3.703e+03                   
#> Path [10] :Best Iter: [44] ELBO (-3697.098601) evaluations: (2573) 
#> Path [11] :Initial log joint density = -481198.330050 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      6.968e-03   1.235e-01    1.000e+00  1.000e+00      2408 -3.701e+03 -3.702e+03                   
#> Path [11] :Best Iter: [42] ELBO (-3700.664079) evaluations: (2408) 
#> Path [12] :Initial log joint density = -482710.342764 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.532e-02   3.123e-01    1.000e+00  1.000e+00      2576 -3.701e+03 -3.703e+03                   
#> Path [12] :Best Iter: [41] ELBO (-3701.363255) evaluations: (2576) 
#> Path [13] :Initial log joint density = -481106.336319 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.623e-02   1.289e-01    1.000e+00  1.000e+00      2344 -3.699e+03 -3.703e+03                   
#> Path [13] :Best Iter: [34] ELBO (-3699.259045) evaluations: (2344) 
#> Path [14] :Initial log joint density = -481614.649944 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.714e-03   2.457e-01    8.841e-01  8.841e-01      3468 -3.696e+03 -3.697e+03                   
#> Path [14] :Best Iter: [56] ELBO (-3696.462971) evaluations: (3468) 
#> Path [15] :Initial log joint density = -484216.421770 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.154e-02   2.217e-01    1.000e+00  1.000e+00      2667 -3.697e+03 -3.698e+03                   
#> Path [15] :Best Iter: [43] ELBO (-3697.016922) evaluations: (2667) 
#> Path [16] :Initial log joint density = -481344.194370 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      7.264e-03   2.074e-01    5.229e-01  5.229e-01      2534 -3.700e+03 -3.700e+03                   
#> Path [16] :Best Iter: [47] ELBO (-3699.620552) evaluations: (2534) 
#> Path [17] :Initial log joint density = -481410.327093 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.060e-02   1.663e-01    1.000e+00  1.000e+00      2622 -3.700e+03 -3.705e+03                   
#> Path [17] :Best Iter: [39] ELBO (-3700.427039) evaluations: (2622) 
#> Path [18] :Initial log joint density = -481197.170760 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.383e-03   1.683e-01    7.920e-01  7.920e-01      3240 -3.695e+03 -3.697e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3694.718951) evaluations: (3240) 
#> Path [19] :Initial log joint density = -481491.608308 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.197e-03   1.218e-01    1.000e+00  1.000e+00      2696 -3.697e+03 -3.705e+03                   
#> Path [19] :Best Iter: [46] ELBO (-3696.548758) evaluations: (2696) 
#> Path [20] :Initial log joint density = -486876.254596 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.212e-02   1.919e-01    1.000e+00  1.000e+00      2621 -3.701e+03 -3.702e+03                   
#> Path [20] :Best Iter: [36] ELBO (-3701.490380) evaluations: (2621) 
#> Path [21] :Initial log joint density = -481303.828341 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      2.885e-02   1.947e-01    1.000e+00  1.000e+00      2338 -3.700e+03 -3.701e+03                   
#> Path [21] :Best Iter: [44] ELBO (-3699.861093) evaluations: (2338) 
#> Path [22] :Initial log joint density = -481336.531249 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.824e-03   1.856e-01    8.961e-01  8.961e-01      2495 -3.700e+03 -3.702e+03                   
#> Path [22] :Best Iter: [34] ELBO (-3700.101296) evaluations: (2495) 
#> Path [23] :Initial log joint density = -485523.540843 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      2.279e-02   2.660e-01    1.000e+00  1.000e+00      2625 -3.701e+03 -3.704e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3700.678132) evaluations: (2625) 
#> Path [24] :Initial log joint density = -482176.463034 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.280e-03   2.592e-01    4.365e-01  1.000e+00      3122 -3.699e+03 -3.702e+03                   
#> Path [24] :Best Iter: [48] ELBO (-3698.876815) evaluations: (3122) 
#> Path [25] :Initial log joint density = -481467.417216 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.020e-02   1.694e-01    1.000e+00  1.000e+00      2592 -3.700e+03 -3.698e+03                   
#> Path [25] :Best Iter: [47] ELBO (-3698.127675) evaluations: (2592) 
#> Path [26] :Initial log joint density = -481315.481786 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.038e-02   1.477e-01    1.000e+00  1.000e+00      3296 -3.695e+03 -3.697e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3694.772840) evaluations: (3296) 
#> Path [27] :Initial log joint density = -481621.750006 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.310e-02   2.262e-01    1.000e+00  1.000e+00      3674 -3.694e+03 -3.695e+03                   
#> Path [27] :Best Iter: [60] ELBO (-3693.906185) evaluations: (3674) 
#> Path [28] :Initial log joint density = -481383.137373 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.061e-02   2.329e-01    1.000e+00  1.000e+00      2998 -3.700e+03 -3.700e+03                   
#> Path [28] :Best Iter: [42] ELBO (-3699.656586) evaluations: (2998) 
#> Path [29] :Initial log joint density = -481118.181688 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      6.706e-03   1.305e-01    1.000e+00  1.000e+00      2432 -3.697e+03 -3.700e+03                   
#> Path [29] :Best Iter: [41] ELBO (-3697.232908) evaluations: (2432) 
#> Path [30] :Initial log joint density = -481305.093639 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.984e-03   1.948e-01    5.059e-01  9.764e-01      2515 -3.698e+03 -3.705e+03                   
#> Path [30] :Best Iter: [44] ELBO (-3698.428968) evaluations: (2515) 
#> Path [31] :Initial log joint density = -481895.570305 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.788e+05      1.371e-02   2.348e-01    1.000e+00  1.000e+00      4300 -3.694e+03 -3.694e+03                   
#> Path [31] :Best Iter: [64] ELBO (-3694.153819) evaluations: (4300) 
#> Path [32] :Initial log joint density = -481736.508736 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.298e-02   3.135e-01    1.000e+00  1.000e+00      3978 -3.695e+03 -3.697e+03                   
#> Path [32] :Best Iter: [60] ELBO (-3694.998239) evaluations: (3978) 
#> Path [33] :Initial log joint density = -481356.971424 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      9.275e-03   8.612e-02    1.000e+00  1.000e+00      2491 -3.701e+03 -3.701e+03                   
#> Path [33] :Best Iter: [46] ELBO (-3700.736933) evaluations: (2491) 
#> Path [34] :Initial log joint density = -481263.422409 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.344e-02   1.569e-01    9.971e-01  9.971e-01      2385 -3.700e+03 -3.699e+03                   
#> Path [34] :Best Iter: [44] ELBO (-3698.665323) evaluations: (2385) 
#> Path [35] :Initial log joint density = -481357.422061 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.319e-02   3.124e-01    1.000e+00  1.000e+00      2520 -3.701e+03 -3.703e+03                   
#> Path [35] :Best Iter: [45] ELBO (-3701.077360) evaluations: (2520) 
#> Path [36] :Initial log joint density = -481530.382233 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      1.414e-02   2.639e-01    1.000e+00  1.000e+00      3983 -3.695e+03 -3.696e+03                   
#> Path [36] :Best Iter: [64] ELBO (-3694.872790) evaluations: (3983) 
#> Path [37] :Initial log joint density = -481237.819084 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.228e-02   1.309e-01    1.000e+00  1.000e+00      2885 -3.698e+03 -3.701e+03                   
#> Path [37] :Best Iter: [47] ELBO (-3697.513891) evaluations: (2885) 
#> Path [38] :Initial log joint density = -481591.813578 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.788e+05      3.332e-03   2.697e-01    5.918e-01  5.918e-01      4478 -3.696e+03 -3.696e+03                   
#> Path [38] :Best Iter: [67] ELBO (-3695.610174) evaluations: (4478) 
#> Path [39] :Initial log joint density = -481401.833557 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      9.353e-03   1.737e-01    1.000e+00  1.000e+00      2565 -3.699e+03 -3.703e+03                   
#> Path [39] :Best Iter: [39] ELBO (-3698.576818) evaluations: (2565) 
#> Path [40] :Initial log joint density = -481789.546559 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.358e-02   2.098e-01    1.000e+00  1.000e+00      4064 -3.695e+03 -3.697e+03                   
#> Path [40] :Best Iter: [58] ELBO (-3694.588073) evaluations: (4064) 
#> Path [41] :Initial log joint density = -485004.336774 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.950e-02   2.375e-01    1.000e+00  1.000e+00      2952 -3.701e+03 -3.701e+03                   
#> Path [41] :Best Iter: [45] ELBO (-3701.017442) evaluations: (2952) 
#> Path [42] :Initial log joint density = -483745.512304 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.311e-02   2.857e-01    1.000e+00  1.000e+00      3084 -3.698e+03 -3.700e+03                   
#> Path [42] :Best Iter: [49] ELBO (-3698.178926) evaluations: (3084) 
#> Path [43] :Initial log joint density = -481457.001859 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.597e-02   1.938e-01    1.000e+00  1.000e+00      2607 -3.698e+03 -3.700e+03                   
#> Path [43] :Best Iter: [40] ELBO (-3697.826437) evaluations: (2607) 
#> Path [44] :Initial log joint density = -483478.139328 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.172e-02   3.592e-01    1.000e+00  1.000e+00      2911 -3.700e+03 -3.700e+03                   
#> Path [44] :Best Iter: [49] ELBO (-3699.503625) evaluations: (2911) 
#> Path [45] :Initial log joint density = -481398.975762 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      2.360e-02   2.558e-01    7.550e-01  7.550e-01      2264 -3.699e+03 -3.709e+03                   
#> Path [45] :Best Iter: [36] ELBO (-3698.561923) evaluations: (2264) 
#> Path [46] :Initial log joint density = -481413.744157 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      8.802e-03   9.211e-02    1.000e+00  1.000e+00      2412 -3.700e+03 -3.699e+03                   
#> Path [46] :Best Iter: [46] ELBO (-3699.221927) evaluations: (2412) 
#> Path [47] :Initial log joint density = -481647.273117 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.358e-03   1.616e-01    8.814e-01  8.814e-01      3593 -3.695e+03 -3.696e+03                   
#> Path [47] :Best Iter: [59] ELBO (-3694.920018) evaluations: (3593) 
#> Path [48] :Initial log joint density = -481429.166742 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.970e-03   1.495e-01    1.000e+00  1.000e+00      3069 -3.700e+03 -3.698e+03                   
#> Path [48] :Best Iter: [54] ELBO (-3698.238071) evaluations: (3069) 
#> Path [49] :Initial log joint density = -481885.804974 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      9.580e-03   1.640e-01    1.000e+00  1.000e+00      3913 -3.694e+03 -3.696e+03                   
#> Path [49] :Best Iter: [62] ELBO (-3693.716282) evaluations: (3913) 
#> Path [50] :Initial log joint density = -481262.282957 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.212e-02   2.345e-01    7.925e-01  7.925e-01      2618 -3.698e+03 -3.703e+03                   
#> Path [50] :Best Iter: [47] ELBO (-3697.722031) evaluations: (2618) 
#> Finished in  12.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
