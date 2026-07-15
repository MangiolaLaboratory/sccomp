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
#> Path [1] :Initial log joint density = -481649.249314 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      1.961e-02   1.600e-01    1.000e+00  1.000e+00      2450 -3.697e+03 -3.699e+03                   
#> Path [1] :Best Iter: [43] ELBO (-3697.409222) evaluations: (2450) 
#> Path [2] :Initial log joint density = -481499.696366 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.788e+05      9.387e-03   1.894e-01    1.000e+00  1.000e+00      4336 -3.696e+03 -3.695e+03                   
#> Path [2] :Best Iter: [66] ELBO (-3695.117599) evaluations: (4336) 
#> Path [3] :Initial log joint density = -481467.482474 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.215e-02   2.244e-01    1.000e+00  1.000e+00      2693 -3.698e+03 -3.701e+03                   
#> Path [3] :Best Iter: [48] ELBO (-3697.957041) evaluations: (2693) 
#> Path [4] :Initial log joint density = -482845.566925 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.481e-03   2.513e-01    1.000e+00  1.000e+00      2931 -3.698e+03 -3.701e+03                   
#> Path [4] :Best Iter: [49] ELBO (-3697.927076) evaluations: (2931) 
#> Path [5] :Initial log joint density = -484455.836607 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.254e-02   2.538e-01    1.000e+00  1.000e+00      2748 -3.700e+03 -3.702e+03                   
#> Path [5] :Best Iter: [46] ELBO (-3700.387663) evaluations: (2748) 
#> Path [6] :Initial log joint density = -481332.035405 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      7.628e-03   2.452e-01    8.864e-01  8.864e-01      2445 -3.699e+03 -3.703e+03                   
#> Path [6] :Best Iter: [42] ELBO (-3698.671179) evaluations: (2445) 
#> Path [7] :Initial log joint density = -482795.520529 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      4.412e-03   1.932e-01    7.544e-01  7.544e-01      2866 -3.701e+03 -3.713e+03                   
#> Path [7] :Best Iter: [47] ELBO (-3701.286629) evaluations: (2866) 
#> Path [8] :Initial log joint density = -481672.199972 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.140e-02   2.213e-01    8.997e-01  8.997e-01      3657 -3.695e+03 -3.701e+03                   
#> Path [8] :Best Iter: [60] ELBO (-3694.957719) evaluations: (3657) 
#> Path [9] :Initial log joint density = -484344.698417 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.516e-02   2.458e-01    1.000e+00  1.000e+00      2583 -3.701e+03 -3.702e+03                   
#> Path [9] :Best Iter: [40] ELBO (-3701.097213) evaluations: (2583) 
#> Path [10] :Initial log joint density = -483083.162162 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.099e-02   1.786e-01    1.000e+00  1.000e+00      3617 -3.692e+03 -3.696e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3692.448245) evaluations: (3617) 
#> Path [11] :Initial log joint density = -481410.737968 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      3.822e-03   2.650e-01    6.340e-01  6.340e-01      2522 -3.700e+03 -3.706e+03                   
#> Path [11] :Best Iter: [46] ELBO (-3700.338122) evaluations: (2522) 
#> Path [12] :Initial log joint density = -482949.150566 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.299e-02   1.819e-01    1.000e+00  1.000e+00      2593 -3.700e+03 -3.700e+03                   
#> Path [12] :Best Iter: [36] ELBO (-3699.621637) evaluations: (2593) 
#> Path [13] :Initial log joint density = -481395.608047 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      8.942e-03   1.363e-01    1.000e+00  1.000e+00      2724 -3.699e+03 -3.701e+03                   
#> Path [13] :Best Iter: [48] ELBO (-3699.472983) evaluations: (2724) 
#> Path [14] :Initial log joint density = -481443.912339 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      5.554e-03   2.285e-01    8.162e-01  8.162e-01      2306 -3.700e+03 -3.703e+03                   
#> Path [14] :Best Iter: [40] ELBO (-3700.434527) evaluations: (2306) 
#> Path [15] :Initial log joint density = -481272.804877 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      8.333e-03   1.873e-01    6.003e-01  6.003e-01      2492 -3.700e+03 -3.703e+03                   
#> Path [15] :Best Iter: [42] ELBO (-3699.945764) evaluations: (2492) 
#> Path [16] :Initial log joint density = -481154.229052 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.220e-02   3.082e-01    9.189e-01  9.189e-01      2563 -3.700e+03 -3.702e+03                   
#> Path [16] :Best Iter: [31] ELBO (-3700.286144) evaluations: (2563) 
#> Path [17] :Initial log joint density = -481795.260240 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.788e+05      5.739e-03   1.259e-01    1.000e+00  1.000e+00      4632 -3.695e+03 -3.696e+03                   
#> Path [17] :Best Iter: [69] ELBO (-3694.906158) evaluations: (4632) 
#> Path [18] :Initial log joint density = -481682.278416 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      1.497e-02   2.245e-01    1.000e+00  1.000e+00      4149 -3.695e+03 -3.694e+03                   
#> Path [18] :Best Iter: [65] ELBO (-3694.286957) evaluations: (4149) 
#> Path [19] :Initial log joint density = -482206.179259 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      8.738e-03   1.897e-01    1.000e+00  1.000e+00      2341 -3.706e+03 -3.711e+03                   
#> Path [19] :Best Iter: [32] ELBO (-3705.908784) evaluations: (2341) 
#> Path [20] :Initial log joint density = -481987.416191 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.031e-03   1.502e-01    1.000e+00  1.000e+00      3621 -3.695e+03 -3.696e+03                   
#> Path [20] :Best Iter: [56] ELBO (-3694.657057) evaluations: (3621) 
#> Path [21] :Initial log joint density = -481425.422244 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.466e-02   2.078e-01    1.000e+00  1.000e+00      2517 -3.700e+03 -3.703e+03                   
#> Path [21] :Best Iter: [44] ELBO (-3700.235086) evaluations: (2517) 
#> Path [22] :Initial log joint density = -481553.807791 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.341e-03   1.498e-01    1.000e+00  1.000e+00      3311 -3.696e+03 -3.696e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3695.925778) evaluations: (3311) 
#> Path [23] :Initial log joint density = -481244.837141 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.726e-02   1.455e-01    1.000e+00  1.000e+00      2593 -3.700e+03 -3.704e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3699.919373) evaluations: (2593) 
#> Path [24] :Initial log joint density = -481632.050541 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.069e-03   1.137e-01    1.000e+00  1.000e+00      3596 -3.696e+03 -3.696e+03                   
#> Path [24] :Best Iter: [57] ELBO (-3695.897306) evaluations: (3596) 
#> Path [25] :Initial log joint density = -481483.974904 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.284e-02   3.289e-01    1.000e+00  1.000e+00      2506 -3.701e+03 -3.705e+03                   
#> Path [25] :Best Iter: [43] ELBO (-3701.354058) evaluations: (2506) 
#> Path [26] :Initial log joint density = -481248.526693 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      5.913e-03   1.443e-01    1.000e+00  1.000e+00      2374 -3.701e+03 -3.702e+03                   
#> Path [26] :Best Iter: [40] ELBO (-3700.705321) evaluations: (2374) 
#> Path [27] :Initial log joint density = -482155.696207 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      9.688e-03   1.955e-01    1.000e+00  1.000e+00      2647 -3.702e+03 -3.701e+03                   
#> Path [27] :Best Iter: [49] ELBO (-3701.308163) evaluations: (2647) 
#> Path [28] :Initial log joint density = -481461.965431 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      9.903e-03   2.731e-01    1.000e+00  1.000e+00      2743 -3.701e+03 -3.705e+03                   
#> Path [28] :Best Iter: [47] ELBO (-3700.930331) evaluations: (2743) 
#> Path [29] :Initial log joint density = -481472.866667 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.568e-02   2.875e-01    1.000e+00  1.000e+00      3098 -3.696e+03 -3.698e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3695.717161) evaluations: (3098) 
#> Path [30] :Initial log joint density = -481262.668540 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.282e-02   2.745e-01    3.729e-01  1.000e+00      2399 -3.701e+03 -3.702e+03                   
#> Path [30] :Best Iter: [45] ELBO (-3701.207095) evaluations: (2399) 
#> Path [31] :Initial log joint density = -482278.768008 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      7.062e-03   1.742e-01    1.000e+00  1.000e+00      2408 -3.703e+03 -3.705e+03                   
#> Path [31] :Best Iter: [35] ELBO (-3703.064657) evaluations: (2408) 
#> Path [32] :Initial log joint density = -481504.943513 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.999e-02   2.061e-01    1.000e+00  1.000e+00      2563 -3.699e+03 -3.701e+03                   
#> Path [32] :Best Iter: [44] ELBO (-3699.451220) evaluations: (2563) 
#> Path [33] :Initial log joint density = -481457.447637 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.440e-02   1.567e-01    1.000e+00  1.000e+00      3801 -3.694e+03 -3.694e+03                   
#> Path [33] :Best Iter: [63] ELBO (-3693.835993) evaluations: (3801) 
#> Path [34] :Initial log joint density = -481292.453367 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.025e-02   2.348e-01    7.680e-01  7.680e-01      2361 -3.702e+03 -3.705e+03                   
#> Path [34] :Best Iter: [33] ELBO (-3702.207341) evaluations: (2361) 
#> Path [35] :Initial log joint density = -481451.029330 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.916e-02   1.929e-01    1.000e+00  1.000e+00      2481 -3.700e+03 -3.703e+03                   
#> Path [35] :Best Iter: [46] ELBO (-3700.491874) evaluations: (2481) 
#> Path [36] :Initial log joint density = -481656.873994 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.788e+05      1.179e-02   1.588e-01    1.000e+00  1.000e+00      4482 -3.695e+03 -3.696e+03                   
#> Path [36] :Best Iter: [68] ELBO (-3694.568561) evaluations: (4482) 
#> Path [37] :Initial log joint density = -481267.543218 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.201e-02   1.694e-01    1.000e+00  1.000e+00      2476 -3.697e+03 -3.698e+03                   
#> Path [37] :Best Iter: [43] ELBO (-3697.395505) evaluations: (2476) 
#> Path [38] :Initial log joint density = -481403.305846 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.058e-02   2.006e-01    1.000e+00  1.000e+00      2589 -3.698e+03 -3.700e+03                   
#> Path [38] :Best Iter: [36] ELBO (-3697.997928) evaluations: (2589) 
#> Path [39] :Initial log joint density = -486961.205236 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      5.441e-03   1.530e-01    1.000e+00  1.000e+00      2345 -3.699e+03 -3.704e+03                   
#> Path [39] :Best Iter: [41] ELBO (-3698.654617) evaluations: (2345) 
#> Path [40] :Initial log joint density = -481385.622422 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.275e-02   1.566e-01    1.000e+00  1.000e+00      2313 -3.700e+03 -3.703e+03                   
#> Path [40] :Best Iter: [42] ELBO (-3700.064791) evaluations: (2313) 
#> Path [41] :Initial log joint density = -483415.030508 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.415e-02   2.055e-01    1.000e+00  1.000e+00      2822 -3.702e+03 -3.701e+03                   
#> Path [41] :Best Iter: [49] ELBO (-3701.168186) evaluations: (2822) 
#> Path [42] :Initial log joint density = -481228.015515 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.580e-02   2.010e-01    1.000e+00  1.000e+00      2786 -3.700e+03 -3.701e+03                   
#> Path [42] :Best Iter: [47] ELBO (-3699.518168) evaluations: (2786) 
#> Path [43] :Initial log joint density = -481301.989475 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      5.242e-03   2.058e-01    6.759e-01  6.759e-01      2294 -3.699e+03 -3.701e+03                   
#> Path [43] :Best Iter: [39] ELBO (-3699.145807) evaluations: (2294) 
#> Path [44] :Initial log joint density = -481370.903520 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.042e-02   1.672e-01    9.440e-01  9.440e-01      2972 -3.699e+03 -3.705e+03                   
#> Path [44] :Best Iter: [50] ELBO (-3698.975701) evaluations: (2972) 
#> Path [45] :Initial log joint density = -481316.392964 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      7.121e-03   1.544e-01    6.886e-01  6.886e-01      2404 -3.697e+03 -3.704e+03                   
#> Path [45] :Best Iter: [44] ELBO (-3696.891977) evaluations: (2404) 
#> Path [46] :Initial log joint density = -485165.467892 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.046e-02   3.816e-01    1.000e+00  1.000e+00      2737 -3.702e+03 -3.707e+03                   
#> Path [46] :Best Iter: [35] ELBO (-3702.174747) evaluations: (2737) 
#> Path [47] :Initial log joint density = -483746.627002 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.316e-02   2.121e-01    1.000e+00  1.000e+00      2790 -3.698e+03 -3.698e+03                   
#> Path [47] :Best Iter: [41] ELBO (-3698.132251) evaluations: (2790) 
#> Path [48] :Initial log joint density = -481159.205780 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      3.693e-03   2.060e-01    7.618e-01  7.618e-01      2536 -3.699e+03 -3.706e+03                   
#> Path [48] :Best Iter: [46] ELBO (-3699.368124) evaluations: (2536) 
#> Path [49] :Initial log joint density = -481261.610785 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      7.211e-03   2.316e-01    1.000e+00  1.000e+00      2489 -3.702e+03 -3.704e+03                   
#> Path [49] :Best Iter: [46] ELBO (-3701.557961) evaluations: (2489) 
#> Path [50] :Initial log joint density = -481116.393883 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      4.406e-02   1.988e-01    1.000e+00  1.000e+00      2403 -3.700e+03 -3.698e+03                   
#> Path [50] :Best Iter: [44] ELBO (-3697.756766) evaluations: (2403) 
#> Finished in  12.0 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
