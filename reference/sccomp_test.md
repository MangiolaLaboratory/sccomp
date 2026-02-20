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
#> Path [1] :Initial log joint density = -481303.530436 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      9.809e-03   1.883e-01    1.000e+00  1.000e+00      2291 -3.711e+03 -3.718e+03                   
#> Path [1] :Best Iter: [43] ELBO (-3710.710309) evaluations: (2291) 
#> Path [2] :Initial log joint density = -481566.115092 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      8.739e-03   2.117e-01    1.000e+00  1.000e+00      2609 -3.711e+03 -3.713e+03                   
#> Path [2] :Best Iter: [47] ELBO (-3711.234545) evaluations: (2609) 
#> Path [3] :Initial log joint density = -481544.454138 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      1.124e-02   2.488e-01    1.000e+00  1.000e+00      2261 -3.712e+03 -3.713e+03                   
#> Path [3] :Best Iter: [42] ELBO (-3712.276695) evaluations: (2261) 
#> Path [4] :Initial log joint density = -481545.658552 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.241e-02   1.763e-01    1.000e+00  1.000e+00      2845 -3.710e+03 -3.713e+03                   
#> Path [4] :Best Iter: [51] ELBO (-3710.395538) evaluations: (2845) 
#> Path [5] :Initial log joint density = -481673.297594 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.038e-02   1.998e-01    1.000e+00  1.000e+00      2776 -3.709e+03 -3.712e+03                   
#> Path [5] :Best Iter: [39] ELBO (-3708.971079) evaluations: (2776) 
#> Path [6] :Initial log joint density = -481504.442156 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      9.437e-03   1.768e-01    9.215e-01  9.215e-01      2410 -3.710e+03 -3.713e+03                   
#> Path [6] :Best Iter: [38] ELBO (-3710.120462) evaluations: (2410) 
#> Path [7] :Initial log joint density = -482389.462721 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.546e-02   2.710e-01    9.756e-01  9.756e-01      2564 -3.712e+03 -3.715e+03                   
#> Path [7] :Best Iter: [45] ELBO (-3711.701259) evaluations: (2564) 
#> Path [8] :Initial log joint density = -482679.552962 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.620e-02   2.303e-01    1.000e+00  1.000e+00      2747 -3.712e+03 -3.731e+03                   
#> Path [8] :Best Iter: [44] ELBO (-3712.026057) evaluations: (2747) 
#> Path [9] :Initial log joint density = -481192.480635 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      1.187e-02   2.548e-01    1.000e+00  1.000e+00      2574 -3.709e+03 -3.711e+03                   
#> Path [9] :Best Iter: [37] ELBO (-3708.716667) evaluations: (2574) 
#> Path [10] :Initial log joint density = -481547.872009 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      8.353e-03   2.079e-01    1.000e+00  1.000e+00      2309 -3.709e+03 -3.712e+03                   
#> Path [10] :Best Iter: [39] ELBO (-3708.868858) evaluations: (2309) 
#> Path [11] :Initial log joint density = -481762.886737 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.092e-02   1.548e-01    1.000e+00  1.000e+00      2753 -3.710e+03 -3.715e+03                   
#> Path [11] :Best Iter: [48] ELBO (-3710.259701) evaluations: (2753) 
#> Path [12] :Initial log joint density = -481640.004775 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.119e-02   2.287e-01    1.000e+00  1.000e+00      2699 -3.711e+03 -3.712e+03                   
#> Path [12] :Best Iter: [50] ELBO (-3711.338625) evaluations: (2699) 
#> Path [13] :Initial log joint density = -481265.547735 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.454e-03   1.731e-01    1.000e+00  1.000e+00      2545 -3.712e+03 -3.715e+03                   
#> Path [13] :Best Iter: [43] ELBO (-3712.212144) evaluations: (2545) 
#> Path [14] :Initial log joint density = -483259.060751 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              45      -4.788e+05      5.792e-03   1.423e-01    7.922e-01  7.922e-01      2374 -3.710e+03 -3.718e+03                   
#> Path [14] :Best Iter: [44] ELBO (-3709.788362) evaluations: (2374) 
#> Path [15] :Initial log joint density = -481336.208675 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      6.643e-03   1.970e-01    8.887e-01  8.887e-01      2398 -3.711e+03 -3.714e+03                   
#> Path [15] :Best Iter: [33] ELBO (-3710.600606) evaluations: (2398) 
#> Path [16] :Initial log joint density = -481405.505811 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      5.336e-03   1.979e-01    8.991e-01  8.991e-01      2608 -3.712e+03 -3.717e+03                   
#> Path [16] :Best Iter: [30] ELBO (-3712.277035) evaluations: (2608) 
#> Path [17] :Initial log joint density = -482076.645066 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      8.471e-03   1.924e-01    1.000e+00  1.000e+00      2619 -3.710e+03 -3.717e+03                   
#> Path [17] :Best Iter: [39] ELBO (-3709.682892) evaluations: (2619) 
#> Path [18] :Initial log joint density = -482205.594986 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.852e-02   3.489e-01    1.000e+00  1.000e+00      2750 -3.711e+03 -3.711e+03                   
#> Path [18] :Best Iter: [35] ELBO (-3710.620356) evaluations: (2750) 
#> Path [19] :Initial log joint density = -483911.134038 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.399e-03   1.998e-01    4.277e-01  1.000e+00      2789 -3.711e+03 -3.715e+03                   
#> Path [19] :Best Iter: [42] ELBO (-3711.120869) evaluations: (2789) 
#> Path [20] :Initial log joint density = -481497.229873 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      2.799e-03   2.544e-01    5.401e-01  5.401e-01      2892 -3.709e+03 -3.714e+03                   
#> Path [20] :Best Iter: [43] ELBO (-3708.644446) evaluations: (2892) 
#> Path [21] :Initial log joint density = -481563.395372 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      7.955e-03   2.008e-01    1.000e+00  1.000e+00      2584 -3.712e+03 -3.717e+03                   
#> Path [21] :Best Iter: [36] ELBO (-3711.785068) evaluations: (2584) 
#> Path [22] :Initial log joint density = -481687.013146 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.075e-03   1.645e-01    1.000e+00  1.000e+00      2825 -3.710e+03 -3.714e+03                   
#> Path [22] :Best Iter: [44] ELBO (-3710.211144) evaluations: (2825) 
#> Path [23] :Initial log joint density = -481903.871266 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.512e-03   1.868e-01    9.959e-01  9.959e-01      2710 -3.712e+03 -3.711e+03                   
#> Path [23] :Best Iter: [51] ELBO (-3711.389995) evaluations: (2710) 
#> Path [24] :Initial log joint density = -481281.586262 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      6.029e-03   1.988e-01    7.676e-01  7.676e-01      2736 -3.711e+03 -3.718e+03                   
#> Path [24] :Best Iter: [38] ELBO (-3710.904851) evaluations: (2736) 
#> Path [25] :Initial log joint density = -481854.442627 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.256e-03   2.291e-01    1.000e+00  1.000e+00      2893 -3.711e+03 -3.712e+03                   
#> Path [25] :Best Iter: [46] ELBO (-3710.742305) evaluations: (2893) 
#> Path [26] :Initial log joint density = -481790.495561 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.268e-02   2.476e-01    1.000e+00  1.000e+00      2674 -3.711e+03 -3.716e+03                   
#> Path [26] :Best Iter: [46] ELBO (-3711.122310) evaluations: (2674) 
#> Path [27] :Initial log joint density = -481393.526342 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.400e-02   3.706e-01    1.000e+00  1.000e+00      2369 -3.710e+03 -3.714e+03                   
#> Path [27] :Best Iter: [34] ELBO (-3709.980558) evaluations: (2369) 
#> Path [28] :Initial log joint density = -481380.382651 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.470e-02   2.428e-01    8.265e-01  8.265e-01      2462 -3.710e+03 -3.715e+03                   
#> Path [28] :Best Iter: [37] ELBO (-3709.965186) evaluations: (2462) 
#> Path [29] :Initial log joint density = -481436.840489 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      1.007e-02   1.900e-01    1.000e+00  1.000e+00      2550 -3.712e+03 -3.713e+03                   
#> Path [29] :Best Iter: [45] ELBO (-3711.624550) evaluations: (2550) 
#> Path [30] :Initial log joint density = -483023.935290 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.186e-03   1.557e-01    1.000e+00  1.000e+00      2936 -3.710e+03 -3.713e+03                   
#> Path [30] :Best Iter: [48] ELBO (-3709.645776) evaluations: (2936) 
#> Path [31] :Initial log joint density = -483940.045618 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.491e-03   1.522e-01    1.000e+00  1.000e+00      2711 -3.710e+03 -3.713e+03                   
#> Path [31] :Best Iter: [43] ELBO (-3710.071356) evaluations: (2711) 
#> Path [32] :Initial log joint density = -481420.701420 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      9.308e-03   2.871e-01    1.000e+00  1.000e+00      2408 -3.711e+03 -3.715e+03                   
#> Path [32] :Best Iter: [43] ELBO (-3710.870258) evaluations: (2408) 
#> Path [33] :Initial log joint density = -481789.072788 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.595e-02   3.247e-01    1.000e+00  1.000e+00      2890 -3.710e+03 -3.713e+03                   
#> Path [33] :Best Iter: [45] ELBO (-3710.047868) evaluations: (2890) 
#> Path [34] :Initial log joint density = -481423.366572 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              42      -4.788e+05      7.808e-03   1.899e-01    1.000e+00  1.000e+00      2071 -3.711e+03 -3.711e+03                   
#> Path [34] :Best Iter: [37] ELBO (-3710.660337) evaluations: (2071) 
#> Path [35] :Initial log joint density = -481418.324533 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.873e-03   1.854e-01    1.000e+00  1.000e+00      2724 -3.710e+03 -3.712e+03                   
#> Path [35] :Best Iter: [45] ELBO (-3709.909431) evaluations: (2724) 
#> Path [36] :Initial log joint density = -482466.485402 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              47      -4.788e+05      6.808e-03   1.042e-01    8.277e-01  8.277e-01      2587 -3.709e+03 -3.713e+03                   
#> Path [36] :Best Iter: [34] ELBO (-3709.485421) evaluations: (2587) 
#> Path [37] :Initial log joint density = -481653.071416 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.444e-02   1.716e-01    1.000e+00  1.000e+00      2773 -3.711e+03 -3.712e+03                   
#> Path [37] :Best Iter: [50] ELBO (-3710.765877) evaluations: (2773) 
#> Path [38] :Initial log joint density = -481466.840303 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.161e-02   1.952e-01    1.000e+00  1.000e+00      2549 -3.711e+03 -3.713e+03                   
#> Path [38] :Best Iter: [39] ELBO (-3710.716774) evaluations: (2549) 
#> Path [39] :Initial log joint density = -484268.077088 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.362e-02   3.291e-01    1.000e+00  1.000e+00      2730 -3.709e+03 -3.713e+03                   
#> Path [39] :Best Iter: [49] ELBO (-3709.209259) evaluations: (2730) 
#> Path [40] :Initial log joint density = -481381.909206 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              44      -4.788e+05      9.861e-03   2.563e-01    9.038e-01  9.038e-01      2253 -3.713e+03 -3.714e+03                   
#> Path [40] :Best Iter: [39] ELBO (-3712.721501) evaluations: (2253) 
#> Path [41] :Initial log joint density = -481561.268138 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.215e-02   2.320e-01    1.000e+00  1.000e+00      2969 -3.709e+03 -3.710e+03                   
#> Path [41] :Best Iter: [50] ELBO (-3709.331225) evaluations: (2969) 
#> Path [42] :Initial log joint density = -481364.751981 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.209e-02   2.419e-01    1.000e+00  1.000e+00      2588 -3.712e+03 -3.711e+03                   
#> Path [42] :Best Iter: [49] ELBO (-3710.691017) evaluations: (2588) 
#> Path [43] :Initial log joint density = -481475.680451 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.607e-02   3.270e-01    1.000e+00  1.000e+00      2873 -3.711e+03 -3.717e+03                   
#> Path [43] :Best Iter: [50] ELBO (-3711.428870) evaluations: (2873) 
#> Path [44] :Initial log joint density = -481440.917209 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              46      -4.788e+05      1.014e-02   1.964e-01    1.000e+00  1.000e+00      2433 -3.712e+03 -3.717e+03                   
#> Path [44] :Best Iter: [33] ELBO (-3711.573676) evaluations: (2433) 
#> Path [45] :Initial log joint density = -481210.885838 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.075e-02   1.768e-01    1.000e+00  1.000e+00      2697 -3.707e+03 -3.712e+03                   
#> Path [45] :Best Iter: [46] ELBO (-3706.555120) evaluations: (2697) 
#> Path [46] :Initial log joint density = -484940.549868 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.004e-03   1.792e-01    1.000e+00  1.000e+00      2929 -3.711e+03 -3.722e+03                   
#> Path [46] :Best Iter: [45] ELBO (-3711.099750) evaluations: (2929) 
#> Path [47] :Initial log joint density = -481418.017771 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      2.035e-02   2.219e-01    1.000e+00  1.000e+00      3059 -3.708e+03 -3.712e+03                   
#> Path [47] :Best Iter: [52] ELBO (-3708.039753) evaluations: (3059) 
#> Path [48] :Initial log joint density = -481515.755614 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.052e-03   1.663e-01    1.000e+00  1.000e+00      2887 -3.709e+03 -3.711e+03                   
#> Path [48] :Best Iter: [49] ELBO (-3709.111802) evaluations: (2887) 
#> Path [49] :Initial log joint density = -481993.181087 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.881e-03   2.140e-01    1.000e+00  1.000e+00      2710 -3.711e+03 -3.712e+03                   
#> Path [49] :Best Iter: [43] ELBO (-3711.377316) evaluations: (2710) 
#> Path [50] :Initial log joint density = -481588.145061 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      4.460e-03   1.943e-01    6.902e-01  6.902e-01      2719 -3.711e+03 -3.720e+03                   
#> Path [50] :Best Iter: [31] ELBO (-3710.997115) evaluations: (2719) 
#> Finished in  11.7 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
