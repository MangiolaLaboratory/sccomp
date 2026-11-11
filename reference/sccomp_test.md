# sccomp_test

This function test contrasts from a sccomp result.

## Usage

``` r
sccomp_test(
  .data,
  contrasts = NULL,
  percent_false_positive = 5,
  test_composition_above_CLR_fold_change = 0.1,
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

- test_composition_above_CLR_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the centered log-ratio
  uncontrained scale.

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
  `test_composition_above_CLR_fold_change` argument.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: typebenign, typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481192.007316 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.787e+05      1.251e-02   3.371e-01    1.000e+00  1.000e+00      2806 -3.720e+03 -3.721e+03                   
#> Path [1] :Best Iter: [46] ELBO (-3720.243301) evaluations: (2806) 
#> Path [2] :Initial log joint density = -483081.806353 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.893e-02   2.442e-01    1.000e+00  1.000e+00      3522 -3.717e+03 -3.718e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3716.642063) evaluations: (3522) 
#> Path [3] :Initial log joint density = -481754.959887 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.272e-02   2.771e-01    5.075e-01  1.000e+00      4233 -3.719e+03 -3.720e+03                   
#> Path [3] :Best Iter: [64] ELBO (-3718.916594) evaluations: (4233) 
#> Path [4] :Initial log joint density = -485707.511888 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.787e+05      1.323e-02   2.382e-01    1.000e+00  1.000e+00      3041 -3.722e+03 -3.722e+03                   
#> Path [4] :Best Iter: [50] ELBO (-3721.993902) evaluations: (3041) 
#> Path [5] :Initial log joint density = -481436.813805 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      7.396e-03   1.902e-01    8.538e-01  8.538e-01      3019 -3.722e+03 -3.722e+03                   
#> Path [5] :Best Iter: [52] ELBO (-3721.590235) evaluations: (3019) 
#> Path [6] :Initial log joint density = -481477.658430 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.340e-02   1.991e-01    1.000e+00  1.000e+00      4203 -3.718e+03 -3.717e+03                   
#> Path [6] :Best Iter: [65] ELBO (-3716.515899) evaluations: (4203) 
#> Path [7] :Initial log joint density = -481346.081258 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.844e-02   9.030e-02    7.636e-01  7.636e-01      4238 -3.718e+03 -3.718e+03                   
#> Path [7] :Best Iter: [67] ELBO (-3718.359641) evaluations: (4238) 
#> Path [8] :Initial log joint density = -481339.966611 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.001e-02   2.063e-01    1.000e+00  1.000e+00      4670 -3.717e+03 -3.718e+03                   
#> Path [8] :Best Iter: [68] ELBO (-3716.805618) evaluations: (4670) 
#> Path [9] :Initial log joint density = -481372.589730 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.513e-02   1.712e-01    1.000e+00  1.000e+00      3098 -3.721e+03 -3.722e+03                   
#> Path [9] :Best Iter: [52] ELBO (-3720.765272) evaluations: (3098) 
#> Path [10] :Initial log joint density = -481510.881740 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.342e-02   2.305e-01    6.224e-01  6.224e-01      4348 -3.717e+03 -3.717e+03                   
#> Path [10] :Best Iter: [65] ELBO (-3717.138432) evaluations: (4348) 
#> Path [11] :Initial log joint density = -481170.269856 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.061e-02   1.461e-01    1.000e+00  1.000e+00      3088 -3.721e+03 -3.722e+03                   
#> Path [11] :Best Iter: [49] ELBO (-3720.771576) evaluations: (3088) 
#> Path [12] :Initial log joint density = -482675.183762 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.370e-02   2.278e-01    1.000e+00  1.000e+00      3549 -3.718e+03 -3.718e+03                   
#> Path [12] :Best Iter: [57] ELBO (-3718.060634) evaluations: (3549) 
#> Path [13] :Initial log joint density = -481076.469160 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.135e-02   7.725e-02    1.000e+00  1.000e+00      3238 -3.720e+03 -3.717e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3717.484489) evaluations: (3238) 
#> Path [14] :Initial log joint density = -481583.087341 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.801e-02   1.873e-01    1.000e+00  1.000e+00      4189 -3.717e+03 -3.718e+03                   
#> Path [14] :Best Iter: [65] ELBO (-3717.250575) evaluations: (4189) 
#> Path [15] :Initial log joint density = -484186.531922 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.787e+05      1.452e-02   2.406e-01    8.900e-01  8.900e-01      3121 -3.721e+03 -3.722e+03                   
#> Path [15] :Best Iter: [46] ELBO (-3720.724554) evaluations: (3121) 
#> Path [16] :Initial log joint density = -481314.398292 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      7.695e-03   1.980e-01    7.946e-01  7.946e-01      2932 -3.721e+03 -3.721e+03                   
#> Path [16] :Best Iter: [52] ELBO (-3720.651560) evaluations: (2932) 
#> Path [17] :Initial log joint density = -481379.490367 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.294e-03   1.062e-01    1.000e+00  1.000e+00      3367 -3.718e+03 -3.717e+03                   
#> Path [17] :Best Iter: [57] ELBO (-3717.492590) evaluations: (3367) 
#> Path [18] :Initial log joint density = -481167.187481 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.526e-02   1.726e-01    1.000e+00  1.000e+00      4120 -3.717e+03 -3.718e+03                   
#> Path [18] :Best Iter: [65] ELBO (-3716.905804) evaluations: (4120) 
#> Path [19] :Initial log joint density = -481458.968105 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      3.860e-03   2.115e-01    3.809e-01  3.809e-01      3465 -3.719e+03 -3.718e+03                   
#> Path [19] :Best Iter: [57] ELBO (-3717.699507) evaluations: (3465) 
#> Path [20] :Initial log joint density = -486845.690434 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      8.501e-03   1.259e-01    4.699e-01  1.000e+00      3446 -3.721e+03 -3.720e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3719.658992) evaluations: (3446) 
#> Path [21] :Initial log joint density = -481272.254415 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.787e+05      2.208e-02   1.916e-01    1.000e+00  1.000e+00      2764 -3.720e+03 -3.721e+03                   
#> Path [21] :Best Iter: [43] ELBO (-3719.802707) evaluations: (2764) 
#> Path [22] :Initial log joint density = -481296.642409 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.044e-02   1.196e-01    1.000e+00  1.000e+00      3151 -3.718e+03 -3.719e+03                   
#> Path [22] :Best Iter: [49] ELBO (-3717.867596) evaluations: (3151) 
#> Path [23] :Initial log joint density = -485496.219010 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      4.539e-03   2.415e-01    3.882e-01  3.882e-01      3666 -3.720e+03 -3.719e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3718.666401) evaluations: (3666) 
#> Path [24] :Initial log joint density = -482150.696941 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.733e-02   1.541e-01    1.000e+00  1.000e+00      3354 -3.720e+03 -3.722e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3720.284182) evaluations: (3354) 
#> Path [25] :Initial log joint density = -481437.551302 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.194e-02   1.639e-01    1.000e+00  1.000e+00      3097 -3.720e+03 -3.724e+03                   
#> Path [25] :Best Iter: [51] ELBO (-3719.713626) evaluations: (3097) 
#> Path [26] :Initial log joint density = -481287.804258 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.460e-02   9.052e-02    1.000e+00  1.000e+00      4219 -3.717e+03 -3.717e+03                   
#> Path [26] :Best Iter: [65] ELBO (-3717.199844) evaluations: (4219) 
#> Path [27] :Initial log joint density = -481587.456470 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.145e-03   1.847e-01    7.295e-01  7.295e-01      4526 -3.717e+03 -3.718e+03                   
#> Path [27] :Best Iter: [65] ELBO (-3716.850766) evaluations: (4526) 
#> Path [28] :Initial log joint density = -481348.115203 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      5.367e-03   1.301e-01    7.361e-01  7.361e-01      3446 -3.716e+03 -3.718e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3716.343377) evaluations: (3446) 
#> Path [29] :Initial log joint density = -481085.950291 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.787e+05      2.162e-02   2.944e-01    1.000e+00  1.000e+00      2853 -3.720e+03 -3.721e+03                   
#> Path [29] :Best Iter: [46] ELBO (-3720.286490) evaluations: (2853) 
#> Path [30] :Initial log joint density = -481273.537475 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.056e-02   7.274e-02    1.000e+00  1.000e+00      3083 -3.720e+03 -3.719e+03                   
#> Path [30] :Best Iter: [54] ELBO (-3718.873787) evaluations: (3083) 
#> Path [31] :Initial log joint density = -481858.170765 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.053e-03   1.340e-01    8.068e-01  8.068e-01      5056 -3.717e+03 -3.718e+03                   
#> Path [31] :Best Iter: [72] ELBO (-3716.943171) evaluations: (5056) 
#> Path [32] :Initial log joint density = -481704.358252 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      9.291e-03   1.910e-01    1.000e+00  1.000e+00      4597 -3.717e+03 -3.718e+03                   
#> Path [32] :Best Iter: [65] ELBO (-3717.250080) evaluations: (4597) 
#> Path [33] :Initial log joint density = -481322.938824 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.787e+05      1.593e-02   2.472e-01    4.514e-01  1.000e+00      2723 -3.722e+03 -3.723e+03                   
#> Path [33] :Best Iter: [45] ELBO (-3721.848949) evaluations: (2723) 
#> Path [34] :Initial log joint density = -481232.689057 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.787e+05      1.290e-02   1.977e-01    1.000e+00  1.000e+00      2956 -3.722e+03 -3.724e+03                   
#> Path [34] :Best Iter: [43] ELBO (-3721.968913) evaluations: (2956) 
#> Path [35] :Initial log joint density = -481325.878293 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.938e-02   2.020e-01    1.000e+00  1.000e+00      3064 -3.718e+03 -3.722e+03                   
#> Path [35] :Best Iter: [53] ELBO (-3718.071540) evaluations: (3064) 
#> Path [36] :Initial log joint density = -481495.628667 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.447e-02   1.333e-01    1.000e+00  1.000e+00      4827 -3.716e+03 -3.717e+03                   
#> Path [36] :Best Iter: [70] ELBO (-3716.379670) evaluations: (4827) 
#> Path [37] :Initial log joint density = -481211.628356 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.437e-03   1.921e-01    7.017e-01  7.017e-01      3943 -3.717e+03 -3.719e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3717.333595) evaluations: (3943) 
#> Path [38] :Initial log joint density = -481564.868767 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.135e-02   2.199e-01    1.000e+00  1.000e+00      4787 -3.718e+03 -3.719e+03                   
#> Path [38] :Best Iter: [68] ELBO (-3718.122965) evaluations: (4787) 
#> Path [39] :Initial log joint density = -481368.454024 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      8.961e-03   9.912e-02    1.000e+00  1.000e+00      4360 -3.717e+03 -3.717e+03                   
#> Path [39] :Best Iter: [67] ELBO (-3717.278643) evaluations: (4360) 
#> Path [40] :Initial log joint density = -481761.580553 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.149e-02   1.530e-01    1.000e+00  1.000e+00      4738 -3.716e+03 -3.717e+03                   
#> Path [40] :Best Iter: [64] ELBO (-3715.636581) evaluations: (4738) 
#> Path [41] :Initial log joint density = -484968.636292 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      7.962e-03   1.797e-01    1.000e+00  1.000e+00      3282 -3.721e+03 -3.719e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3718.678711) evaluations: (3282) 
#> Path [42] :Initial log joint density = -483719.587481 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      9.962e-03   1.618e-01    1.000e+00  1.000e+00      3482 -3.718e+03 -3.718e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3717.813417) evaluations: (3482) 
#> Path [43] :Initial log joint density = -481431.142817 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      1.146e-02   1.673e-01    1.000e+00  1.000e+00      3072 -3.721e+03 -3.721e+03                   
#> Path [43] :Best Iter: [53] ELBO (-3721.395231) evaluations: (3072) 
#> Path [44] :Initial log joint density = -483445.217752 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      8.978e-03   1.649e-01    9.350e-01  9.350e-01      3316 -3.722e+03 -3.724e+03                   
#> Path [44] :Best Iter: [50] ELBO (-3721.691937) evaluations: (3316) 
#> Path [45] :Initial log joint density = -481363.698139 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      3.042e-03   2.260e-01    5.914e-01  5.914e-01      2884 -3.721e+03 -3.725e+03                   
#> Path [45] :Best Iter: [48] ELBO (-3721.100751) evaluations: (2884) 
#> Path [46] :Initial log joint density = -481382.310075 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      4.753e-03   1.530e-01    7.406e-01  7.406e-01      2842 -3.722e+03 -3.725e+03                   
#> Path [46] :Best Iter: [51] ELBO (-3722.453434) evaluations: (2842) 
#> Path [47] :Initial log joint density = -481618.116536 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.482e-02   1.911e-01    1.000e+00  1.000e+00      4198 -3.717e+03 -3.719e+03                   
#> Path [47] :Best Iter: [62] ELBO (-3717.084340) evaluations: (4198) 
#> Path [48] :Initial log joint density = -481398.370753 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.075e-02   2.225e-01    1.000e+00  1.000e+00      3516 -3.717e+03 -3.718e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3717.445836) evaluations: (3516) 
#> Path [49] :Initial log joint density = -481855.885815 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      6.350e-03   1.738e-01    9.541e-01  9.541e-01      4347 -3.718e+03 -3.717e+03                   
#> Path [49] :Best Iter: [67] ELBO (-3716.816912) evaluations: (4347) 
#> Path [50] :Initial log joint density = -481236.639469 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.994e-02   1.223e-01    1.000e+00  1.000e+00      3189 -3.722e+03 -3.723e+03                   
#> Path [50] :Best Iter: [54] ELBO (-3721.836583) evaluations: (3189) 
#> Finished in  12.8 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
