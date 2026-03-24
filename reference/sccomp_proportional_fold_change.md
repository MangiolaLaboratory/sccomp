# Calculate Proportional Fold Change from sccomp Estimated Effects

This function calculates the proportional fold change between two
conditions using the estimated effects from a sccomp model. The fold
changes are derived from the model's posterior predictions rather than
raw counts, providing a more robust estimate that accounts for the
model's uncertainty and covariate effects.

Note! This statistic is descriptive and should not be used to define
significance - use sccomp_test() for that. While fold changes in
proportions are easier to interpret than changes in logit space, they
are not linear (the same proportional change has different meaning for
rare vs abundant cell types). In contrast, the logit scale used
internally by sccomp provides linear effects that are more appropriate
for statistical inference.

## Usage

``` r
sccomp_proportional_fold_change(.data, formula_composition, from, to)
```

## Arguments

- .data:

  A sccomp estimate object (of class 'sccomp_tbl') obtained from running
  sccomp_estimate(). This object contains the fitted model and estimated
  effects.

- formula_composition:

  The formula specifying which model effects to use for calculating fold
  changes. This should match or be a subset of the formula used in the
  original sccomp_estimate() call.

- from:

  Character string specifying the reference/control condition (e.g.,
  "benign").

- to:

  Character string specifying the comparison condition (e.g., "cancer").

## Value

A tibble with the following columns:

- cell_group - The cell group identifier

- proportion_fold_change - The estimated fold change in proportions
  between conditions. Positive values indicate increases, negative
  values indicate decreases.

- average_uncertainty - The average uncertainty in the fold change
  estimate, derived from the credible intervals

- statement - A text description of the fold change, including the
  direction and the estimated proportions

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))


# \donttest{
if (instantiate::stan_cmdstan_exists()) {
  # Load example data
  data("counts_obj")

  # First estimate the composition effects
  estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
  )
 
  # Calculate proportional fold changes from the estimated effects
  estimate |> 
  sccomp_proportional_fold_change(
    formula_composition = ~  type, 
    from = "benign", 
    to = "cancer"
  ) 
}
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481713.119262 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.432e-02   3.232e-01    1.000e+00  1.000e+00      3050 -3.691e+03 -3.695e+03                   
#> Path [1] :Best Iter: [37] ELBO (-3691.406733) evaluations: (3050) 
#> Path [2] :Initial log joint density = -484112.904491 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      9.037e-03   1.714e-01    9.083e-01  9.083e-01      3840 -3.681e+03 -3.688e+03                   
#> Path [2] :Best Iter: [61] ELBO (-3680.879612) evaluations: (3840) 
#> Path [3] :Initial log joint density = -481931.858749 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.656e-03   2.320e-01    1.000e+00  1.000e+00      3373 -3.684e+03 -3.684e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3683.567078) evaluations: (3373) 
#> Path [4] :Initial log joint density = -481390.497564 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.368e-02   2.192e-01    1.000e+00  1.000e+00      3379 -3.680e+03 -3.685e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3679.650069) evaluations: (3379) 
#> Path [5] :Initial log joint density = -481625.660098 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.718e-03   2.633e-01    1.000e+00  1.000e+00      3105 -3.692e+03 -3.687e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3686.891355) evaluations: (3105) 
#> Path [6] :Initial log joint density = -481474.757013 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.047e-03   2.530e-01    8.675e-01  8.675e-01      2890 -3.693e+03 -3.714e+03                   
#> Path [6] :Best Iter: [52] ELBO (-3692.625485) evaluations: (2890) 
#> Path [7] :Initial log joint density = -481546.824663 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.051e-03   2.089e-01    9.565e-01  9.565e-01      3582 -3.682e+03 -3.691e+03                   
#> Path [7] :Best Iter: [59] ELBO (-3682.341199) evaluations: (3582) 
#> Path [8] :Initial log joint density = -481724.350590 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.049e-02   2.852e-01    1.000e+00  1.000e+00      3343 -3.681e+03 -3.688e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3681.487490) evaluations: (3343) 
#> Path [9] :Initial log joint density = -481687.868872 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.356e-02   2.608e-01    1.000e+00  1.000e+00      2796 -3.689e+03 -3.692e+03                   
#> Path [9] :Best Iter: [47] ELBO (-3688.990083) evaluations: (2796) 
#> Path [10] :Initial log joint density = -481904.061169 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.579e-03   2.268e-01    1.000e+00  1.000e+00      3391 -3.684e+03 -3.687e+03                   
#> Path [10] :Best Iter: [59] ELBO (-3683.864874) evaluations: (3391) 
#> Path [11] :Initial log joint density = -481845.720101 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.959e-03   2.187e-01    7.514e-01  7.514e-01      3026 -3.689e+03 -3.698e+03                   
#> Path [11] :Best Iter: [53] ELBO (-3689.178020) evaluations: (3026) 
#> Path [12] :Initial log joint density = -481656.507610 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.557e-03   2.392e-01    9.149e-01  9.149e-01      2890 -3.686e+03 -3.698e+03                   
#> Path [12] :Best Iter: [52] ELBO (-3686.079207) evaluations: (2890) 
#> Path [13] :Initial log joint density = -481473.222749 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.023e-03   1.754e-01    1.000e+00  1.000e+00      3156 -3.692e+03 -3.692e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3691.755074) evaluations: (3156) 
#> Path [14] :Initial log joint density = -481655.716999 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.166e-03   2.014e-01    1.000e+00  1.000e+00      3396 -3.684e+03 -3.694e+03                   
#> Path [14] :Best Iter: [56] ELBO (-3683.771029) evaluations: (3396) 
#> Path [15] :Initial log joint density = -482977.226193 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.027e-03   2.620e-01    1.000e+00  1.000e+00      3605 -3.683e+03 -3.694e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3682.708356) evaluations: (3605) 
#> Path [16] :Initial log joint density = -481615.748475 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.616e-03   2.796e-01    1.000e+00  1.000e+00      3241 -3.690e+03 -3.692e+03                   
#> Path [16] :Best Iter: [54] ELBO (-3689.853311) evaluations: (3241) 
#> Path [17] :Initial log joint density = -481693.967367 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.266e-02   3.594e-01    1.000e+00  1.000e+00      3422 -3.684e+03 -3.695e+03                   
#> Path [17] :Best Iter: [58] ELBO (-3684.205255) evaluations: (3422) 
#> Path [18] :Initial log joint density = -481644.281513 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.588e-02   4.154e-01    9.545e-01  9.545e-01      3096 -3.684e+03 -3.695e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3684.023660) evaluations: (3096) 
#> Path [19] :Initial log joint density = -481811.168877 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.206e-02   2.096e-01    1.000e+00  1.000e+00      3633 -3.687e+03 -3.686e+03                   
#> Path [19] :Best Iter: [62] ELBO (-3685.864858) evaluations: (3633) 
#> Path [20] :Initial log joint density = -481975.970484 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.999e-03   2.916e-01    8.363e-01  8.363e-01      3363 -3.686e+03 -3.694e+03                   
#> Path [20] :Best Iter: [56] ELBO (-3685.791613) evaluations: (3363) 
#> Path [21] :Initial log joint density = -481971.133214 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.315e-03   2.136e-01    1.000e+00  1.000e+00      2956 -3.686e+03 -3.693e+03                   
#> Path [21] :Best Iter: [50] ELBO (-3686.460825) evaluations: (2956) 
#> Path [22] :Initial log joint density = -481644.760932 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.173e-02   2.775e-01    1.000e+00  1.000e+00      3277 -3.686e+03 -3.691e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3685.899847) evaluations: (3277) 
#> Path [23] :Initial log joint density = -481673.617389 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.132e-03   2.114e-01    1.000e+00  1.000e+00      3437 -3.685e+03 -3.691e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3685.328782) evaluations: (3437) 
#> Path [24] :Initial log joint density = -481797.017186 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.396e-02   2.558e-01    8.972e-01  8.972e-01      3549 -3.683e+03 -3.689e+03                   
#> Path [24] :Best Iter: [57] ELBO (-3683.428770) evaluations: (3549) 
#> Path [25] :Initial log joint density = -481625.611766 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.764e-03   2.239e-01    1.000e+00  1.000e+00      3428 -3.679e+03 -3.690e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3679.112213) evaluations: (3428) 
#> Path [26] :Initial log joint density = -481568.340551 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.739e-03   2.150e-01    1.000e+00  1.000e+00      3186 -3.690e+03 -3.695e+03                   
#> Path [26] :Best Iter: [49] ELBO (-3690.377981) evaluations: (3186) 
#> Path [27] :Initial log joint density = -481459.558344 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      4.134e-03   2.197e-01    7.763e-01  7.763e-01      2704 -3.691e+03 -3.713e+03                   
#> Path [27] :Best Iter: [48] ELBO (-3691.125509) evaluations: (2704) 
#> Path [28] :Initial log joint density = -482397.035136 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      5.665e-03   3.372e-01    6.197e-01  6.197e-01      3751 -3.683e+03 -3.691e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3682.662579) evaluations: (3751) 
#> Path [29] :Initial log joint density = -481997.574898 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.053e-02   2.385e-01    1.000e+00  1.000e+00      3751 -3.682e+03 -3.682e+03                   
#> Path [29] :Best Iter: [63] ELBO (-3681.711498) evaluations: (3751) 
#> Path [30] :Initial log joint density = -484617.935750 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.829e-03   2.441e-01    1.000e+00  1.000e+00      3653 -3.685e+03 -3.687e+03                   
#> Path [30] :Best Iter: [58] ELBO (-3685.209725) evaluations: (3653) 
#> Path [31] :Initial log joint density = -481643.151956 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.239e-02   2.146e-01    9.717e-01  9.717e-01      3278 -3.682e+03 -3.692e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3681.518328) evaluations: (3278) 
#> Path [32] :Initial log joint density = -481854.024498 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.556e-02   1.902e-01    1.000e+00  1.000e+00      3216 -3.685e+03 -3.686e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3685.098166) evaluations: (3216) 
#> Path [33] :Initial log joint density = -481658.591859 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.887e-02   3.118e-01    1.000e+00  1.000e+00      3047 -3.687e+03 -3.698e+03                   
#> Path [33] :Best Iter: [52] ELBO (-3686.995006) evaluations: (3047) 
#> Path [34] :Initial log joint density = -481311.388049 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.807e-03   2.075e-01    1.000e+00  1.000e+00      2890 -3.691e+03 -3.695e+03                   
#> Path [34] :Best Iter: [43] ELBO (-3690.645522) evaluations: (2890) 
#> Path [35] :Initial log joint density = -481543.874104 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.974e-03   2.523e-01    1.000e+00  1.000e+00      2833 -3.692e+03 -3.694e+03                   
#> Path [35] :Best Iter: [44] ELBO (-3691.893515) evaluations: (2833) 
#> Path [36] :Initial log joint density = -481374.409409 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      7.249e-03   2.277e-01    1.000e+00  1.000e+00      2749 -3.689e+03 -3.694e+03                   
#> Path [36] :Best Iter: [47] ELBO (-3689.481200) evaluations: (2749) 
#> Path [37] :Initial log joint density = -481484.728746 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.552e-02   3.512e-01    1.000e+00  1.000e+00      3302 -3.684e+03 -3.686e+03                   
#> Path [37] :Best Iter: [57] ELBO (-3684.216286) evaluations: (3302) 
#> Path [38] :Initial log joint density = -482200.940492 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.179e-02   2.006e-01    1.000e+00  1.000e+00      3453 -3.683e+03 -3.685e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3682.984902) evaluations: (3453) 
#> Path [39] :Initial log joint density = -481260.632701 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.075e-03   2.088e-01    1.000e+00  1.000e+00      2877 -3.693e+03 -3.694e+03                   
#> Path [39] :Best Iter: [47] ELBO (-3693.330937) evaluations: (2877) 
#> Path [40] :Initial log joint density = -481938.082448 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.078e-03   3.005e-01    9.401e-01  9.401e-01      3391 -3.682e+03 -3.694e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3681.985840) evaluations: (3391) 
#> Path [41] :Initial log joint density = -482457.258874 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.530e-02   3.791e-01    1.000e+00  1.000e+00      3339 -3.684e+03 -3.690e+03                   
#> Path [41] :Best Iter: [56] ELBO (-3684.213852) evaluations: (3339) 
#> Path [42] :Initial log joint density = -481694.594644 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.545e-02   2.576e-01    1.000e+00  1.000e+00      3242 -3.679e+03 -3.685e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3679.407980) evaluations: (3242) 
#> Path [43] :Initial log joint density = -481623.924788 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.027e-02   2.412e-01    1.000e+00  1.000e+00      3742 -3.685e+03 -3.684e+03                   
#> Path [43] :Best Iter: [62] ELBO (-3683.589670) evaluations: (3742) 
#> Path [44] :Initial log joint density = -481745.860016 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.480e-02   3.283e-01    4.142e-01  1.000e+00      3734 -3.680e+03 -3.691e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3679.980206) evaluations: (3734) 
#> Path [45] :Initial log joint density = -484606.476085 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.978e-03   1.839e-01    9.329e-01  9.329e-01      3440 -3.683e+03 -3.687e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3683.450685) evaluations: (3440) 
#> Path [46] :Initial log joint density = -481609.124429 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.747e-03   2.758e-01    7.949e-01  7.949e-01      3061 -3.688e+03 -3.702e+03                   
#> Path [46] :Best Iter: [52] ELBO (-3687.807413) evaluations: (3061) 
#> Path [47] :Initial log joint density = -482354.276597 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.619e-03   2.199e-01    1.000e+00  1.000e+00      3440 -3.683e+03 -3.688e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3683.365370) evaluations: (3440) 
#> Path [48] :Initial log joint density = -481762.094639 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.901e-03   2.238e-01    1.000e+00  1.000e+00      3410 -3.686e+03 -3.686e+03                   
#> Path [48] :Best Iter: [59] ELBO (-3685.589706) evaluations: (3410) 
#> Path [49] :Initial log joint density = -481555.175152 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.257e-02   2.101e-01    1.000e+00  1.000e+00      3695 -3.680e+03 -3.686e+03                   
#> Path [49] :Best Iter: [58] ELBO (-3680.285806) evaluations: (3695) 
#> Path [50] :Initial log joint density = -483625.381771 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.114e-03   2.056e-01    1.000e+00  1.000e+00      3647 -3.680e+03 -3.682e+03                   
#> Path [50] :Best Iter: [57] ELBO (-3680.485444) evaluations: (3647) 
#> Finished in  13.7 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 4
#>    cell_group proportion_fold_change average_uncertainty statement              
#>    <chr>                       <dbl>               <dbl> <glue>                 
#>  1 B1                          -1.94             0.244   1.9-fold decrease (fro…
#>  2 B2                          -2.05             0.304   2-fold decrease (from …
#>  3 B3                          -1.39             0.155   1.4-fold decrease (fro…
#>  4 BM                          -1.39             0.228   1.4-fold decrease (fro…
#>  5 CD4 1                        1.16             0.00175 1.2-fold increase (fro…
#>  6 CD4 2                        1.43             0.0412  1.4-fold increase (fro…
#>  7 CD4 3                       -2.34             0.393   2.3-fold decrease (fro…
#>  8 CD4 4                       -1.07             0.0226  1.1-fold decrease (fro…
#>  9 CD4 5                        1.04             0.0984  1-fold increase (from …
#> 10 CD8 1                        1.10             0.0265  1.1-fold increase (fro…
#> # ℹ 26 more rows
# }
```
