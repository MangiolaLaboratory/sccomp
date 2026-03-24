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
#> Path [1] :Initial log joint density = -481714.173107 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.396e-03   2.266e-01    7.881e-01  7.881e-01      3065 -3.709e+03 -3.721e+03                   
#> Path [1] :Best Iter: [37] ELBO (-3708.715395) evaluations: (3065) 
#> Path [2] :Initial log joint density = -484114.918764 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.495e-03   1.683e-01    1.000e+00  1.000e+00      3473 -3.702e+03 -3.706e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3701.959152) evaluations: (3473) 
#> Path [3] :Initial log joint density = -481931.863895 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.419e-03   2.601e-01    1.000e+00  1.000e+00      3286 -3.701e+03 -3.708e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3701.109324) evaluations: (3286) 
#> Path [4] :Initial log joint density = -481391.006217 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.434e-03   2.093e-01    1.000e+00  1.000e+00      3040 -3.708e+03 -3.714e+03                   
#> Path [4] :Best Iter: [37] ELBO (-3708.247713) evaluations: (3040) 
#> Path [5] :Initial log joint density = -481626.445925 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.609e-03   1.947e-01    8.258e-01  8.258e-01      3624 -3.698e+03 -3.710e+03                   
#> Path [5] :Best Iter: [60] ELBO (-3698.258714) evaluations: (3624) 
#> Path [6] :Initial log joint density = -481474.835386 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.235e-03   1.946e-01    1.000e+00  1.000e+00      3124 -3.703e+03 -3.699e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3699.414142) evaluations: (3124) 
#> Path [7] :Initial log joint density = -481548.013394 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.108e-02   2.093e-01    1.000e+00  1.000e+00      3493 -3.697e+03 -3.702e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3696.656345) evaluations: (3493) 
#> Path [8] :Initial log joint density = -481726.217954 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.901e-03   2.844e-01    7.234e-01  7.234e-01      3479 -3.698e+03 -3.713e+03                   
#> Path [8] :Best Iter: [57] ELBO (-3698.415169) evaluations: (3479) 
#> Path [9] :Initial log joint density = -481688.900100 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.098e-03   2.303e-01    1.000e+00  1.000e+00      2758 -3.708e+03 -3.708e+03                   
#> Path [9] :Best Iter: [42] ELBO (-3708.070011) evaluations: (2758) 
#> Path [10] :Initial log joint density = -481906.883693 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.491e-03   2.213e-01    7.910e-01  7.910e-01      3347 -3.701e+03 -3.710e+03                   
#> Path [10] :Best Iter: [56] ELBO (-3701.194730) evaluations: (3347) 
#> Path [11] :Initial log joint density = -481846.461400 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.248e-03   1.704e-01    1.000e+00  1.000e+00      3026 -3.708e+03 -3.713e+03                   
#> Path [11] :Best Iter: [46] ELBO (-3707.752784) evaluations: (3026) 
#> Path [12] :Initial log joint density = -481658.767616 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.812e-03   1.701e-01    7.548e-01  7.548e-01      3220 -3.703e+03 -3.713e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3703.108212) evaluations: (3220) 
#> Path [13] :Initial log joint density = -481473.787144 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.560e-03   2.767e-01    1.000e+00  1.000e+00      2783 -3.707e+03 -3.717e+03                   
#> Path [13] :Best Iter: [45] ELBO (-3707.154792) evaluations: (2783) 
#> Path [14] :Initial log joint density = -481656.443238 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.347e-02   2.408e-01    1.000e+00  1.000e+00      3312 -3.703e+03 -3.701e+03                   
#> Path [14] :Best Iter: [58] ELBO (-3700.694880) evaluations: (3312) 
#> Path [15] :Initial log joint density = -482977.773738 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.604e-02   3.170e-01    1.000e+00  1.000e+00      3426 -3.703e+03 -3.703e+03                   
#> Path [15] :Best Iter: [58] ELBO (-3702.534859) evaluations: (3426) 
#> Path [16] :Initial log joint density = -481618.429086 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.506e-03   1.776e-01    1.000e+00  1.000e+00      3168 -3.703e+03 -3.707e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3703.440031) evaluations: (3168) 
#> Path [17] :Initial log joint density = -481694.445157 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.323e-03   1.775e-01    1.000e+00  1.000e+00      3598 -3.699e+03 -3.704e+03                   
#> Path [17] :Best Iter: [58] ELBO (-3699.361230) evaluations: (3598) 
#> Path [18] :Initial log joint density = -481646.573372 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.549e-02   3.909e-01    9.097e-01  9.097e-01      3059 -3.702e+03 -3.714e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3702.362117) evaluations: (3059) 
#> Path [19] :Initial log joint density = -481812.336195 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.147e-03   2.210e-01    1.000e+00  1.000e+00      3657 -3.702e+03 -3.700e+03                   
#> Path [19] :Best Iter: [62] ELBO (-3699.728351) evaluations: (3657) 
#> Path [20] :Initial log joint density = -481978.313663 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.872e-03   1.955e-01    1.000e+00  1.000e+00      3218 -3.707e+03 -3.712e+03                   
#> Path [20] :Best Iter: [40] ELBO (-3706.645184) evaluations: (3218) 
#> Path [21] :Initial log joint density = -481972.933075 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.270e-02   2.423e-01    1.000e+00  1.000e+00      3361 -3.699e+03 -3.702e+03                   
#> Path [21] :Best Iter: [56] ELBO (-3699.463348) evaluations: (3361) 
#> Path [22] :Initial log joint density = -481644.767241 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.276e-02   3.239e-01    1.000e+00  1.000e+00      3117 -3.702e+03 -3.704e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3702.339028) evaluations: (3117) 
#> Path [23] :Initial log joint density = -481674.865290 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.301e-03   2.241e-01    1.000e+00  1.000e+00      3165 -3.707e+03 -3.710e+03                   
#> Path [23] :Best Iter: [52] ELBO (-3706.756348) evaluations: (3165) 
#> Path [24] :Initial log joint density = -481798.414537 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.728e-02   2.342e-01    1.000e+00  1.000e+00      3277 -3.702e+03 -3.700e+03                   
#> Path [24] :Best Iter: [56] ELBO (-3700.486050) evaluations: (3277) 
#> Path [25] :Initial log joint density = -481627.537262 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.591e-02   2.564e-01    1.000e+00  1.000e+00      3697 -3.698e+03 -3.701e+03                   
#> Path [25] :Best Iter: [61] ELBO (-3698.153620) evaluations: (3697) 
#> Path [26] :Initial log joint density = -481569.774387 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.125e-02   2.817e-01    1.000e+00  1.000e+00      3357 -3.701e+03 -3.703e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3700.658497) evaluations: (3357) 
#> Path [27] :Initial log joint density = -481460.711691 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.275e-03   2.329e-01    1.000e+00  1.000e+00      2741 -3.707e+03 -3.718e+03                   
#> Path [27] :Best Iter: [48] ELBO (-3706.623768) evaluations: (2741) 
#> Path [28] :Initial log joint density = -482398.833253 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.829e-03   1.972e-01    1.000e+00  1.000e+00      3729 -3.699e+03 -3.707e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3699.045263) evaluations: (3729) 
#> Path [29] :Initial log joint density = -481999.357837 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.529e-03   2.058e-01    1.000e+00  1.000e+00      3451 -3.697e+03 -3.700e+03                   
#> Path [29] :Best Iter: [56] ELBO (-3696.569382) evaluations: (3451) 
#> Path [30] :Initial log joint density = -484619.196827 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.014e-02   2.717e-01    1.000e+00  1.000e+00      3532 -3.702e+03 -3.699e+03                   
#> Path [30] :Best Iter: [60] ELBO (-3699.397487) evaluations: (3532) 
#> Path [31] :Initial log joint density = -481646.422118 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.354e-02   2.790e-01    9.601e-01  9.601e-01      3327 -3.699e+03 -3.712e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3699.350327) evaluations: (3327) 
#> Path [32] :Initial log joint density = -481854.032859 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.914e-03   2.337e-01    9.033e-01  9.033e-01      3128 -3.709e+03 -3.709e+03                   
#> Path [32] :Best Iter: [42] ELBO (-3708.970790) evaluations: (3128) 
#> Path [33] :Initial log joint density = -481660.228656 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.980e-02   2.439e-01    1.000e+00  1.000e+00      3309 -3.699e+03 -3.701e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3698.930320) evaluations: (3309) 
#> Path [34] :Initial log joint density = -481311.733916 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.306e-03   2.132e-01    7.913e-01  7.913e-01      3017 -3.705e+03 -3.709e+03                   
#> Path [34] :Best Iter: [48] ELBO (-3704.817667) evaluations: (3017) 
#> Path [35] :Initial log joint density = -481545.177934 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.507e-02   2.769e-01    1.000e+00  1.000e+00      2922 -3.709e+03 -3.709e+03                   
#> Path [35] :Best Iter: [50] ELBO (-3708.943652) evaluations: (2922) 
#> Path [36] :Initial log joint density = -481374.960072 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.287e-03   1.978e-01    1.000e+00  1.000e+00      2798 -3.707e+03 -3.709e+03                   
#> Path [36] :Best Iter: [47] ELBO (-3707.305482) evaluations: (2798) 
#> Path [37] :Initial log joint density = -481484.907479 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.183e-03   2.742e-01    1.000e+00  1.000e+00      3026 -3.706e+03 -3.704e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3703.736631) evaluations: (3026) 
#> Path [38] :Initial log joint density = -482201.085189 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.945e-03   2.216e-01    1.000e+00  1.000e+00      3493 -3.701e+03 -3.704e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3700.775098) evaluations: (3493) 
#> Path [39] :Initial log joint density = -481262.155689 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.153e-03   3.245e-01    4.688e-01  1.000e+00      2918 -3.708e+03 -3.715e+03                   
#> Path [39] :Best Iter: [41] ELBO (-3708.384001) evaluations: (2918) 
#> Path [40] :Initial log joint density = -481939.071420 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.996e-03   2.804e-01    9.963e-01  9.963e-01      3385 -3.699e+03 -3.706e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3698.796534) evaluations: (3385) 
#> Path [41] :Initial log joint density = -482457.725013 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.409e-03   2.190e-01    1.000e+00  1.000e+00      3339 -3.703e+03 -3.700e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3699.511272) evaluations: (3339) 
#> Path [42] :Initial log joint density = -481696.143196 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.290e-02   2.207e-01    1.000e+00  1.000e+00      3278 -3.696e+03 -3.700e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3696.180338) evaluations: (3278) 
#> Path [43] :Initial log joint density = -481624.289050 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.836e-03   2.081e-01    1.000e+00  1.000e+00      3197 -3.710e+03 -3.710e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3709.549138) evaluations: (3197) 
#> Path [44] :Initial log joint density = -481747.502224 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.908e-03   2.065e-01    1.000e+00  1.000e+00      3165 -3.710e+03 -3.712e+03                   
#> Path [44] :Best Iter: [42] ELBO (-3709.715416) evaluations: (3165) 
#> Path [45] :Initial log joint density = -484608.207402 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.925e-03   1.889e-01    8.602e-01  8.602e-01      3353 -3.703e+03 -3.717e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3702.794160) evaluations: (3353) 
#> Path [46] :Initial log joint density = -481609.474828 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.957e-03   2.907e-01    8.338e-01  8.338e-01      3045 -3.706e+03 -3.721e+03                   
#> Path [46] :Best Iter: [51] ELBO (-3706.097196) evaluations: (3045) 
#> Path [47] :Initial log joint density = -482354.345044 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.433e-02   2.981e-01    1.000e+00  1.000e+00      3327 -3.700e+03 -3.702e+03                   
#> Path [47] :Best Iter: [57] ELBO (-3699.601425) evaluations: (3327) 
#> Path [48] :Initial log joint density = -481763.513025 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.705e-03   2.583e-01    1.000e+00  1.000e+00      3318 -3.705e+03 -3.701e+03                   
#> Path [48] :Best Iter: [57] ELBO (-3700.595986) evaluations: (3318) 
#> Path [49] :Initial log joint density = -481555.970194 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.714e-03   2.487e-01    1.000e+00  1.000e+00      3126 -3.709e+03 -3.708e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3708.116936) evaluations: (3126) 
#> Path [50] :Initial log joint density = -483627.253501 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.606e-02   2.995e-01    1.000e+00  1.000e+00      3207 -3.699e+03 -3.707e+03                   
#> Path [50] :Best Iter: [55] ELBO (-3698.939996) evaluations: (3207) 
#> Finished in  13.5 seconds.
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
#>  1 B1                          -1.96             0.271   2-fold decrease (from …
#>  2 B2                          -2.06             0.234   2.1-fold decrease (fro…
#>  3 B3                          -1.36             0.193   1.4-fold decrease (fro…
#>  4 BM                          -1.39             0.170   1.4-fold decrease (fro…
#>  5 CD4 1                        1.16             0.00675 1.2-fold increase (fro…
#>  6 CD4 2                        1.43             0.0390  1.4-fold increase (fro…
#>  7 CD4 3                       -2.32             0.344   2.3-fold decrease (fro…
#>  8 CD4 4                       -1.07             0.0806  1.1-fold decrease (fro…
#>  9 CD4 5                        1.05             0.0842  1-fold increase (from …
#> 10 CD8 1                        1.10             0.0202  1.1-fold increase (fro…
#> # ℹ 26 more rows
# }
```
