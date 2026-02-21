# Remove unwanted effects from a sccomp_tbl object

This method removes unwanted effects from a dataset using the model
estimates. For example, if you fit your data with the formula
`~ factor_1 + factor_2` and use the formula `~ factor_1` to remove
unwanted variation, the `factor_2` effect will be factored out.

## Usage

``` r
sccomp_remove_unwanted_effects(
  .data,
  formula_composition_keep = NULL,
  formula_composition = NULL,
  cores = detectCores()
)
```

## Arguments

- .data:

  A `sccomp_tbl` object. The result of `sccomp_estimate`.

- formula_composition_keep:

  A formula. The formula describing the model for differential
  abundance, for example `~type`. In this case, only the effect of the
  `type` factor will be preserved, while all other factors will be
  factored out.

- formula_composition:

  DEPRECATED. Use `formula_composition_keep` instead.

- cores:

  Integer, the number of cores to be used for parallel calculations.

## Value

A tibble (`tbl`) with the following columns:

- **sample** - A character column representing the sample name for which
  data was adjusted.

- **cell_group** - A character column representing the cell group being
  tested.

- **adjusted_proportion** - A numeric column representing the adjusted
  proportion after removing unwanted variation.

- **adjusted_counts** - A numeric column representing the adjusted
  counts after removing unwanted variation.

- **logit_residuals** - A numeric column representing the logit
  residuals calculated after adjustment.

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
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_remove_unwanted_effects()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482669.393161 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.332e-02   3.228e-01    9.490e-01  9.490e-01      3073 -3.707e+03 -3.713e+03                   
#> Path [1] :Best Iter: [49] ELBO (-3707.344667) evaluations: (3073) 
#> Path [2] :Initial log joint density = -481552.415132 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.617e-03   2.165e-01    1.000e+00  1.000e+00      2863 -3.708e+03 -3.708e+03                   
#> Path [2] :Best Iter: [52] ELBO (-3707.908708) evaluations: (2863) 
#> Path [3] :Initial log joint density = -481678.075023 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.351e-02   2.778e-01    1.000e+00  1.000e+00      3495 -3.699e+03 -3.700e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3698.988322) evaluations: (3495) 
#> Path [4] :Initial log joint density = -481355.071676 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.238e-03   1.574e-01    1.000e+00  1.000e+00      3372 -3.700e+03 -3.708e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3700.250550) evaluations: (3372) 
#> Path [5] :Initial log joint density = -481827.993582 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.284e-02   3.398e-01    1.000e+00  1.000e+00      2996 -3.707e+03 -3.724e+03                   
#> Path [5] :Best Iter: [47] ELBO (-3707.000669) evaluations: (2996) 
#> Path [6] :Initial log joint density = -483001.333934 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.061e-03   3.730e-01    4.910e-01  1.000e+00      3344 -3.700e+03 -3.710e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3700.129494) evaluations: (3344) 
#> Path [7] :Initial log joint density = -481490.782579 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.886e-02   3.430e-01    9.937e-01  9.937e-01      3162 -3.697e+03 -3.707e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3696.942511) evaluations: (3162) 
#> Path [8] :Initial log joint density = -481489.040280 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.866e-03   2.879e-01    7.269e-01  7.269e-01      3037 -3.709e+03 -3.726e+03                   
#> Path [8] :Best Iter: [44] ELBO (-3708.885109) evaluations: (3037) 
#> Path [9] :Initial log joint density = -481455.140199 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.635e-03   2.340e-01    1.000e+00  1.000e+00      2842 -3.707e+03 -3.709e+03                   
#> Path [9] :Best Iter: [46] ELBO (-3706.653499) evaluations: (2842) 
#> Path [10] :Initial log joint density = -484199.918640 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.843e-03   1.774e-01    7.740e-01  7.740e-01      3575 -3.698e+03 -3.710e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3698.329796) evaluations: (3575) 
#> Path [11] :Initial log joint density = -481645.794233 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.494e-03   1.882e-01    1.000e+00  1.000e+00      3524 -3.700e+03 -3.701e+03                   
#> Path [11] :Best Iter: [56] ELBO (-3699.963620) evaluations: (3524) 
#> Path [12] :Initial log joint density = -481630.109781 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.673e-02   3.356e-01    1.000e+00  1.000e+00      3079 -3.709e+03 -3.703e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3703.140683) evaluations: (3079) 
#> Path [13] :Initial log joint density = -481999.636630 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.788e-03   1.785e-01    7.785e-01  7.785e-01      3305 -3.700e+03 -3.715e+03                   
#> Path [13] :Best Iter: [57] ELBO (-3699.750638) evaluations: (3305) 
#> Path [14] :Initial log joint density = -481641.087199 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.069e-03   2.105e-01    9.580e-01  9.580e-01      2998 -3.707e+03 -3.716e+03                   
#> Path [14] :Best Iter: [50] ELBO (-3706.718999) evaluations: (2998) 
#> Path [15] :Initial log joint density = -482254.841793 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.889e-03   2.593e-01    1.000e+00  1.000e+00      3394 -3.700e+03 -3.702e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3700.468993) evaluations: (3394) 
#> Path [16] :Initial log joint density = -481497.055254 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.915e-03   2.677e-01    1.000e+00  1.000e+00      3131 -3.707e+03 -3.714e+03                   
#> Path [16] :Best Iter: [45] ELBO (-3707.182107) evaluations: (3131) 
#> Path [17] :Initial log joint density = -481901.044395 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.008e-02   1.821e-01    1.000e+00  1.000e+00      3270 -3.701e+03 -3.703e+03                   
#> Path [17] :Best Iter: [55] ELBO (-3701.139222) evaluations: (3270) 
#> Path [18] :Initial log joint density = -482132.504563 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.015e-03   2.285e-01    1.000e+00  1.000e+00      3174 -3.706e+03 -3.706e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3705.692441) evaluations: (3174) 
#> Path [19] :Initial log joint density = -481391.685368 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.313e-03   3.126e-01    8.639e-01  8.639e-01      2993 -3.707e+03 -3.722e+03                   
#> Path [19] :Best Iter: [47] ELBO (-3707.456615) evaluations: (2993) 
#> Path [20] :Initial log joint density = -481617.243305 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.103e-02   2.279e-01    1.000e+00  1.000e+00      3199 -3.700e+03 -3.700e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3699.538112) evaluations: (3199) 
#> Path [21] :Initial log joint density = -481677.920701 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.652e-03   1.999e-01    1.000e+00  1.000e+00      2794 -3.707e+03 -3.723e+03                   
#> Path [21] :Best Iter: [47] ELBO (-3707.334656) evaluations: (2794) 
#> Path [22] :Initial log joint density = -481761.028332 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.201e-02   1.525e-01    1.000e+00  1.000e+00      3088 -3.708e+03 -3.717e+03                   
#> Path [22] :Best Iter: [38] ELBO (-3707.665095) evaluations: (3088) 
#> Path [23] :Initial log joint density = -483119.232802 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.598e-02   3.278e-01    1.000e+00  1.000e+00      3525 -3.699e+03 -3.702e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3698.504941) evaluations: (3525) 
#> Path [24] :Initial log joint density = -481621.586441 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.837e-02   2.306e-01    1.000e+00  1.000e+00      3515 -3.701e+03 -3.702e+03                   
#> Path [24] :Best Iter: [58] ELBO (-3701.019452) evaluations: (3515) 
#> Path [25] :Initial log joint density = -481723.673284 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.899e-03   2.305e-01    1.000e+00  1.000e+00      3207 -3.705e+03 -3.710e+03                   
#> Path [25] :Best Iter: [48] ELBO (-3705.181854) evaluations: (3207) 
#> Path [26] :Initial log joint density = -481772.093997 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.479e-02   3.118e-01    1.000e+00  1.000e+00      3478 -3.698e+03 -3.701e+03                   
#> Path [26] :Best Iter: [60] ELBO (-3697.607715) evaluations: (3478) 
#> Path [27] :Initial log joint density = -481602.072464 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.267e-03   2.606e-01    1.000e+00  1.000e+00      3101 -3.708e+03 -3.704e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3703.578275) evaluations: (3101) 
#> Path [28] :Initial log joint density = -482197.893221 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.369e-02   2.923e-01    1.000e+00  1.000e+00      3446 -3.699e+03 -3.707e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3699.095461) evaluations: (3446) 
#> Path [29] :Initial log joint density = -482203.853817 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      3.667e-03   3.025e-01    5.550e-01  5.550e-01      3628 -3.699e+03 -3.715e+03                   
#> Path [29] :Best Iter: [59] ELBO (-3698.659369) evaluations: (3628) 
#> Path [30] :Initial log joint density = -481558.003458 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.692e-03   2.368e-01    1.000e+00  1.000e+00      2926 -3.706e+03 -3.713e+03                   
#> Path [30] :Best Iter: [50] ELBO (-3705.714760) evaluations: (2926) 
#> Path [31] :Initial log joint density = -481787.635446 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.031e-02   2.475e-01    8.639e-01  8.639e-01      3517 -3.703e+03 -3.716e+03                   
#> Path [31] :Best Iter: [60] ELBO (-3703.374577) evaluations: (3517) 
#> Path [32] :Initial log joint density = -481761.537049 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.394e-02   2.838e-01    8.963e-01  8.963e-01      3158 -3.704e+03 -3.713e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3704.029984) evaluations: (3158) 
#> Path [33] :Initial log joint density = -481554.694675 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.932e-03   2.132e-01    1.000e+00  1.000e+00      3136 -3.705e+03 -3.702e+03                   
#> Path [33] :Best Iter: [57] ELBO (-3701.993579) evaluations: (3136) 
#> Path [34] :Initial log joint density = -481963.799370 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.377e-03   2.332e-01    1.000e+00  1.000e+00      2984 -3.706e+03 -3.713e+03                   
#> Path [34] :Best Iter: [50] ELBO (-3705.867559) evaluations: (2984) 
#> Path [35] :Initial log joint density = -482960.293528 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.154e-02   2.818e-01    1.000e+00  1.000e+00      3526 -3.697e+03 -3.702e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3697.086207) evaluations: (3526) 
#> Path [36] :Initial log joint density = -481656.756750 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.788e-03   3.025e-01    6.878e-01  6.878e-01      3009 -3.707e+03 -3.712e+03                   
#> Path [36] :Best Iter: [53] ELBO (-3707.039803) evaluations: (3009) 
#> Path [37] :Initial log joint density = -481939.012720 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.886e-03   2.658e-01    3.790e-01  1.000e+00      3194 -3.700e+03 -3.710e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3700.056662) evaluations: (3194) 
#> Path [38] :Initial log joint density = -481822.324303 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.754e-03   2.633e-01    1.000e+00  1.000e+00      2863 -3.708e+03 -3.713e+03                   
#> Path [38] :Best Iter: [39] ELBO (-3708.178629) evaluations: (2863) 
#> Path [39] :Initial log joint density = -481739.106475 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.189e-02   2.290e-01    1.000e+00  1.000e+00      2843 -3.708e+03 -3.714e+03                   
#> Path [39] :Best Iter: [38] ELBO (-3707.636182) evaluations: (2843) 
#> Path [40] :Initial log joint density = -482800.539786 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.300e-03   1.928e-01    1.000e+00  1.000e+00      3274 -3.709e+03 -3.711e+03                   
#> Path [40] :Best Iter: [46] ELBO (-3709.411645) evaluations: (3274) 
#> Path [41] :Initial log joint density = -483385.661844 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.836e-03   2.210e-01    1.000e+00  1.000e+00      3000 -3.707e+03 -3.715e+03                   
#> Path [41] :Best Iter: [50] ELBO (-3707.167832) evaluations: (3000) 
#> Path [42] :Initial log joint density = -483825.465883 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.685e-03   2.274e-01    1.000e+00  1.000e+00      3339 -3.704e+03 -3.701e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3700.840420) evaluations: (3339) 
#> Path [43] :Initial log joint density = -481728.987534 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.318e-02   3.417e-01    1.000e+00  1.000e+00      3192 -3.699e+03 -3.710e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3699.200214) evaluations: (3192) 
#> Path [44] :Initial log joint density = -481467.210701 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.414e-02   2.725e-01    1.000e+00  1.000e+00      3388 -3.702e+03 -3.700e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3699.660046) evaluations: (3388) 
#> Path [45] :Initial log joint density = -481645.531630 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      5.013e-03   2.613e-01    5.816e-01  5.816e-01      3699 -3.702e+03 -3.708e+03                   
#> Path [45] :Best Iter: [61] ELBO (-3701.514624) evaluations: (3699) 
#> Path [46] :Initial log joint density = -482283.050293 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.269e-03   2.756e-01    7.214e-01  7.214e-01      3242 -3.703e+03 -3.710e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3702.814929) evaluations: (3242) 
#> Path [47] :Initial log joint density = -481681.147632 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.959e-03   2.583e-01    8.825e-01  8.825e-01      2841 -3.710e+03 -3.722e+03                   
#> Path [47] :Best Iter: [43] ELBO (-3709.520545) evaluations: (2841) 
#> Path [48] :Initial log joint density = -481507.753719 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.682e-03   2.902e-01    1.000e+00  1.000e+00      2783 -3.711e+03 -3.716e+03                   
#> Path [48] :Best Iter: [48] ELBO (-3710.964563) evaluations: (2783) 
#> Path [49] :Initial log joint density = -481603.454721 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.728e-03   2.482e-01    1.000e+00  1.000e+00      3098 -3.707e+03 -3.710e+03                   
#> Path [49] :Best Iter: [52] ELBO (-3707.484961) evaluations: (3098) 
#> Path [50] :Initial log joint density = -484470.912936 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.285e-03   2.172e-01    1.000e+00  1.000e+00      3230 -3.708e+03 -3.707e+03                   
#> Path [50] :Best Iter: [57] ELBO (-3707.318398) evaluations: (3230) 
#> Finished in  13.8 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
# }
```
