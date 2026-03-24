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
#> Path [1] :Initial log joint density = -481474.687108 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.302e-03   2.606e-01    1.000e+00  1.000e+00      2987 -3.707e+03 -3.717e+03                   
#> Path [1] :Best Iter: [46] ELBO (-3706.550778) evaluations: (2987) 
#> Path [2] :Initial log joint density = -483032.030842 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.738e-03   3.577e-01    6.395e-01  6.395e-01      3069 -3.705e+03 -3.713e+03                   
#> Path [2] :Best Iter: [47] ELBO (-3704.904630) evaluations: (3069) 
#> Path [3] :Initial log joint density = -481574.081704 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.001e-03   1.869e-01    8.429e-01  8.429e-01      3403 -3.702e+03 -3.708e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3701.568913) evaluations: (3403) 
#> Path [4] :Initial log joint density = -483624.590945 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.694e-03   2.190e-01    1.000e+00  1.000e+00      3386 -3.705e+03 -3.702e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3701.641753) evaluations: (3386) 
#> Path [5] :Initial log joint density = -481626.755381 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.691e-02   3.403e-01    1.000e+00  1.000e+00      3356 -3.700e+03 -3.702e+03                   
#> Path [5] :Best Iter: [58] ELBO (-3700.228694) evaluations: (3356) 
#> Path [6] :Initial log joint density = -481648.628928 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.038e-02   2.242e-01    1.000e+00  1.000e+00      3137 -3.707e+03 -3.701e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3701.231707) evaluations: (3137) 
#> Path [7] :Initial log joint density = -481620.678814 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.603e-03   2.082e-01    7.555e-01  7.555e-01      3329 -3.699e+03 -3.712e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3698.502110) evaluations: (3329) 
#> Path [8] :Initial log joint density = -482314.615695 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.179e-03   2.312e-01    1.000e+00  1.000e+00      3438 -3.703e+03 -3.701e+03                   
#> Path [8] :Best Iter: [58] ELBO (-3700.975177) evaluations: (3438) 
#> Path [9] :Initial log joint density = -482062.123630 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.403e-03   1.730e-01    8.677e-01  8.677e-01      3220 -3.703e+03 -3.714e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3702.691637) evaluations: (3220) 
#> Path [10] :Initial log joint density = -481785.331376 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.706e-02   2.950e-01    1.000e+00  1.000e+00      2788 -3.706e+03 -3.709e+03                   
#> Path [10] :Best Iter: [50] ELBO (-3706.471417) evaluations: (2788) 
#> Path [11] :Initial log joint density = -481577.022472 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.422e-03   1.672e-01    1.000e+00  1.000e+00      2993 -3.706e+03 -3.716e+03                   
#> Path [11] :Best Iter: [50] ELBO (-3705.513156) evaluations: (2993) 
#> Path [12] :Initial log joint density = -481838.520181 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.958e-03   1.513e-01    8.648e-01  8.648e-01      3416 -3.700e+03 -3.709e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3699.766945) evaluations: (3416) 
#> Path [13] :Initial log joint density = -481505.785887 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.826e-03   2.347e-01    1.000e+00  1.000e+00      2783 -3.707e+03 -3.708e+03                   
#> Path [13] :Best Iter: [49] ELBO (-3706.754839) evaluations: (2783) 
#> Path [14] :Initial log joint density = -481931.668621 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.010e-02   2.608e-01    1.000e+00  1.000e+00      3149 -3.705e+03 -3.706e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3704.624229) evaluations: (3149) 
#> Path [15] :Initial log joint density = -481565.674598 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.790e-03   2.495e-01    1.000e+00  1.000e+00      3367 -3.701e+03 -3.711e+03                   
#> Path [15] :Best Iter: [55] ELBO (-3701.420936) evaluations: (3367) 
#> Path [16] :Initial log joint density = -481407.206403 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.045e-02   2.834e-01    9.640e-01  9.640e-01      2944 -3.707e+03 -3.719e+03                   
#> Path [16] :Best Iter: [51] ELBO (-3706.701579) evaluations: (2944) 
#> Path [17] :Initial log joint density = -481447.599622 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.222e-02   4.304e-01    1.000e+00  1.000e+00      3157 -3.700e+03 -3.714e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3699.666266) evaluations: (3157) 
#> Path [18] :Initial log joint density = -481701.318688 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.224e-02   1.283e-01    1.000e+00  1.000e+00      3302 -3.703e+03 -3.701e+03                   
#> Path [18] :Best Iter: [57] ELBO (-3700.636692) evaluations: (3302) 
#> Path [19] :Initial log joint density = -482535.092188 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.400e-03   2.936e-01    5.135e-01  1.000e+00      3343 -3.703e+03 -3.710e+03                   
#> Path [19] :Best Iter: [55] ELBO (-3703.261429) evaluations: (3343) 
#> Path [20] :Initial log joint density = -481584.819260 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.633e-03   2.615e-01    1.000e+00  1.000e+00      3391 -3.702e+03 -3.702e+03                   
#> Path [20] :Best Iter: [60] ELBO (-3701.586012) evaluations: (3391) 
#> Path [21] :Initial log joint density = -482624.230943 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.498e-03   1.556e-01    1.000e+00  1.000e+00      3382 -3.705e+03 -3.708e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3705.013251) evaluations: (3382) 
#> Path [22] :Initial log joint density = -482819.444027 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.210e-03   1.960e-01    9.069e-01  9.069e-01      3592 -3.701e+03 -3.710e+03                   
#> Path [22] :Best Iter: [59] ELBO (-3700.528395) evaluations: (3592) 
#> Path [23] :Initial log joint density = -481660.815843 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.139e-03   2.338e-01    1.000e+00  1.000e+00      2748 -3.706e+03 -3.712e+03                   
#> Path [23] :Best Iter: [47] ELBO (-3705.981993) evaluations: (2748) 
#> Path [24] :Initial log joint density = -482091.411052 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      6.922e-03   2.241e-01    1.000e+00  1.000e+00      3744 -3.700e+03 -3.707e+03                   
#> Path [24] :Best Iter: [60] ELBO (-3699.502830) evaluations: (3744) 
#> Path [25] :Initial log joint density = -481961.146370 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.835e-03   2.117e-01    9.269e-01  9.269e-01      3257 -3.700e+03 -3.709e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3699.916739) evaluations: (3257) 
#> Path [26] :Initial log joint density = -482112.057980 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.081e-02   3.412e-01    1.000e+00  1.000e+00      3560 -3.698e+03 -3.703e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3697.534909) evaluations: (3560) 
#> Path [27] :Initial log joint density = -481378.391622 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.436e-03   2.226e-01    1.000e+00  1.000e+00      3405 -3.704e+03 -3.703e+03                   
#> Path [27] :Best Iter: [59] ELBO (-3703.439378) evaluations: (3405) 
#> Path [28] :Initial log joint density = -481866.862562 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.747e-03   2.103e-01    1.000e+00  1.000e+00      3185 -3.705e+03 -3.707e+03                   
#> Path [28] :Best Iter: [54] ELBO (-3705.359868) evaluations: (3185) 
#> Path [29] :Initial log joint density = -482156.737821 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.797e-03   2.354e-01    1.000e+00  1.000e+00      3363 -3.703e+03 -3.704e+03                   
#> Path [29] :Best Iter: [56] ELBO (-3702.921946) evaluations: (3363) 
#> Path [30] :Initial log joint density = -481799.400440 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.171e-03   4.036e-01    1.000e+00  1.000e+00      3288 -3.699e+03 -3.710e+03                   
#> Path [30] :Best Iter: [56] ELBO (-3699.014799) evaluations: (3288) 
#> Path [31] :Initial log joint density = -481644.907273 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.282e-02   2.014e-01    1.000e+00  1.000e+00      3192 -3.701e+03 -3.706e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3700.860603) evaluations: (3192) 
#> Path [32] :Initial log joint density = -481702.909847 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.838e-03   2.258e-01    1.000e+00  1.000e+00      2831 -3.708e+03 -3.711e+03                   
#> Path [32] :Best Iter: [38] ELBO (-3707.862580) evaluations: (2831) 
#> Path [33] :Initial log joint density = -481489.691187 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.081e-02   2.022e-01    1.000e+00  1.000e+00      3029 -3.708e+03 -3.702e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3702.356931) evaluations: (3029) 
#> Path [34] :Initial log joint density = -481843.274119 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.630e-02   3.509e-01    1.000e+00  1.000e+00      3531 -3.700e+03 -3.704e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3700.441077) evaluations: (3531) 
#> Path [35] :Initial log joint density = -481816.050483 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.207e-02   3.096e-01    1.000e+00  1.000e+00      3489 -3.699e+03 -3.707e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3699.154691) evaluations: (3489) 
#> Path [36] :Initial log joint density = -481865.707859 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.787e-03   2.699e-01    1.000e+00  1.000e+00      3136 -3.706e+03 -3.705e+03                   
#> Path [36] :Best Iter: [57] ELBO (-3704.712162) evaluations: (3136) 
#> Path [37] :Initial log joint density = -482011.400149 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.179e-02   3.393e-01    1.000e+00  1.000e+00      3248 -3.701e+03 -3.711e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3700.542289) evaluations: (3248) 
#> Path [38] :Initial log joint density = -481799.835106 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.323e-02   2.583e-01    9.709e-01  9.709e-01      3794 -3.700e+03 -3.705e+03                   
#> Path [38] :Best Iter: [62] ELBO (-3699.698401) evaluations: (3794) 
#> Path [39] :Initial log joint density = -481897.678420 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.002e-03   1.845e-01    7.874e-01  7.874e-01      3620 -3.699e+03 -3.711e+03                   
#> Path [39] :Best Iter: [59] ELBO (-3698.585167) evaluations: (3620) 
#> Path [40] :Initial log joint density = -482758.547467 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.192e-02   3.057e-01    1.000e+00  1.000e+00      3471 -3.699e+03 -3.705e+03                   
#> Path [40] :Best Iter: [58] ELBO (-3698.789389) evaluations: (3471) 
#> Path [41] :Initial log joint density = -481763.095700 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.103e-02   2.349e-01    1.000e+00  1.000e+00      3597 -3.695e+03 -3.701e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3695.074497) evaluations: (3597) 
#> Path [42] :Initial log joint density = -481630.894666 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.742e-03   2.102e-01    4.366e-01  1.000e+00      3191 -3.701e+03 -3.714e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3700.767782) evaluations: (3191) 
#> Path [43] :Initial log joint density = -481508.925573 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.207e-02   3.418e-01    1.000e+00  1.000e+00      2926 -3.708e+03 -3.709e+03                   
#> Path [43] :Best Iter: [42] ELBO (-3707.965904) evaluations: (2926) 
#> Path [44] :Initial log joint density = -481383.137932 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.147e-03   2.745e-01    1.000e+00  1.000e+00      2870 -3.706e+03 -3.719e+03                   
#> Path [44] :Best Iter: [40] ELBO (-3705.775300) evaluations: (2870) 
#> Path [45] :Initial log joint density = -481585.386590 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.972e-03   2.042e-01    1.000e+00  1.000e+00      3242 -3.704e+03 -3.703e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3703.364734) evaluations: (3242) 
#> Path [46] :Initial log joint density = -482352.320709 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.696e-03   2.514e-01    8.760e-01  8.760e-01      3327 -3.704e+03 -3.710e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3704.141719) evaluations: (3327) 
#> Path [47] :Initial log joint density = -481680.763159 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.823e-03   2.265e-01    8.814e-01  8.814e-01      3517 -3.700e+03 -3.708e+03                   
#> Path [47] :Best Iter: [56] ELBO (-3699.644224) evaluations: (3517) 
#> Path [48] :Initial log joint density = -483358.619094 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.609e-03   2.417e-01    5.207e-01  5.207e-01      3415 -3.701e+03 -3.705e+03                   
#> Path [48] :Best Iter: [56] ELBO (-3700.733899) evaluations: (3415) 
#> Path [49] :Initial log joint density = -481705.019570 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.311e-03   2.067e-01    8.604e-01  8.604e-01      3017 -3.707e+03 -3.713e+03                   
#> Path [49] :Best Iter: [53] ELBO (-3707.482759) evaluations: (3017) 
#> Path [50] :Initial log joint density = -481577.041417 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.949e-03   2.449e-01    8.510e-01  8.510e-01      2941 -3.710e+03 -3.720e+03                   
#> Path [50] :Best Iter: [51] ELBO (-3709.887273) evaluations: (2941) 
#> Finished in  13.6 seconds.
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
