# DEPRECATED: Remove Unwanted Variation from sccomp Estimates

This function is DEPRECATED. Please use
[`sccomp_remove_unwanted_effects`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_remove_unwanted_effects.md)
instead. This function uses the model to remove unwanted variation from
a dataset using the estimates of the model. For example, if you fit your
data with the formula `~ factor_1 + factor_2` and use the formula
`~ factor_1` to remove unwanted variation, the `factor_2` effect will be
factored out.

## Usage

``` r
sccomp_remove_unwanted_variation(
  .data,
  formula_composition_keep = NULL,
  formula_composition = NULL,
  formula_variability = NULL,
  cores = detectCores()
)
```

## Arguments

- .data:

  A tibble. The result of `sccomp_estimate`.

- formula_composition_keep:

  A formula. The formula describing the model for differential
  abundance, for example `~type`. In this case, only the effect of the
  `type` factor will be preserved, while all other factors will be
  factored out.

- formula_composition:

  DEPRECATED. Use `formula_composition_keep` instead.

- formula_variability:

  DEPRECATED. Use `formula_variability_keep` instead.

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
    sccomp_remove_unwanted_variation()
  }
#> Warning: `sccomp_remove_unwanted_variation()` was deprecated in sccomp 1.99.20.
#> ℹ sccomp says: sccomp_remove_unwanted_variation is deprecated. Please use
#>   sccomp_remove_unwanted_effects() instead.
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481765.635924 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.128e-03   2.784e-01    8.518e-01  8.518e-01      3115 -3.703e+03 -3.717e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3703.387542) evaluations: (3115) 
#> Path [2] :Initial log joint density = -481687.186896 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.127e-02   3.558e-01    1.000e+00  1.000e+00      3257 -3.702e+03 -3.705e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3702.414434) evaluations: (3257) 
#> Path [3] :Initial log joint density = -481594.534771 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.433e-03   1.794e-01    7.306e-01  7.306e-01      3274 -3.701e+03 -3.713e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3701.336626) evaluations: (3274) 
#> Path [4] :Initial log joint density = -482536.763458 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.578e-03   3.600e-01    4.160e-01  1.000e+00      3318 -3.704e+03 -3.713e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3704.494041) evaluations: (3318) 
#> Path [5] :Initial log joint density = -482502.190773 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.348e-03   3.178e-01    4.532e-01  1.000e+00      3289 -3.701e+03 -3.712e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3701.141503) evaluations: (3289) 
#> Path [6] :Initial log joint density = -481616.067916 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.248e-03   2.172e-01    1.000e+00  1.000e+00      3053 -3.706e+03 -3.705e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3704.585453) evaluations: (3053) 
#> Path [7] :Initial log joint density = -481572.076490 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.875e-02   3.493e-01    1.000e+00  1.000e+00      3328 -3.699e+03 -3.706e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3699.321772) evaluations: (3328) 
#> Path [8] :Initial log joint density = -485794.574387 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.943e-03   2.202e-01    1.000e+00  1.000e+00      3508 -3.703e+03 -3.708e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3703.311926) evaluations: (3508) 
#> Path [9] :Initial log joint density = -481516.542224 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.667e-02   2.716e-01    1.000e+00  1.000e+00      2912 -3.710e+03 -3.712e+03                   
#> Path [9] :Best Iter: [52] ELBO (-3709.688943) evaluations: (2912) 
#> Path [10] :Initial log joint density = -481621.275578 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.652e-03   2.499e-01    8.466e-01  8.466e-01      3212 -3.701e+03 -3.713e+03                   
#> Path [10] :Best Iter: [55] ELBO (-3700.845160) evaluations: (3212) 
#> Path [11] :Initial log joint density = -481508.361792 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.165e-02   4.491e-01    1.000e+00  1.000e+00      2746 -3.709e+03 -3.722e+03                   
#> Path [11] :Best Iter: [44] ELBO (-3709.293736) evaluations: (2746) 
#> Path [12] :Initial log joint density = -481871.389016 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.893e-03   2.116e-01    1.000e+00  1.000e+00      3182 -3.708e+03 -3.707e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3707.213917) evaluations: (3182) 
#> Path [13] :Initial log joint density = -484357.691314 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.808e-03   2.016e-01    1.000e+00  1.000e+00      3155 -3.706e+03 -3.699e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3699.481548) evaluations: (3155) 
#> Path [14] :Initial log joint density = -481629.367799 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.452e-03   1.791e-01    1.000e+00  1.000e+00      2894 -3.708e+03 -3.705e+03                   
#> Path [14] :Best Iter: [53] ELBO (-3705.150464) evaluations: (2894) 
#> Path [15] :Initial log joint density = -482054.225378 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.336e-02   3.248e-01    1.000e+00  1.000e+00      3243 -3.697e+03 -3.702e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3697.228335) evaluations: (3243) 
#> Path [16] :Initial log joint density = -485457.909861 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.173e-02   3.161e-01    1.000e+00  1.000e+00      3523 -3.699e+03 -3.707e+03                   
#> Path [16] :Best Iter: [58] ELBO (-3698.885659) evaluations: (3523) 
#> Path [17] :Initial log joint density = -482363.277135 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.917e-03   2.340e-01    1.000e+00  1.000e+00      2892 -3.709e+03 -3.711e+03                   
#> Path [17] :Best Iter: [39] ELBO (-3708.924526) evaluations: (2892) 
#> Path [18] :Initial log joint density = -480932.341129 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      4.668e-03   2.739e-01    8.074e-01  8.074e-01      2752 -3.706e+03 -3.721e+03                   
#> Path [18] :Best Iter: [45] ELBO (-3705.807103) evaluations: (2752) 
#> Path [19] :Initial log joint density = -482766.920848 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.055e-02   3.613e-01    9.975e-01  9.975e-01      3300 -3.699e+03 -3.714e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3698.642290) evaluations: (3300) 
#> Path [20] :Initial log joint density = -481571.942977 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.225e-02   2.816e-01    1.000e+00  1.000e+00      3232 -3.699e+03 -3.703e+03                   
#> Path [20] :Best Iter: [56] ELBO (-3698.884059) evaluations: (3232) 
#> Path [21] :Initial log joint density = -482117.833040 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.796e-03   2.075e-01    1.000e+00  1.000e+00      3511 -3.701e+03 -3.702e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3700.905017) evaluations: (3511) 
#> Path [22] :Initial log joint density = -481346.615986 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.153e-03   2.445e-01    7.865e-01  7.865e-01      2918 -3.708e+03 -3.725e+03                   
#> Path [22] :Best Iter: [47] ELBO (-3708.247044) evaluations: (2918) 
#> Path [23] :Initial log joint density = -481552.143388 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.461e-02   1.988e-01    1.000e+00  1.000e+00      3171 -3.706e+03 -3.703e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3703.385387) evaluations: (3171) 
#> Path [24] :Initial log joint density = -481762.992777 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.441e-02   2.320e-01    1.000e+00  1.000e+00      3566 -3.700e+03 -3.698e+03                   
#> Path [24] :Best Iter: [62] ELBO (-3698.071968) evaluations: (3566) 
#> Path [25] :Initial log joint density = -481640.646995 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.864e-02   3.713e-01    1.000e+00  1.000e+00      3542 -3.699e+03 -3.707e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3698.865527) evaluations: (3542) 
#> Path [26] :Initial log joint density = -481537.029577 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.722e-03   2.346e-01    1.000e+00  1.000e+00      3238 -3.698e+03 -3.705e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3697.588792) evaluations: (3238) 
#> Path [27] :Initial log joint density = -482082.330368 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.058e-02   2.445e-01    1.000e+00  1.000e+00      3286 -3.706e+03 -3.702e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3702.403225) evaluations: (3286) 
#> Path [28] :Initial log joint density = -481874.743790 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.075e-02   2.094e-01    1.000e+00  1.000e+00      3136 -3.705e+03 -3.701e+03                   
#> Path [28] :Best Iter: [57] ELBO (-3701.410251) evaluations: (3136) 
#> Path [29] :Initial log joint density = -483427.241035 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.265e-03   2.506e-01    4.773e-01  1.000e+00      3420 -3.698e+03 -3.709e+03                   
#> Path [29] :Best Iter: [58] ELBO (-3698.356165) evaluations: (3420) 
#> Path [30] :Initial log joint density = -481843.883836 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.061e-02   2.671e-01    1.000e+00  1.000e+00      2985 -3.705e+03 -3.715e+03                   
#> Path [30] :Best Iter: [52] ELBO (-3705.493415) evaluations: (2985) 
#> Path [31] :Initial log joint density = -481641.597126 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.019e-02   2.333e-01    1.000e+00  1.000e+00      2913 -3.709e+03 -3.710e+03                   
#> Path [31] :Best Iter: [52] ELBO (-3709.269964) evaluations: (2913) 
#> Path [32] :Initial log joint density = -481866.081271 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.607e-02   4.177e-01    1.000e+00  1.000e+00      3364 -3.701e+03 -3.705e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3700.801914) evaluations: (3364) 
#> Path [33] :Initial log joint density = -481481.567464 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.321e-03   2.723e-01    4.636e-01  1.000e+00      3190 -3.705e+03 -3.712e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3704.895984) evaluations: (3190) 
#> Path [34] :Initial log joint density = -483289.791725 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.886e-02   3.182e-01    1.000e+00  1.000e+00      3072 -3.706e+03 -3.703e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3702.683512) evaluations: (3072) 
#> Path [35] :Initial log joint density = -481830.877049 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.807e-03   3.142e-01    1.000e+00  1.000e+00      3391 -3.703e+03 -3.705e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3702.617607) evaluations: (3391) 
#> Path [36] :Initial log joint density = -481415.496164 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.720e-03   2.132e-01    1.000e+00  1.000e+00      3092 -3.709e+03 -3.699e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3699.262279) evaluations: (3092) 
#> Path [37] :Initial log joint density = -481463.914666 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.939e-03   1.373e-01    1.000e+00  1.000e+00      3244 -3.706e+03 -3.703e+03                   
#> Path [37] :Best Iter: [57] ELBO (-3703.379606) evaluations: (3244) 
#> Path [38] :Initial log joint density = -483015.783797 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.056e-02   1.936e-01    8.403e-01  8.403e-01      3420 -3.700e+03 -3.711e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3699.791203) evaluations: (3420) 
#> Path [39] :Initial log joint density = -481808.876488 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.416e-02   3.156e-01    1.000e+00  1.000e+00      3064 -3.709e+03 -3.703e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3702.531224) evaluations: (3064) 
#> Path [40] :Initial log joint density = -481663.955652 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.258e-02   3.055e-01    9.219e-01  9.219e-01      2626 -3.707e+03 -3.720e+03                   
#> Path [40] :Best Iter: [44] ELBO (-3706.671271) evaluations: (2626) 
#> Path [41] :Initial log joint density = -481573.307452 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.979e-03   2.123e-01    8.575e-01  8.575e-01      3149 -3.706e+03 -3.711e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3705.580433) evaluations: (3149) 
#> Path [42] :Initial log joint density = -481607.927776 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.483e-02   3.316e-01    1.000e+00  1.000e+00      3447 -3.702e+03 -3.703e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3702.111787) evaluations: (3447) 
#> Path [43] :Initial log joint density = -481656.989012 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.741e-03   2.438e-01    1.000e+00  1.000e+00      3476 -3.701e+03 -3.712e+03                   
#> Path [43] :Best Iter: [58] ELBO (-3701.464300) evaluations: (3476) 
#> Path [44] :Initial log joint density = -481420.742833 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.760e-03   2.396e-01    1.000e+00  1.000e+00      3133 -3.703e+03 -3.710e+03                   
#> Path [44] :Best Iter: [55] ELBO (-3703.335016) evaluations: (3133) 
#> Path [45] :Initial log joint density = -484293.801891 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.835e-03   2.637e-01    1.000e+00  1.000e+00      3204 -3.704e+03 -3.706e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3704.186240) evaluations: (3204) 
#> Path [46] :Initial log joint density = -481905.217635 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.446e-03   2.259e-01    9.981e-01  9.981e-01      3478 -3.701e+03 -3.708e+03                   
#> Path [46] :Best Iter: [59] ELBO (-3701.053807) evaluations: (3478) 
#> Path [47] :Initial log joint density = -483843.145240 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.046e-02   1.648e-01    7.911e-01  7.911e-01      3754 -3.701e+03 -3.711e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3700.882420) evaluations: (3754) 
#> Path [48] :Initial log joint density = -481582.853450 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.377e-03   2.957e-01    1.000e+00  1.000e+00      2990 -3.705e+03 -3.727e+03                   
#> Path [48] :Best Iter: [40] ELBO (-3705.179201) evaluations: (2990) 
#> Path [49] :Initial log joint density = -481750.788377 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.957e-03   1.988e-01    9.060e-01  9.060e-01      3158 -3.703e+03 -3.711e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3703.279820) evaluations: (3158) 
#> Path [50] :Initial log joint density = -487105.400458 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.026e-02   2.747e-01    9.953e-01  9.953e-01      3543 -3.702e+03 -3.709e+03                   
#> Path [50] :Best Iter: [57] ELBO (-3701.947628) evaluations: (3543) 
#> Finished in  13.5 seconds.
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
