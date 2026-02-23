# Plot 1D Intervals for Cell-group Effects

This function creates a series of 1D interval plots for cell-group
effects, highlighting significant differences based on a given
significance threshold.

## Usage

``` r
plot_1D_intervals(
  .data,
  significance_threshold = 0.05,
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
  show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR")
)
```

## Arguments

- .data:

  Data frame containing the main data.

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences.

- test_composition_above_logit_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the logit uncontrained
  scale.

- show_fdr_message:

  Logical. Whether to show the Bayesian FDR interpretation message on
  the plot. Default is TRUE.

- significance_statistic:

  Character vector indicating which statistic to highlight. Default is
  "pH0".

## Value

A combined plot of 1D interval plots.

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."

# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
    ) |> 
    sccomp_test()
    
  # Example usage:
  my_plot = plot_1D_intervals(estimate)
    
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481754.568746 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.221e-02   2.968e-01    9.889e-01  9.889e-01      2991 -3.704e+03 -3.715e+03                   
#> Path [1] :Best Iter: [51] ELBO (-3704.416173) evaluations: (2991) 
#> Path [2] :Initial log joint density = -481785.846975 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.502e-02   2.313e-01    1.000e+00  1.000e+00      3332 -3.702e+03 -3.700e+03                   
#> Path [2] :Best Iter: [57] ELBO (-3699.890726) evaluations: (3332) 
#> Path [3] :Initial log joint density = -481491.975251 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.817e-03   2.698e-01    7.864e-01  7.864e-01      3122 -3.709e+03 -3.711e+03                   
#> Path [3] :Best Iter: [45] ELBO (-3709.147393) evaluations: (3122) 
#> Path [4] :Initial log joint density = -481830.294515 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.086e-02   2.888e-01    1.000e+00  1.000e+00      3135 -3.708e+03 -3.704e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3704.214864) evaluations: (3135) 
#> Path [5] :Initial log joint density = -481741.991847 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.293e-02   2.075e-01    1.000e+00  1.000e+00      2942 -3.710e+03 -3.708e+03                   
#> Path [5] :Best Iter: [54] ELBO (-3707.831771) evaluations: (2942) 
#> Path [6] :Initial log joint density = -482471.373633 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.710e-03   2.361e-01    1.000e+00  1.000e+00      3544 -3.702e+03 -3.705e+03                   
#> Path [6] :Best Iter: [58] ELBO (-3702.052228) evaluations: (3544) 
#> Path [7] :Initial log joint density = -481965.897185 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.830e-03   2.418e-01    8.937e-01  8.937e-01      3499 -3.700e+03 -3.711e+03                   
#> Path [7] :Best Iter: [58] ELBO (-3699.751546) evaluations: (3499) 
#> Path [8] :Initial log joint density = -481755.415030 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      2.782e-03   2.070e-01    7.133e-01  7.133e-01      2922 -3.708e+03 -3.720e+03                   
#> Path [8] :Best Iter: [46] ELBO (-3707.728772) evaluations: (2922) 
#> Path [9] :Initial log joint density = -481673.952912 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.543e-03   2.682e-01    1.000e+00  1.000e+00      3197 -3.701e+03 -3.701e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3700.611106) evaluations: (3197) 
#> Path [10] :Initial log joint density = -481560.123427 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.473e-02   2.081e-01    1.000e+00  1.000e+00      3333 -3.704e+03 -3.703e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3702.589930) evaluations: (3333) 
#> Path [11] :Initial log joint density = -483488.375822 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      3.900e-03   3.206e-01    6.492e-01  6.492e-01      3346 -3.698e+03 -3.713e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3698.098798) evaluations: (3346) 
#> Path [12] :Initial log joint density = -481897.797231 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.282e-03   1.805e-01    1.000e+00  1.000e+00      3531 -3.699e+03 -3.709e+03                   
#> Path [12] :Best Iter: [59] ELBO (-3699.488629) evaluations: (3531) 
#> Path [13] :Initial log joint density = -481768.406761 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.674e-03   2.944e-01    1.000e+00  1.000e+00      3531 -3.702e+03 -3.706e+03                   
#> Path [13] :Best Iter: [59] ELBO (-3701.794902) evaluations: (3531) 
#> Path [14] :Initial log joint density = -481964.449515 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.503e-02   4.198e-01    1.000e+00  1.000e+00      3429 -3.699e+03 -3.710e+03                   
#> Path [14] :Best Iter: [57] ELBO (-3698.518181) evaluations: (3429) 
#> Path [15] :Initial log joint density = -481574.867409 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.462e-03   2.434e-01    8.940e-01  8.940e-01      2918 -3.709e+03 -3.722e+03                   
#> Path [15] :Best Iter: [49] ELBO (-3708.738281) evaluations: (2918) 
#> Path [16] :Initial log joint density = -485557.292072 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.777e-03   2.188e-01    1.000e+00  1.000e+00      3303 -3.705e+03 -3.701e+03                   
#> Path [16] :Best Iter: [57] ELBO (-3701.183045) evaluations: (3303) 
#> Path [17] :Initial log joint density = -481955.566221 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.715e-03   3.094e-01    7.556e-01  7.556e-01      3140 -3.710e+03 -3.713e+03                   
#> Path [17] :Best Iter: [47] ELBO (-3710.040816) evaluations: (3140) 
#> Path [18] :Initial log joint density = -481774.854036 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.343e-03   1.714e-01    1.000e+00  1.000e+00      3077 -3.710e+03 -3.710e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3710.030897) evaluations: (3077) 
#> Path [19] :Initial log joint density = -481759.861169 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.816e-03   2.247e-01    1.000e+00  1.000e+00      2911 -3.705e+03 -3.710e+03                   
#> Path [19] :Best Iter: [51] ELBO (-3704.852771) evaluations: (2911) 
#> Path [20] :Initial log joint density = -481843.536235 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.945e-03   2.040e-01    1.000e+00  1.000e+00      3084 -3.709e+03 -3.705e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3704.712641) evaluations: (3084) 
#> Path [21] :Initial log joint density = -482796.878201 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.179e-03   1.922e-01    9.577e-01  9.577e-01      3437 -3.701e+03 -3.711e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3700.925643) evaluations: (3437) 
#> Path [22] :Initial log joint density = -481605.077876 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.219e-03   2.943e-01    3.962e-01  1.000e+00      3042 -3.710e+03 -3.716e+03                   
#> Path [22] :Best Iter: [48] ELBO (-3709.967213) evaluations: (3042) 
#> Path [23] :Initial log joint density = -482396.649373 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.189e-02   3.608e-01    1.000e+00  1.000e+00      3538 -3.701e+03 -3.704e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3701.447496) evaluations: (3538) 
#> Path [24] :Initial log joint density = -481645.737650 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.141e-02   1.913e-01    1.000e+00  1.000e+00      3274 -3.699e+03 -3.700e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3699.195064) evaluations: (3274) 
#> Path [25] :Initial log joint density = -481504.576848 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.985e-03   2.963e-01    6.498e-01  6.498e-01      3077 -3.706e+03 -3.712e+03                   
#> Path [25] :Best Iter: [44] ELBO (-3706.110057) evaluations: (3077) 
#> Path [26] :Initial log joint density = -481957.866553 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.009e-02   2.608e-01    1.000e+00  1.000e+00      3136 -3.704e+03 -3.706e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3703.612800) evaluations: (3136) 
#> Path [27] :Initial log joint density = -481590.724357 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.073e-02   1.249e-01    1.000e+00  1.000e+00      3281 -3.704e+03 -3.705e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3703.675036) evaluations: (3281) 
#> Path [28] :Initial log joint density = -481650.839364 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.662e-03   2.790e-01    1.000e+00  1.000e+00      3033 -3.708e+03 -3.716e+03                   
#> Path [28] :Best Iter: [52] ELBO (-3707.859214) evaluations: (3033) 
#> Path [29] :Initial log joint density = -481749.025949 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.011e-02   2.426e-01    1.000e+00  1.000e+00      3561 -3.704e+03 -3.705e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3703.981574) evaluations: (3561) 
#> Path [30] :Initial log joint density = -482384.514328 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.952e-03   2.462e-01    1.000e+00  1.000e+00      3499 -3.701e+03 -3.700e+03                   
#> Path [30] :Best Iter: [60] ELBO (-3699.952920) evaluations: (3499) 
#> Path [31] :Initial log joint density = -481781.166601 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.257e-02   4.080e-01    5.102e-01  1.000e+00      3027 -3.708e+03 -3.711e+03                   
#> Path [31] :Best Iter: [43] ELBO (-3707.717634) evaluations: (3027) 
#> Path [32] :Initial log joint density = -481753.197809 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.741e-03   2.321e-01    9.357e-01  9.357e-01      3189 -3.698e+03 -3.712e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3698.087840) evaluations: (3189) 
#> Path [33] :Initial log joint density = -484011.144950 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.835e-03   3.011e-01    5.920e-01  5.920e-01      3022 -3.709e+03 -3.714e+03                   
#> Path [33] :Best Iter: [52] ELBO (-3709.457992) evaluations: (3022) 
#> Path [34] :Initial log joint density = -481587.778842 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.147e-03   2.346e-01    1.000e+00  1.000e+00      3554 -3.697e+03 -3.711e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3696.654199) evaluations: (3554) 
#> Path [35] :Initial log joint density = -481482.835385 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.883e-03   2.152e-01    1.000e+00  1.000e+00      2994 -3.707e+03 -3.718e+03                   
#> Path [35] :Best Iter: [45] ELBO (-3706.760889) evaluations: (2994) 
#> Path [36] :Initial log joint density = -481700.226466 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.223e-02   2.241e-01    1.000e+00  1.000e+00      3469 -3.699e+03 -3.700e+03                   
#> Path [36] :Best Iter: [58] ELBO (-3698.861147) evaluations: (3469) 
#> Path [37] :Initial log joint density = -482822.070281 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.977e-03   2.915e-01    4.269e-01  1.000e+00      3400 -3.701e+03 -3.713e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3701.087611) evaluations: (3400) 
#> Path [38] :Initial log joint density = -481704.376004 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.540e-02   1.682e-01    1.000e+00  1.000e+00      3646 -3.701e+03 -3.702e+03                   
#> Path [38] :Best Iter: [57] ELBO (-3700.542339) evaluations: (3646) 
#> Path [39] :Initial log joint density = -483670.112624 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.569e-03   2.583e-01    4.069e-01  1.000e+00      3038 -3.709e+03 -3.719e+03                   
#> Path [39] :Best Iter: [42] ELBO (-3708.842960) evaluations: (3038) 
#> Path [40] :Initial log joint density = -484189.224009 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.209e-02   2.022e-01    1.000e+00  1.000e+00      3267 -3.702e+03 -3.706e+03                   
#> Path [40] :Best Iter: [55] ELBO (-3701.902614) evaluations: (3267) 
#> Path [41] :Initial log joint density = -481473.585143 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.571e-03   2.176e-01    1.000e+00  1.000e+00      3210 -3.707e+03 -3.707e+03                   
#> Path [41] :Best Iter: [56] ELBO (-3706.681160) evaluations: (3210) 
#> Path [42] :Initial log joint density = -481627.394501 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.161e-03   2.191e-01    9.131e-01  9.131e-01      3103 -3.706e+03 -3.707e+03                   
#> Path [42] :Best Iter: [49] ELBO (-3706.213051) evaluations: (3103) 
#> Path [43] :Initial log joint density = -481888.897573 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.040e-02   2.831e-01    1.000e+00  1.000e+00      3630 -3.701e+03 -3.708e+03                   
#> Path [43] :Best Iter: [59] ELBO (-3701.429252) evaluations: (3630) 
#> Path [44] :Initial log joint density = -481773.031118 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.118e-02   4.167e-01    9.590e-01  9.590e-01      3102 -3.709e+03 -3.715e+03                   
#> Path [44] :Best Iter: [43] ELBO (-3709.292487) evaluations: (3102) 
#> Path [45] :Initial log joint density = -482019.909857 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.927e-03   1.610e-01    7.737e-01  7.737e-01      3578 -3.701e+03 -3.714e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3700.666997) evaluations: (3578) 
#> Path [46] :Initial log joint density = -481719.045977 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.554e-02   2.896e-01    9.671e-01  9.671e-01      3200 -3.702e+03 -3.708e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3701.867338) evaluations: (3200) 
#> Path [47] :Initial log joint density = -482722.245715 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.827e-03   2.278e-01    8.892e-01  8.892e-01      3510 -3.704e+03 -3.710e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3704.016008) evaluations: (3510) 
#> Path [48] :Initial log joint density = -482070.276128 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.527e-03   2.646e-01    8.945e-01  8.945e-01      3670 -3.700e+03 -3.711e+03                   
#> Path [48] :Best Iter: [60] ELBO (-3699.923818) evaluations: (3670) 
#> Path [49] :Initial log joint density = -481587.849431 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.041e-02   2.023e-01    1.000e+00  1.000e+00      3563 -3.700e+03 -3.700e+03                   
#> Path [49] :Best Iter: [61] ELBO (-3699.908241) evaluations: (3563) 
#> Path [50] :Initial log joint density = -481573.299649 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.819e-03   3.220e-01    8.951e-01  8.951e-01      2903 -3.710e+03 -3.724e+03                   
#> Path [50] :Best Iter: [42] ELBO (-3709.584588) evaluations: (2903) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }

```
