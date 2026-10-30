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

- **logit_residuals** - A numeric column representing the centered
  log-ratio residuals calculated after adjustment.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482518.273125 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.252e-02   2.809e-01    1.000e+00  1.000e+00      3886 -3.702e+03 -3.708e+03                   
#> Path [1] :Best Iter: [61] ELBO (-3702.426406) evaluations: (3886) 
#> Path [2] :Initial log joint density = -481727.935616 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.235e-02   1.903e-01    1.000e+00  1.000e+00      5185 -3.703e+03 -3.709e+03                   
#> Path [2] :Best Iter: [68] ELBO (-3702.632885) evaluations: (5185) 
#> Path [3] :Initial log joint density = -482701.433017 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.872e-03   2.029e-01    8.876e-01  8.876e-01      3867 -3.699e+03 -3.710e+03                   
#> Path [3] :Best Iter: [61] ELBO (-3699.245174) evaluations: (3867) 
#> Path [4] :Initial log joint density = -481821.242091 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.190e-03   2.464e-01    6.347e-01  6.347e-01      4621 -3.702e+03 -3.708e+03                   
#> Path [4] :Best Iter: [66] ELBO (-3702.300563) evaluations: (4621) 
#> Path [5] :Initial log joint density = -481756.820973 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      4.049e-03   2.209e-01    6.971e-01  6.971e-01      4044 -3.698e+03 -3.713e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3698.267365) evaluations: (4044) 
#> Path [6] :Initial log joint density = -483187.820426 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      4.880e-03   2.035e-01    7.652e-01  7.652e-01      4348 -3.699e+03 -3.711e+03                   
#> Path [6] :Best Iter: [65] ELBO (-3699.147640) evaluations: (4348) 
#> Path [7] :Initial log joint density = -481446.791238 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.039e-02   2.128e-01    1.000e+00  1.000e+00      5398 -3.701e+03 -3.711e+03                   
#> Path [7] :Best Iter: [78] ELBO (-3700.995662) evaluations: (5398) 
#> Path [8] :Initial log joint density = -482055.083633 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      2.118e-02   1.885e-01    1.000e+00  1.000e+00      5187 -3.702e+03 -3.701e+03                   
#> Path [8] :Best Iter: [77] ELBO (-3700.864665) evaluations: (5187) 
#> Path [9] :Initial log joint density = -481477.608383 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.846e-02   3.425e-01    9.121e-01  9.121e-01      3686 -3.702e+03 -3.712e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3701.625046) evaluations: (3686) 
#> Path [10] :Initial log joint density = -481768.768513 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.410e-02   2.812e-01    9.579e-01  9.579e-01      4221 -3.699e+03 -3.716e+03                   
#> Path [10] :Best Iter: [66] ELBO (-3699.286256) evaluations: (4221) 
#> Path [11] :Initial log joint density = -481569.094866 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      5.652e-03   2.871e-01    5.734e-01  5.734e-01      3882 -3.704e+03 -3.717e+03                   
#> Path [11] :Best Iter: [57] ELBO (-3703.964489) evaluations: (3882) 
#> Path [12] :Initial log joint density = -482163.983829 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.344e-02   2.554e-01    1.000e+00  1.000e+00      4774 -3.701e+03 -3.704e+03                   
#> Path [12] :Best Iter: [70] ELBO (-3700.728505) evaluations: (4774) 
#> Path [13] :Initial log joint density = -481679.980549 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.809e-02   1.857e-01    1.000e+00  1.000e+00      5018 -3.702e+03 -3.700e+03                   
#> Path [13] :Best Iter: [75] ELBO (-3699.850084) evaluations: (5018) 
#> Path [14] :Initial log joint density = -481464.078047 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.511e-02   1.577e-01    1.000e+00  1.000e+00      4151 -3.700e+03 -3.705e+03                   
#> Path [14] :Best Iter: [60] ELBO (-3699.892144) evaluations: (4151) 
#> Path [15] :Initial log joint density = -483895.905240 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              88      -4.787e+05      5.477e-03   2.669e-01    7.206e-01  7.206e-01      6502 -3.699e+03 -3.710e+03                   
#> Path [15] :Best Iter: [83] ELBO (-3698.716561) evaluations: (6502) 
#> Path [16] :Initial log joint density = -481768.837640 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.605e-02   2.762e-01    9.711e-01  9.711e-01      5289 -3.701e+03 -3.711e+03                   
#> Path [16] :Best Iter: [76] ELBO (-3701.220015) evaluations: (5289) 
#> Path [17] :Initial log joint density = -481850.709518 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      1.149e-02   2.293e-01    8.216e-01  8.216e-01      6082 -3.695e+03 -3.707e+03                   
#> Path [17] :Best Iter: [75] ELBO (-3695.348313) evaluations: (6082) 
#> Path [18] :Initial log joint density = -481871.562123 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      5.825e-03   1.824e-01    1.000e+00  1.000e+00      4426 -3.703e+03 -3.711e+03                   
#> Path [18] :Best Iter: [67] ELBO (-3702.850095) evaluations: (4426) 
#> Path [19] :Initial log joint density = -481550.291208 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.680e-03   2.136e-01    1.000e+00  1.000e+00      3246 -3.705e+03 -3.709e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3704.834979) evaluations: (3246) 
#> Path [20] :Initial log joint density = -481491.853270 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.980e-02   2.411e-01    1.000e+00  1.000e+00      4027 -3.701e+03 -3.702e+03                   
#> Path [20] :Best Iter: [62] ELBO (-3700.885960) evaluations: (4027) 
#> Path [21] :Initial log joint density = -481670.831122 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      9.161e-03   2.519e-01    5.070e-01  5.070e-01      5111 -3.702e+03 -3.703e+03                   
#> Path [21] :Best Iter: [74] ELBO (-3701.727654) evaluations: (5111) 
#> Path [22] :Initial log joint density = -481798.250172 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.727e-02   3.262e-01    1.000e+00  1.000e+00      4767 -3.702e+03 -3.706e+03                   
#> Path [22] :Best Iter: [71] ELBO (-3702.389935) evaluations: (4767) 
#> Path [23] :Initial log joint density = -482332.292778 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      6.223e-03   2.547e-01    1.000e+00  1.000e+00      3503 -3.704e+03 -3.717e+03                   
#> Path [23] :Best Iter: [57] ELBO (-3704.066275) evaluations: (3503) 
#> Path [24] :Initial log joint density = -481518.698036 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.969e-03   1.819e-01    1.000e+00  1.000e+00      4115 -3.700e+03 -3.702e+03                   
#> Path [24] :Best Iter: [64] ELBO (-3700.088568) evaluations: (4115) 
#> Path [25] :Initial log joint density = -487304.361559 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.618e-02   2.053e-01    1.000e+00  1.000e+00      4565 -3.700e+03 -3.701e+03                   
#> Path [25] :Best Iter: [68] ELBO (-3700.022328) evaluations: (4565) 
#> Path [26] :Initial log joint density = -481706.443325 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.696e-02   2.123e-01    1.000e+00  1.000e+00      4005 -3.704e+03 -3.704e+03                   
#> Path [26] :Best Iter: [64] ELBO (-3703.709561) evaluations: (4005) 
#> Path [27] :Initial log joint density = -481587.974838 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.772e-02   3.092e-01    1.000e+00  1.000e+00      3784 -3.700e+03 -3.712e+03                   
#> Path [27] :Best Iter: [62] ELBO (-3700.495496) evaluations: (3784) 
#> Path [28] :Initial log joint density = -481577.813020 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.458e-02   2.190e-01    1.000e+00  1.000e+00      3616 -3.703e+03 -3.704e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3703.447960) evaluations: (3616) 
#> Path [29] :Initial log joint density = -481348.454582 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.052e-02   2.147e-01    1.000e+00  1.000e+00      3743 -3.702e+03 -3.707e+03                   
#> Path [29] :Best Iter: [58] ELBO (-3701.956881) evaluations: (3743) 
#> Path [30] :Initial log joint density = -483975.724387 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      8.036e-03   1.785e-01    1.000e+00  1.000e+00      4108 -3.700e+03 -3.708e+03                   
#> Path [30] :Best Iter: [63] ELBO (-3699.814664) evaluations: (4108) 
#> Path [31] :Initial log joint density = -482021.883127 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      6.587e-03   2.587e-01    1.000e+00  1.000e+00      4210 -3.705e+03 -3.710e+03                   
#> Path [31] :Best Iter: [65] ELBO (-3704.815476) evaluations: (4210) 
#> Path [32] :Initial log joint density = -481543.174381 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      9.664e-03   2.096e-01    1.000e+00  1.000e+00      4244 -3.700e+03 -3.706e+03                   
#> Path [32] :Best Iter: [62] ELBO (-3699.794423) evaluations: (4244) 
#> Path [33] :Initial log joint density = -481945.634980 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.138e-02   2.188e-01    1.000e+00  1.000e+00      3642 -3.703e+03 -3.706e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3703.166944) evaluations: (3642) 
#> Path [34] :Initial log joint density = -481975.606532 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.050e-02   2.032e-01    1.000e+00  1.000e+00      5118 -3.701e+03 -3.700e+03                   
#> Path [34] :Best Iter: [75] ELBO (-3699.633615) evaluations: (5118) 
#> Path [35] :Initial log joint density = -481551.552675 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.604e-02   2.122e-01    1.000e+00  1.000e+00      3788 -3.702e+03 -3.706e+03                   
#> Path [35] :Best Iter: [62] ELBO (-3702.023500) evaluations: (3788) 
#> Path [36] :Initial log joint density = -481642.765114 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.764e-02   2.264e-01    1.000e+00  1.000e+00      4753 -3.702e+03 -3.709e+03                   
#> Path [36] :Best Iter: [68] ELBO (-3701.851730) evaluations: (4753) 
#> Path [37] :Initial log joint density = -481830.458664 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.796e-02   3.170e-01    1.000e+00  1.000e+00      5507 -3.699e+03 -3.710e+03                   
#> Path [37] :Best Iter: [78] ELBO (-3699.329393) evaluations: (5507) 
#> Path [38] :Initial log joint density = -481684.339804 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.492e-03   1.155e-01    1.000e+00  1.000e+00      3836 -3.702e+03 -3.713e+03                   
#> Path [38] :Best Iter: [59] ELBO (-3701.958824) evaluations: (3836) 
#> Path [39] :Initial log joint density = -482314.171347 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.003e-02   2.328e-01    8.455e-01  8.455e-01      3917 -3.698e+03 -3.711e+03                   
#> Path [39] :Best Iter: [62] ELBO (-3697.678602) evaluations: (3917) 
#> Path [40] :Initial log joint density = -482067.573487 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      8.757e-03   1.969e-01    5.234e-01  5.234e-01      5406 -3.700e+03 -3.705e+03                   
#> Path [40] :Best Iter: [77] ELBO (-3700.046596) evaluations: (5406) 
#> Path [41] :Initial log joint density = -481416.107437 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.425e-02   2.126e-01    1.000e+00  1.000e+00      4214 -3.703e+03 -3.702e+03                   
#> Path [41] :Best Iter: [67] ELBO (-3702.482762) evaluations: (4214) 
#> Path [42] :Initial log joint density = -481820.549925 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      8.705e-03   2.064e-01    1.000e+00  1.000e+00      4130 -3.700e+03 -3.705e+03                   
#> Path [42] :Best Iter: [62] ELBO (-3700.423862) evaluations: (4130) 
#> Path [43] :Initial log joint density = -481454.507716 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      2.372e-02   2.724e-01    9.942e-01  9.942e-01      4036 -3.702e+03 -3.711e+03                   
#> Path [43] :Best Iter: [59] ELBO (-3702.465481) evaluations: (4036) 
#> Path [44] :Initial log joint density = -481459.296150 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.782e-03   2.977e-01    6.621e-01  6.621e-01      3867 -3.702e+03 -3.714e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3701.795258) evaluations: (3867) 
#> Path [45] :Initial log joint density = -481280.379033 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.435e-02   2.255e-01    1.000e+00  1.000e+00      4687 -3.702e+03 -3.700e+03                   
#> Path [45] :Best Iter: [72] ELBO (-3699.517415) evaluations: (4687) 
#> Path [46] :Initial log joint density = -484340.710538 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.031e-02   1.780e-01    1.000e+00  1.000e+00      4966 -3.699e+03 -3.700e+03                   
#> Path [46] :Best Iter: [71] ELBO (-3698.802222) evaluations: (4966) 
#> Path [47] :Initial log joint density = -481699.587056 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.567e-02   3.922e-01    1.000e+00  1.000e+00      4300 -3.702e+03 -3.711e+03                   
#> Path [47] :Best Iter: [67] ELBO (-3701.506841) evaluations: (4300) 
#> Path [48] :Initial log joint density = -482598.418757 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.403e-03   2.278e-01    1.000e+00  1.000e+00      3972 -3.702e+03 -3.705e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3701.677137) evaluations: (3972) 
#> Path [49] :Initial log joint density = -481859.173657 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.464e-02   1.290e-01    1.000e+00  1.000e+00      5213 -3.699e+03 -3.703e+03                   
#> Path [49] :Best Iter: [75] ELBO (-3699.242267) evaluations: (5213) 
#> Path [50] :Initial log joint density = -481443.539221 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.030e-02   3.508e-01    1.000e+00  1.000e+00      3075 -3.711e+03 -3.719e+03                   
#> Path [50] :Best Iter: [54] ELBO (-3710.665500) evaluations: (3075) 
#> Finished in  15.0 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.519 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.52 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
