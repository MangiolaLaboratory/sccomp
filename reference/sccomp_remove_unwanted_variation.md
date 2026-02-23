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
#> Path [1] :Initial log joint density = -481906.161217 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.562e-03   2.484e-01    1.000e+00  1.000e+00      3294 -3.707e+03 -3.703e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3702.602492) evaluations: (3294) 
#> Path [2] :Initial log joint density = -484098.341258 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.778e-03   2.713e-01    8.686e-01  8.686e-01      2940 -3.708e+03 -3.718e+03                   
#> Path [2] :Best Iter: [44] ELBO (-3708.092441) evaluations: (2940) 
#> Path [3] :Initial log joint density = -482276.976929 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.453e-03   1.509e-01    1.000e+00  1.000e+00      3627 -3.701e+03 -3.710e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3701.206528) evaluations: (3627) 
#> Path [4] :Initial log joint density = -481454.682013 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.851e-03   2.083e-01    1.000e+00  1.000e+00      2776 -3.708e+03 -3.710e+03                   
#> Path [4] :Best Iter: [50] ELBO (-3708.296816) evaluations: (2776) 
#> Path [5] :Initial log joint density = -481243.288960 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.534e-03   1.992e-01    1.000e+00  1.000e+00      2798 -3.707e+03 -3.715e+03                   
#> Path [5] :Best Iter: [47] ELBO (-3707.417779) evaluations: (2798) 
#> Path [6] :Initial log joint density = -486651.479622 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.255e-03   2.415e-01    1.000e+00  1.000e+00      3364 -3.702e+03 -3.702e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3702.250691) evaluations: (3364) 
#> Path [7] :Initial log joint density = -481626.469133 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.114e-03   1.857e-01    9.258e-01  9.258e-01      3533 -3.698e+03 -3.709e+03                   
#> Path [7] :Best Iter: [59] ELBO (-3697.642291) evaluations: (3533) 
#> Path [8] :Initial log joint density = -481873.033489 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.098e-03   2.112e-01    1.000e+00  1.000e+00      3068 -3.709e+03 -3.703e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3703.021570) evaluations: (3068) 
#> Path [9] :Initial log joint density = -481714.868405 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.257e-03   2.101e-01    1.000e+00  1.000e+00      3251 -3.707e+03 -3.713e+03                   
#> Path [9] :Best Iter: [53] ELBO (-3706.920381) evaluations: (3251) 
#> Path [10] :Initial log joint density = -481834.244627 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.244e-03   2.458e-01    1.000e+00  1.000e+00      3053 -3.709e+03 -3.708e+03                   
#> Path [10] :Best Iter: [56] ELBO (-3707.610559) evaluations: (3053) 
#> Path [11] :Initial log joint density = -481881.983585 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.431e-03   2.743e-01    1.000e+00  1.000e+00      3215 -3.706e+03 -3.712e+03                   
#> Path [11] :Best Iter: [53] ELBO (-3706.324859) evaluations: (3215) 
#> Path [12] :Initial log joint density = -481528.501027 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.006e-02   2.023e-01    1.000e+00  1.000e+00      3238 -3.702e+03 -3.703e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3702.058273) evaluations: (3238) 
#> Path [13] :Initial log joint density = -481650.421273 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.183e-03   1.570e-01    7.086e-01  7.086e-01      3584 -3.700e+03 -3.715e+03                   
#> Path [13] :Best Iter: [59] ELBO (-3699.517733) evaluations: (3584) 
#> Path [14] :Initial log joint density = -481600.832143 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.775e-03   2.206e-01    1.000e+00  1.000e+00      3440 -3.701e+03 -3.703e+03                   
#> Path [14] :Best Iter: [58] ELBO (-3700.700434) evaluations: (3440) 
#> Path [15] :Initial log joint density = -482074.882928 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      8.848e-03   2.698e-01    4.538e-01  1.000e+00      3901 -3.699e+03 -3.709e+03                   
#> Path [15] :Best Iter: [64] ELBO (-3699.202306) evaluations: (3901) 
#> Path [16] :Initial log joint density = -481627.182239 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      6.911e-03   2.022e-01    1.000e+00  1.000e+00      3769 -3.699e+03 -3.705e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3698.764343) evaluations: (3769) 
#> Path [17] :Initial log joint density = -485594.612149 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.096e-02   2.724e-01    9.589e-01  9.589e-01      3656 -3.701e+03 -3.710e+03                   
#> Path [17] :Best Iter: [60] ELBO (-3700.683571) evaluations: (3656) 
#> Path [18] :Initial log joint density = -481612.348749 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.437e-02   2.514e-01    1.000e+00  1.000e+00      3825 -3.696e+03 -3.705e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3695.666436) evaluations: (3825) 
#> Path [19] :Initial log joint density = -481575.489154 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.790e-03   1.924e-01    1.000e+00  1.000e+00      3013 -3.706e+03 -3.711e+03                   
#> Path [19] :Best Iter: [52] ELBO (-3705.751978) evaluations: (3013) 
#> Path [20] :Initial log joint density = -481735.824434 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.120e-02   3.114e-01    1.000e+00  1.000e+00      2844 -3.711e+03 -3.712e+03                   
#> Path [20] :Best Iter: [39] ELBO (-3710.750126) evaluations: (2844) 
#> Path [21] :Initial log joint density = -483270.734682 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.694e-03   2.365e-01    1.000e+00  1.000e+00      3374 -3.703e+03 -3.705e+03                   
#> Path [21] :Best Iter: [56] ELBO (-3703.451520) evaluations: (3374) 
#> Path [22] :Initial log joint density = -485169.552171 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.146e-03   3.034e-01    5.913e-01  5.913e-01      3477 -3.703e+03 -3.718e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3703.089370) evaluations: (3477) 
#> Path [23] :Initial log joint density = -481723.851120 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.504e-02   2.457e-01    8.213e-01  8.213e-01      3422 -3.701e+03 -3.704e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3701.141613) evaluations: (3422) 
#> Path [24] :Initial log joint density = -482368.596519 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.199e-03   1.851e-01    1.000e+00  1.000e+00      3344 -3.706e+03 -3.704e+03                   
#> Path [24] :Best Iter: [57] ELBO (-3703.788527) evaluations: (3344) 
#> Path [25] :Initial log joint density = -483501.409755 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.932e-03   2.538e-01    9.943e-01  9.943e-01      3373 -3.697e+03 -3.706e+03                   
#> Path [25] :Best Iter: [56] ELBO (-3696.557105) evaluations: (3373) 
#> Path [26] :Initial log joint density = -481654.926918 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.171e-02   2.235e-01    1.000e+00  1.000e+00      3391 -3.701e+03 -3.704e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3700.885647) evaluations: (3391) 
#> Path [27] :Initial log joint density = -483015.505479 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.689e-03   2.727e-01    1.000e+00  1.000e+00      2910 -3.709e+03 -3.721e+03                   
#> Path [27] :Best Iter: [47] ELBO (-3709.226986) evaluations: (2910) 
#> Path [28] :Initial log joint density = -481717.659059 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.087e-02   3.249e-01    8.796e-01  8.796e-01      2869 -3.707e+03 -3.721e+03                   
#> Path [28] :Best Iter: [50] ELBO (-3706.640062) evaluations: (2869) 
#> Path [29] :Initial log joint density = -482474.937570 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.701e-03   2.364e-01    6.160e-01  6.160e-01      3478 -3.703e+03 -3.710e+03                   
#> Path [29] :Best Iter: [56] ELBO (-3702.601270) evaluations: (3478) 
#> Path [30] :Initial log joint density = -484801.993153 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.259e-03   2.707e-01    7.233e-01  7.233e-01      3187 -3.700e+03 -3.714e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3700.238412) evaluations: (3187) 
#> Path [31] :Initial log joint density = -481932.006402 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.600e-02   2.528e-01    1.000e+00  1.000e+00      3122 -3.705e+03 -3.712e+03                   
#> Path [31] :Best Iter: [47] ELBO (-3704.507582) evaluations: (3122) 
#> Path [32] :Initial log joint density = -483725.520998 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.539e-03   3.448e-01    3.649e-01  1.000e+00      3333 -3.700e+03 -3.710e+03                   
#> Path [32] :Best Iter: [57] ELBO (-3700.047236) evaluations: (3333) 
#> Path [33] :Initial log joint density = -481949.457732 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.210e-03   3.383e-01    5.668e-01  5.668e-01      3450 -3.702e+03 -3.710e+03                   
#> Path [33] :Best Iter: [58] ELBO (-3701.594214) evaluations: (3450) 
#> Path [34] :Initial log joint density = -482072.537571 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      5.517e-03   1.381e-01    1.000e+00  1.000e+00      3970 -3.702e+03 -3.709e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3701.931801) evaluations: (3970) 
#> Path [35] :Initial log joint density = -481945.853146 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.160e-03   2.308e-01    9.493e-01  9.493e-01      3519 -3.701e+03 -3.707e+03                   
#> Path [35] :Best Iter: [59] ELBO (-3700.728299) evaluations: (3519) 
#> Path [36] :Initial log joint density = -481703.255672 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.485e-03   1.680e-01    1.000e+00  1.000e+00      3102 -3.707e+03 -3.717e+03                   
#> Path [36] :Best Iter: [45] ELBO (-3706.579302) evaluations: (3102) 
#> Path [37] :Initial log joint density = -482806.970879 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.435e-03   3.826e-01    5.103e-01  1.000e+00      3211 -3.707e+03 -3.708e+03                   
#> Path [37] :Best Iter: [45] ELBO (-3707.151635) evaluations: (3211) 
#> Path [38] :Initial log joint density = -482020.748849 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.403e-03   2.646e-01    1.000e+00  1.000e+00      3322 -3.701e+03 -3.708e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3701.462230) evaluations: (3322) 
#> Path [39] :Initial log joint density = -481584.191255 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.195e-02   3.600e-01    1.000e+00  1.000e+00      3109 -3.700e+03 -3.707e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3699.935281) evaluations: (3109) 
#> Path [40] :Initial log joint density = -481986.122779 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.667e-03   2.958e-01    1.000e+00  1.000e+00      3371 -3.700e+03 -3.708e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3699.811210) evaluations: (3371) 
#> Path [41] :Initial log joint density = -481667.392480 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.503e-03   2.471e-01    1.000e+00  1.000e+00      2890 -3.709e+03 -3.717e+03                   
#> Path [41] :Best Iter: [41] ELBO (-3709.207836) evaluations: (2890) 
#> Path [42] :Initial log joint density = -481712.727730 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.028e-03   2.922e-01    1.000e+00  1.000e+00      3352 -3.697e+03 -3.709e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3697.246763) evaluations: (3352) 
#> Path [43] :Initial log joint density = -481476.995909 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.238e-02   2.913e-01    1.000e+00  1.000e+00      3148 -3.705e+03 -3.704e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3704.277374) evaluations: (3148) 
#> Path [44] :Initial log joint density = -481846.253158 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.136e-03   2.213e-01    7.975e-01  7.975e-01      3439 -3.700e+03 -3.709e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3699.911809) evaluations: (3439) 
#> Path [45] :Initial log joint density = -484331.085093 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.495e-03   3.019e-01    6.102e-01  6.102e-01      3453 -3.701e+03 -3.716e+03                   
#> Path [45] :Best Iter: [56] ELBO (-3701.002240) evaluations: (3453) 
#> Path [46] :Initial log joint density = -481475.604204 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.280e-02   2.882e-01    1.000e+00  1.000e+00      2748 -3.708e+03 -3.718e+03                   
#> Path [46] :Best Iter: [40] ELBO (-3707.825072) evaluations: (2748) 
#> Path [47] :Initial log joint density = -481290.039507 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.688e-03   1.885e-01    1.000e+00  1.000e+00      2909 -3.710e+03 -3.721e+03                   
#> Path [47] :Best Iter: [39] ELBO (-3709.964485) evaluations: (2909) 
#> Path [48] :Initial log joint density = -481453.840799 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.847e-03   2.046e-01    1.000e+00  1.000e+00      3570 -3.701e+03 -3.704e+03                   
#> Path [48] :Best Iter: [56] ELBO (-3701.425825) evaluations: (3570) 
#> Path [49] :Initial log joint density = -481688.501072 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.642e-03   2.070e-01    1.000e+00  1.000e+00      3212 -3.704e+03 -3.699e+03                   
#> Path [49] :Best Iter: [56] ELBO (-3698.753639) evaluations: (3212) 
#> Path [50] :Initial log joint density = -481880.497602 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.253e-02   4.072e-01    1.000e+00  1.000e+00      3248 -3.701e+03 -3.710e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3700.642528) evaluations: (3248) 
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
