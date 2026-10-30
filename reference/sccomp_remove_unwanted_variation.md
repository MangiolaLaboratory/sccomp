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
    sccomp_remove_unwanted_variation()
  }
#> Warning: `sccomp_remove_unwanted_variation()` was deprecated in sccomp 1.99.20.
#> ℹ sccomp says: sccomp_remove_unwanted_variation is deprecated. Please use
#>   sccomp_remove_unwanted_effects() instead.
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -483859.353615 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.134e-03   2.403e-01    1.000e+00  1.000e+00      3617 -3.705e+03 -3.706e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3704.553929) evaluations: (3617) 
#> Path [2] :Initial log joint density = -481803.865839 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      3.779e-03   2.173e-01    6.910e-01  6.910e-01      4392 -3.701e+03 -3.718e+03                   
#> Path [2] :Best Iter: [62] ELBO (-3701.167505) evaluations: (4392) 
#> Path [3] :Initial log joint density = -481700.828998 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.289e-02   2.009e-01    1.000e+00  1.000e+00      4634 -3.700e+03 -3.708e+03                   
#> Path [3] :Best Iter: [69] ELBO (-3700.480816) evaluations: (4634) 
#> Path [4] :Initial log joint density = -481770.289114 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.034e-02   2.155e-01    1.000e+00  1.000e+00      3652 -3.705e+03 -3.709e+03                   
#> Path [4] :Best Iter: [58] ELBO (-3704.635949) evaluations: (3652) 
#> Path [5] :Initial log joint density = -481785.917695 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.435e-02   1.860e-01    1.000e+00  1.000e+00      4530 -3.700e+03 -3.700e+03                   
#> Path [5] :Best Iter: [71] ELBO (-3699.968798) evaluations: (4530) 
#> Path [6] :Initial log joint density = -482268.080943 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      2.715e-02   1.713e-01    1.000e+00  1.000e+00      5172 -3.700e+03 -3.702e+03                   
#> Path [6] :Best Iter: [71] ELBO (-3699.535238) evaluations: (5172) 
#> Path [7] :Initial log joint density = -481358.332723 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      8.065e-03   2.588e-01    1.000e+00  1.000e+00      3729 -3.704e+03 -3.716e+03                   
#> Path [7] :Best Iter: [58] ELBO (-3704.137494) evaluations: (3729) 
#> Path [8] :Initial log joint density = -481395.600367 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      5.501e-03   1.941e-01    1.000e+00  1.000e+00      3357 -3.705e+03 -3.709e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3705.325724) evaluations: (3357) 
#> Path [9] :Initial log joint density = -481546.256869 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      8.255e-03   1.685e-01    8.878e-01  8.878e-01      5556 -3.700e+03 -3.712e+03                   
#> Path [9] :Best Iter: [79] ELBO (-3700.485690) evaluations: (5556) 
#> Path [10] :Initial log joint density = -483784.392181 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.086e-02   1.685e-01    1.000e+00  1.000e+00      4122 -3.700e+03 -3.708e+03                   
#> Path [10] :Best Iter: [64] ELBO (-3699.521887) evaluations: (4122) 
#> Path [11] :Initial log joint density = -481554.472407 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      3.909e-03   1.854e-01    7.673e-01  7.673e-01      3357 -3.705e+03 -3.722e+03                   
#> Path [11] :Best Iter: [57] ELBO (-3705.175437) evaluations: (3357) 
#> Path [12] :Initial log joint density = -481367.459344 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.810e-02   2.317e-01    1.000e+00  1.000e+00      4211 -3.703e+03 -3.703e+03                   
#> Path [12] :Best Iter: [58] ELBO (-3702.921514) evaluations: (4211) 
#> Path [13] :Initial log joint density = -481825.008622 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      8.073e-03   1.926e-01    1.000e+00  1.000e+00      5650 -3.697e+03 -3.709e+03                   
#> Path [13] :Best Iter: [80] ELBO (-3697.389043) evaluations: (5650) 
#> Path [14] :Initial log joint density = -482229.366677 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.766e-02   2.524e-01    1.000e+00  1.000e+00      4942 -3.697e+03 -3.701e+03                   
#> Path [14] :Best Iter: [73] ELBO (-3697.158354) evaluations: (4942) 
#> Path [15] :Initial log joint density = -481479.330336 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.705e-03   1.641e-01    8.743e-01  8.743e-01      3821 -3.702e+03 -3.712e+03                   
#> Path [15] :Best Iter: [60] ELBO (-3702.182552) evaluations: (3821) 
#> Path [16] :Initial log joint density = -481401.636871 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      4.938e-03   1.560e-01    8.159e-01  8.159e-01      3644 -3.702e+03 -3.712e+03                   
#> Path [16] :Best Iter: [61] ELBO (-3702.319322) evaluations: (3644) 
#> Path [17] :Initial log joint density = -481561.392241 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      2.316e-02   1.480e-01    1.000e+00  1.000e+00      5144 -3.699e+03 -3.701e+03                   
#> Path [17] :Best Iter: [73] ELBO (-3699.265200) evaluations: (5144) 
#> Path [18] :Initial log joint density = -481431.002315 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      6.475e-03   1.906e-01    8.082e-01  8.082e-01      3739 -3.703e+03 -3.715e+03                   
#> Path [18] :Best Iter: [60] ELBO (-3702.703693) evaluations: (3739) 
#> Path [19] :Initial log joint density = -483166.832878 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.513e-02   2.230e-01    9.221e-01  9.221e-01      4277 -3.702e+03 -3.710e+03                   
#> Path [19] :Best Iter: [66] ELBO (-3702.042927) evaluations: (4277) 
#> Path [20] :Initial log joint density = -482308.013449 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.690e-02   1.838e-01    9.186e-01  9.186e-01      4624 -3.698e+03 -3.708e+03                   
#> Path [20] :Best Iter: [69] ELBO (-3697.813642) evaluations: (4624) 
#> Path [21] :Initial log joint density = -481872.032222 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.051e-02   1.874e-01    1.000e+00  1.000e+00      4612 -3.704e+03 -3.701e+03                   
#> Path [21] :Best Iter: [71] ELBO (-3700.634052) evaluations: (4612) 
#> Path [22] :Initial log joint density = -481285.104304 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      4.705e-03   2.045e-01    7.859e-01  7.859e-01      4957 -3.700e+03 -3.705e+03                   
#> Path [22] :Best Iter: [73] ELBO (-3700.007615) evaluations: (4957) 
#> Path [23] :Initial log joint density = -481566.875315 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.495e-02   3.195e-01    1.000e+00  1.000e+00      3590 -3.701e+03 -3.710e+03                   
#> Path [23] :Best Iter: [59] ELBO (-3700.963222) evaluations: (3590) 
#> Path [24] :Initial log joint density = -484343.301320 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.358e-02   2.312e-01    1.000e+00  1.000e+00      4602 -3.699e+03 -3.704e+03                   
#> Path [24] :Best Iter: [69] ELBO (-3698.755869) evaluations: (4602) 
#> Path [25] :Initial log joint density = -481461.792992 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      9.223e-03   2.394e-01    9.269e-01  9.269e-01      3880 -3.699e+03 -3.713e+03                   
#> Path [25] :Best Iter: [61] ELBO (-3699.336516) evaluations: (3880) 
#> Path [26] :Initial log joint density = -481825.086213 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      7.097e-03   1.809e-01    1.000e+00  1.000e+00      3956 -3.704e+03 -3.710e+03                   
#> Path [26] :Best Iter: [62] ELBO (-3703.735143) evaluations: (3956) 
#> Path [27] :Initial log joint density = -481384.647944 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.547e-02   2.675e-01    1.000e+00  1.000e+00      4279 -3.699e+03 -3.705e+03                   
#> Path [27] :Best Iter: [68] ELBO (-3698.616801) evaluations: (4279) 
#> Path [28] :Initial log joint density = -481435.055990 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.424e-02   1.715e-01    1.000e+00  1.000e+00      5467 -3.697e+03 -3.705e+03                   
#> Path [28] :Best Iter: [78] ELBO (-3697.292047) evaluations: (5467) 
#> Path [29] :Initial log joint density = -481663.903316 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      8.914e-03   1.773e-01    1.000e+00  1.000e+00      4697 -3.698e+03 -3.701e+03                   
#> Path [29] :Best Iter: [69] ELBO (-3698.458294) evaluations: (4697) 
#> Path [30] :Initial log joint density = -484557.316914 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      4.701e-03   2.456e-01    8.907e-01  8.907e-01      3778 -3.706e+03 -3.718e+03                   
#> Path [30] :Best Iter: [59] ELBO (-3705.724501) evaluations: (3778) 
#> Path [31] :Initial log joint density = -482191.071798 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.437e-02   2.765e-01    4.769e-01  1.000e+00      4638 -3.702e+03 -3.709e+03                   
#> Path [31] :Best Iter: [66] ELBO (-3702.267108) evaluations: (4638) 
#> Path [32] :Initial log joint density = -481519.303578 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.370e-02   1.990e-01    1.000e+00  1.000e+00      5119 -3.701e+03 -3.705e+03                   
#> Path [32] :Best Iter: [74] ELBO (-3700.892029) evaluations: (5119) 
#> Path [33] :Initial log joint density = -481330.580057 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.526e-02   2.229e-01    1.000e+00  1.000e+00      3384 -3.704e+03 -3.713e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3704.013671) evaluations: (3384) 
#> Path [34] :Initial log joint density = -481577.300438 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.858e-02   2.964e-01    1.000e+00  1.000e+00      3646 -3.701e+03 -3.707e+03                   
#> Path [34] :Best Iter: [59] ELBO (-3701.163762) evaluations: (3646) 
#> Path [35] :Initial log joint density = -482161.258004 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.407e-02   1.956e-01    1.000e+00  1.000e+00      5541 -3.701e+03 -3.701e+03                   
#> Path [35] :Best Iter: [80] ELBO (-3700.914931) evaluations: (5541) 
#> Path [36] :Initial log joint density = -481280.595763 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      7.833e-03   2.366e-01    1.000e+00  1.000e+00      3736 -3.706e+03 -3.713e+03                   
#> Path [36] :Best Iter: [58] ELBO (-3706.348122) evaluations: (3736) 
#> Path [37] :Initial log joint density = -485054.493552 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.063e-02   3.419e-01    1.000e+00  1.000e+00      4058 -3.701e+03 -3.707e+03                   
#> Path [37] :Best Iter: [64] ELBO (-3700.762731) evaluations: (4058) 
#> Path [38] :Initial log joint density = -481932.417454 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.375e-02   1.463e-01    1.000e+00  1.000e+00      5524 -3.701e+03 -3.702e+03                   
#> Path [38] :Best Iter: [73] ELBO (-3701.112213) evaluations: (5524) 
#> Path [39] :Initial log joint density = -481610.511162 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.545e-03   2.490e-01    1.000e+00  1.000e+00      3402 -3.704e+03 -3.715e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3704.477809) evaluations: (3402) 
#> Path [40] :Initial log joint density = -481621.605657 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      7.250e-03   1.921e-01    9.185e-01  9.185e-01      4924 -3.697e+03 -3.706e+03                   
#> Path [40] :Best Iter: [73] ELBO (-3697.082773) evaluations: (4924) 
#> Path [41] :Initial log joint density = -482574.412602 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.116e-02   1.919e-01    8.640e-01  8.640e-01      4592 -3.701e+03 -3.713e+03                   
#> Path [41] :Best Iter: [67] ELBO (-3700.685234) evaluations: (4592) 
#> Path [42] :Initial log joint density = -481624.925636 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.301e-02   2.467e-01    9.706e-01  9.706e-01      4501 -3.700e+03 -3.707e+03                   
#> Path [42] :Best Iter: [68] ELBO (-3699.636056) evaluations: (4501) 
#> Path [43] :Initial log joint density = -485666.962652 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      7.956e-03   2.372e-01    1.000e+00  1.000e+00      4845 -3.703e+03 -3.709e+03                   
#> Path [43] :Best Iter: [70] ELBO (-3703.395732) evaluations: (4845) 
#> Path [44] :Initial log joint density = -481908.261271 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.654e-02   2.228e-01    1.000e+00  1.000e+00      4581 -3.700e+03 -3.701e+03                   
#> Path [44] :Best Iter: [68] ELBO (-3699.500996) evaluations: (4581) 
#> Path [45] :Initial log joint density = -481169.751074 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.149e-02   2.317e-01    9.403e-01  9.403e-01      3559 -3.702e+03 -3.711e+03                   
#> Path [45] :Best Iter: [59] ELBO (-3702.072987) evaluations: (3559) 
#> Path [46] :Initial log joint density = -481451.925177 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.696e-03   2.844e-01    5.882e-01  5.882e-01      4056 -3.700e+03 -3.715e+03                   
#> Path [46] :Best Iter: [64] ELBO (-3699.875744) evaluations: (4056) 
#> Path [47] :Initial log joint density = -481492.659552 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      9.070e-03   2.296e-01    8.029e-01  8.029e-01      4326 -3.700e+03 -3.711e+03                   
#> Path [47] :Best Iter: [66] ELBO (-3699.733115) evaluations: (4326) 
#> Path [48] :Initial log joint density = -481748.425265 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.320e-02   2.319e-01    1.000e+00  1.000e+00      3744 -3.705e+03 -3.706e+03                   
#> Path [48] :Best Iter: [59] ELBO (-3704.934832) evaluations: (3744) 
#> Path [49] :Initial log joint density = -481708.492702 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.297e-02   2.729e-01    1.000e+00  1.000e+00      3545 -3.702e+03 -3.710e+03                   
#> Path [49] :Best Iter: [60] ELBO (-3702.392220) evaluations: (3545) 
#> Path [50] :Initial log joint density = -482961.146901 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      8.107e-03   2.502e-01    1.000e+00  1.000e+00      3734 -3.704e+03 -3.712e+03                   
#> Path [50] :Best Iter: [58] ELBO (-3704.039987) evaluations: (3734) 
#> Finished in  14.7 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.517 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.52 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
