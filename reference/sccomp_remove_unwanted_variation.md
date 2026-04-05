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
#> Path [1] :Initial log joint density = -481765.057867 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.442e-02   2.595e-01    1.000e+00  1.000e+00      2994 -3.691e+03 -3.692e+03                   
#> Path [1] :Best Iter: [51] ELBO (-3691.038719) evaluations: (2994) 
#> Path [2] :Initial log joint density = -481685.463002 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.447e-03   2.152e-01    1.000e+00  1.000e+00      3095 -3.690e+03 -3.694e+03                   
#> Path [2] :Best Iter: [48] ELBO (-3689.688748) evaluations: (3095) 
#> Path [3] :Initial log joint density = -481593.964316 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.721e-03   2.323e-01    1.000e+00  1.000e+00      3274 -3.684e+03 -3.691e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3683.660769) evaluations: (3274) 
#> Path [4] :Initial log joint density = -482536.464759 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.100e-02   3.063e-01    1.000e+00  1.000e+00      3485 -3.682e+03 -3.686e+03                   
#> Path [4] :Best Iter: [58] ELBO (-3681.889278) evaluations: (3485) 
#> Path [5] :Initial log joint density = -482500.074690 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.953e-03   2.761e-01    6.917e-01  6.917e-01      3415 -3.685e+03 -3.701e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3685.081259) evaluations: (3415) 
#> Path [6] :Initial log joint density = -481615.893042 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.589e-03   2.163e-01    9.549e-01  9.549e-01      3493 -3.684e+03 -3.692e+03                   
#> Path [6] :Best Iter: [57] ELBO (-3684.149387) evaluations: (3493) 
#> Path [7] :Initial log joint density = -481572.012556 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.979e-02   3.225e-01    1.000e+00  1.000e+00      3435 -3.683e+03 -3.689e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3682.708317) evaluations: (3435) 
#> Path [8] :Initial log joint density = -485791.724450 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.998e-03   1.445e-01    1.000e+00  1.000e+00      3418 -3.691e+03 -3.695e+03                   
#> Path [8] :Best Iter: [53] ELBO (-3691.186000) evaluations: (3418) 
#> Path [9] :Initial log joint density = -481514.937431 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.084e-03   2.315e-01    9.039e-01  9.039e-01      2907 -3.689e+03 -3.704e+03                   
#> Path [9] :Best Iter: [41] ELBO (-3689.181986) evaluations: (2907) 
#> Path [10] :Initial log joint density = -481621.153115 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.911e-03   2.270e-01    8.292e-01  8.292e-01      3147 -3.688e+03 -3.692e+03                   
#> Path [10] :Best Iter: [48] ELBO (-3687.687483) evaluations: (3147) 
#> Path [11] :Initial log joint density = -481506.369872 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.676e-03   2.385e-01    7.585e-01  7.585e-01      2962 -3.691e+03 -3.705e+03                   
#> Path [11] :Best Iter: [50] ELBO (-3691.458202) evaluations: (2962) 
#> Path [12] :Initial log joint density = -481869.173530 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.498e-02   3.048e-01    1.000e+00  1.000e+00      3079 -3.692e+03 -3.687e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3687.104582) evaluations: (3079) 
#> Path [13] :Initial log joint density = -484355.440014 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.097e-02   3.290e-01    1.000e+00  1.000e+00      3276 -3.682e+03 -3.693e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3682.124264) evaluations: (3276) 
#> Path [14] :Initial log joint density = -481628.769914 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.373e-02   2.749e-01    8.848e-01  8.848e-01      2894 -3.692e+03 -3.700e+03                   
#> Path [14] :Best Iter: [46] ELBO (-3692.120490) evaluations: (2894) 
#> Path [15] :Initial log joint density = -482052.416554 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.197e-02   2.027e-01    1.000e+00  1.000e+00      3329 -3.682e+03 -3.682e+03                   
#> Path [15] :Best Iter: [58] ELBO (-3681.706783) evaluations: (3329) 
#> Path [16] :Initial log joint density = -485455.505455 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.001e-02   3.033e-01    1.000e+00  1.000e+00      3434 -3.683e+03 -3.686e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3682.819167) evaluations: (3434) 
#> Path [17] :Initial log joint density = -482361.772552 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.131e-03   2.154e-01    8.190e-01  8.190e-01      3102 -3.690e+03 -3.698e+03                   
#> Path [17] :Best Iter: [53] ELBO (-3689.783287) evaluations: (3102) 
#> Path [18] :Initial log joint density = -480930.847723 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      6.849e-03   2.324e-01    1.000e+00  1.000e+00      2672 -3.690e+03 -3.693e+03                   
#> Path [18] :Best Iter: [42] ELBO (-3690.458861) evaluations: (2672) 
#> Path [19] :Initial log joint density = -482764.630209 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.176e-02   3.193e-01    1.000e+00  1.000e+00      3325 -3.684e+03 -3.690e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3683.598508) evaluations: (3325) 
#> Path [20] :Initial log joint density = -481570.875174 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.224e-03   2.185e-01    1.000e+00  1.000e+00      3273 -3.682e+03 -3.695e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3682.101106) evaluations: (3273) 
#> Path [21] :Initial log joint density = -482117.709853 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.694e-03   2.981e-01    1.000e+00  1.000e+00      3503 -3.685e+03 -3.686e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3684.567952) evaluations: (3503) 
#> Path [22] :Initial log joint density = -481345.834300 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.977e-03   2.998e-01    1.000e+00  1.000e+00      3100 -3.690e+03 -3.692e+03                   
#> Path [22] :Best Iter: [54] ELBO (-3690.201707) evaluations: (3100) 
#> Path [23] :Initial log joint density = -481551.501125 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.166e-02   2.012e-01    1.000e+00  1.000e+00      3169 -3.691e+03 -3.688e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3687.953262) evaluations: (3169) 
#> Path [24] :Initial log joint density = -481762.352955 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      2.371e-03   2.620e-01    6.618e-01  6.618e-01      3854 -3.685e+03 -3.697e+03                   
#> Path [24] :Best Iter: [62] ELBO (-3684.593868) evaluations: (3854) 
#> Path [25] :Initial log joint density = -481638.937653 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.091e-03   2.059e-01    7.974e-01  7.974e-01      3590 -3.682e+03 -3.694e+03                   
#> Path [25] :Best Iter: [58] ELBO (-3681.751348) evaluations: (3590) 
#> Path [26] :Initial log joint density = -481536.476555 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.418e-02   2.702e-01    1.000e+00  1.000e+00      2982 -3.690e+03 -3.685e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3685.228293) evaluations: (2982) 
#> Path [27] :Initial log joint density = -482081.544962 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.717e-03   2.352e-01    1.000e+00  1.000e+00      3538 -3.685e+03 -3.686e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3685.330156) evaluations: (3538) 
#> Path [28] :Initial log joint density = -481872.932284 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.748e-03   2.294e-01    7.212e-01  7.212e-01      3270 -3.683e+03 -3.696e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3682.827285) evaluations: (3270) 
#> Path [29] :Initial log joint density = -483425.873062 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.045e-02   3.706e-01    1.000e+00  1.000e+00      3338 -3.682e+03 -3.691e+03                   
#> Path [29] :Best Iter: [57] ELBO (-3682.386724) evaluations: (3338) 
#> Path [30] :Initial log joint density = -481843.014569 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.452e-02   2.502e-01    1.000e+00  1.000e+00      2991 -3.690e+03 -3.694e+03                   
#> Path [30] :Best Iter: [52] ELBO (-3689.938700) evaluations: (2991) 
#> Path [31] :Initial log joint density = -481641.206361 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.809e-02   2.129e-01    1.000e+00  1.000e+00      2940 -3.693e+03 -3.694e+03                   
#> Path [31] :Best Iter: [52] ELBO (-3693.006372) evaluations: (2940) 
#> Path [32] :Initial log joint density = -481864.392271 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.005e-02   2.435e-01    1.000e+00  1.000e+00      3115 -3.692e+03 -3.699e+03                   
#> Path [32] :Best Iter: [43] ELBO (-3692.450213) evaluations: (3115) 
#> Path [33] :Initial log joint density = -481480.128804 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.136e-03   1.857e-01    9.029e-01  9.029e-01      3077 -3.689e+03 -3.692e+03                   
#> Path [33] :Best Iter: [52] ELBO (-3688.562373) evaluations: (3077) 
#> Path [34] :Initial log joint density = -483287.367826 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.543e-03   2.240e-01    6.567e-01  6.567e-01      3327 -3.684e+03 -3.698e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3683.688634) evaluations: (3327) 
#> Path [35] :Initial log joint density = -481829.639080 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.771e-03   1.979e-01    1.000e+00  1.000e+00      3426 -3.682e+03 -3.692e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3682.157161) evaluations: (3426) 
#> Path [36] :Initial log joint density = -481413.581204 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.407e-03   2.961e-01    8.279e-01  8.279e-01      2902 -3.693e+03 -3.710e+03                   
#> Path [36] :Best Iter: [41] ELBO (-3692.558995) evaluations: (2902) 
#> Path [37] :Initial log joint density = -481462.205053 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.344e-03   1.888e-01    1.000e+00  1.000e+00      3129 -3.691e+03 -3.690e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3690.389534) evaluations: (3129) 
#> Path [38] :Initial log joint density = -483012.707860 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.400e-03   2.874e-01    8.237e-01  8.237e-01      3335 -3.682e+03 -3.695e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3682.290586) evaluations: (3335) 
#> Path [39] :Initial log joint density = -481808.554745 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.249e-03   1.833e-01    1.000e+00  1.000e+00      3576 -3.688e+03 -3.696e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3687.872318) evaluations: (3576) 
#> Path [40] :Initial log joint density = -481662.209275 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      5.103e-03   2.094e-01    1.000e+00  1.000e+00      2761 -3.691e+03 -3.695e+03                   
#> Path [40] :Best Iter: [47] ELBO (-3690.584168) evaluations: (2761) 
#> Path [41] :Initial log joint density = -481572.261717 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.352e-02   2.432e-01    1.000e+00  1.000e+00      3065 -3.692e+03 -3.690e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3690.341451) evaluations: (3065) 
#> Path [42] :Initial log joint density = -481606.383392 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.011e-03   1.916e-01    1.000e+00  1.000e+00      3284 -3.685e+03 -3.684e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3684.005011) evaluations: (3284) 
#> Path [43] :Initial log joint density = -481655.815610 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.294e-02   2.660e-01    1.000e+00  1.000e+00      3478 -3.684e+03 -3.685e+03                   
#> Path [43] :Best Iter: [60] ELBO (-3683.651113) evaluations: (3478) 
#> Path [44] :Initial log joint density = -481420.687200 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.728e-03   2.023e-01    8.726e-01  8.726e-01      3180 -3.686e+03 -3.691e+03                   
#> Path [44] :Best Iter: [55] ELBO (-3686.160760) evaluations: (3180) 
#> Path [45] :Initial log joint density = -484292.762473 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.929e-02   2.978e-01    9.575e-01  9.575e-01      3578 -3.682e+03 -3.690e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3682.180075) evaluations: (3578) 
#> Path [46] :Initial log joint density = -481903.575735 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.802e-03   2.234e-01    1.000e+00  1.000e+00      3341 -3.685e+03 -3.688e+03                   
#> Path [46] :Best Iter: [57] ELBO (-3684.671760) evaluations: (3341) 
#> Path [47] :Initial log joint density = -483841.546837 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.670e-03   2.810e-01    9.450e-01  9.450e-01      3662 -3.683e+03 -3.693e+03                   
#> Path [47] :Best Iter: [59] ELBO (-3683.416251) evaluations: (3662) 
#> Path [48] :Initial log joint density = -481581.736831 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.116e-03   2.727e-01    1.000e+00  1.000e+00      2908 -3.689e+03 -3.699e+03                   
#> Path [48] :Best Iter: [48] ELBO (-3688.532794) evaluations: (2908) 
#> Path [49] :Initial log joint density = -481750.581890 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.590e-03   1.930e-01    1.000e+00  1.000e+00      3082 -3.691e+03 -3.687e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3686.544096) evaluations: (3082) 
#> Path [50] :Initial log joint density = -487102.844626 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.137e-02   1.921e-01    8.032e-01  8.032e-01      3510 -3.687e+03 -3.691e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3687.061160) evaluations: (3510) 
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
