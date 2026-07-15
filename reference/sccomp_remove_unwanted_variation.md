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
#> Path [1] :Initial log joint density = -483880.423727 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.853e-03   2.349e-01    1.000e+00  1.000e+00      3350 -3.695e+03 -3.700e+03                   
#> Path [1] :Best Iter: [49] ELBO (-3694.874125) evaluations: (3350) 
#> Path [2] :Initial log joint density = -481821.668042 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      8.580e-03   1.687e-01    9.557e-01  9.557e-01      4437 -3.683e+03 -3.690e+03                   
#> Path [2] :Best Iter: [67] ELBO (-3683.411454) evaluations: (4437) 
#> Path [3] :Initial log joint density = -481717.551194 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.185e-02   2.172e-01    1.000e+00  1.000e+00      4689 -3.683e+03 -3.691e+03                   
#> Path [3] :Best Iter: [69] ELBO (-3683.379572) evaluations: (4689) 
#> Path [4] :Initial log joint density = -481788.102709 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      9.841e-03   2.237e-01    1.000e+00  1.000e+00      3332 -3.688e+03 -3.693e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3688.464000) evaluations: (3332) 
#> Path [5] :Initial log joint density = -481798.809288 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.287e-02   2.157e-01    1.000e+00  1.000e+00      5092 -3.684e+03 -3.685e+03                   
#> Path [5] :Best Iter: [75] ELBO (-3683.728772) evaluations: (5092) 
#> Path [6] :Initial log joint density = -482287.192260 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      7.527e-03   2.124e-01    1.000e+00  1.000e+00      4278 -3.685e+03 -3.693e+03                   
#> Path [6] :Best Iter: [66] ELBO (-3685.176442) evaluations: (4278) 
#> Path [7] :Initial log joint density = -481374.033864 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      9.633e-03   2.248e-01    8.514e-01  8.514e-01      3159 -3.693e+03 -3.701e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3692.987350) evaluations: (3159) 
#> Path [8] :Initial log joint density = -481411.417774 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      8.325e-03   3.066e-01    9.813e-01  9.813e-01      3163 -3.692e+03 -3.698e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3691.889636) evaluations: (3163) 
#> Path [9] :Initial log joint density = -481564.616306 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.453e-02   2.067e-01    1.000e+00  1.000e+00      4509 -3.688e+03 -3.690e+03                   
#> Path [9] :Best Iter: [64] ELBO (-3688.302519) evaluations: (4509) 
#> Path [10] :Initial log joint density = -483800.482002 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.061e-03   2.242e-01    1.000e+00  1.000e+00      3542 -3.690e+03 -3.689e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3688.505663) evaluations: (3542) 
#> Path [11] :Initial log joint density = -481571.167192 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.176e-03   2.486e-01    1.000e+00  1.000e+00      3271 -3.691e+03 -3.698e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3690.557394) evaluations: (3271) 
#> Path [12] :Initial log joint density = -481380.383799 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      5.505e-03   2.328e-01    7.516e-01  7.516e-01      3293 -3.687e+03 -3.698e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3687.204820) evaluations: (3293) 
#> Path [13] :Initial log joint density = -481841.073873 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.170e-02   2.241e-01    1.000e+00  1.000e+00      4855 -3.686e+03 -3.687e+03                   
#> Path [13] :Best Iter: [73] ELBO (-3686.008573) evaluations: (4855) 
#> Path [14] :Initial log joint density = -482242.410501 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.604e-02   1.863e-01    1.000e+00  1.000e+00      5044 -3.681e+03 -3.678e+03                   
#> Path [14] :Best Iter: [75] ELBO (-3677.577661) evaluations: (5044) 
#> Path [15] :Initial log joint density = -481499.785612 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      9.374e-03   3.664e-01    8.647e-01  8.647e-01      2872 -3.695e+03 -3.709e+03                   
#> Path [15] :Best Iter: [50] ELBO (-3694.773830) evaluations: (2872) 
#> Path [16] :Initial log joint density = -481413.590938 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.405e-02   3.592e-01    1.000e+00  1.000e+00      3207 -3.694e+03 -3.692e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3691.966397) evaluations: (3207) 
#> Path [17] :Initial log joint density = -481576.747545 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      8.389e-03   2.237e-01    1.000e+00  1.000e+00      4969 -3.682e+03 -3.688e+03                   
#> Path [17] :Best Iter: [71] ELBO (-3681.739795) evaluations: (4969) 
#> Path [18] :Initial log joint density = -481443.364995 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.659e-03   2.688e-01    1.000e+00  1.000e+00      3167 -3.693e+03 -3.694e+03                   
#> Path [18] :Best Iter: [54] ELBO (-3692.664949) evaluations: (3167) 
#> Path [19] :Initial log joint density = -483181.897979 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.387e-02   2.273e-01    8.717e-01  8.717e-01      4371 -3.684e+03 -3.696e+03                   
#> Path [19] :Best Iter: [62] ELBO (-3683.942096) evaluations: (4371) 
#> Path [20] :Initial log joint density = -482325.067190 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.157e-02   2.077e-01    8.394e-01  8.394e-01      4210 -3.684e+03 -3.695e+03                   
#> Path [20] :Best Iter: [65] ELBO (-3684.338670) evaluations: (4210) 
#> Path [21] :Initial log joint density = -481887.872066 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.308e-02   2.019e-01    1.000e+00  1.000e+00      4252 -3.689e+03 -3.689e+03                   
#> Path [21] :Best Iter: [67] ELBO (-3688.669407) evaluations: (4252) 
#> Path [22] :Initial log joint density = -481301.774741 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.521e-02   2.215e-01    1.000e+00  1.000e+00      4178 -3.684e+03 -3.685e+03                   
#> Path [22] :Best Iter: [65] ELBO (-3684.240964) evaluations: (4178) 
#> Path [23] :Initial log joint density = -481582.017109 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.057e-02   2.114e-01    1.000e+00  1.000e+00      3158 -3.693e+03 -3.689e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3688.549385) evaluations: (3158) 
#> Path [24] :Initial log joint density = -484360.665999 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      8.124e-03   1.416e-01    1.000e+00  1.000e+00      4910 -3.683e+03 -3.693e+03                   
#> Path [24] :Best Iter: [66] ELBO (-3683.461641) evaluations: (4910) 
#> Path [25] :Initial log joint density = -481476.465906 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      4.703e-03   1.922e-01    1.000e+00  1.000e+00      3307 -3.695e+03 -3.694e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3693.885419) evaluations: (3307) 
#> Path [26] :Initial log joint density = -481843.948460 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      2.975e-03   3.211e-01    5.847e-01  5.847e-01      3326 -3.690e+03 -3.705e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3689.584693) evaluations: (3326) 
#> Path [27] :Initial log joint density = -481400.114779 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      2.164e-02   2.686e-01    1.000e+00  1.000e+00      3960 -3.684e+03 -3.687e+03                   
#> Path [27] :Best Iter: [59] ELBO (-3684.467792) evaluations: (3960) 
#> Path [28] :Initial log joint density = -481446.803548 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.347e-03   2.197e-01    1.000e+00  1.000e+00      3275 -3.692e+03 -3.697e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3691.741352) evaluations: (3275) 
#> Path [29] :Initial log joint density = -481680.327844 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.025e-02   1.998e-01    1.000e+00  1.000e+00      4289 -3.687e+03 -3.694e+03                   
#> Path [29] :Best Iter: [62] ELBO (-3687.498939) evaluations: (4289) 
#> Path [30] :Initial log joint density = -484568.437984 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      8.210e-03   1.860e-01    1.000e+00  1.000e+00      3627 -3.686e+03 -3.696e+03                   
#> Path [30] :Best Iter: [60] ELBO (-3685.658918) evaluations: (3627) 
#> Path [31] :Initial log joint density = -482205.992522 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.406e-03   2.427e-01    8.678e-01  8.678e-01      3798 -3.688e+03 -3.699e+03                   
#> Path [31] :Best Iter: [61] ELBO (-3688.384975) evaluations: (3798) 
#> Path [32] :Initial log joint density = -481532.477318 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      3.734e-03   1.677e-01    6.539e-01  6.539e-01      4813 -3.684e+03 -3.698e+03                   
#> Path [32] :Best Iter: [68] ELBO (-3684.425771) evaluations: (4813) 
#> Path [33] :Initial log joint density = -481340.608043 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.787e+05      1.249e-02   3.900e-01    1.000e+00  1.000e+00      2716 -3.694e+03 -3.698e+03                   
#> Path [33] :Best Iter: [41] ELBO (-3694.477210) evaluations: (2716) 
#> Path [34] :Initial log joint density = -481595.592409 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      4.430e-03   1.912e-01    1.000e+00  1.000e+00      3036 -3.694e+03 -3.706e+03                   
#> Path [34] :Best Iter: [50] ELBO (-3693.801245) evaluations: (3036) 
#> Path [35] :Initial log joint density = -482180.870139 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      8.588e-03   2.279e-01    1.000e+00  1.000e+00      4947 -3.688e+03 -3.695e+03                   
#> Path [35] :Best Iter: [74] ELBO (-3687.669421) evaluations: (4947) 
#> Path [36] :Initial log joint density = -481292.884196 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.064e-02   2.069e-01    1.000e+00  1.000e+00      4173 -3.684e+03 -3.692e+03                   
#> Path [36] :Best Iter: [64] ELBO (-3683.871216) evaluations: (4173) 
#> Path [37] :Initial log joint density = -485072.432955 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      8.140e-03   3.346e-01    4.492e-01  1.000e+00      3301 -3.688e+03 -3.699e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3687.804865) evaluations: (3301) 
#> Path [38] :Initial log joint density = -481943.407331 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      3.998e-03   1.941e-01    7.138e-01  7.138e-01      4982 -3.685e+03 -3.703e+03                   
#> Path [38] :Best Iter: [72] ELBO (-3685.078769) evaluations: (4982) 
#> Path [39] :Initial log joint density = -481629.308280 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.324e-02   2.308e-01    1.000e+00  1.000e+00      2978 -3.695e+03 -3.701e+03                   
#> Path [39] :Best Iter: [53] ELBO (-3694.752905) evaluations: (2978) 
#> Path [40] :Initial log joint density = -481637.214247 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.974e-02   2.830e-01    9.738e-01  9.738e-01      3364 -3.689e+03 -3.697e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3688.854521) evaluations: (3364) 
#> Path [41] :Initial log joint density = -482588.274572 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.539e-02   2.243e-01    1.000e+00  1.000e+00      4340 -3.684e+03 -3.683e+03                   
#> Path [41] :Best Iter: [67] ELBO (-3683.364591) evaluations: (4340) 
#> Path [42] :Initial log joint density = -481641.615888 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.495e-02   1.908e-01    1.000e+00  1.000e+00      4590 -3.683e+03 -3.688e+03                   
#> Path [42] :Best Iter: [65] ELBO (-3683.316717) evaluations: (4590) 
#> Path [43] :Initial log joint density = -485685.500299 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.990e-02   2.474e-01    9.422e-01  9.422e-01      4909 -3.683e+03 -3.690e+03                   
#> Path [43] :Best Iter: [72] ELBO (-3682.671499) evaluations: (4909) 
#> Path [44] :Initial log joint density = -481920.640105 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.100e-02   3.673e-01    4.043e-01  1.000e+00      4376 -3.685e+03 -3.695e+03                   
#> Path [44] :Best Iter: [67] ELBO (-3685.186082) evaluations: (4376) 
#> Path [45] :Initial log joint density = -481182.736116 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.960e-03   2.141e-01    9.748e-01  9.748e-01      3209 -3.693e+03 -3.702e+03                   
#> Path [45] :Best Iter: [54] ELBO (-3692.812398) evaluations: (3209) 
#> Path [46] :Initial log joint density = -481465.873751 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      2.394e-02   1.627e-01    1.000e+00  1.000e+00      4250 -3.684e+03 -3.686e+03                   
#> Path [46] :Best Iter: [65] ELBO (-3683.694665) evaluations: (4250) 
#> Path [47] :Initial log joint density = -481515.563332 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.356e-02   1.704e-01    1.000e+00  1.000e+00      4311 -3.685e+03 -3.685e+03                   
#> Path [47] :Best Iter: [65] ELBO (-3684.733455) evaluations: (4311) 
#> Path [48] :Initial log joint density = -481761.839262 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      6.632e-03   1.953e-01    4.555e-01  1.000e+00      3061 -3.695e+03 -3.705e+03                   
#> Path [48] :Best Iter: [42] ELBO (-3695.194277) evaluations: (3061) 
#> Path [49] :Initial log joint density = -481727.131929 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.256e-02   3.057e-01    9.866e-01  9.866e-01      3075 -3.694e+03 -3.700e+03                   
#> Path [49] :Best Iter: [53] ELBO (-3693.840781) evaluations: (3075) 
#> Path [50] :Initial log joint density = -482980.819244 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      2.667e-03   2.269e-01    6.204e-01  6.204e-01      3577 -3.687e+03 -3.701e+03                   
#> Path [50] :Best Iter: [58] ELBO (-3687.314445) evaluations: (3577) 
#> Finished in  15.2 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.556 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.558 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
