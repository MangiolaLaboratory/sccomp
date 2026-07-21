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
#> Path [1] :Initial log joint density = -482535.830703 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.051e-02   2.491e-01    1.000e+00  1.000e+00      3513 -3.687e+03 -3.698e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3686.962295) evaluations: (3513) 
#> Path [2] :Initial log joint density = -481743.379542 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      8.299e-03   1.785e-01    8.961e-01  8.961e-01      5623 -3.681e+03 -3.691e+03                   
#> Path [2] :Best Iter: [76] ELBO (-3680.888414) evaluations: (5623) 
#> Path [3] :Initial log joint density = -482715.935554 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.341e-02   2.624e-01    1.000e+00  1.000e+00      3247 -3.690e+03 -3.687e+03                   
#> Path [3] :Best Iter: [57] ELBO (-3687.299458) evaluations: (3247) 
#> Path [4] :Initial log joint density = -481839.052831 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.070e-02   1.423e-01    1.000e+00  1.000e+00      4410 -3.684e+03 -3.685e+03                   
#> Path [4] :Best Iter: [65] ELBO (-3683.851556) evaluations: (4410) 
#> Path [5] :Initial log joint density = -481771.742582 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      8.671e-03   2.228e-01    1.000e+00  1.000e+00      3975 -3.684e+03 -3.692e+03                   
#> Path [5] :Best Iter: [61] ELBO (-3683.793735) evaluations: (3975) 
#> Path [6] :Initial log joint density = -483204.973713 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.722e-03   2.539e-01    6.343e-01  6.343e-01      4461 -3.681e+03 -3.697e+03                   
#> Path [6] :Best Iter: [68] ELBO (-3681.437366) evaluations: (4461) 
#> Path [7] :Initial log joint density = -481460.540810 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      7.631e-03   1.920e-01    1.000e+00  1.000e+00      4902 -3.687e+03 -3.700e+03                   
#> Path [7] :Best Iter: [72] ELBO (-3687.426325) evaluations: (4902) 
#> Path [8] :Initial log joint density = -482072.191225 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      9.180e-03   2.097e-01    1.000e+00  1.000e+00      5078 -3.688e+03 -3.688e+03                   
#> Path [8] :Best Iter: [74] ELBO (-3687.895817) evaluations: (5078) 
#> Path [9] :Initial log joint density = -481494.219662 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.183e-02   2.224e-01    1.000e+00  1.000e+00      4176 -3.686e+03 -3.686e+03                   
#> Path [9] :Best Iter: [66] ELBO (-3685.887973) evaluations: (4176) 
#> Path [10] :Initial log joint density = -481789.349024 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.097e-02   2.772e-01    1.000e+00  1.000e+00      4390 -3.686e+03 -3.694e+03                   
#> Path [10] :Best Iter: [64] ELBO (-3686.244319) evaluations: (4390) 
#> Path [11] :Initial log joint density = -481585.016950 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      9.462e-03   1.841e-01    1.000e+00  1.000e+00      3083 -3.694e+03 -3.691e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3690.925428) evaluations: (3083) 
#> Path [12] :Initial log joint density = -482179.063554 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.348e-03   1.862e-01    1.000e+00  1.000e+00      4316 -3.688e+03 -3.692e+03                   
#> Path [12] :Best Iter: [65] ELBO (-3688.364448) evaluations: (4316) 
#> Path [13] :Initial log joint density = -481695.596803 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      6.144e-03   2.538e-01    4.840e-01  4.840e-01      4526 -3.686e+03 -3.690e+03                   
#> Path [13] :Best Iter: [68] ELBO (-3685.535451) evaluations: (4526) 
#> Path [14] :Initial log joint density = -481481.518519 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      9.348e-03   2.345e-01    1.000e+00  1.000e+00      4349 -3.687e+03 -3.687e+03                   
#> Path [14] :Best Iter: [69] ELBO (-3686.544764) evaluations: (4349) 
#> Path [15] :Initial log joint density = -483911.629168 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.197e-02   2.393e-01    1.000e+00  1.000e+00      3831 -3.685e+03 -3.692e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3685.404213) evaluations: (3831) 
#> Path [16] :Initial log joint density = -481787.791318 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      2.264e-03   2.285e-01    5.734e-01  5.734e-01      4905 -3.685e+03 -3.697e+03                   
#> Path [16] :Best Iter: [73] ELBO (-3684.555504) evaluations: (4905) 
#> Path [17] :Initial log joint density = -481870.260794 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      2.011e-02   1.916e-01    1.000e+00  1.000e+00      4975 -3.682e+03 -3.684e+03                   
#> Path [17] :Best Iter: [75] ELBO (-3682.368738) evaluations: (4975) 
#> Path [18] :Initial log joint density = -481886.452891 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.048e-03   1.881e-01    1.000e+00  1.000e+00      4569 -3.686e+03 -3.695e+03                   
#> Path [18] :Best Iter: [68] ELBO (-3686.172301) evaluations: (4569) 
#> Path [19] :Initial log joint density = -481571.275632 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      8.484e-03   1.711e-01    9.078e-01  9.078e-01      3469 -3.693e+03 -3.698e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3692.853986) evaluations: (3469) 
#> Path [20] :Initial log joint density = -481508.200005 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.110e-03   2.214e-01    6.923e-01  6.923e-01      4269 -3.682e+03 -3.699e+03                   
#> Path [20] :Best Iter: [62] ELBO (-3681.861006) evaluations: (4269) 
#> Path [21] :Initial log joint density = -481688.347249 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.728e-03   2.071e-01    1.000e+00  1.000e+00      4135 -3.690e+03 -3.689e+03                   
#> Path [21] :Best Iter: [67] ELBO (-3689.379584) evaluations: (4135) 
#> Path [22] :Initial log joint density = -481812.512243 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.070e-02   2.524e-01    1.000e+00  1.000e+00      4048 -3.689e+03 -3.690e+03                   
#> Path [22] :Best Iter: [64] ELBO (-3689.102884) evaluations: (4048) 
#> Path [23] :Initial log joint density = -482348.671170 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.278e-03   2.233e-01    1.000e+00  1.000e+00      3717 -3.690e+03 -3.693e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3689.558643) evaluations: (3717) 
#> Path [24] :Initial log joint density = -481534.535348 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.164e-03   2.339e-01    9.239e-01  9.239e-01      3738 -3.686e+03 -3.696e+03                   
#> Path [24] :Best Iter: [61] ELBO (-3686.469193) evaluations: (3738) 
#> Path [25] :Initial log joint density = -487318.726233 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.645e-03   1.938e-01    1.000e+00  1.000e+00      4666 -3.683e+03 -3.691e+03                   
#> Path [25] :Best Iter: [68] ELBO (-3682.696619) evaluations: (4666) 
#> Path [26] :Initial log joint density = -481721.134191 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      5.154e-03   1.866e-01    1.000e+00  1.000e+00      3376 -3.687e+03 -3.698e+03                   
#> Path [26] :Best Iter: [55] ELBO (-3686.556899) evaluations: (3376) 
#> Path [27] :Initial log joint density = -481605.594279 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.691e-02   1.871e-01    1.000e+00  1.000e+00      3901 -3.684e+03 -3.685e+03                   
#> Path [27] :Best Iter: [64] ELBO (-3684.355277) evaluations: (3901) 
#> Path [28] :Initial log joint density = -481596.014489 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      5.611e-03   2.673e-01    1.000e+00  1.000e+00      3212 -3.699e+03 -3.701e+03                   
#> Path [28] :Best Iter: [38] ELBO (-3699.080917) evaluations: (3212) 
#> Path [29] :Initial log joint density = -481365.089186 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.923e-03   2.476e-01    5.699e-01  5.699e-01      4474 -3.685e+03 -3.698e+03                   
#> Path [29] :Best Iter: [68] ELBO (-3685.256894) evaluations: (4474) 
#> Path [30] :Initial log joint density = -483992.906998 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      3.884e-03   2.168e-01    5.328e-01  5.328e-01      4473 -3.683e+03 -3.698e+03                   
#> Path [30] :Best Iter: [64] ELBO (-3683.137043) evaluations: (4473) 
#> Path [31] :Initial log joint density = -482038.540459 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.327e-02   4.614e-01    1.000e+00  1.000e+00      4031 -3.690e+03 -3.696e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3689.507301) evaluations: (4031) 
#> Path [32] :Initial log joint density = -481560.927653 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.258e-03   2.162e-01    8.494e-01  8.494e-01      4499 -3.684e+03 -3.696e+03                   
#> Path [32] :Best Iter: [68] ELBO (-3683.734290) evaluations: (4499) 
#> Path [33] :Initial log joint density = -481961.249379 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.060e-02   3.155e-01    1.000e+00  1.000e+00      3745 -3.687e+03 -3.692e+03                   
#> Path [33] :Best Iter: [63] ELBO (-3686.916849) evaluations: (3745) 
#> Path [34] :Initial log joint density = -481991.054571 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.121e-02   3.002e-01    1.000e+00  1.000e+00      4579 -3.686e+03 -3.697e+03                   
#> Path [34] :Best Iter: [70] ELBO (-3685.900810) evaluations: (4579) 
#> Path [35] :Initial log joint density = -481567.744580 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      6.627e-03   2.193e-01    8.289e-01  8.289e-01      3379 -3.691e+03 -3.703e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3691.128150) evaluations: (3379) 
#> Path [36] :Initial log joint density = -481658.280279 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.969e-03   2.674e-01    6.415e-01  6.415e-01      3328 -3.688e+03 -3.704e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3687.983002) evaluations: (3328) 
#> Path [37] :Initial log joint density = -481844.991680 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      4.090e-03   1.389e-01    7.234e-01  7.234e-01      4879 -3.685e+03 -3.698e+03                   
#> Path [37] :Best Iter: [72] ELBO (-3684.845650) evaluations: (4879) 
#> Path [38] :Initial log joint density = -481701.949360 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      6.837e-03   2.206e-01    9.877e-01  9.877e-01      3078 -3.695e+03 -3.689e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3688.883638) evaluations: (3078) 
#> Path [39] :Initial log joint density = -482331.299827 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.308e-02   2.204e-01    1.000e+00  1.000e+00      4343 -3.684e+03 -3.688e+03                   
#> Path [39] :Best Iter: [62] ELBO (-3684.104128) evaluations: (4343) 
#> Path [40] :Initial log joint density = -482085.579113 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      9.536e-03   2.389e-01    1.000e+00  1.000e+00      4602 -3.687e+03 -3.696e+03                   
#> Path [40] :Best Iter: [69] ELBO (-3686.811312) evaluations: (4602) 
#> Path [41] :Initial log joint density = -481434.956192 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      5.176e-03   2.715e-01    7.617e-01  7.617e-01      4151 -3.684e+03 -3.693e+03                   
#> Path [41] :Best Iter: [63] ELBO (-3684.162302) evaluations: (4151) 
#> Path [42] :Initial log joint density = -481843.404000 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.905e-02   3.058e-01    1.000e+00  1.000e+00      4554 -3.682e+03 -3.689e+03                   
#> Path [42] :Best Iter: [68] ELBO (-3682.332155) evaluations: (4554) 
#> Path [43] :Initial log joint density = -481473.089461 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      3.793e-02   1.843e-01    1.000e+00  1.000e+00      4471 -3.685e+03 -3.683e+03                   
#> Path [43] :Best Iter: [69] ELBO (-3683.153485) evaluations: (4471) 
#> Path [44] :Initial log joint density = -481477.789040 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      5.565e-03   2.118e-01    7.795e-01  7.795e-01      4149 -3.686e+03 -3.698e+03                   
#> Path [44] :Best Iter: [64] ELBO (-3685.657959) evaluations: (4149) 
#> Path [45] :Initial log joint density = -481294.946179 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.122e-02   1.671e-01    1.000e+00  1.000e+00      5400 -3.680e+03 -3.689e+03                   
#> Path [45] :Best Iter: [75] ELBO (-3680.119238) evaluations: (5400) 
#> Path [46] :Initial log joint density = -484360.810234 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      8.158e-03   1.615e-01    1.000e+00  1.000e+00      4205 -3.689e+03 -3.692e+03                   
#> Path [46] :Best Iter: [64] ELBO (-3688.887836) evaluations: (4205) 
#> Path [47] :Initial log joint density = -481715.126034 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      5.601e-03   1.832e-01    1.000e+00  1.000e+00      3375 -3.687e+03 -3.699e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3687.398261) evaluations: (3375) 
#> Path [48] :Initial log joint density = -482610.451813 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      1.313e-02   2.770e-01    1.000e+00  1.000e+00      3260 -3.692e+03 -3.689e+03                   
#> Path [48] :Best Iter: [57] ELBO (-3689.493413) evaluations: (3260) 
#> Path [49] :Initial log joint density = -481876.028277 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.205e-02   2.726e-01    1.000e+00  1.000e+00      4438 -3.686e+03 -3.690e+03                   
#> Path [49] :Best Iter: [63] ELBO (-3686.368946) evaluations: (4438) 
#> Path [50] :Initial log joint density = -481459.841401 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.319e-03   2.118e-01    1.000e+00  1.000e+00      3417 -3.691e+03 -3.693e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3691.219029) evaluations: (3417) 
#> Finished in  15.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.518 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.512 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
