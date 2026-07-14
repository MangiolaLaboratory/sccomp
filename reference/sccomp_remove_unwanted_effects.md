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
#> Path [1] :Initial log joint density = -481616.311559 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.359e-02   2.529e-01    1.000e+00  1.000e+00      3220 -3.685e+03 -3.683e+03                   
#> Path [1] :Best Iter: [58] ELBO (-3683.468338) evaluations: (3220) 
#> Path [2] :Initial log joint density = -481743.918103 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.325e-03   1.989e-01    1.000e+00  1.000e+00      3422 -3.683e+03 -3.686e+03                   
#> Path [2] :Best Iter: [58] ELBO (-3683.291583) evaluations: (3422) 
#> Path [3] :Initial log joint density = -481893.559388 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.592e-02   2.757e-01    1.000e+00  1.000e+00      3532 -3.683e+03 -3.686e+03                   
#> Path [3] :Best Iter: [59] ELBO (-3682.722391) evaluations: (3532) 
#> Path [4] :Initial log joint density = -481637.091079 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.820e-03   1.990e-01    7.928e-01  7.928e-01      3563 -3.682e+03 -3.690e+03                   
#> Path [4] :Best Iter: [58] ELBO (-3682.098283) evaluations: (3563) 
#> Path [5] :Initial log joint density = -481657.492370 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.901e-03   2.904e-01    9.336e-01  9.336e-01      3366 -3.685e+03 -3.694e+03                   
#> Path [5] :Best Iter: [58] ELBO (-3684.681327) evaluations: (3366) 
#> Path [6] :Initial log joint density = -482371.019663 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.183e-03   1.998e-01    8.994e-01  8.994e-01      3304 -3.683e+03 -3.699e+03                   
#> Path [6] :Best Iter: [56] ELBO (-3682.575101) evaluations: (3304) 
#> Path [7] :Initial log joint density = -481817.369292 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.602e-03   1.855e-01    1.000e+00  1.000e+00      3486 -3.681e+03 -3.681e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3680.512322) evaluations: (3486) 
#> Path [8] :Initial log joint density = -481957.539467 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.326e-03   2.001e-01    1.000e+00  1.000e+00      3640 -3.685e+03 -3.685e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3684.976580) evaluations: (3640) 
#> Path [9] :Initial log joint density = -481560.577923 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.926e-03   1.264e-01    1.000e+00  1.000e+00      3331 -3.687e+03 -3.688e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3686.633133) evaluations: (3331) 
#> Path [10] :Initial log joint density = -481464.932411 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.694e-03   2.000e-01    7.063e-01  7.063e-01      3649 -3.679e+03 -3.695e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3678.720205) evaluations: (3649) 
#> Path [11] :Initial log joint density = -481572.655030 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.435e-03   2.229e-01    1.000e+00  1.000e+00      2783 -3.691e+03 -3.693e+03                   
#> Path [11] :Best Iter: [40] ELBO (-3691.365721) evaluations: (2783) 
#> Path [12] :Initial log joint density = -481678.938098 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.826e-03   2.580e-01    1.000e+00  1.000e+00      3112 -3.691e+03 -3.694e+03                   
#> Path [12] :Best Iter: [44] ELBO (-3690.609704) evaluations: (3112) 
#> Path [13] :Initial log joint density = -482179.849017 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.689e-03   1.764e-01    1.000e+00  1.000e+00      3685 -3.683e+03 -3.682e+03                   
#> Path [13] :Best Iter: [62] ELBO (-3682.435577) evaluations: (3685) 
#> Path [14] :Initial log joint density = -481750.089564 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.653e-03   2.512e-01    6.427e-01  6.427e-01      3392 -3.686e+03 -3.694e+03                   
#> Path [14] :Best Iter: [56] ELBO (-3685.995961) evaluations: (3392) 
#> Path [15] :Initial log joint density = -482005.277867 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.098e-02   2.551e-01    1.000e+00  1.000e+00      3466 -3.683e+03 -3.685e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3683.131731) evaluations: (3466) 
#> Path [16] :Initial log joint density = -481550.811852 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.191e-02   2.020e-01    1.000e+00  1.000e+00      3332 -3.684e+03 -3.688e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3684.265058) evaluations: (3332) 
#> Path [17] :Initial log joint density = -481756.751558 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.034e-02   2.297e-01    7.847e-01  7.847e-01      3406 -3.686e+03 -3.693e+03                   
#> Path [17] :Best Iter: [57] ELBO (-3686.007634) evaluations: (3406) 
#> Path [18] :Initial log joint density = -481638.318561 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      8.048e-03   2.902e-01    1.000e+00  1.000e+00      2758 -3.692e+03 -3.704e+03                   
#> Path [18] :Best Iter: [37] ELBO (-3692.370129) evaluations: (2758) 
#> Path [19] :Initial log joint density = -481679.297111 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.513e-02   2.533e-01    1.000e+00  1.000e+00      3923 -3.682e+03 -3.685e+03                   
#> Path [19] :Best Iter: [59] ELBO (-3682.370916) evaluations: (3923) 
#> Path [20] :Initial log joint density = -481917.330300 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.933e-03   2.071e-01    1.000e+00  1.000e+00      3493 -3.681e+03 -3.687e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3681.484821) evaluations: (3493) 
#> Path [21] :Initial log joint density = -481749.309441 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.347e-02   2.329e-01    1.000e+00  1.000e+00      3278 -3.685e+03 -3.690e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3684.970421) evaluations: (3278) 
#> Path [22] :Initial log joint density = -481453.729917 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.219e-02   2.001e-01    1.000e+00  1.000e+00      3686 -3.680e+03 -3.681e+03                   
#> Path [22] :Best Iter: [59] ELBO (-3680.319030) evaluations: (3686) 
#> Path [23] :Initial log joint density = -481751.509370 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.193e-02   3.328e-01    1.000e+00  1.000e+00      3427 -3.682e+03 -3.689e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3682.495493) evaluations: (3427) 
#> Path [24] :Initial log joint density = -481554.397554 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.160e-02   2.808e-01    1.000e+00  1.000e+00      3128 -3.688e+03 -3.684e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3684.495448) evaluations: (3128) 
#> Path [25] :Initial log joint density = -484424.649722 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.323e-03   2.078e-01    1.000e+00  1.000e+00      3458 -3.684e+03 -3.693e+03                   
#> Path [25] :Best Iter: [56] ELBO (-3683.798032) evaluations: (3458) 
#> Path [26] :Initial log joint density = -482697.678493 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.450e-02   3.770e-01    1.000e+00  1.000e+00      3431 -3.682e+03 -3.689e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3681.948632) evaluations: (3431) 
#> Path [27] :Initial log joint density = -481382.537477 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.551e-02   2.748e-01    1.000e+00  1.000e+00      3189 -3.686e+03 -3.684e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3684.122623) evaluations: (3189) 
#> Path [28] :Initial log joint density = -482035.793874 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.152e-02   1.777e-01    1.000e+00  1.000e+00      3346 -3.685e+03 -3.683e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3683.114850) evaluations: (3346) 
#> Path [29] :Initial log joint density = -482540.399977 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.238e-03   1.965e-01    1.000e+00  1.000e+00      3192 -3.680e+03 -3.684e+03                   
#> Path [29] :Best Iter: [56] ELBO (-3680.260618) evaluations: (3192) 
#> Path [30] :Initial log joint density = -481630.636974 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.491e-03   2.325e-01    1.000e+00  1.000e+00      3271 -3.685e+03 -3.688e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3685.214210) evaluations: (3271) 
#> Path [31] :Initial log joint density = -481848.825453 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.077e-02   2.679e-01    1.000e+00  1.000e+00      3044 -3.689e+03 -3.693e+03                   
#> Path [31] :Best Iter: [46] ELBO (-3688.573003) evaluations: (3044) 
#> Path [32] :Initial log joint density = -482401.166348 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.468e-03   2.071e-01    1.000e+00  1.000e+00      2870 -3.694e+03 -3.701e+03                   
#> Path [32] :Best Iter: [35] ELBO (-3694.033214) evaluations: (2870) 
#> Path [33] :Initial log joint density = -481809.063620 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.426e-02   2.739e-01    1.000e+00  1.000e+00      3323 -3.686e+03 -3.684e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3683.879961) evaluations: (3323) 
#> Path [34] :Initial log joint density = -481793.810447 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.411e-03   2.022e-01    1.000e+00  1.000e+00      3492 -3.687e+03 -3.688e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3686.613284) evaluations: (3492) 
#> Path [35] :Initial log joint density = -481673.125300 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.401e-03   2.200e-01    7.764e-01  7.764e-01      3077 -3.691e+03 -3.697e+03                   
#> Path [35] :Best Iter: [40] ELBO (-3690.846875) evaluations: (3077) 
#> Path [36] :Initial log joint density = -481611.546740 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.278e-02   3.843e-01    1.000e+00  1.000e+00      3166 -3.684e+03 -3.692e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3683.626586) evaluations: (3166) 
#> Path [37] :Initial log joint density = -481860.915239 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.080e-03   2.545e-01    1.000e+00  1.000e+00      3339 -3.683e+03 -3.688e+03                   
#> Path [37] :Best Iter: [58] ELBO (-3683.120056) evaluations: (3339) 
#> Path [38] :Initial log joint density = -481638.835171 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.841e-03   1.604e-01    1.000e+00  1.000e+00      3321 -3.685e+03 -3.689e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3684.786464) evaluations: (3321) 
#> Path [39] :Initial log joint density = -481488.604529 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.058e-02   4.226e-01    9.426e-01  9.426e-01      2880 -3.693e+03 -3.704e+03                   
#> Path [39] :Best Iter: [51] ELBO (-3692.537577) evaluations: (2880) 
#> Path [40] :Initial log joint density = -481984.507708 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.250e-02   3.383e-01    1.000e+00  1.000e+00      3453 -3.682e+03 -3.689e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3682.273366) evaluations: (3453) 
#> Path [41] :Initial log joint density = -489207.467148 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.311e-03   1.962e-01    1.000e+00  1.000e+00      3135 -3.690e+03 -3.695e+03                   
#> Path [41] :Best Iter: [43] ELBO (-3689.661126) evaluations: (3135) 
#> Path [42] :Initial log joint density = -481475.933637 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      3.591e-03   3.245e-01    5.680e-01  5.680e-01      3558 -3.682e+03 -3.694e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3681.625644) evaluations: (3558) 
#> Path [43] :Initial log joint density = -481627.797447 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.350e-03   1.744e-01    1.000e+00  1.000e+00      3020 -3.692e+03 -3.695e+03                   
#> Path [43] :Best Iter: [51] ELBO (-3692.313464) evaluations: (3020) 
#> Path [44] :Initial log joint density = -481645.079340 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.238e-03   1.747e-01    1.000e+00  1.000e+00      2910 -3.691e+03 -3.697e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3690.832490) evaluations: (2910) 
#> Path [45] :Initial log joint density = -481686.570970 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.038e-02   3.279e-01    8.303e-01  8.303e-01      3739 -3.685e+03 -3.694e+03                   
#> Path [45] :Best Iter: [56] ELBO (-3684.557529) evaluations: (3739) 
#> Path [46] :Initial log joint density = -482011.548668 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.858e-03   3.461e-01    5.994e-01  5.994e-01      3620 -3.683e+03 -3.697e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3682.821957) evaluations: (3620) 
#> Path [47] :Initial log joint density = -481572.787679 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      3.893e-03   3.298e-01    6.457e-01  6.457e-01      3624 -3.685e+03 -3.695e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3684.776421) evaluations: (3624) 
#> Path [48] :Initial log joint density = -481679.972667 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.540e-03   2.322e-01    1.000e+00  1.000e+00      3054 -3.688e+03 -3.685e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3684.587766) evaluations: (3054) 
#> Path [49] :Initial log joint density = -484454.418863 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.167e-02   3.263e-01    1.000e+00  1.000e+00      3405 -3.683e+03 -3.698e+03                   
#> Path [49] :Best Iter: [58] ELBO (-3682.790858) evaluations: (3405) 
#> Path [50] :Initial log joint density = -481642.227406 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.982e-03   1.708e-01    1.000e+00  1.000e+00      3319 -3.681e+03 -3.690e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3681.418723) evaluations: (3319) 
#> Finished in  13.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: calculating residuals
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.567 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.567 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
