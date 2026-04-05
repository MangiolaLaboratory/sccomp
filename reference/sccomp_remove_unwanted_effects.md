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
#> Path [1] :Initial log joint density = -481474.464514 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.350e-02   2.365e-01    1.000e+00  1.000e+00      3250 -3.685e+03 -3.687e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3685.235909) evaluations: (3250) 
#> Path [2] :Initial log joint density = -483030.338695 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.001e-02   3.094e-01    9.779e-01  9.779e-01      3338 -3.684e+03 -3.695e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3684.087209) evaluations: (3338) 
#> Path [3] :Initial log joint density = -481572.339319 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.474e-02   2.594e-01    1.000e+00  1.000e+00      3186 -3.689e+03 -3.686e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3685.757967) evaluations: (3186) 
#> Path [4] :Initial log joint density = -483623.572425 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.819e-03   2.353e-01    7.233e-01  7.233e-01      3298 -3.690e+03 -3.696e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3690.297873) evaluations: (3298) 
#> Path [5] :Initial log joint density = -481626.669834 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.693e-03   2.301e-01    6.674e-01  6.674e-01      3547 -3.682e+03 -3.695e+03                   
#> Path [5] :Best Iter: [59] ELBO (-3682.401626) evaluations: (3547) 
#> Path [6] :Initial log joint density = -481645.731927 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.073e-02   3.712e-01    1.000e+00  1.000e+00      2961 -3.690e+03 -3.697e+03                   
#> Path [6] :Best Iter: [51] ELBO (-3690.291877) evaluations: (2961) 
#> Path [7] :Initial log joint density = -481619.437690 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.841e-03   2.355e-01    8.532e-01  8.532e-01      3381 -3.681e+03 -3.696e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3681.468766) evaluations: (3381) 
#> Path [8] :Initial log joint density = -482314.052579 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.525e-03   2.166e-01    7.799e-01  7.799e-01      3331 -3.686e+03 -3.695e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3686.133323) evaluations: (3331) 
#> Path [9] :Initial log joint density = -482060.227865 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.430e-03   2.131e-01    1.000e+00  1.000e+00      3334 -3.683e+03 -3.694e+03                   
#> Path [9] :Best Iter: [56] ELBO (-3683.491327) evaluations: (3334) 
#> Path [10] :Initial log joint density = -481784.365869 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      3.366e-03   2.476e-01    5.826e-01  5.826e-01      2810 -3.687e+03 -3.705e+03                   
#> Path [10] :Best Iter: [45] ELBO (-3687.290867) evaluations: (2810) 
#> Path [11] :Initial log joint density = -481575.919999 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.225e-02   2.181e-01    1.000e+00  1.000e+00      2935 -3.689e+03 -3.688e+03                   
#> Path [11] :Best Iter: [53] ELBO (-3688.480765) evaluations: (2935) 
#> Path [12] :Initial log joint density = -481836.998957 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.608e-03   2.474e-01    1.000e+00  1.000e+00      3278 -3.681e+03 -3.682e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3680.838431) evaluations: (3278) 
#> Path [13] :Initial log joint density = -481504.379147 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.709e-03   2.044e-01    1.000e+00  1.000e+00      3184 -3.692e+03 -3.690e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3689.776047) evaluations: (3184) 
#> Path [14] :Initial log joint density = -481931.247017 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.208e-03   2.060e-01    1.000e+00  1.000e+00      3149 -3.687e+03 -3.688e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3687.463212) evaluations: (3149) 
#> Path [15] :Initial log joint density = -481563.888888 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.334e-03   2.832e-01    6.920e-01  6.920e-01      3272 -3.690e+03 -3.694e+03                   
#> Path [15] :Best Iter: [54] ELBO (-3690.293266) evaluations: (3272) 
#> Path [16] :Initial log joint density = -481405.698776 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.617e-03   2.126e-01    9.486e-01  9.486e-01      3071 -3.690e+03 -3.692e+03                   
#> Path [16] :Best Iter: [53] ELBO (-3690.053215) evaluations: (3071) 
#> Path [17] :Initial log joint density = -481446.897890 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.738e-03   1.891e-01    1.000e+00  1.000e+00      3305 -3.681e+03 -3.683e+03                   
#> Path [17] :Best Iter: [58] ELBO (-3681.486896) evaluations: (3305) 
#> Path [18] :Initial log joint density = -481701.018483 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.307e-02   1.644e-01    1.000e+00  1.000e+00      3425 -3.684e+03 -3.682e+03                   
#> Path [18] :Best Iter: [58] ELBO (-3682.226687) evaluations: (3425) 
#> Path [19] :Initial log joint density = -482534.399378 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.612e-03   3.135e-01    4.934e-01  1.000e+00      3470 -3.686e+03 -3.692e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3685.506467) evaluations: (3470) 
#> Path [20] :Initial log joint density = -481582.959641 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.721e-03   1.529e-01    1.000e+00  1.000e+00      3391 -3.683e+03 -3.693e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3682.765216) evaluations: (3391) 
#> Path [21] :Initial log joint density = -482621.289916 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.941e-03   1.761e-01    1.000e+00  1.000e+00      3295 -3.689e+03 -3.684e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3684.278695) evaluations: (3295) 
#> Path [22] :Initial log joint density = -482816.210165 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.683e-02   3.654e-01    1.000e+00  1.000e+00      3682 -3.683e+03 -3.690e+03                   
#> Path [22] :Best Iter: [61] ELBO (-3683.245935) evaluations: (3682) 
#> Path [23] :Initial log joint density = -481659.525817 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      6.940e-03   2.464e-01    1.000e+00  1.000e+00      2632 -3.693e+03 -3.702e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3692.829744) evaluations: (2632) 
#> Path [24] :Initial log joint density = -482089.066629 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.511e-02   2.867e-01    9.497e-01  9.497e-01      3744 -3.683e+03 -3.691e+03                   
#> Path [24] :Best Iter: [61] ELBO (-3683.174209) evaluations: (3744) 
#> Path [25] :Initial log joint density = -481960.866784 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.902e-03   2.647e-01    8.720e-01  8.720e-01      3239 -3.685e+03 -3.696e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3685.127571) evaluations: (3239) 
#> Path [26] :Initial log joint density = -482111.825831 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.747e-03   1.831e-01    1.000e+00  1.000e+00      3560 -3.683e+03 -3.684e+03                   
#> Path [26] :Best Iter: [57] ELBO (-3682.787648) evaluations: (3560) 
#> Path [27] :Initial log joint density = -481377.543819 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.339e-03   1.614e-01    1.000e+00  1.000e+00      3405 -3.686e+03 -3.694e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3686.457509) evaluations: (3405) 
#> Path [28] :Initial log joint density = -481865.836229 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.829e-03   1.805e-01    1.000e+00  1.000e+00      3477 -3.685e+03 -3.685e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3684.524910) evaluations: (3477) 
#> Path [29] :Initial log joint density = -482156.219188 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.101e-02   3.062e-01    1.000e+00  1.000e+00      3400 -3.684e+03 -3.689e+03                   
#> Path [29] :Best Iter: [58] ELBO (-3684.099798) evaluations: (3400) 
#> Path [30] :Initial log joint density = -481798.954504 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.106e-02   2.710e-01    1.000e+00  1.000e+00      3372 -3.684e+03 -3.687e+03                   
#> Path [30] :Best Iter: [56] ELBO (-3683.790730) evaluations: (3372) 
#> Path [31] :Initial log joint density = -481642.726839 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.309e-02   3.009e-01    1.000e+00  1.000e+00      3272 -3.683e+03 -3.687e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3683.337287) evaluations: (3272) 
#> Path [32] :Initial log joint density = -481702.857655 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.577e-03   2.039e-01    1.000e+00  1.000e+00      2912 -3.690e+03 -3.704e+03                   
#> Path [32] :Best Iter: [50] ELBO (-3689.786522) evaluations: (2912) 
#> Path [33] :Initial log joint density = -481489.636748 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.473e-03   1.671e-01    8.246e-01  8.246e-01      3109 -3.685e+03 -3.699e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3685.461101) evaluations: (3109) 
#> Path [34] :Initial log joint density = -481842.031273 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.432e-02   2.707e-01    1.000e+00  1.000e+00      3821 -3.683e+03 -3.682e+03                   
#> Path [34] :Best Iter: [63] ELBO (-3681.881983) evaluations: (3821) 
#> Path [35] :Initial log joint density = -481815.816776 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.126e-02   2.274e-01    1.000e+00  1.000e+00      3533 -3.683e+03 -3.686e+03                   
#> Path [35] :Best Iter: [59] ELBO (-3683.388269) evaluations: (3533) 
#> Path [36] :Initial log joint density = -481864.709645 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.335e-03   1.752e-01    1.000e+00  1.000e+00      3345 -3.682e+03 -3.689e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3681.932133) evaluations: (3345) 
#> Path [37] :Initial log joint density = -482009.672163 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.231e-02   3.599e-01    1.000e+00  1.000e+00      3163 -3.686e+03 -3.690e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3685.791725) evaluations: (3163) 
#> Path [38] :Initial log joint density = -481799.104237 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      4.734e-03   1.496e-01    1.000e+00  1.000e+00      3804 -3.683e+03 -3.690e+03                   
#> Path [38] :Best Iter: [60] ELBO (-3682.793698) evaluations: (3804) 
#> Path [39] :Initial log joint density = -481895.772797 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      3.418e-03   2.876e-01    5.724e-01  5.724e-01      3710 -3.680e+03 -3.688e+03                   
#> Path [39] :Best Iter: [61] ELBO (-3679.934671) evaluations: (3710) 
#> Path [40] :Initial log joint density = -482756.553597 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.784e-02   2.037e-01    1.000e+00  1.000e+00      3666 -3.684e+03 -3.684e+03                   
#> Path [40] :Best Iter: [61] ELBO (-3683.869243) evaluations: (3666) 
#> Path [41] :Initial log joint density = -481759.341554 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.113e-02   2.398e-01    1.000e+00  1.000e+00      3700 -3.680e+03 -3.686e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3680.119223) evaluations: (3700) 
#> Path [42] :Initial log joint density = -481630.131972 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.597e-02   3.691e-01    1.000e+00  1.000e+00      3106 -3.684e+03 -3.691e+03                   
#> Path [42] :Best Iter: [55] ELBO (-3684.490967) evaluations: (3106) 
#> Path [43] :Initial log joint density = -481506.367504 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.408e-02   4.004e-01    9.329e-01  9.329e-01      3224 -3.690e+03 -3.694e+03                   
#> Path [43] :Best Iter: [40] ELBO (-3690.351929) evaluations: (3224) 
#> Path [44] :Initial log joint density = -481382.447436 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.531e-02   3.067e-01    1.000e+00  1.000e+00      2829 -3.690e+03 -3.691e+03                   
#> Path [44] :Best Iter: [43] ELBO (-3690.458259) evaluations: (2829) 
#> Path [45] :Initial log joint density = -481585.319720 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.274e-02   2.590e-01    1.000e+00  1.000e+00      3503 -3.684e+03 -3.689e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3683.511716) evaluations: (3503) 
#> Path [46] :Initial log joint density = -482350.731077 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.125e-02   1.949e-01    1.000e+00  1.000e+00      3330 -3.684e+03 -3.686e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3683.740902) evaluations: (3330) 
#> Path [47] :Initial log joint density = -481680.247611 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.612e-03   2.227e-01    1.000e+00  1.000e+00      3227 -3.692e+03 -3.689e+03                   
#> Path [47] :Best Iter: [57] ELBO (-3689.306750) evaluations: (3227) 
#> Path [48] :Initial log joint density = -483356.942326 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.857e-03   2.054e-01    9.075e-01  9.075e-01      3155 -3.691e+03 -3.693e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3691.078297) evaluations: (3155) 
#> Path [49] :Initial log joint density = -481703.821138 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.614e-03   2.095e-01    9.684e-01  9.684e-01      3100 -3.689e+03 -3.693e+03                   
#> Path [49] :Best Iter: [53] ELBO (-3689.449000) evaluations: (3100) 
#> Path [50] :Initial log joint density = -481575.347293 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.791e-03   2.585e-01    1.000e+00  1.000e+00      2750 -3.690e+03 -3.697e+03                   
#> Path [50] :Best Iter: [46] ELBO (-3690.495387) evaluations: (2750) 
#> Finished in  13.7 seconds.
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
