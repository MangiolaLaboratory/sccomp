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
#> Path [1] :Initial log joint density = -481412.729645 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      4.143e-03   3.119e-01    6.320e-01  6.320e-01      2937 -3.689e+03 -3.710e+03                   
#> Path [1] :Best Iter: [44] ELBO (-3689.211347) evaluations: (2937) 
#> Path [2] :Initial log joint density = -483370.292828 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.458e-03   2.252e-01    7.945e-01  7.945e-01      3161 -3.687e+03 -3.693e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3686.582483) evaluations: (3161) 
#> Path [3] :Initial log joint density = -481353.835049 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.161e-03   2.910e-01    1.000e+00  1.000e+00      3080 -3.691e+03 -3.690e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3690.162428) evaluations: (3080) 
#> Path [4] :Initial log joint density = -481590.319138 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.277e-02   2.228e-01    1.000e+00  1.000e+00      3290 -3.689e+03 -3.692e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3689.182809) evaluations: (3290) 
#> Path [5] :Initial log joint density = -481316.538921 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.943e-03   2.289e-01    6.647e-01  6.647e-01      3169 -3.691e+03 -3.698e+03                   
#> Path [5] :Best Iter: [48] ELBO (-3690.753217) evaluations: (3169) 
#> Path [6] :Initial log joint density = -481823.126763 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.817e-03   3.472e-01    6.939e-01  6.939e-01      3069 -3.691e+03 -3.696e+03                   
#> Path [6] :Best Iter: [47] ELBO (-3691.195802) evaluations: (3069) 
#> Path [7] :Initial log joint density = -481551.143973 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.583e-03   2.405e-01    7.498e-01  7.498e-01      3088 -3.689e+03 -3.705e+03                   
#> Path [7] :Best Iter: [40] ELBO (-3688.947957) evaluations: (3088) 
#> Path [8] :Initial log joint density = -481971.431973 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.230e-02   2.275e-01    1.000e+00  1.000e+00      3193 -3.681e+03 -3.684e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3681.092070) evaluations: (3193) 
#> Path [9] :Initial log joint density = -481583.200574 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.248e-02   2.076e-01    1.000e+00  1.000e+00      3105 -3.692e+03 -3.686e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3686.099476) evaluations: (3105) 
#> Path [10] :Initial log joint density = -481709.977564 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.252e-03   2.828e-01    1.000e+00  1.000e+00      2992 -3.691e+03 -3.693e+03                   
#> Path [10] :Best Iter: [52] ELBO (-3691.056310) evaluations: (2992) 
#> Path [11] :Initial log joint density = -481642.178775 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.054e-02   1.703e-01    1.000e+00  1.000e+00      3539 -3.685e+03 -3.683e+03                   
#> Path [11] :Best Iter: [59] ELBO (-3682.618949) evaluations: (3539) 
#> Path [12] :Initial log joint density = -481832.041547 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      7.384e-03   1.735e-01    1.000e+00  1.000e+00      3710 -3.683e+03 -3.685e+03                   
#> Path [12] :Best Iter: [60] ELBO (-3683.229695) evaluations: (3710) 
#> Path [13] :Initial log joint density = -481515.980977 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.274e-02   3.377e-01    1.000e+00  1.000e+00      3374 -3.686e+03 -3.688e+03                   
#> Path [13] :Best Iter: [57] ELBO (-3686.483783) evaluations: (3374) 
#> Path [14] :Initial log joint density = -481691.147554 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.744e-03   2.330e-01    1.000e+00  1.000e+00      2987 -3.690e+03 -3.704e+03                   
#> Path [14] :Best Iter: [45] ELBO (-3690.383122) evaluations: (2987) 
#> Path [15] :Initial log joint density = -481819.312273 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.774e-03   1.655e-01    1.000e+00  1.000e+00      3075 -3.692e+03 -3.692e+03                   
#> Path [15] :Best Iter: [39] ELBO (-3692.248182) evaluations: (3075) 
#> Path [16] :Initial log joint density = -481539.607416 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.002e-02   2.158e-01    1.000e+00  1.000e+00      3123 -3.688e+03 -3.683e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3683.455291) evaluations: (3123) 
#> Path [17] :Initial log joint density = -481198.927795 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.283e-02   1.885e-01    1.000e+00  1.000e+00      3651 -3.680e+03 -3.684e+03                   
#> Path [17] :Best Iter: [58] ELBO (-3679.738115) evaluations: (3651) 
#> Path [18] :Initial log joint density = -482012.611215 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.488e-02   2.968e-01    1.000e+00  1.000e+00      3827 -3.682e+03 -3.686e+03                   
#> Path [18] :Best Iter: [62] ELBO (-3681.812867) evaluations: (3827) 
#> Path [19] :Initial log joint density = -481942.406685 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.316e-02   2.627e-01    1.000e+00  1.000e+00      3220 -3.684e+03 -3.687e+03                   
#> Path [19] :Best Iter: [57] ELBO (-3684.485100) evaluations: (3220) 
#> Path [20] :Initial log joint density = -481463.339337 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.038e-02   2.487e-01    1.000e+00  1.000e+00      3157 -3.685e+03 -3.687e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3684.702717) evaluations: (3157) 
#> Path [21] :Initial log joint density = -481575.559666 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.468e-03   1.928e-01    1.000e+00  1.000e+00      3049 -3.692e+03 -3.700e+03                   
#> Path [21] :Best Iter: [52] ELBO (-3691.994049) evaluations: (3049) 
#> Path [22] :Initial log joint density = -481785.376998 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.039e-03   2.221e-01    1.000e+00  1.000e+00      3073 -3.690e+03 -3.693e+03                   
#> Path [22] :Best Iter: [41] ELBO (-3689.946908) evaluations: (3073) 
#> Path [23] :Initial log joint density = -481367.476538 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.994e-03   2.845e-01    1.000e+00  1.000e+00      2987 -3.692e+03 -3.696e+03                   
#> Path [23] :Best Iter: [53] ELBO (-3692.202592) evaluations: (2987) 
#> Path [24] :Initial log joint density = -482498.560765 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.123e-02   3.641e-01    1.000e+00  1.000e+00      3349 -3.684e+03 -3.685e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3684.431832) evaluations: (3349) 
#> Path [25] :Initial log joint density = -481506.046978 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      5.608e-03   2.074e-01    1.000e+00  1.000e+00      2595 -3.689e+03 -3.699e+03                   
#> Path [25] :Best Iter: [41] ELBO (-3688.918355) evaluations: (2595) 
#> Path [26] :Initial log joint density = -481485.072787 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.186e-02   4.015e-01    1.000e+00  1.000e+00      3021 -3.691e+03 -3.692e+03                   
#> Path [26] :Best Iter: [54] ELBO (-3691.324688) evaluations: (3021) 
#> Path [27] :Initial log joint density = -481916.399681 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      5.665e-03   1.472e-01    7.848e-01  7.848e-01      3818 -3.684e+03 -3.696e+03                   
#> Path [27] :Best Iter: [62] ELBO (-3684.493906) evaluations: (3818) 
#> Path [28] :Initial log joint density = -481716.915588 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.306e-02   2.948e-01    1.000e+00  1.000e+00      3407 -3.682e+03 -3.689e+03                   
#> Path [28] :Best Iter: [56] ELBO (-3681.679990) evaluations: (3407) 
#> Path [29] :Initial log joint density = -483015.170450 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.995e-03   2.024e-01    1.000e+00  1.000e+00      3255 -3.689e+03 -3.686e+03                   
#> Path [29] :Best Iter: [57] ELBO (-3686.247044) evaluations: (3255) 
#> Path [30] :Initial log joint density = -481739.638557 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.044e-03   1.777e-01    6.730e-01  6.730e-01      3488 -3.683e+03 -3.696e+03                   
#> Path [30] :Best Iter: [58] ELBO (-3682.881835) evaluations: (3488) 
#> Path [31] :Initial log joint density = -481959.237031 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.910e-03   2.004e-01    1.000e+00  1.000e+00      3443 -3.685e+03 -3.688e+03                   
#> Path [31] :Best Iter: [59] ELBO (-3685.232639) evaluations: (3443) 
#> Path [32] :Initial log joint density = -483039.929643 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.715e-03   2.997e-01    1.000e+00  1.000e+00      3025 -3.691e+03 -3.689e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3688.521788) evaluations: (3025) 
#> Path [33] :Initial log joint density = -486426.947044 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.057e-02   2.355e-01    8.313e-01  8.313e-01      3382 -3.686e+03 -3.696e+03                   
#> Path [33] :Best Iter: [57] ELBO (-3686.266769) evaluations: (3382) 
#> Path [34] :Initial log joint density = -481517.510697 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.689e-03   2.050e-01    1.000e+00  1.000e+00      2890 -3.691e+03 -3.697e+03                   
#> Path [34] :Best Iter: [52] ELBO (-3691.482087) evaluations: (2890) 
#> Path [35] :Initial log joint density = -482417.030228 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.894e-03   3.277e-01    1.000e+00  1.000e+00      2995 -3.691e+03 -3.703e+03                   
#> Path [35] :Best Iter: [51] ELBO (-3690.562887) evaluations: (2995) 
#> Path [36] :Initial log joint density = -483561.969625 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      4.266e-03   2.886e-01    4.933e-01  4.933e-01      3794 -3.683e+03 -3.689e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3682.795820) evaluations: (3794) 
#> Path [37] :Initial log joint density = -481618.814035 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.788e-03   1.795e-01    1.000e+00  1.000e+00      3597 -3.688e+03 -3.686e+03                   
#> Path [37] :Best Iter: [61] ELBO (-3686.281309) evaluations: (3597) 
#> Path [38] :Initial log joint density = -481436.379278 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.280e-03   1.741e-01    1.000e+00  1.000e+00      3158 -3.691e+03 -3.695e+03                   
#> Path [38] :Best Iter: [54] ELBO (-3691.443780) evaluations: (3158) 
#> Path [39] :Initial log joint density = -481644.568109 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.281e-03   2.035e-01    1.000e+00  1.000e+00      3163 -3.689e+03 -3.689e+03                   
#> Path [39] :Best Iter: [48] ELBO (-3688.619991) evaluations: (3163) 
#> Path [40] :Initial log joint density = -481655.104461 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.774e-03   2.530e-01    9.592e-01  9.592e-01      3271 -3.684e+03 -3.697e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3684.231731) evaluations: (3271) 
#> Path [41] :Initial log joint density = -483097.811793 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.956e-03   3.010e-01    2.975e-01  1.000e+00      3242 -3.683e+03 -3.695e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3683.439851) evaluations: (3242) 
#> Path [42] :Initial log joint density = -481634.336921 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.107e-03   1.788e-01    8.826e-01  8.826e-01      3555 -3.680e+03 -3.689e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3680.498556) evaluations: (3555) 
#> Path [43] :Initial log joint density = -481370.578117 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.219e-02   2.846e-01    1.000e+00  1.000e+00      2686 -3.689e+03 -3.697e+03                   
#> Path [43] :Best Iter: [47] ELBO (-3689.304222) evaluations: (2686) 
#> Path [44] :Initial log joint density = -481639.579924 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.117e-03   2.214e-01    7.686e-01  7.686e-01      2863 -3.689e+03 -3.707e+03                   
#> Path [44] :Best Iter: [47] ELBO (-3689.266360) evaluations: (2863) 
#> Path [45] :Initial log joint density = -481743.967156 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.070e-03   2.870e-01    1.000e+00  1.000e+00      3077 -3.692e+03 -3.693e+03                   
#> Path [45] :Best Iter: [38] ELBO (-3691.790869) evaluations: (3077) 
#> Path [46] :Initial log joint density = -481793.072904 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.712e-03   3.068e-01    1.000e+00  1.000e+00      3134 -3.689e+03 -3.686e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3685.866711) evaluations: (3134) 
#> Path [47] :Initial log joint density = -481532.363900 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.041e-03   2.617e-01    7.192e-01  7.192e-01      3118 -3.688e+03 -3.698e+03                   
#> Path [47] :Best Iter: [49] ELBO (-3688.208339) evaluations: (3118) 
#> Path [48] :Initial log joint density = -482322.966973 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.661e-03   2.778e-01    1.000e+00  1.000e+00      3590 -3.686e+03 -3.694e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3685.590771) evaluations: (3590) 
#> Path [49] :Initial log joint density = -483263.531424 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      2.278e-02   3.380e-01    1.000e+00  1.000e+00      4014 -3.680e+03 -3.684e+03                   
#> Path [49] :Best Iter: [61] ELBO (-3679.794383) evaluations: (4014) 
#> Path [50] :Initial log joint density = -481974.852138 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.026e-03   2.881e-01    9.925e-01  9.925e-01      3408 -3.684e+03 -3.696e+03                   
#> Path [50] :Best Iter: [57] ELBO (-3684.163153) evaluations: (3408) 
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
#> Chain 1  Elapsed Time: 0.566 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: regressing out unwanted factors
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.571 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
# }
```
