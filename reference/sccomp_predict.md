# sccomp_predict

This function replicates counts from a real-world dataset.

## Usage

``` r
sccomp_predict(
  fit,
  formula_composition = NULL,
  new_data = NULL,
  number_of_draws = 500,
  mcmc_seed = sample_seed(),
  summary_instead_of_draws = TRUE,
  robust = FALSE
)
```

## Arguments

- fit:

  The result of sccomp_estimate.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment. This formula can be a sub-formula
  of your estimated model; in this case all other factor will be
  factored out.

- new_data:

  A sample-wise data frame including the column that represent the
  factors in your formula. If you want to predict proportions for 10
  samples, there should be 10 rows. T

- number_of_draws:

  An integer. How may copies of the data you want to draw from the model
  joint posterior distribution.

- mcmc_seed:

  An integer. Used for Markov-chain Monte Carlo reproducibility. By
  default a random number is sampled from 1 to 999999. This itself can
  be controlled by set.seed()

- summary_instead_of_draws:

  Return the summary values (i.e. mean and quantiles) of the predicted
  proportions, or return single draws. Single draws can be helful to
  better analyse the uncertainty of the prediction.

- robust:

  A logical. If TRUE, use robust statistics (median and median absolute
  deviation) instead of classical statistics (mean and standard
  deviation) for the summary calculations.

## Value

A tibble (`tbl`) with the following columns:

- **cell_group** - A character column representing the cell group being
  tested.

- **sample** - A factor column representing the sample name for which
  the predictions are made.

- **proportion_mean** - A numeric column representing the predicted mean
  (or median when robust=TRUE) proportions from the model.

- **proportion_lower** - A numeric column representing the lower bound
  (2.5%) of the 95% credible interval for the predicted proportions.

- **proportion_upper** - A numeric column representing the upper bound
  (97.5%) of the 95% credible interval for the predicted proportions.

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
  if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
    data("counts_obj")

    sccomp_estimate(
      counts_obj,
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_predict()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481572.345045 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.965e-03   2.179e-01    1.000e+00  1.000e+00      2911 -3.706e+03 -3.713e+03                   
#> Path [1] :Best Iter: [42] ELBO (-3705.770457) evaluations: (2911) 
#> Path [2] :Initial log joint density = -482163.812312 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.227e-03   1.792e-01    1.000e+00  1.000e+00      3327 -3.703e+03 -3.706e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3702.900783) evaluations: (3327) 
#> Path [3] :Initial log joint density = -481584.869724 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.359e-03   2.094e-01    8.768e-01  8.768e-01      3393 -3.703e+03 -3.712e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3703.130462) evaluations: (3393) 
#> Path [4] :Initial log joint density = -481744.069203 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.695e-03   2.538e-01    1.000e+00  1.000e+00      3080 -3.706e+03 -3.699e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3699.163097) evaluations: (3080) 
#> Path [5] :Initial log joint density = -481579.928171 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.247e-03   2.239e-01    1.000e+00  1.000e+00      3079 -3.709e+03 -3.706e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3705.537905) evaluations: (3079) 
#> Path [6] :Initial log joint density = -481820.864967 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.757e-03   2.351e-01    7.193e-01  7.193e-01      3443 -3.698e+03 -3.716e+03                   
#> Path [6] :Best Iter: [58] ELBO (-3698.258585) evaluations: (3443) 
#> Path [7] :Initial log joint density = -486660.743002 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.375e-02   4.224e-01    1.000e+00  1.000e+00      3149 -3.701e+03 -3.709e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3700.674039) evaluations: (3149) 
#> Path [8] :Initial log joint density = -481613.729385 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.278e-03   2.976e-01    1.000e+00  1.000e+00      3207 -3.703e+03 -3.710e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3703.103574) evaluations: (3207) 
#> Path [9] :Initial log joint density = -481379.325218 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.529e-03   2.251e-01    1.000e+00  1.000e+00      2971 -3.706e+03 -3.700e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3699.689382) evaluations: (2971) 
#> Path [10] :Initial log joint density = -483265.715914 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.937e-03   2.177e-01    4.708e-01  1.000e+00      3222 -3.699e+03 -3.711e+03                   
#> Path [10] :Best Iter: [55] ELBO (-3699.239659) evaluations: (3222) 
#> Path [11] :Initial log joint density = -481703.986761 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.309e-03   2.440e-01    1.000e+00  1.000e+00      3397 -3.701e+03 -3.704e+03                   
#> Path [11] :Best Iter: [57] ELBO (-3700.790994) evaluations: (3397) 
#> Path [12] :Initial log joint density = -481419.473390 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.173e-03   2.747e-01    8.305e-01  8.305e-01      2725 -3.708e+03 -3.721e+03                   
#> Path [12] :Best Iter: [46] ELBO (-3708.139005) evaluations: (2725) 
#> Path [13] :Initial log joint density = -481446.332676 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.848e-03   2.585e-01    1.000e+00  1.000e+00      2965 -3.707e+03 -3.711e+03                   
#> Path [13] :Best Iter: [49] ELBO (-3707.436596) evaluations: (2965) 
#> Path [14] :Initial log joint density = -483423.085044 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.786e-03   1.995e-01    1.000e+00  1.000e+00      3332 -3.699e+03 -3.702e+03                   
#> Path [14] :Best Iter: [56] ELBO (-3699.003075) evaluations: (3332) 
#> Path [15] :Initial log joint density = -481862.504843 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.459e-03   2.746e-01    6.815e-01  6.815e-01      3439 -3.701e+03 -3.714e+03                   
#> Path [15] :Best Iter: [57] ELBO (-3700.517591) evaluations: (3439) 
#> Path [16] :Initial log joint density = -481708.765653 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.386e-03   2.066e-01    1.000e+00  1.000e+00      2971 -3.710e+03 -3.710e+03                   
#> Path [16] :Best Iter: [51] ELBO (-3709.507804) evaluations: (2971) 
#> Path [17] :Initial log joint density = -481880.431975 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.492e-02   2.793e-01    1.000e+00  1.000e+00      3290 -3.700e+03 -3.702e+03                   
#> Path [17] :Best Iter: [55] ELBO (-3699.962185) evaluations: (3290) 
#> Path [18] :Initial log joint density = -481707.095036 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.819e-03   1.722e-01    9.394e-01  9.394e-01      3021 -3.711e+03 -3.725e+03                   
#> Path [18] :Best Iter: [43] ELBO (-3710.694627) evaluations: (3021) 
#> Path [19] :Initial log joint density = -481823.237548 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      9.307e-03   1.858e-01    1.000e+00  1.000e+00      3747 -3.699e+03 -3.699e+03                   
#> Path [19] :Best Iter: [60] ELBO (-3698.714421) evaluations: (3747) 
#> Path [20] :Initial log joint density = -481747.961449 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.225e-02   1.763e-01    1.000e+00  1.000e+00      3744 -3.698e+03 -3.698e+03                   
#> Path [20] :Best Iter: [59] ELBO (-3697.837811) evaluations: (3744) 
#> Path [21] :Initial log joint density = -481727.197861 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.394e-03   2.078e-01    7.674e-01  7.674e-01      2863 -3.709e+03 -3.724e+03                   
#> Path [21] :Best Iter: [52] ELBO (-3708.846214) evaluations: (2863) 
#> Path [22] :Initial log joint density = -481613.060041 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.112e-02   2.624e-01    1.000e+00  1.000e+00      3375 -3.701e+03 -3.704e+03                   
#> Path [22] :Best Iter: [57] ELBO (-3700.872766) evaluations: (3375) 
#> Path [23] :Initial log joint density = -482008.349734 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.251e-02   3.497e-01    1.000e+00  1.000e+00      3844 -3.699e+03 -3.711e+03                   
#> Path [23] :Best Iter: [63] ELBO (-3698.613880) evaluations: (3844) 
#> Path [24] :Initial log joint density = -482035.623761 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.449e-03   2.160e-01    7.006e-01  7.006e-01      3249 -3.702e+03 -3.720e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3702.023401) evaluations: (3249) 
#> Path [25] :Initial log joint density = -481681.150638 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.300e-02   4.696e-01    1.000e+00  1.000e+00      2866 -3.706e+03 -3.725e+03                   
#> Path [25] :Best Iter: [51] ELBO (-3706.233380) evaluations: (2866) 
#> Path [26] :Initial log joint density = -482495.331163 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.409e-03   2.089e-01    1.000e+00  1.000e+00      3662 -3.701e+03 -3.708e+03                   
#> Path [26] :Best Iter: [57] ELBO (-3700.871267) evaluations: (3662) 
#> Path [27] :Initial log joint density = -481997.545751 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.556e-03   2.096e-01    8.969e-01  8.969e-01      3478 -3.696e+03 -3.710e+03                   
#> Path [27] :Best Iter: [59] ELBO (-3695.983125) evaluations: (3478) 
#> Path [28] :Initial log joint density = -481511.496566 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.347e-02   1.836e-01    1.000e+00  1.000e+00      3203 -3.700e+03 -3.703e+03                   
#> Path [28] :Best Iter: [56] ELBO (-3699.592340) evaluations: (3203) 
#> Path [29] :Initial log joint density = -484682.698319 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.111e-03   2.472e-01    6.568e-01  6.568e-01      3678 -3.703e+03 -3.714e+03                   
#> Path [29] :Best Iter: [61] ELBO (-3703.344670) evaluations: (3678) 
#> Path [30] :Initial log joint density = -481822.271280 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.941e-03   1.551e-01    7.803e-01  7.803e-01      3076 -3.709e+03 -3.725e+03                   
#> Path [30] :Best Iter: [45] ELBO (-3708.766086) evaluations: (3076) 
#> Path [31] :Initial log joint density = -481615.799590 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.635e-02   2.864e-01    1.000e+00  1.000e+00      3328 -3.698e+03 -3.702e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3697.516012) evaluations: (3328) 
#> Path [32] :Initial log joint density = -481551.276677 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.575e-03   1.690e-01    1.000e+00  1.000e+00      3082 -3.707e+03 -3.709e+03                   
#> Path [32] :Best Iter: [52] ELBO (-3706.872637) evaluations: (3082) 
#> Path [33] :Initial log joint density = -481745.865107 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.007e-02   2.185e-01    1.000e+00  1.000e+00      3273 -3.706e+03 -3.701e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3700.698378) evaluations: (3273) 
#> Path [34] :Initial log joint density = -482402.517522 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.676e-03   1.881e-01    1.000e+00  1.000e+00      3417 -3.699e+03 -3.708e+03                   
#> Path [34] :Best Iter: [55] ELBO (-3698.835785) evaluations: (3417) 
#> Path [35] :Initial log joint density = -481406.249447 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.173e-02   1.823e-01    1.000e+00  1.000e+00      3295 -3.701e+03 -3.699e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3699.266704) evaluations: (3295) 
#> Path [36] :Initial log joint density = -481680.013285 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.120e-02   1.973e-01    1.000e+00  1.000e+00      3192 -3.705e+03 -3.699e+03                   
#> Path [36] :Best Iter: [57] ELBO (-3699.482462) evaluations: (3192) 
#> Path [37] :Initial log joint density = -481922.925608 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.553e-03   2.397e-01    1.000e+00  1.000e+00      3047 -3.709e+03 -3.701e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3701.486226) evaluations: (3047) 
#> Path [38] :Initial log joint density = -481928.748559 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.492e-03   2.525e-01    8.607e-01  8.607e-01      3478 -3.702e+03 -3.711e+03                   
#> Path [38] :Best Iter: [59] ELBO (-3702.215112) evaluations: (3478) 
#> Path [39] :Initial log joint density = -484068.486920 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.144e-02   2.160e-01    1.000e+00  1.000e+00      3506 -3.703e+03 -3.711e+03                   
#> Path [39] :Best Iter: [58] ELBO (-3703.287172) evaluations: (3506) 
#> Path [40] :Initial log joint density = -481584.754644 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.054e-03   1.823e-01    1.000e+00  1.000e+00      3347 -3.707e+03 -3.705e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3704.836455) evaluations: (3347) 
#> Path [41] :Initial log joint density = -487494.357628 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.450e-03   2.454e-01    1.000e+00  1.000e+00      3451 -3.706e+03 -3.710e+03                   
#> Path [41] :Best Iter: [56] ELBO (-3705.647378) evaluations: (3451) 
#> Path [42] :Initial log joint density = -481886.876588 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.005e-03   2.476e-01    1.000e+00  1.000e+00      3479 -3.700e+03 -3.700e+03                   
#> Path [42] :Best Iter: [60] ELBO (-3699.532545) evaluations: (3479) 
#> Path [43] :Initial log joint density = -481884.123088 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.127e-02   2.778e-01    9.313e-01  9.313e-01      3035 -3.710e+03 -3.722e+03                   
#> Path [43] :Best Iter: [45] ELBO (-3710.131764) evaluations: (3035) 
#> Path [44] :Initial log joint density = -481532.821379 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.299e-02   3.690e-01    1.000e+00  1.000e+00      3508 -3.697e+03 -3.701e+03                   
#> Path [44] :Best Iter: [59] ELBO (-3697.250301) evaluations: (3508) 
#> Path [45] :Initial log joint density = -481533.493056 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.082e-03   1.654e-01    1.000e+00  1.000e+00      3005 -3.711e+03 -3.711e+03                   
#> Path [45] :Best Iter: [51] ELBO (-3711.287981) evaluations: (3005) 
#> Path [46] :Initial log joint density = -481592.864295 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      7.482e-03   2.559e-01    1.000e+00  1.000e+00      2581 -3.709e+03 -3.720e+03                   
#> Path [46] :Best Iter: [48] ELBO (-3709.130132) evaluations: (2581) 
#> Path [47] :Initial log joint density = -481333.357196 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.182e-02   2.434e-01    1.000e+00  1.000e+00      3102 -3.704e+03 -3.706e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3703.622163) evaluations: (3102) 
#> Path [48] :Initial log joint density = -481622.406052 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.310e-03   2.465e-01    6.974e-01  6.974e-01      3579 -3.700e+03 -3.711e+03                   
#> Path [48] :Best Iter: [59] ELBO (-3699.964106) evaluations: (3579) 
#> Path [49] :Initial log joint density = -481593.618794 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.957e-03   2.009e-01    1.000e+00  1.000e+00      2994 -3.709e+03 -3.713e+03                   
#> Path [49] :Best Iter: [45] ELBO (-3708.871498) evaluations: (2994) 
#> Path [50] :Initial log joint density = -481928.421695 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.417e-02   4.070e-01    1.000e+00  1.000e+00      3026 -3.709e+03 -3.706e+03                   
#> Path [50] :Best Iter: [55] ELBO (-3706.044293) evaluations: (3026) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 6
#>    sample type   cell_group proportion_mean proportion_lower proportion_upper
#>    <fct>  <fct>  <chr>                <dbl>            <dbl>            <dbl>
#>  1 10x_6K benign B1                 0.0587           0.0461           0.0720 
#>  2 10x_6K benign B2                 0.0376           0.0291           0.0469 
#>  3 10x_6K benign B3                 0.0127           0.00962          0.0160 
#>  4 10x_6K benign BM                 0.00674          0.00512          0.00859
#>  5 10x_6K benign CD4 1              0.0254           0.0215           0.0297 
#>  6 10x_6K benign CD4 2              0.0509           0.0426           0.0608 
#>  7 10x_6K benign CD4 3              0.0825           0.0631           0.103  
#>  8 10x_6K benign CD4 4              0.00171          0.00111          0.00242
#>  9 10x_6K benign CD4 5              0.0304           0.0237           0.0373 
#> 10 10x_6K benign CD8 1              0.111            0.0949           0.129  
#> # ℹ 710 more rows
# }

```
