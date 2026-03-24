# sccomp_remove_outliers main

The `sccomp_remove_outliers` function takes as input a table of cell
counts with columns for cell-group identifier, sample identifier,
integer count, and factors (continuous or discrete). The user can define
a linear model using an input R formula, where the first factor is the
factor of interest. Alternatively, `sccomp` accepts single-cell data
containers (e.g., Seurat, SingleCellExperiment, cell metadata, or
group-size) and derives the count data from cell metadata.

## Usage

``` r
sccomp_remove_outliers(
  .estimate,
  percent_false_positive = 5,
  cores = detectCores(),
  inference_method = attr(.estimate, "inference_method"),
  output_directory = "sccomp_draws_files",
  verbose = TRUE,
  mcmc_seed = sample_seed(),
  max_sampling_iterations = 20000,
  enable_loo = FALSE,
  sig_figs = 9,
  cache_stan_model = sccomp_stan_models_cache_dir,
  cleanup_draw_files = TRUE,
  approximate_posterior_inference = NULL,
  variational_inference = NULL,
  ...
)
```

## Arguments

- .estimate:

  A tibble including a cell_group name column, sample name column, read
  counts column (optional depending on the input class), and factor
  columns.

- percent_false_positive:

  A real number between 0 and 100 (not inclusive), used to identify
  outliers with a specific false positive rate.

- cores:

  Integer, the number of cores to be used for parallel calculations.

- inference_method:

  Character string specifying the inference method to use ('pathfinder',
  'hmc', or 'variational').

- output_directory:

  A character string specifying the output directory for Stan draws.

- verbose:

  Logical, whether to print progression details.

- mcmc_seed:

  Integer, used for Markov-chain Monte Carlo reproducibility. By
  default, a random number is sampled from 1 to 999999.

- max_sampling_iterations:

  Integer, limits the maximum number of iterations in case a large
  dataset is used, to limit computation time.

- enable_loo:

  Logical, whether to enable model comparison using the R package LOO.
  This is useful for comparing fits between models, similar to ANOVA.

- sig_figs:

  Number of significant figures to use for Stan model output. Default is
  9.

- cache_stan_model:

  A character string specifying the cache directory for compiled Stan
  models. The sccomp version will be automatically appended to ensure
  version isolation. Default is `sccomp_stan_models_cache_dir` which
  points to `~/.sccomp_models`.

- cleanup_draw_files:

  Logical, whether to automatically delete Stan draw CSV files after
  extracting results. These files can be large (MBs to GBs) and are
  typically only needed during the analysis session. Default is TRUE to
  save disk space. Set to FALSE if you need to inspect draw files for
  debugging.

- approximate_posterior_inference:

  DEPRECATED, use the `variational_inference` argument.

- variational_inference:

  DEPRECATED Logical, whether to use variational Bayes for posterior
  inference. It is faster and convenient. Setting this argument to
  `FALSE` runs full Bayesian (Hamiltonian Monte Carlo) inference, which
  is slower but the gold standard.

- ...:

  Additional arguments passed to the
  [`cmdstanr::sample`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)
  function.

## Value

A tibble (`tbl`), with the following columns:

- cell_group - The cell groups being tested.

- parameter - The parameter being estimated from the design matrix
  described by the input formula_composition and formula_variability.

- factor - The covariate factor in the formula, if applicable (e.g., not
  present for Intercept or contrasts).

- c_lower - Lower (2.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_effect - Mean of the posterior distribution for a composition (c)
  parameter.

- c_upper - Upper (97.5%) quantile of the posterior distribution for a
  composition (c) parameter.

- c_pH0 - Probability of the c_effect being smaller or bigger than the
  `test_composition_above_logit_fold_change` argument.

- c_FDR - False discovery rate of the c_effect being smaller or bigger
  than the `test_composition_above_logit_fold_change` argument. False
  discovery rate for Bayesian models is calculated differently from
  frequentists models, as detailed in Mangiola et al, PNAS 2023.

- c_n_eff - Effective sample size, the number of independent draws in
  the sample. The higher, the better.

- c_R_k_hat - R statistic, a measure of chain equilibrium, should be
  within 0.05 of 1.0.

- v_lower - Lower (2.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_effect - Mean of the posterior distribution for a variability (v)
  parameter.

- v_upper - Upper (97.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_pH0 - Probability of the v_effect being smaller or bigger than the
  `test_composition_above_logit_fold_change` argument.

- v_FDR - False discovery rate of the v_effect being smaller or bigger
  than the `test_composition_above_logit_fold_change` argument. False
  discovery rate for Bayesian models is calculated differently from
  frequentists models, as detailed in Mangiola et al, PNAS 2023.

- v_n_eff - Effective sample size for a variability (v) parameter.

- v_R_k_hat - R statistic for a variability (v) parameter, a measure of
  chain equilibrium.

The function also attaches several attributes to the result:

- count_data - The original count data used in the analysis, stored as
  an attribute for efficient access.

- model_input - The model input data used for fitting.

- formula_composition - The formula used for composition modeling.

- formula_variability - The formula used for variability modeling.

- fit - The Stan fit object (if pass_fit = TRUE).

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
    
    estimate = sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
    ) |>
    sccomp_remove_outliers(cores = 1)
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481680.232296 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.540e-03   2.761e-01    1.000e+00  1.000e+00      2944 -3.704e+03 -3.711e+03                   
#> Path [1] :Best Iter: [50] ELBO (-3704.092856) evaluations: (2944) 
#> Path [2] :Initial log joint density = -481542.591647 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.887e-03   2.088e-01    1.000e+00  1.000e+00      3169 -3.709e+03 -3.709e+03                   
#> Path [2] :Best Iter: [42] ELBO (-3708.693020) evaluations: (3169) 
#> Path [3] :Initial log joint density = -481564.188856 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.381e-03   2.349e-01    1.000e+00  1.000e+00      2991 -3.708e+03 -3.717e+03                   
#> Path [3] :Best Iter: [51] ELBO (-3708.499393) evaluations: (2991) 
#> Path [4] :Initial log joint density = -481103.869819 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.177e-02   2.549e-01    1.000e+00  1.000e+00      2749 -3.708e+03 -3.717e+03                   
#> Path [4] :Best Iter: [43] ELBO (-3707.991849) evaluations: (2749) 
#> Path [5] :Initial log joint density = -481773.054526 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.961e-03   1.372e-01    1.000e+00  1.000e+00      3560 -3.701e+03 -3.711e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3700.968991) evaluations: (3560) 
#> Path [6] :Initial log joint density = -481923.659340 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      4.546e-03   2.197e-01    8.655e-01  8.655e-01      3763 -3.699e+03 -3.713e+03                   
#> Path [6] :Best Iter: [61] ELBO (-3698.686262) evaluations: (3763) 
#> Path [7] :Initial log joint density = -481876.677054 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.468e-02   2.256e-01    1.000e+00  1.000e+00      2947 -3.709e+03 -3.710e+03                   
#> Path [7] :Best Iter: [50] ELBO (-3708.894268) evaluations: (2947) 
#> Path [8] :Initial log joint density = -481690.677026 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.316e-02   2.458e-01    1.000e+00  1.000e+00      2944 -3.707e+03 -3.709e+03                   
#> Path [8] :Best Iter: [43] ELBO (-3707.030016) evaluations: (2944) 
#> Path [9] :Initial log joint density = -481844.644896 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.145e-02   3.017e-01    1.000e+00  1.000e+00      2953 -3.707e+03 -3.723e+03                   
#> Path [9] :Best Iter: [51] ELBO (-3706.811202) evaluations: (2953) 
#> Path [10] :Initial log joint density = -482217.535647 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.289e-03   1.361e-01    6.556e-01  6.556e-01      3633 -3.702e+03 -3.718e+03                   
#> Path [10] :Best Iter: [60] ELBO (-3701.654024) evaluations: (3633) 
#> Path [11] :Initial log joint density = -481494.623315 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.437e-03   2.055e-01    1.000e+00  1.000e+00      3054 -3.710e+03 -3.720e+03                   
#> Path [11] :Best Iter: [45] ELBO (-3710.084562) evaluations: (3054) 
#> Path [12] :Initial log joint density = -481680.659676 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.237e-03   1.976e-01    1.000e+00  1.000e+00      3443 -3.700e+03 -3.700e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3699.634202) evaluations: (3443) 
#> Path [13] :Initial log joint density = -481540.579172 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.410e-03   2.170e-01    1.000e+00  1.000e+00      2910 -3.706e+03 -3.709e+03                   
#> Path [13] :Best Iter: [37] ELBO (-3705.631002) evaluations: (2910) 
#> Path [14] :Initial log joint density = -481633.961566 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.153e-03   2.745e-01    6.170e-01  6.170e-01      3586 -3.699e+03 -3.714e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3699.118978) evaluations: (3586) 
#> Path [15] :Initial log joint density = -482486.284222 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.026e-03   2.776e-01    1.000e+00  1.000e+00      3207 -3.712e+03 -3.709e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3709.280606) evaluations: (3207) 
#> Path [16] :Initial log joint density = -481499.980728 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.736e-03   2.482e-01    9.618e-01  9.618e-01      2941 -3.708e+03 -3.719e+03                   
#> Path [16] :Best Iter: [45] ELBO (-3708.042270) evaluations: (2941) 
#> Path [17] :Initial log joint density = -484937.595106 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.313e-03   2.842e-01    1.000e+00  1.000e+00      3247 -3.704e+03 -3.707e+03                   
#> Path [17] :Best Iter: [48] ELBO (-3704.494238) evaluations: (3247) 
#> Path [18] :Initial log joint density = -481538.125972 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.304e-02   2.300e-01    1.000e+00  1.000e+00      3271 -3.704e+03 -3.704e+03                   
#> Path [18] :Best Iter: [52] ELBO (-3704.010967) evaluations: (3271) 
#> Path [19] :Initial log joint density = -481827.799516 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.347e-03   2.148e-01    1.000e+00  1.000e+00      3320 -3.704e+03 -3.701e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3700.687039) evaluations: (3320) 
#> Path [20] :Initial log joint density = -481961.221529 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.244e-03   2.106e-01    7.016e-01  7.016e-01      3120 -3.706e+03 -3.718e+03                   
#> Path [20] :Best Iter: [43] ELBO (-3706.059399) evaluations: (3120) 
#> Path [21] :Initial log joint density = -481796.435256 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.849e-03   2.828e-01    5.860e-01  5.860e-01      3263 -3.705e+03 -3.711e+03                   
#> Path [21] :Best Iter: [54] ELBO (-3705.264430) evaluations: (3263) 
#> Path [22] :Initial log joint density = -482422.806484 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.232e-02   1.425e-01    1.000e+00  1.000e+00      3472 -3.702e+03 -3.705e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3701.593143) evaluations: (3472) 
#> Path [23] :Initial log joint density = -483302.705764 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.610e-03   2.191e-01    1.000e+00  1.000e+00      3040 -3.705e+03 -3.712e+03                   
#> Path [23] :Best Iter: [40] ELBO (-3704.930650) evaluations: (3040) 
#> Path [24] :Initial log joint density = -481479.050570 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.322e-02   2.330e-01    1.000e+00  1.000e+00      2998 -3.706e+03 -3.714e+03                   
#> Path [24] :Best Iter: [42] ELBO (-3706.405345) evaluations: (2998) 
#> Path [25] :Initial log joint density = -481928.631132 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.337e-03   2.287e-01    7.401e-01  7.401e-01      3237 -3.704e+03 -3.702e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3701.505431) evaluations: (3237) 
#> Path [26] :Initial log joint density = -481635.395886 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.535e-03   1.818e-01    1.000e+00  1.000e+00      3220 -3.703e+03 -3.703e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3702.871254) evaluations: (3220) 
#> Path [27] :Initial log joint density = -481479.988574 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.435e-03   1.808e-01    9.453e-01  9.453e-01      3471 -3.697e+03 -3.704e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3696.815580) evaluations: (3471) 
#> Path [28] :Initial log joint density = -481538.198825 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.220e-02   2.802e-01    1.000e+00  1.000e+00      3321 -3.702e+03 -3.699e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3699.156470) evaluations: (3321) 
#> Path [29] :Initial log joint density = -481487.691881 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      2.593e-03   2.649e-01    7.312e-01  7.312e-01      3253 -3.702e+03 -3.711e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3702.158187) evaluations: (3253) 
#> Path [30] :Initial log joint density = -481768.300300 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.255e-02   3.692e-01    1.000e+00  1.000e+00      2980 -3.707e+03 -3.714e+03                   
#> Path [30] :Best Iter: [39] ELBO (-3707.006795) evaluations: (2980) 
#> Path [31] :Initial log joint density = -481849.864078 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.566e-03   2.167e-01    1.000e+00  1.000e+00      3192 -3.709e+03 -3.708e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3707.925640) evaluations: (3192) 
#> Path [32] :Initial log joint density = -481658.505034 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      5.057e-03   2.042e-01    9.988e-01  9.988e-01      3678 -3.703e+03 -3.705e+03                   
#> Path [32] :Best Iter: [58] ELBO (-3703.416752) evaluations: (3678) 
#> Path [33] :Initial log joint density = -482149.041325 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.010e-02   2.227e-01    1.000e+00  1.000e+00      3466 -3.700e+03 -3.707e+03                   
#> Path [33] :Best Iter: [57] ELBO (-3700.178251) evaluations: (3466) 
#> Path [34] :Initial log joint density = -481848.719421 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.900e-03   1.761e-01    1.000e+00  1.000e+00      3136 -3.707e+03 -3.699e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3698.917458) evaluations: (3136) 
#> Path [35] :Initial log joint density = -483905.726286 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.033e-02   3.590e-01    1.000e+00  1.000e+00      3186 -3.702e+03 -3.707e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3701.785236) evaluations: (3186) 
#> Path [36] :Initial log joint density = -483400.862807 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.067e-03   1.645e-01    7.627e-01  7.627e-01      3300 -3.701e+03 -3.714e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3701.282616) evaluations: (3300) 
#> Path [37] :Initial log joint density = -481475.597588 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.824e-02   2.927e-01    1.000e+00  1.000e+00      3245 -3.697e+03 -3.701e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3696.740763) evaluations: (3245) 
#> Path [38] :Initial log joint density = -481560.913976 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.911e-03   1.282e-01    7.995e-01  7.995e-01      3328 -3.700e+03 -3.712e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3699.768217) evaluations: (3328) 
#> Path [39] :Initial log joint density = -481557.870273 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.216e-02   3.442e-01    1.000e+00  1.000e+00      3391 -3.699e+03 -3.705e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3698.872637) evaluations: (3391) 
#> Path [40] :Initial log joint density = -483190.950891 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.652e-03   2.224e-01    1.000e+00  1.000e+00      3108 -3.708e+03 -3.706e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3705.682084) evaluations: (3108) 
#> Path [41] :Initial log joint density = -481730.207334 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.300e-03   2.375e-01    6.491e-01  6.491e-01      3158 -3.708e+03 -3.715e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3707.965913) evaluations: (3158) 
#> Path [42] :Initial log joint density = -481749.320833 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.822e-03   2.412e-01    1.000e+00  1.000e+00      3033 -3.710e+03 -3.714e+03                   
#> Path [42] :Best Iter: [50] ELBO (-3710.416705) evaluations: (3033) 
#> Path [43] :Initial log joint density = -483009.064425 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.619e-03   3.149e-01    8.615e-01  8.615e-01      3470 -3.701e+03 -3.710e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3700.831508) evaluations: (3470) 
#> Path [44] :Initial log joint density = -481682.196995 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.331e-03   2.162e-01    1.000e+00  1.000e+00      2751 -3.707e+03 -3.711e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3707.026530) evaluations: (2751) 
#> Path [45] :Initial log joint density = -482066.760406 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.727e-03   1.826e-01    1.000e+00  1.000e+00      3491 -3.701e+03 -3.707e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3700.950536) evaluations: (3491) 
#> Path [46] :Initial log joint density = -481364.621882 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.671e-03   3.409e-01    4.418e-01  1.000e+00      3082 -3.709e+03 -3.721e+03                   
#> Path [46] :Best Iter: [34] ELBO (-3709.487833) evaluations: (3082) 
#> Path [47] :Initial log joint density = -482177.412466 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.216e-03   2.165e-01    4.049e-01  1.000e+00      3644 -3.700e+03 -3.704e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3700.107744) evaluations: (3644) 
#> Path [48] :Initial log joint density = -482080.464177 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      6.407e-03   2.487e-01    9.078e-01  9.078e-01      3821 -3.698e+03 -3.711e+03                   
#> Path [48] :Best Iter: [61] ELBO (-3698.201674) evaluations: (3821) 
#> Path [49] :Initial log joint density = -481561.704521 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.515e-03   2.930e-01    5.735e-01  5.735e-01      3202 -3.708e+03 -3.711e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3707.631251) evaluations: (3202) 
#> Path [50] :Initial log joint density = -481687.519296 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.842e-03   2.522e-01    1.000e+00  1.000e+00      3189 -3.704e+03 -3.706e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3704.444422) evaluations: (3189) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier identification - step 1/2
#> Loading model from cache...
#> Path [1] :Initial log joint density = -429349.632479 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.964e-02   1.510e+03    2.532e-02  4.570e-02      8355 -3.319e+03 -1.780e+05                   
#> Path [1] :Best Iter: [53] ELBO (-3318.960850) evaluations: (8355) 
#> Path [2] :Initial log joint density = -429118.487965 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.027e-01   1.418e+03    3.140e-02  3.140e-02      8126 -3.317e+03 -3.632e+05                   
#> Path [2] :Best Iter: [48] ELBO (-3316.880479) evaluations: (8126) 
#> Path [3] :Initial log joint density = -429125.782370 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.413e-02   5.104e+03    1.627e-02  1.627e-02      8174 -3.319e+03 -9.010e+03                   
#> Path [3] :Best Iter: [40] ELBO (-3319.304228) evaluations: (8174) 
#> Path [4] :Initial log joint density = -429142.382963 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.942e-02   1.387e+03    4.474e-02  4.474e-02      8184 -3.317e+03 -2.969e+04                   
#> Path [4] :Best Iter: [44] ELBO (-3317.055101) evaluations: (8184) 
#> Path [5] :Initial log joint density = -429285.921097 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.054e-02   2.080e+03    3.130e-02  3.130e-02      8329 -3.319e+03 -7.656e+03                   
#> Path [5] :Best Iter: [41] ELBO (-3318.721473) evaluations: (8329) 
#> Path [6] :Initial log joint density = -429507.664724 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.624e-02   3.744e+03    2.727e-02  2.727e-02      8374 -3.321e+03 -9.860e+05                   
#> Path [6] :Best Iter: [51] ELBO (-3320.617363) evaluations: (8374) 
#> Path [7] :Initial log joint density = -430290.926839 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.859e-02   1.251e+03    1.949e-02  1.949e-02      8359 -3.314e+03 -1.038e+04                   
#> Path [7] :Best Iter: [55] ELBO (-3314.001620) evaluations: (8359) 
#> Path [8] :Initial log joint density = -429166.644964 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.006e-02   1.802e+03    1.435e-02  4.982e-02      8387 -3.318e+03 -2.397e+05                   
#> Path [8] :Best Iter: [45] ELBO (-3317.609926) evaluations: (8387) 
#> Path [9] :Initial log joint density = -429169.120485 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.958e-02   1.620e+03    2.389e-02  2.389e-02      8051 -3.323e+03 -6.035e+04                   
#> Path [9] :Best Iter: [44] ELBO (-3322.544193) evaluations: (8051) 
#> Path [10] :Initial log joint density = -429012.253677 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.633e-02   2.262e+03    1.867e-02  4.538e-02      8150 -3.320e+03 -6.983e+03                   
#> Path [10] :Best Iter: [43] ELBO (-3319.890155) evaluations: (8150) 
#> Path [11] :Initial log joint density = -428931.318195 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.568e-02   2.380e+03    2.106e-02  2.106e-02      8141 -3.320e+03 -7.593e+03                   
#> Path [11] :Best Iter: [41] ELBO (-3320.237203) evaluations: (8141) 
#> Path [12] :Initial log joint density = -429184.153311 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.414e-02   2.339e+03    2.933e-02  2.933e-02      8010 -3.321e+03 -1.144e+05                   
#> Path [12] :Best Iter: [45] ELBO (-3320.895763) evaluations: (8010) 
#> Path [13] :Initial log joint density = -429484.582430 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.949e-02   3.936e+03    1.454e-02  6.162e-02      8272 -3.321e+03 -1.094e+05                   
#> Path [13] :Best Iter: [47] ELBO (-3320.775241) evaluations: (8272) 
#> Path [14] :Initial log joint density = -431943.406643 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.564e-02   1.537e+03    2.486e-02  2.486e-02      8120 -3.318e+03 -9.493e+03                   
#> Path [14] :Best Iter: [46] ELBO (-3317.646825) evaluations: (8120) 
#> Path [15] :Initial log joint density = -430089.443863 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.323e-02   9.445e+02    3.643e-02  3.643e-02      8081 -3.311e+03 -1.494e+04                   
#> Path [15] :Best Iter: [55] ELBO (-3310.790762) evaluations: (8081) 
#> Path [16] :Initial log joint density = -429472.725377 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.517e-02   1.550e+03    2.557e-02  2.557e-02      8606 -3.318e+03 -2.600e+04                   
#> Path [16] :Best Iter: [60] ELBO (-3318.130487) evaluations: (8606) 
#> Path [17] :Initial log joint density = -429663.586130 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.871e-02   9.956e+02    3.863e-02  3.863e-02      8159 -3.310e+03 -1.261e+04                   
#> Path [17] :Best Iter: [55] ELBO (-3309.771536) evaluations: (8159) 
#> Path [18] :Initial log joint density = -429442.161702 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.358e-02   2.236e+03    2.840e-02  2.840e-02      8271 -3.320e+03 -1.441e+04                   
#> Path [18] :Best Iter: [45] ELBO (-3320.123731) evaluations: (8271) 
#> Path [19] :Initial log joint density = -434040.196518 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.102e-02   2.244e+03    4.225e-02  4.225e-02      8156 -3.319e+03 -4.615e+03                   
#> Path [19] :Best Iter: [53] ELBO (-3318.782446) evaluations: (8156) 
#> Path [20] :Initial log joint density = -429588.954497 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.945e-02   1.198e+03    3.499e-02  3.499e-02      7984 -3.312e+03 -9.192e+03                   
#> Path [20] :Best Iter: [60] ELBO (-3311.522470) evaluations: (7984) 
#> Path [21] :Initial log joint density = -429606.600175 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.050e-02   1.785e+03    2.240e-02  2.240e-02      8308 -3.319e+03 -1.185e+05                   
#> Path [21] :Best Iter: [41] ELBO (-3318.687529) evaluations: (8308) 
#> Path [22] :Initial log joint density = -429427.659350 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.245e-02   4.061e+03    2.512e-02  2.512e-02      8067 -3.324e+03 -1.515e+06                   
#> Path [22] :Best Iter: [50] ELBO (-3323.636911) evaluations: (8067) 
#> Path [23] :Initial log joint density = -430457.305145 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.263e-02   3.183e+03    2.112e-02  2.112e-02      8162 -3.317e+03 -1.575e+06                   
#> Path [23] :Best Iter: [56] ELBO (-3316.667058) evaluations: (8162) 
#> Path [24] :Initial log joint density = -429057.818246 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.297e-02   1.766e+03    1.461e-02  1.461e-02      8152 -3.320e+03 -1.290e+07                   
#> Path [24] :Best Iter: [48] ELBO (-3320.032393) evaluations: (8152) 
#> Path [25] :Initial log joint density = -429118.436234 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.285e-02   1.648e+03    2.485e-02  2.485e-02      8261 -3.319e+03 -7.169e+06                   
#> Path [25] :Best Iter: [48] ELBO (-3318.667224) evaluations: (8261) 
#> Path [26] :Initial log joint density = -429509.106869 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.730e-02   1.943e+03    2.427e-02  2.427e-02      8336 -3.321e+03 -7.745e+04                   
#> Path [26] :Best Iter: [46] ELBO (-3320.523653) evaluations: (8336) 
#> Path [27] :Initial log joint density = -429447.426919 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.163e-02   2.906e+03    2.388e-02  2.388e-02      7913 -3.320e+03 -2.046e+04                   
#> Path [27] :Best Iter: [45] ELBO (-3319.560583) evaluations: (7913) 
#> Path [28] :Initial log joint density = -429903.367236 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.069e-02   6.804e+02    6.874e-02  6.874e-02      8251 -3.317e+03 -6.143e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3317.367314) evaluations: (8251) 
#> Path [29] :Initial log joint density = -429478.327319 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.291e-02   2.315e+03    1.326e-02  2.583e-02      8600 -3.318e+03 -2.581e+04                   
#> Path [29] :Best Iter: [44] ELBO (-3317.830038) evaluations: (8600) 
#> Path [30] :Initial log joint density = -430698.595508 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.103e-02   1.341e+03    4.278e-02  4.278e-02      8167 -3.318e+03 -1.523e+05                   
#> Path [30] :Best Iter: [56] ELBO (-3317.549422) evaluations: (8167) 
#> Path [31] :Initial log joint density = -429233.394967 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.263e-02   3.000e+03    2.532e-02  2.532e-02      8152 -3.318e+03 -1.686e+04                   
#> Path [31] :Best Iter: [51] ELBO (-3318.012906) evaluations: (8152) 
#> Path [32] :Initial log joint density = -429152.789248 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.798e-02   2.203e+03    4.004e-02  4.004e-02      8269 -3.321e+03 -1.822e+05                   
#> Path [32] :Best Iter: [46] ELBO (-3321.163169) evaluations: (8269) 
#> Path [33] :Initial log joint density = -429471.797509 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.069e-01   9.063e+02    5.410e-02  5.410e-02      8038 -3.314e+03 -1.524e+05                   
#> Path [33] :Best Iter: [55] ELBO (-3314.183861) evaluations: (8038) 
#> Path [34] :Initial log joint density = -430431.604586 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.508e-02   1.240e+03    3.796e-02  7.320e-02      8131 -3.319e+03 -1.831e+06                   
#> Path [34] :Best Iter: [55] ELBO (-3319.179315) evaluations: (8131) 
#> Path [35] :Initial log joint density = -429796.713446 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.188e-01   2.996e+03    4.328e-02  4.328e-02      8221 -3.311e+03 -1.803e+06                   
#> Path [35] :Best Iter: [58] ELBO (-3311.022020) evaluations: (8221) 
#> Path [36] :Initial log joint density = -429166.250360 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.260e-01   4.338e+03    1.402e-02  4.529e-02      8057 -3.314e+03 -7.782e+09                   
#> Path [36] :Best Iter: [44] ELBO (-3314.382335) evaluations: (8057) 
#> Path [37] :Initial log joint density = -429144.991625 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.685e-02   1.648e+03    3.654e-02  3.654e-02      8259 -3.319e+03 -8.277e+05                   
#> Path [37] :Best Iter: [48] ELBO (-3318.961688) evaluations: (8259) 
#> Path [38] :Initial log joint density = -429224.180586 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.626e-02   4.282e+03    1.700e-02  1.700e-02      8181 -3.320e+03 -1.609e+08                   
#> Path [38] :Best Iter: [46] ELBO (-3319.892983) evaluations: (8181) 
#> Path [39] :Initial log joint density = -429228.228753 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.421e-02   3.381e+03    1.896e-02  1.896e-02      8303 -3.318e+03 -5.323e+08                   
#> Path [39] :Best Iter: [50] ELBO (-3318.361180) evaluations: (8303) 
#> Path [40] :Initial log joint density = -430241.658497 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.902e-02   2.131e+03    2.191e-02  4.738e-02      8347 -3.314e+03 -2.612e+04                   
#> Path [40] :Best Iter: [47] ELBO (-3313.823163) evaluations: (8347) 
#> Path [41] :Initial log joint density = -429292.563606 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.010e-01   1.173e+03    2.665e-02  5.475e-02      7835 -3.312e+03 -6.955e+07                   
#> Path [41] :Best Iter: [61] ELBO (-3311.618365) evaluations: (7835) 
#> Path [42] :Initial log joint density = -429437.747973 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.962e-02   1.817e+03    2.019e-02  5.744e-02      8351 -3.323e+03 -1.360e+07                   
#> Path [42] :Best Iter: [45] ELBO (-3322.551684) evaluations: (8351) 
#> Path [43] :Initial log joint density = -429404.420228 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.901e-02   1.506e+03    5.221e-02  5.221e-02      8073 -3.314e+03 -1.651e+04                   
#> Path [43] :Best Iter: [55] ELBO (-3313.848011) evaluations: (8073) 
#> Path [44] :Initial log joint density = -432723.509579 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.428e-02   1.010e+03    3.585e-02  3.585e-02      8209 -3.322e+03 -6.987e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3321.662808) evaluations: (8209) 
#> Path [45] :Initial log joint density = -429369.264751 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.718e-02   1.292e+03    4.681e-02  4.681e-02      8349 -3.318e+03 -1.202e+04                   
#> Path [45] :Best Iter: [48] ELBO (-3318.349427) evaluations: (8349) 
#> Path [46] :Initial log joint density = -429313.221164 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.776e-02   4.526e+03    1.254e-02  1.254e-02      8522 -3.316e+03 -7.509e+05                   
#> Path [46] :Best Iter: [55] ELBO (-3316.246835) evaluations: (8522) 
#> Path [47] :Initial log joint density = -429249.559302 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.159e-02   2.608e+03    3.189e-02  3.189e-02      8367 -3.319e+03 -2.696e+05                   
#> Path [47] :Best Iter: [48] ELBO (-3319.060306) evaluations: (8367) 
#> Path [48] :Initial log joint density = -431590.898810 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.846e-02   1.464e+03    2.475e-02  2.475e-02      8071 -3.313e+03 -1.345e+04                   
#> Path [48] :Best Iter: [58] ELBO (-3312.977577) evaluations: (8071) 
#> Path [49] :Initial log joint density = -429264.550487 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.804e-02   1.431e+03    3.442e-02  3.442e-02      8202 -3.322e+03 -2.768e+07                   
#> Path [49] :Best Iter: [51] ELBO (-3321.622896) evaluations: (8202) 
#> Path [50] :Initial log joint density = -429906.425362 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.115e-02   6.752e+02    3.781e-02  3.781e-02      8206 -3.311e+03 -4.209e+04                   
#> Path [50] :Best Iter: [55] ELBO (-3311.161607) evaluations: (8206) 
#> Finished in  30.6 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -429349.632479 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.964e-02   1.510e+03    2.532e-02  4.570e-02      8355 -3.319e+03 -1.780e+05                   
#> Path [1] :Best Iter: [53] ELBO (-3318.960850) evaluations: (8355) 
#> Path [2] :Initial log joint density = -429118.487965 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.027e-01   1.418e+03    3.140e-02  3.140e-02      8126 -3.317e+03 -3.632e+05                   
#> Path [2] :Best Iter: [48] ELBO (-3316.880479) evaluations: (8126) 
#> Path [3] :Initial log joint density = -429125.782370 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.413e-02   5.104e+03    1.627e-02  1.627e-02      8174 -3.319e+03 -9.010e+03                   
#> Path [3] :Best Iter: [40] ELBO (-3319.304228) evaluations: (8174) 
#> Path [4] :Initial log joint density = -429142.382963 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.942e-02   1.387e+03    4.474e-02  4.474e-02      8184 -3.317e+03 -2.969e+04                   
#> Path [4] :Best Iter: [44] ELBO (-3317.055101) evaluations: (8184) 
#> Path [5] :Initial log joint density = -429285.921097 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.054e-02   2.080e+03    3.130e-02  3.130e-02      8329 -3.319e+03 -7.656e+03                   
#> Path [5] :Best Iter: [41] ELBO (-3318.721473) evaluations: (8329) 
#> Path [6] :Initial log joint density = -429507.664724 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.624e-02   3.744e+03    2.727e-02  2.727e-02      8374 -3.321e+03 -9.860e+05                   
#> Path [6] :Best Iter: [51] ELBO (-3320.617363) evaluations: (8374) 
#> Path [7] :Initial log joint density = -430290.926839 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.859e-02   1.251e+03    1.949e-02  1.949e-02      8359 -3.314e+03 -1.038e+04                   
#> Path [7] :Best Iter: [55] ELBO (-3314.001620) evaluations: (8359) 
#> Path [8] :Initial log joint density = -429166.644964 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.006e-02   1.802e+03    1.435e-02  4.982e-02      8387 -3.318e+03 -2.397e+05                   
#> Path [8] :Best Iter: [45] ELBO (-3317.609926) evaluations: (8387) 
#> Path [9] :Initial log joint density = -429169.120485 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.958e-02   1.620e+03    2.389e-02  2.389e-02      8051 -3.323e+03 -6.035e+04                   
#> Path [9] :Best Iter: [44] ELBO (-3322.544193) evaluations: (8051) 
#> Path [10] :Initial log joint density = -429012.253677 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.633e-02   2.262e+03    1.867e-02  4.538e-02      8150 -3.320e+03 -6.983e+03                   
#> Path [10] :Best Iter: [43] ELBO (-3319.890155) evaluations: (8150) 
#> Path [11] :Initial log joint density = -428931.318195 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.568e-02   2.380e+03    2.106e-02  2.106e-02      8141 -3.320e+03 -7.593e+03                   
#> Path [11] :Best Iter: [41] ELBO (-3320.237203) evaluations: (8141) 
#> Path [12] :Initial log joint density = -429184.153311 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.414e-02   2.339e+03    2.933e-02  2.933e-02      8010 -3.321e+03 -1.144e+05                   
#> Path [12] :Best Iter: [45] ELBO (-3320.895763) evaluations: (8010) 
#> Path [13] :Initial log joint density = -429484.582430 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.949e-02   3.936e+03    1.454e-02  6.162e-02      8272 -3.321e+03 -1.094e+05                   
#> Path [13] :Best Iter: [47] ELBO (-3320.775241) evaluations: (8272) 
#> Path [14] :Initial log joint density = -431943.406643 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.564e-02   1.537e+03    2.486e-02  2.486e-02      8120 -3.318e+03 -9.493e+03                   
#> Path [14] :Best Iter: [46] ELBO (-3317.646825) evaluations: (8120) 
#> Path [15] :Initial log joint density = -430089.443863 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.323e-02   9.445e+02    3.643e-02  3.643e-02      8081 -3.311e+03 -1.494e+04                   
#> Path [15] :Best Iter: [55] ELBO (-3310.790762) evaluations: (8081) 
#> Path [16] :Initial log joint density = -429472.725377 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.517e-02   1.550e+03    2.557e-02  2.557e-02      8606 -3.318e+03 -2.600e+04                   
#> Path [16] :Best Iter: [60] ELBO (-3318.130487) evaluations: (8606) 
#> Path [17] :Initial log joint density = -429663.586130 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.871e-02   9.956e+02    3.863e-02  3.863e-02      8159 -3.310e+03 -1.261e+04                   
#> Path [17] :Best Iter: [55] ELBO (-3309.771536) evaluations: (8159) 
#> Path [18] :Initial log joint density = -429442.161702 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.358e-02   2.236e+03    2.840e-02  2.840e-02      8271 -3.320e+03 -1.441e+04                   
#> Path [18] :Best Iter: [45] ELBO (-3320.123731) evaluations: (8271) 
#> Path [19] :Initial log joint density = -434040.196518 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.102e-02   2.244e+03    4.225e-02  4.225e-02      8156 -3.319e+03 -4.615e+03                   
#> Path [19] :Best Iter: [53] ELBO (-3318.782446) evaluations: (8156) 
#> Path [20] :Initial log joint density = -429588.954497 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.945e-02   1.198e+03    3.499e-02  3.499e-02      7984 -3.312e+03 -9.192e+03                   
#> Path [20] :Best Iter: [60] ELBO (-3311.522470) evaluations: (7984) 
#> Path [21] :Initial log joint density = -429606.600175 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.050e-02   1.785e+03    2.240e-02  2.240e-02      8308 -3.319e+03 -1.185e+05                   
#> Path [21] :Best Iter: [41] ELBO (-3318.687529) evaluations: (8308) 
#> Path [22] :Initial log joint density = -429427.659350 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.245e-02   4.061e+03    2.512e-02  2.512e-02      8067 -3.324e+03 -1.515e+06                   
#> Path [22] :Best Iter: [50] ELBO (-3323.636911) evaluations: (8067) 
#> Path [23] :Initial log joint density = -430457.305145 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.263e-02   3.183e+03    2.112e-02  2.112e-02      8162 -3.317e+03 -1.575e+06                   
#> Path [23] :Best Iter: [56] ELBO (-3316.667058) evaluations: (8162) 
#> Path [24] :Initial log joint density = -429057.818246 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.297e-02   1.766e+03    1.461e-02  1.461e-02      8152 -3.320e+03 -1.290e+07                   
#> Path [24] :Best Iter: [48] ELBO (-3320.032393) evaluations: (8152) 
#> Path [25] :Initial log joint density = -429118.436234 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.285e-02   1.648e+03    2.485e-02  2.485e-02      8261 -3.319e+03 -7.169e+06                   
#> Path [25] :Best Iter: [48] ELBO (-3318.667224) evaluations: (8261) 
#> Path [26] :Initial log joint density = -429509.106869 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.730e-02   1.943e+03    2.427e-02  2.427e-02      8336 -3.321e+03 -7.745e+04                   
#> Path [26] :Best Iter: [46] ELBO (-3320.523653) evaluations: (8336) 
#> Path [27] :Initial log joint density = -429447.426919 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.163e-02   2.906e+03    2.388e-02  2.388e-02      7913 -3.320e+03 -2.046e+04                   
#> Path [27] :Best Iter: [45] ELBO (-3319.560583) evaluations: (7913) 
#> Path [28] :Initial log joint density = -429903.367236 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.069e-02   6.804e+02    6.874e-02  6.874e-02      8251 -3.317e+03 -6.143e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3317.367314) evaluations: (8251) 
#> Path [29] :Initial log joint density = -429478.327319 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.291e-02   2.315e+03    1.326e-02  2.583e-02      8600 -3.318e+03 -2.581e+04                   
#> Path [29] :Best Iter: [44] ELBO (-3317.830038) evaluations: (8600) 
#> Path [30] :Initial log joint density = -430698.595508 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.103e-02   1.341e+03    4.278e-02  4.278e-02      8167 -3.318e+03 -1.523e+05                   
#> Path [30] :Best Iter: [56] ELBO (-3317.549422) evaluations: (8167) 
#> Path [31] :Initial log joint density = -429233.394967 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.263e-02   3.000e+03    2.532e-02  2.532e-02      8152 -3.318e+03 -1.686e+04                   
#> Path [31] :Best Iter: [51] ELBO (-3318.012906) evaluations: (8152) 
#> Path [32] :Initial log joint density = -429152.789248 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.798e-02   2.203e+03    4.004e-02  4.004e-02      8269 -3.321e+03 -1.822e+05                   
#> Path [32] :Best Iter: [46] ELBO (-3321.163169) evaluations: (8269) 
#> Path [33] :Initial log joint density = -429471.797509 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.069e-01   9.063e+02    5.410e-02  5.410e-02      8038 -3.314e+03 -1.524e+05                   
#> Path [33] :Best Iter: [55] ELBO (-3314.183861) evaluations: (8038) 
#> Path [34] :Initial log joint density = -430431.604586 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.508e-02   1.240e+03    3.796e-02  7.320e-02      8131 -3.319e+03 -1.831e+06                   
#> Path [34] :Best Iter: [55] ELBO (-3319.179315) evaluations: (8131) 
#> Path [35] :Initial log joint density = -429796.713446 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.188e-01   2.996e+03    4.328e-02  4.328e-02      8221 -3.311e+03 -1.803e+06                   
#> Path [35] :Best Iter: [58] ELBO (-3311.022020) evaluations: (8221) 
#> Path [36] :Initial log joint density = -429166.250360 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.260e-01   4.338e+03    1.402e-02  4.529e-02      8057 -3.314e+03 -7.782e+09                   
#> Path [36] :Best Iter: [44] ELBO (-3314.382335) evaluations: (8057) 
#> Path [37] :Initial log joint density = -429144.991625 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.685e-02   1.648e+03    3.654e-02  3.654e-02      8259 -3.319e+03 -8.277e+05                   
#> Path [37] :Best Iter: [48] ELBO (-3318.961688) evaluations: (8259) 
#> Path [38] :Initial log joint density = -429224.180586 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      9.626e-02   4.282e+03    1.700e-02  1.700e-02      8181 -3.320e+03 -1.609e+08                   
#> Path [38] :Best Iter: [46] ELBO (-3319.892983) evaluations: (8181) 
#> Path [39] :Initial log joint density = -429228.228753 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.421e-02   3.381e+03    1.896e-02  1.896e-02      8303 -3.318e+03 -5.323e+08                   
#> Path [39] :Best Iter: [50] ELBO (-3318.361180) evaluations: (8303) 
#> Path [40] :Initial log joint density = -430241.658497 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.902e-02   2.131e+03    2.191e-02  4.738e-02      8347 -3.314e+03 -2.612e+04                   
#> Path [40] :Best Iter: [47] ELBO (-3313.823163) evaluations: (8347) 
#> Path [41] :Initial log joint density = -429292.563606 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      1.010e-01   1.173e+03    2.665e-02  5.475e-02      7835 -3.312e+03 -6.955e+07                   
#> Path [41] :Best Iter: [61] ELBO (-3311.618365) evaluations: (7835) 
#> Path [42] :Initial log joint density = -429437.747973 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      8.962e-02   1.817e+03    2.019e-02  5.744e-02      8351 -3.323e+03 -1.360e+07                   
#> Path [42] :Best Iter: [45] ELBO (-3322.551684) evaluations: (8351) 
#> Path [43] :Initial log joint density = -429404.420228 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      6.901e-02   1.506e+03    5.221e-02  5.221e-02      8073 -3.314e+03 -1.651e+04                   
#> Path [43] :Best Iter: [55] ELBO (-3313.848011) evaluations: (8073) 
#> Path [44] :Initial log joint density = -432723.509579 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.428e-02   1.010e+03    3.585e-02  3.585e-02      8209 -3.322e+03 -6.987e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3321.662808) evaluations: (8209) 
#> Path [45] :Initial log joint density = -429369.264751 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.718e-02   1.292e+03    4.681e-02  4.681e-02      8349 -3.318e+03 -1.202e+04                   
#> Path [45] :Best Iter: [48] ELBO (-3318.349427) evaluations: (8349) 
#> Path [46] :Initial log joint density = -429313.221164 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      2.776e-02   4.526e+03    1.254e-02  1.254e-02      8522 -3.316e+03 -7.509e+05                   
#> Path [46] :Best Iter: [55] ELBO (-3316.246835) evaluations: (8522) 
#> Path [47] :Initial log joint density = -429249.559302 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      5.159e-02   2.608e+03    3.189e-02  3.189e-02      8367 -3.319e+03 -2.696e+05                   
#> Path [47] :Best Iter: [48] ELBO (-3319.060306) evaluations: (8367) 
#> Path [48] :Initial log joint density = -431590.898810 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      3.846e-02   1.464e+03    2.475e-02  2.475e-02      8071 -3.313e+03 -1.345e+04                   
#> Path [48] :Best Iter: [58] ELBO (-3312.977577) evaluations: (8071) 
#> Path [49] :Initial log joint density = -429264.550487 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      7.804e-02   1.431e+03    3.442e-02  3.442e-02      8202 -3.322e+03 -2.768e+07                   
#> Path [49] :Best Iter: [51] ELBO (-3321.622896) evaluations: (8202) 
#> Path [50] :Initial log joint density = -429906.425362 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.264e+05      4.115e-02   6.752e+02    3.781e-02  3.781e-02      8206 -3.311e+03 -4.209e+04                   
#> Path [50] :Best Iter: [55] ELBO (-3311.161607) evaluations: (8206) 
#> Finished in  27.3 seconds.
#> sccomp says: auto-cleanup removed 2 draw files from 'sccomp_draws_files'
# }
```
