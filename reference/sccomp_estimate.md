# Main Function for SCCOMP Estimate

The `sccomp_estimate` function performs linear modeling on a table of
cell counts or proportions, which includes a cell-group identifier,
sample identifier, abundance (counts or proportions), and factors
(continuous or discrete). The user can define a linear model using an R
formula, where the first factor is the factor of interest.
Alternatively, `sccomp` accepts single-cell data containers (e.g.,
Seurat, SingleCellExperiment, cell metadata, or group-size) and derives
the count data from cell metadata.

## Usage

``` r
sccomp_estimate(
  .data,
  formula_composition = ~1,
  formula_variability = ~1,
  sample,
  cell_group,
  abundance = NULL,
  cores = detectCores(),
  bimodal_mean_variability_association = FALSE,
  percent_false_positive = 5,
  inference_method = "pathfinder",
  prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
  prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0, 0.6),
    standard_deviation = c(10, 20)),
  .sample_cell_group_pairs_to_exclude = NULL,
  output_directory = "sccomp_draws_files",
  verbose = TRUE,
  enable_loo = FALSE,
  noise_model = "multi_beta_binomial",
  exclude_priors = FALSE,
  use_data = TRUE,
  mcmc_seed = sample_seed(),
  max_sampling_iterations = 20000,
  pass_fit = TRUE,
  sig_figs = 9,
  cache_stan_model = sccomp_stan_models_cache_dir,
  cleanup_draw_files = TRUE,
  ...,
  .count = NULL,
  approximate_posterior_inference = NULL,
  variational_inference = NULL,
  .sample = NULL,
  .cell_group = NULL,
  .abundance = NULL
)
```

## Arguments

- .data:

  A tibble including cell_group name column, sample name column,
  abundance column (counts or proportions), and factor columns.

- formula_composition:

  A formula describing the model for differential abundance.

- formula_variability:

  A formula describing the model for differential variability.

- sample:

  A column name as a character string for the sample identifier.
  Replaces the deprecated `.sample`.

- cell_group:

  A column name as a character string for the cell-group identifier.
  Replaces the deprecated `.cell_group`.

- abundance:

  A column name as a character string for the cell-group abundance,
  which can be counts (\> 0) or proportions (between 0 and 1, summing to
  1 across `cell_group`). Replaces the deprecated `.abundance` and
  `.count`.

- cores:

  Number of cores to use for parallel calculations.

- bimodal_mean_variability_association:

  Logical, whether to model mean-variability as bimodal.

- percent_false_positive:

  A real number between 0 and 100 for outlier identification.

- inference_method:

  Character string specifying the inference method to use ('pathfinder',
  'hmc', or 'variational'). Replaces the deprecated
  `approximate_posterior_inference` and `variational_inference`.

- prior_mean:

  A list specifying prior knowledge about the mean distribution,
  including intercept and coefficients.

- prior_overdispersion_mean_association:

  A list specifying prior knowledge about mean/variability association.

- .sample_cell_group_pairs_to_exclude:

  A column name indicating sample/cell-group pairs to exclude.

- output_directory:

  A character string specifying the output directory for Stan draws.

- verbose:

  Logical, whether to print progression details.

- enable_loo:

  Logical, whether to enable model comparison using the LOO package.

- noise_model:

  A character string specifying the noise model (e.g.,
  'multi_beta_binomial').

- exclude_priors:

  Logical, whether to run a prior-free model.

- use_data:

  Logical, whether to run the model data-free.

- mcmc_seed:

  An integer seed for MCMC reproducibility.

- max_sampling_iterations:

  Integer to limit the maximum number of iterations for large datasets.

- pass_fit:

  Logical, whether to include the Stan fit as an attribute in the
  output.

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

- ...:

  Additional arguments passed to the
  [`cmdstanr::sample`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)
  function.

- .count:

  **DEPRECATED**. Use `abundance` instead.

- approximate_posterior_inference:

  **DEPRECATED**. Use `inference_method` instead.

- variational_inference:

  **DEPRECATED**. Use `inference_method` instead.

- .sample:

  **DEPRECATED**. Use `sample` instead.

- .cell_group:

  **DEPRECATED**. Use `cell_group` instead.

- .abundance:

  **DEPRECATED**. Use `abundance` instead.

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

    estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
    )
    
   # Note! 
   # If counts are available, do not use proportion.
   # Using proportion ignores the high uncertainty of low counts
   
   estimate_proportion <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "proportion",
      cores = 1
    )
    
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481752.956065 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.199e-03   1.548e-01    1.000e+00  1.000e+00      3350 -3.703e+03 -3.708e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3703.095196) evaluations: (3350) 
#> Path [2] :Initial log joint density = -481687.292454 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.349e-02   3.234e-01    1.000e+00  1.000e+00      3326 -3.698e+03 -3.706e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3698.183162) evaluations: (3326) 
#> Path [3] :Initial log joint density = -481466.714907 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.856e-03   3.143e-01    5.869e-01  5.869e-01      3262 -3.708e+03 -3.712e+03                   
#> Path [3] :Best Iter: [53] ELBO (-3707.734992) evaluations: (3262) 
#> Path [4] :Initial log joint density = -483811.829603 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.137e-03   2.202e-01    1.000e+00  1.000e+00      3123 -3.707e+03 -3.701e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3700.827387) evaluations: (3123) 
#> Path [5] :Initial log joint density = -482057.866596 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.189e-02   2.014e-01    1.000e+00  1.000e+00      3797 -3.700e+03 -3.704e+03                   
#> Path [5] :Best Iter: [59] ELBO (-3699.817403) evaluations: (3797) 
#> Path [6] :Initial log joint density = -481656.698464 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.546e-02   3.888e-01    1.000e+00  1.000e+00      2991 -3.709e+03 -3.719e+03                   
#> Path [6] :Best Iter: [43] ELBO (-3708.812269) evaluations: (2991) 
#> Path [7] :Initial log joint density = -487951.667798 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.531e-03   2.068e-01    1.000e+00  1.000e+00      3481 -3.709e+03 -3.704e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3703.958306) evaluations: (3481) 
#> Path [8] :Initial log joint density = -481682.098005 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.084e-02   3.279e-01    1.000e+00  1.000e+00      3270 -3.702e+03 -3.705e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3701.894358) evaluations: (3270) 
#> Path [9] :Initial log joint density = -482093.445617 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.674e-03   2.155e-01    1.000e+00  1.000e+00      3405 -3.700e+03 -3.709e+03                   
#> Path [9] :Best Iter: [58] ELBO (-3700.390252) evaluations: (3405) 
#> Path [10] :Initial log joint density = -481751.467677 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.560e-03   3.044e-01    7.491e-01  7.491e-01      3079 -3.708e+03 -3.717e+03                   
#> Path [10] :Best Iter: [39] ELBO (-3707.884450) evaluations: (3079) 
#> Path [11] :Initial log joint density = -481921.748885 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.869e-03   1.884e-01    4.860e-01  1.000e+00      3488 -3.701e+03 -3.711e+03                   
#> Path [11] :Best Iter: [58] ELBO (-3700.829818) evaluations: (3488) 
#> Path [12] :Initial log joint density = -483816.265680 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.028e-03   2.469e-01    1.000e+00  1.000e+00      3392 -3.707e+03 -3.707e+03                   
#> Path [12] :Best Iter: [54] ELBO (-3706.935020) evaluations: (3392) 
#> Path [13] :Initial log joint density = -481995.174796 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.383e-03   2.841e-01    7.336e-01  7.336e-01      3457 -3.698e+03 -3.709e+03                   
#> Path [13] :Best Iter: [57] ELBO (-3697.896705) evaluations: (3457) 
#> Path [14] :Initial log joint density = -481548.018394 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.690e-03   3.305e-01    7.378e-01  7.378e-01      2910 -3.709e+03 -3.725e+03                   
#> Path [14] :Best Iter: [47] ELBO (-3708.547600) evaluations: (2910) 
#> Path [15] :Initial log joint density = -483588.691777 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.996e-03   2.743e-01    7.951e-01  7.951e-01      3016 -3.707e+03 -3.709e+03                   
#> Path [15] :Best Iter: [43] ELBO (-3706.542751) evaluations: (3016) 
#> Path [16] :Initial log joint density = -481336.972558 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.137e-03   3.041e-01    6.777e-01  6.777e-01      2789 -3.709e+03 -3.716e+03                   
#> Path [16] :Best Iter: [50] ELBO (-3708.810959) evaluations: (2789) 
#> Path [17] :Initial log joint density = -484267.219000 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.251e-02   3.191e-01    1.000e+00  1.000e+00      3417 -3.700e+03 -3.705e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3700.469463) evaluations: (3417) 
#> Path [18] :Initial log joint density = -485077.588983 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.073e-02   2.830e-01    1.000e+00  1.000e+00      3386 -3.701e+03 -3.703e+03                   
#> Path [18] :Best Iter: [58] ELBO (-3700.745645) evaluations: (3386) 
#> Path [19] :Initial log joint density = -481723.035880 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.734e-02   3.375e-01    1.000e+00  1.000e+00      3612 -3.700e+03 -3.705e+03                   
#> Path [19] :Best Iter: [59] ELBO (-3700.387664) evaluations: (3612) 
#> Path [20] :Initial log joint density = -481704.920175 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.084e-02   2.961e-01    4.550e-01  1.000e+00      3099 -3.705e+03 -3.711e+03                   
#> Path [20] :Best Iter: [48] ELBO (-3704.987147) evaluations: (3099) 
#> Path [21] :Initial log joint density = -481569.658474 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.874e-02   2.877e-01    9.174e-01  9.174e-01      3451 -3.702e+03 -3.704e+03                   
#> Path [21] :Best Iter: [59] ELBO (-3702.179177) evaluations: (3451) 
#> Path [22] :Initial log joint density = -481629.409149 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.570e-02   2.506e-01    1.000e+00  1.000e+00      3595 -3.699e+03 -3.702e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3699.435902) evaluations: (3595) 
#> Path [23] :Initial log joint density = -483420.656945 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.455e-03   2.326e-01    1.000e+00  1.000e+00      3269 -3.701e+03 -3.702e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3700.851420) evaluations: (3269) 
#> Path [24] :Initial log joint density = -481758.653202 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.017e-03   2.189e-01    1.000e+00  1.000e+00      3270 -3.701e+03 -3.704e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3700.848014) evaluations: (3270) 
#> Path [25] :Initial log joint density = -481811.580663 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.933e-03   2.391e-01    1.000e+00  1.000e+00      3296 -3.706e+03 -3.706e+03                   
#> Path [25] :Best Iter: [54] ELBO (-3706.004403) evaluations: (3296) 
#> Path [26] :Initial log joint density = -481892.975849 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.443e-03   2.472e-01    9.797e-01  9.797e-01      2905 -3.704e+03 -3.717e+03                   
#> Path [26] :Best Iter: [51] ELBO (-3703.700661) evaluations: (2905) 
#> Path [27] :Initial log joint density = -481921.234871 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.487e-03   2.578e-01    1.000e+00  1.000e+00      2754 -3.706e+03 -3.718e+03                   
#> Path [27] :Best Iter: [48] ELBO (-3705.894802) evaluations: (2754) 
#> Path [28] :Initial log joint density = -481569.664149 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.116e-03   1.995e-01    9.298e-01  9.298e-01      3466 -3.696e+03 -3.705e+03                   
#> Path [28] :Best Iter: [57] ELBO (-3696.456038) evaluations: (3466) 
#> Path [29] :Initial log joint density = -481879.722822 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      3.897e-03   2.149e-01    8.343e-01  8.343e-01      3136 -3.701e+03 -3.711e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3700.638184) evaluations: (3136) 
#> Path [30] :Initial log joint density = -481471.177111 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.231e-03   3.326e-01    5.363e-01  5.363e-01      2991 -3.709e+03 -3.732e+03                   
#> Path [30] :Best Iter: [42] ELBO (-3708.563162) evaluations: (2991) 
#> Path [31] :Initial log joint density = -483937.438520 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.233e-02   3.342e-01    4.680e-01  1.000e+00      3062 -3.706e+03 -3.709e+03                   
#> Path [31] :Best Iter: [49] ELBO (-3706.480702) evaluations: (3062) 
#> Path [32] :Initial log joint density = -481824.320620 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.171e-02   2.385e-01    1.000e+00  1.000e+00      3491 -3.700e+03 -3.707e+03                   
#> Path [32] :Best Iter: [58] ELBO (-3699.984956) evaluations: (3491) 
#> Path [33] :Initial log joint density = -481528.457661 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.317e-02   2.479e-01    1.000e+00  1.000e+00      2906 -3.708e+03 -3.712e+03                   
#> Path [33] :Best Iter: [39] ELBO (-3707.721575) evaluations: (2906) 
#> Path [34] :Initial log joint density = -481613.328858 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.765e-03   1.607e-01    5.080e-01  1.000e+00      3274 -3.702e+03 -3.712e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3702.282800) evaluations: (3274) 
#> Path [35] :Initial log joint density = -482073.488758 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      3.912e-03   2.474e-01    6.966e-01  6.966e-01      3192 -3.703e+03 -3.711e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3703.153629) evaluations: (3192) 
#> Path [36] :Initial log joint density = -482288.125204 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.720e-03   3.072e-01    1.000e+00  1.000e+00      3025 -3.709e+03 -3.703e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3703.276959) evaluations: (3025) 
#> Path [37] :Initial log joint density = -482511.888327 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.935e-03   1.812e-01    8.389e-01  8.389e-01      3396 -3.701e+03 -3.712e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3700.922313) evaluations: (3396) 
#> Path [38] :Initial log joint density = -481713.261929 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.655e-03   1.942e-01    1.000e+00  1.000e+00      2970 -3.708e+03 -3.721e+03                   
#> Path [38] :Best Iter: [50] ELBO (-3708.470435) evaluations: (2970) 
#> Path [39] :Initial log joint density = -482529.361542 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.007e-03   2.550e-01    7.902e-01  7.902e-01      3085 -3.706e+03 -3.721e+03                   
#> Path [39] :Best Iter: [52] ELBO (-3705.993776) evaluations: (3085) 
#> Path [40] :Initial log joint density = -481375.765488 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.765e-03   2.102e-01    8.273e-01  8.273e-01      3220 -3.702e+03 -3.714e+03                   
#> Path [40] :Best Iter: [55] ELBO (-3701.903709) evaluations: (3220) 
#> Path [41] :Initial log joint density = -481681.098004 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.874e-03   2.194e-01    1.000e+00  1.000e+00      3492 -3.700e+03 -3.702e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3700.387863) evaluations: (3492) 
#> Path [42] :Initial log joint density = -481384.900643 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.976e-03   1.678e-01    1.000e+00  1.000e+00      2990 -3.706e+03 -3.717e+03                   
#> Path [42] :Best Iter: [51] ELBO (-3706.038734) evaluations: (2990) 
#> Path [43] :Initial log joint density = -481949.425692 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.025e-03   2.216e-01    1.000e+00  1.000e+00      3109 -3.710e+03 -3.704e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3704.001584) evaluations: (3109) 
#> Path [44] :Initial log joint density = -481758.182971 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.071e-02   2.647e-01    1.000e+00  1.000e+00      3223 -3.698e+03 -3.703e+03                   
#> Path [44] :Best Iter: [56] ELBO (-3698.023444) evaluations: (3223) 
#> Path [45] :Initial log joint density = -481663.902959 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.338e-03   1.815e-01    1.000e+00  1.000e+00      3044 -3.710e+03 -3.714e+03                   
#> Path [45] :Best Iter: [50] ELBO (-3710.279677) evaluations: (3044) 
#> Path [46] :Initial log joint density = -481756.605441 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.040e-03   2.411e-01    1.000e+00  1.000e+00      3419 -3.698e+03 -3.710e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3698.419379) evaluations: (3419) 
#> Path [47] :Initial log joint density = -484995.239271 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.879e-03   2.007e-01    1.000e+00  1.000e+00      3148 -3.706e+03 -3.700e+03                   
#> Path [47] :Best Iter: [56] ELBO (-3700.342988) evaluations: (3148) 
#> Path [48] :Initial log joint density = -481640.357646 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.272e-02   3.244e-01    1.000e+00  1.000e+00      2909 -3.708e+03 -3.717e+03                   
#> Path [48] :Best Iter: [52] ELBO (-3707.556156) evaluations: (2909) 
#> Path [49] :Initial log joint density = -481592.058339 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.181e-03   1.525e-01    7.681e-01  7.681e-01      3305 -3.704e+03 -3.715e+03                   
#> Path [49] :Best Iter: [57] ELBO (-3704.063688) evaluations: (3305) 
#> Path [50] :Initial log joint density = -482087.454580 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.038e-03   1.261e-01    1.000e+00  1.000e+00      3737 -3.700e+03 -3.710e+03                   
#> Path [50] :Best Iter: [58] ELBO (-3700.314195) evaluations: (3737) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: proportion column is a proportion. The sum-constrained beta model will be used. When possible using counts is preferred as the binomial noise component is often dominating for rare groups (e.g. rare cell types).
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Warning: sccomp says: your proportion values include 0. Assuming that 0s derive from a precision threshold (e.g. deconvolution), 0s are converted to the smaller non 0 proportion value.
#> Path [1] :Initial log joint density = -282.340164 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.543e+03      3.122e-02   2.438e+03    2.387e-02  2.387e-02      7932  2.279e+03 -5.270e+04                   
#> Path [1] :Best Iter: [34] ELBO (2279.332813) evaluations: (7932) 
#> Path [2] :Initial log joint density = -1239.949493 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.530e+03      3.136e-02   2.503e+03    1.939e-02  1.939e-02      9001  2.284e+03 -3.199e+04                   
#> Path [2] :Best Iter: [34] ELBO (2283.623412) evaluations: (9001) 
#> Path [3] :Initial log joint density = -476.016082 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.552e+03      7.224e-02   3.786e+03    3.387e-02  3.387e-02      8410  2.283e+03 -3.124e+04                   
#> Path [3] :Best Iter: [31] ELBO (2283.029464) evaluations: (8410) 
#> Path [4] :Initial log joint density = -548.181501 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.561e+03      1.315e-01   7.775e+03    3.878e-02  3.878e-02      8945  2.280e+03 -4.275e+05                   
#> Path [4] :Best Iter: [27] ELBO (2279.674235) evaluations: (8945) 
#> Path [5] :Initial log joint density = -1145.278184 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.545e+03      3.196e-02   2.305e+03    2.013e-02  2.013e-02      8381  2.275e+03 -7.915e+03                   
#> Path [5] :Best Iter: [38] ELBO (2274.664051) evaluations: (8381) 
#> Path [6] :Initial log joint density = -4661.056874 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.546e+03      1.324e-01   5.624e+03    4.727e-02  4.727e-02      8547  2.276e+03 -3.757e+04                   
#> Path [6] :Best Iter: [34] ELBO (2275.564047) evaluations: (8547) 
#> Path [7] :Initial log joint density = -628.575689 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.548e+03      2.661e-02   4.275e+03    2.207e-02  2.207e-02      8513  2.280e+03 -4.027e+04                   
#> Path [7] :Best Iter: [35] ELBO (2279.812530) evaluations: (8513) 
#> Path [8] :Initial log joint density = -604.277225 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.540e+03      4.785e-02   2.407e+03    1.323e-02  3.675e-02      8335  2.282e+03  8.004e+02                   
#> Path [8] :Best Iter: [32] ELBO (2281.565705) evaluations: (8335) 
#> Path [9] :Initial log joint density = -555.645514 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      4.875e-02   5.572e+03    4.052e-02  4.052e-02      8640  2.283e+03 -2.662e+04                   
#> Path [9] :Best Iter: [25] ELBO (2283.390071) evaluations: (8640) 
#> Path [10] :Initial log joint density = -1856.679996 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.546e+03      4.659e-02   4.860e+03    1.512e-02  1.512e-02      8551  2.286e+03 -2.725e+05                   
#> Path [10] :Best Iter: [36] ELBO (2286.140269) evaluations: (8551) 
#> Path [11] :Initial log joint density = -859.898507 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.541e+03      3.240e-02   3.824e+03    4.845e-02  4.845e-02      8456  2.282e+03 -1.691e+03                   
#> Path [11] :Best Iter: [34] ELBO (2281.656241) evaluations: (8456) 
#> Path [12] :Initial log joint density = -120.040639 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.572e+03      4.321e-02   1.130e+04    1.999e-02  1.999e-02      9143  2.279e+03 -7.221e+05                   
#> Path [12] :Best Iter: [30] ELBO (2279.016064) evaluations: (9143) 
#> Path [13] :Initial log joint density = -2612.748510 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.550e+03      4.781e-02   3.053e+03    1.906e-02  1.906e-02      8595  2.283e+03 -2.266e+05                   
#> Path [13] :Best Iter: [36] ELBO (2283.349388) evaluations: (8595) 
#> Path [14] :Initial log joint density = -1937.275896 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.547e+03      2.841e-02   7.552e+03    5.582e-02  5.582e-02      8300  2.278e+03  9.950e+02                   
#> Path [14] :Best Iter: [32] ELBO (2278.371178) evaluations: (8300) 
#> Path [15] :Initial log joint density = -589.415891 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.555e+03      3.456e-02   2.968e+03    4.087e-02  4.087e-02      8616  2.278e+03 -4.840e+03                   
#> Path [15] :Best Iter: [26] ELBO (2278.286309) evaluations: (8616) 
#> Path [16] :Initial log joint density = -631.602320 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      1.217e-01   8.545e+03    3.292e-02  3.292e-02      8961  2.280e+03 -2.437e+07                   
#> Path [16] :Best Iter: [33] ELBO (2279.732295) evaluations: (8961) 
#> Path [17] :Initial log joint density = -487.117403 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      1.513e-01   1.070e+04    3.613e-02  1.101e-01      8877  2.279e+03 -1.169e+06                   
#> Path [17] :Best Iter: [27] ELBO (2278.761430) evaluations: (8877) 
#> Path [18] :Initial log joint density = -1780.257556 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      7.207e-02   4.100e+03    2.581e-02  2.581e-02      8557  2.281e+03 -3.165e+07                   
#> Path [18] :Best Iter: [30] ELBO (2281.102874) evaluations: (8557) 
#> Path [19] :Initial log joint density = -730.377254 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      8.858e-02   3.387e+03    2.635e-02  5.438e-02      8747  2.279e+03 -1.751e+04                   
#> Path [19] :Best Iter: [29] ELBO (2278.755906) evaluations: (8747) 
#> Path [20] :Initial log joint density = -482.431394 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.544e+03      2.383e-02   3.592e+03    2.971e-02  2.971e-02      8271  2.282e+03 -1.001e+03                   
#> Path [20] :Best Iter: [30] ELBO (2281.715422) evaluations: (8271) 
#> Path [21] :Initial log joint density = -220.736133 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      5.733e-02   9.462e+03    2.666e-02  2.666e-02      8701  2.278e+03 -1.072e+06                   
#> Path [21] :Best Iter: [32] ELBO (2277.885606) evaluations: (8701) 
#> Path [22] :Initial log joint density = -507.105489 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.549e+03      6.224e-02   4.979e+03    3.085e-02  3.085e-02      8799  2.282e+03 -7.937e+05                   
#> Path [22] :Best Iter: [36] ELBO (2282.418688) evaluations: (8799) 
#> Path [23] :Initial log joint density = -510.126646 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.562e+03      6.477e-02   1.050e+04    1.823e-02  1.823e-02      8672  2.276e+03 -3.631e+06                   
#> Path [23] :Best Iter: [30] ELBO (2275.576373) evaluations: (8672) 
#> Path [24] :Initial log joint density = -514.990090 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      3.498e-02   3.516e+03    5.350e-02  5.350e-02      8628  2.281e+03 -1.271e+03                   
#> Path [24] :Best Iter: [29] ELBO (2280.548126) evaluations: (8628) 
#> Path [25] :Initial log joint density = -533.682894 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.567e+03      4.517e-02   1.669e+04    2.544e-02  2.544e-02      8597  2.278e+03 -9.337e+05                   
#> Path [25] :Best Iter: [34] ELBO (2278.434335) evaluations: (8597) 
#> Path [26] :Initial log joint density = -693.089247 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.567e+03      4.500e-02   1.025e+04    1.697e-02  1.697e-02      8792  2.280e+03 -1.011e+04                   
#> Path [26] :Best Iter: [25] ELBO (2279.779427) evaluations: (8792) 
#> Path [27] :Initial log joint density = -848.101108 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.548e+03      8.616e-02   3.291e+03    3.101e-02  3.101e-02      8407  2.281e+03  1.091e+03                   
#> Path [27] :Best Iter: [37] ELBO (2281.138889) evaluations: (8407) 
#> Path [28] :Initial log joint density = -1564.995320 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.541e+03      1.102e-01   2.100e+03    2.912e-02  2.912e-02      8661  2.280e+03 -1.441e+08                   
#> Path [28] :Best Iter: [30] ELBO (2279.658598) evaluations: (8661) 
#> Path [29] :Initial log joint density = -93.844030 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.571e+03      7.331e-02   8.262e+03    3.338e-02  3.338e-02      8254  2.284e+03 -3.227e+04                   
#> Path [29] :Best Iter: [29] ELBO (2283.515007) evaluations: (8254) 
#> Path [30] :Initial log joint density = -322.891211 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.559e+03      6.754e-02   8.021e+03    2.798e-02  2.798e-02      8390  2.283e+03 -3.310e+04                   
#> Path [30] :Best Iter: [37] ELBO (2283.113994) evaluations: (8390) 
#> Path [31] :Initial log joint density = -566.195822 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      3.979e-02   7.521e+03    2.845e-02  2.845e-02      8518  2.282e+03 -1.991e+04                   
#> Path [31] :Best Iter: [29] ELBO (2281.954121) evaluations: (8518) 
#> Path [32] :Initial log joint density = -401.129778 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.538e+03      4.664e-02   2.393e+03    5.067e-02  5.067e-02      9034  2.277e+03 -2.449e+03                   
#> Path [32] :Best Iter: [27] ELBO (2277.145679) evaluations: (9034) 
#> Path [33] :Initial log joint density = -1365.804682 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.546e+03      2.829e-02   4.724e+03    2.483e-02  2.483e-02      9048  2.283e+03 -4.397e+04                   
#> Path [33] :Best Iter: [36] ELBO (2283.023152) evaluations: (9048) 
#> Path [34] :Initial log joint density = -654.929040 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.558e+03      4.744e-02   9.147e+03    2.481e-02  5.835e-02      8821  2.278e+03 -4.962e+06                   
#> Path [34] :Best Iter: [29] ELBO (2278.328609) evaluations: (8821) 
#> Path [35] :Initial log joint density = -279.441597 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      5.755e-02   1.074e+04    3.055e-02  3.055e-02      8667  2.283e+03 -2.763e+05                   
#> Path [35] :Best Iter: [29] ELBO (2282.633919) evaluations: (8667) 
#> Path [36] :Initial log joint density = -421.786197 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.574e+03      3.900e-02   6.733e+03    1.029e-02  3.045e-02      8479  2.283e+03 -2.130e+05                   
#> Path [36] :Best Iter: [29] ELBO (2283.396039) evaluations: (8479) 
#> Path [37] :Initial log joint density = -775.585550 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.537e+03      1.102e-02   3.722e+03    2.268e-02  2.268e-02      8425  2.276e+03  5.017e+02                   
#> Path [37] :Best Iter: [34] ELBO (2275.655071) evaluations: (8425) 
#> Path [38] :Initial log joint density = -420.656977 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.564e+03      3.948e-02   6.262e+03    1.180e-02  2.361e-02      8828  2.283e+03 -6.110e+05                   
#> Path [38] :Best Iter: [31] ELBO (2282.671392) evaluations: (8828) 
#> Path [39] :Initial log joint density = -585.469236 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.550e+03      5.453e-02   4.619e+03    2.768e-02  2.768e-02      8528  2.281e+03 -1.690e+06                   
#> Path [39] :Best Iter: [32] ELBO (2280.520347) evaluations: (8528) 
#> Path [40] :Initial log joint density = -292.023761 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.543e+03      5.257e-02   4.290e+03    2.528e-02  2.528e-02      8682  2.279e+03 -3.311e+03                   
#> Path [40] :Best Iter: [32] ELBO (2278.885536) evaluations: (8682) 
#> Path [41] :Initial log joint density = -752.816932 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.565e+03      7.154e-02   8.488e+03    1.304e-02  4.934e-02      8735  2.281e+03 -4.870e+04                   
#> Path [41] :Best Iter: [28] ELBO (2280.868159) evaluations: (8735) 
#> Path [42] :Initial log joint density = -331.309214 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.560e+03      5.135e-02   1.033e+04    3.292e-02  3.292e-02      8608  2.282e+03 -3.522e+04                   
#> Path [42] :Best Iter: [29] ELBO (2282.189515) evaluations: (8608) 
#> Path [43] :Initial log joint density = -680.815253 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.548e+03      3.841e-02   4.999e+03    2.922e-02  2.922e-02      8731  2.280e+03 -6.525e+03                   
#> Path [43] :Best Iter: [29] ELBO (2280.220362) evaluations: (8731) 
#> Path [44] :Initial log joint density = -458.522086 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.548e+03      2.133e-02   4.095e+03    1.516e-02  3.090e-02      8887  2.281e+03 -3.888e+03                   
#> Path [44] :Best Iter: [34] ELBO (2281.403283) evaluations: (8887) 
#> Path [45] :Initial log joint density = -3155.304473 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.542e+03      6.538e-02   3.592e+03    2.987e-02  2.987e-02      8602  2.281e+03 -9.926e+04                   
#> Path [45] :Best Iter: [38] ELBO (2280.816327) evaluations: (8602) 
#> Path [46] :Initial log joint density = -7396.985493 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.543e+03      5.429e-02   7.821e+03    1.160e-02  1.160e-02      9017  2.285e+03 -6.825e+05                   
#> Path [46] :Best Iter: [36] ELBO (2284.589223) evaluations: (9017) 
#> Path [47] :Initial log joint density = -647.844981 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.552e+03      6.876e-02   3.997e+03    3.086e-02  3.086e-02      8635  2.278e+03 -6.738e+05                   
#> Path [47] :Best Iter: [29] ELBO (2278.430191) evaluations: (8635) 
#> Path [48] :Initial log joint density = -912.671957 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.540e+03      1.395e-02   2.865e+03    1.435e-02  1.435e-02      8161  2.280e+03 -2.047e+03                   
#> Path [48] :Best Iter: [34] ELBO (2279.515821) evaluations: (8161) 
#> Path [49] :Initial log joint density = -898.601656 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.547e+03      1.118e-01   5.865e+03    2.434e-02  4.053e-02      8604  2.279e+03 -6.385e+11                   
#> Path [49] :Best Iter: [39] ELBO (2279.093673) evaluations: (8604) 
#> Path [50] :Initial log joint density = -578.027456 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.563e+03      4.102e-02   3.777e+03    4.283e-02  4.283e-02      8905  2.281e+03 -7.426e+04                   
#> Path [50] :Best Iter: [25] ELBO (2280.646391) evaluations: (8905) 
#> Finished in  16.0 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
