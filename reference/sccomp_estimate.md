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
#> Path [1] :Initial log joint density = -481654.381653 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.167e-03   2.029e-01    6.600e-01  6.600e-01      3424 -3.697e+03 -3.716e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3697.434692) evaluations: (3424) 
#> Path [2] :Initial log joint density = -481500.229411 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.412e-03   2.123e-01    1.000e+00  1.000e+00      2971 -3.706e+03 -3.707e+03                   
#> Path [2] :Best Iter: [39] ELBO (-3705.954588) evaluations: (2971) 
#> Path [3] :Initial log joint density = -481822.102920 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.304e-02   3.860e-01    1.000e+00  1.000e+00      3916 -3.697e+03 -3.708e+03                   
#> Path [3] :Best Iter: [63] ELBO (-3697.272401) evaluations: (3916) 
#> Path [4] :Initial log joint density = -481551.788506 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      8.631e-03   2.809e-01    9.750e-01  9.750e-01      2716 -3.707e+03 -3.719e+03                   
#> Path [4] :Best Iter: [42] ELBO (-3707.097311) evaluations: (2716) 
#> Path [5] :Initial log joint density = -484419.961750 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.471e-03   2.344e-01    9.691e-01  9.691e-01      3279 -3.700e+03 -3.711e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3700.030866) evaluations: (3279) 
#> Path [6] :Initial log joint density = -481631.827885 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.808e-03   2.497e-01    1.000e+00  1.000e+00      2866 -3.709e+03 -3.709e+03                   
#> Path [6] :Best Iter: [45] ELBO (-3708.540375) evaluations: (2866) 
#> Path [7] :Initial log joint density = -481710.030492 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.104e-02   3.024e-01    1.000e+00  1.000e+00      3539 -3.697e+03 -3.701e+03                   
#> Path [7] :Best Iter: [60] ELBO (-3697.158723) evaluations: (3539) 
#> Path [8] :Initial log joint density = -481330.877215 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.186e-02   3.337e-01    1.000e+00  1.000e+00      3147 -3.699e+03 -3.710e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3699.348605) evaluations: (3147) 
#> Path [9] :Initial log joint density = -481611.041983 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.365e-02   2.373e-01    1.000e+00  1.000e+00      3006 -3.707e+03 -3.716e+03                   
#> Path [9] :Best Iter: [49] ELBO (-3707.409962) evaluations: (3006) 
#> Path [10] :Initial log joint density = -481745.971644 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.569e-03   1.859e-01    1.000e+00  1.000e+00      3566 -3.700e+03 -3.701e+03                   
#> Path [10] :Best Iter: [58] ELBO (-3699.962136) evaluations: (3566) 
#> Path [11] :Initial log joint density = -481872.465900 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.808e-03   2.991e-01    1.000e+00  1.000e+00      3019 -3.708e+03 -3.722e+03                   
#> Path [11] :Best Iter: [42] ELBO (-3707.865249) evaluations: (3019) 
#> Path [12] :Initial log joint density = -481559.476654 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.563e-03   2.112e-01    7.010e-01  7.010e-01      3364 -3.699e+03 -3.714e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3698.876768) evaluations: (3364) 
#> Path [13] :Initial log joint density = -481696.894307 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.038e-02   1.895e-01    1.000e+00  1.000e+00      3810 -3.700e+03 -3.707e+03                   
#> Path [13] :Best Iter: [61] ELBO (-3700.025273) evaluations: (3810) 
#> Path [14] :Initial log joint density = -483510.973696 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.308e-02   2.432e-01    1.000e+00  1.000e+00      3301 -3.699e+03 -3.705e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3698.823192) evaluations: (3301) 
#> Path [15] :Initial log joint density = -481957.307460 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.788e+05      6.586e-03   1.940e-01    5.008e-01  1.000e+00      4142 -3.702e+03 -3.711e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3701.545769) evaluations: (4142) 
#> Path [16] :Initial log joint density = -482282.152351 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.781e-03   2.929e-01    1.000e+00  1.000e+00      3503 -3.702e+03 -3.701e+03                   
#> Path [16] :Best Iter: [60] ELBO (-3701.203964) evaluations: (3503) 
#> Path [17] :Initial log joint density = -482230.918879 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.252e-02   2.596e-01    1.000e+00  1.000e+00      3535 -3.701e+03 -3.701e+03                   
#> Path [17] :Best Iter: [59] ELBO (-3700.816082) evaluations: (3535) 
#> Path [18] :Initial log joint density = -482065.406703 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.765e-03   1.403e-01    1.000e+00  1.000e+00      3192 -3.706e+03 -3.710e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3705.670744) evaluations: (3192) 
#> Path [19] :Initial log joint density = -481648.673557 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.943e-03   2.073e-01    8.565e-01  8.565e-01      3427 -3.702e+03 -3.716e+03                   
#> Path [19] :Best Iter: [57] ELBO (-3702.467771) evaluations: (3427) 
#> Path [20] :Initial log joint density = -481614.512691 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.290e-02   3.085e-01    1.000e+00  1.000e+00      2994 -3.708e+03 -3.706e+03                   
#> Path [20] :Best Iter: [54] ELBO (-3706.201107) evaluations: (2994) 
#> Path [21] :Initial log joint density = -481743.073235 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.344e-03   1.823e-01    9.147e-01  9.147e-01      3621 -3.698e+03 -3.710e+03                   
#> Path [21] :Best Iter: [59] ELBO (-3698.046204) evaluations: (3621) 
#> Path [22] :Initial log joint density = -481521.666538 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.577e-03   2.271e-01    1.000e+00  1.000e+00      3312 -3.702e+03 -3.701e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3700.529150) evaluations: (3312) 
#> Path [23] :Initial log joint density = -481608.246511 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.073e-03   2.285e-01    1.000e+00  1.000e+00      3273 -3.701e+03 -3.702e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3701.494855) evaluations: (3273) 
#> Path [24] :Initial log joint density = -481590.757812 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.357e-03   3.162e-01    1.000e+00  1.000e+00      3109 -3.699e+03 -3.711e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3698.947510) evaluations: (3109) 
#> Path [25] :Initial log joint density = -481883.127110 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.923e-03   2.298e-01    1.000e+00  1.000e+00      3265 -3.702e+03 -3.701e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3701.212455) evaluations: (3265) 
#> Path [26] :Initial log joint density = -481610.802137 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      2.457e-03   1.994e-01    6.760e-01  6.760e-01      2923 -3.706e+03 -3.721e+03                   
#> Path [26] :Best Iter: [39] ELBO (-3706.333994) evaluations: (2923) 
#> Path [27] :Initial log joint density = -483043.100237 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.976e-03   2.161e-01    7.877e-01  7.877e-01      3162 -3.708e+03 -3.708e+03                   
#> Path [27] :Best Iter: [52] ELBO (-3707.686934) evaluations: (3162) 
#> Path [28] :Initial log joint density = -481568.764842 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.104e-03   1.708e-01    9.489e-01  9.489e-01      3270 -3.701e+03 -3.709e+03                   
#> Path [28] :Best Iter: [56] ELBO (-3700.609982) evaluations: (3270) 
#> Path [29] :Initial log joint density = -481709.139730 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.489e-03   2.384e-01    6.879e-01  6.879e-01      3115 -3.710e+03 -3.716e+03                   
#> Path [29] :Best Iter: [49] ELBO (-3710.466699) evaluations: (3115) 
#> Path [30] :Initial log joint density = -483078.070446 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.035e-02   2.263e-01    1.000e+00  1.000e+00      3049 -3.709e+03 -3.709e+03                   
#> Path [30] :Best Iter: [51] ELBO (-3708.530262) evaluations: (3049) 
#> Path [31] :Initial log joint density = -481626.358297 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.094e-02   2.166e-01    1.000e+00  1.000e+00      3103 -3.701e+03 -3.699e+03                   
#> Path [31] :Best Iter: [56] ELBO (-3699.436563) evaluations: (3103) 
#> Path [32] :Initial log joint density = -481527.834220 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.462e-03   2.762e-01    8.202e-01  8.202e-01      2876 -3.710e+03 -3.726e+03                   
#> Path [32] :Best Iter: [51] ELBO (-3709.758603) evaluations: (2876) 
#> Path [33] :Initial log joint density = -481526.836811 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.923e-03   1.674e-01    8.020e-01  8.020e-01      3292 -3.701e+03 -3.712e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3700.516092) evaluations: (3292) 
#> Path [34] :Initial log joint density = -481442.115160 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.005e-02   2.795e-01    1.000e+00  1.000e+00      3083 -3.704e+03 -3.712e+03                   
#> Path [34] :Best Iter: [52] ELBO (-3704.168963) evaluations: (3083) 
#> Path [35] :Initial log joint density = -482817.992860 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.217e-03   1.964e-01    1.000e+00  1.000e+00      3383 -3.699e+03 -3.706e+03                   
#> Path [35] :Best Iter: [56] ELBO (-3699.012032) evaluations: (3383) 
#> Path [36] :Initial log joint density = -482371.212869 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.694e-03   1.702e-01    1.000e+00  1.000e+00      3463 -3.706e+03 -3.702e+03                   
#> Path [36] :Best Iter: [59] ELBO (-3702.118871) evaluations: (3463) 
#> Path [37] :Initial log joint density = -481669.243196 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.186e-03   2.434e-01    8.575e-01  8.575e-01      3382 -3.702e+03 -3.708e+03                   
#> Path [37] :Best Iter: [57] ELBO (-3702.403673) evaluations: (3382) 
#> Path [38] :Initial log joint density = -481726.485038 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.826e-03   1.837e-01    1.000e+00  1.000e+00      3478 -3.701e+03 -3.706e+03                   
#> Path [38] :Best Iter: [58] ELBO (-3700.519387) evaluations: (3478) 
#> Path [39] :Initial log joint density = -486797.702276 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.542e-02   3.756e-01    1.000e+00  1.000e+00      3504 -3.700e+03 -3.707e+03                   
#> Path [39] :Best Iter: [59] ELBO (-3700.368903) evaluations: (3504) 
#> Path [40] :Initial log joint density = -482343.229078 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.091e-02   2.573e-01    1.000e+00  1.000e+00      3428 -3.702e+03 -3.704e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3702.001015) evaluations: (3428) 
#> Path [41] :Initial log joint density = -482326.036142 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.060e-02   3.173e-01    4.180e-01  1.000e+00      3465 -3.700e+03 -3.708e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3700.247131) evaluations: (3465) 
#> Path [42] :Initial log joint density = -484233.788559 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.209e-02   3.638e-01    1.000e+00  1.000e+00      3519 -3.702e+03 -3.710e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3702.496905) evaluations: (3519) 
#> Path [43] :Initial log joint density = -481499.079984 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.917e-03   2.701e-01    1.000e+00  1.000e+00      2921 -3.709e+03 -3.722e+03                   
#> Path [43] :Best Iter: [48] ELBO (-3708.901039) evaluations: (2921) 
#> Path [44] :Initial log joint density = -481690.725613 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.690e-02   2.872e-01    1.000e+00  1.000e+00      2791 -3.708e+03 -3.712e+03                   
#> Path [44] :Best Iter: [38] ELBO (-3708.309108) evaluations: (2791) 
#> Path [45] :Initial log joint density = -481978.836665 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.012e-02   3.314e-01    1.000e+00  1.000e+00      3356 -3.700e+03 -3.709e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3699.877352) evaluations: (3356) 
#> Path [46] :Initial log joint density = -482134.522441 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.966e-03   1.688e-01    1.000e+00  1.000e+00      3243 -3.701e+03 -3.702e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3700.875823) evaluations: (3243) 
#> Path [47] :Initial log joint density = -481964.238753 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      9.695e-03   2.962e-01    9.379e-01  9.379e-01      3655 -3.701e+03 -3.709e+03                   
#> Path [47] :Best Iter: [61] ELBO (-3700.571289) evaluations: (3655) 
#> Path [48] :Initial log joint density = -482907.260269 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.646e-03   2.683e-01    8.935e-01  8.935e-01      3277 -3.698e+03 -3.713e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3698.324182) evaluations: (3277) 
#> Path [49] :Initial log joint density = -481811.295994 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.082e-02   2.068e-01    1.000e+00  1.000e+00      3161 -3.702e+03 -3.700e+03                   
#> Path [49] :Best Iter: [56] ELBO (-3699.722007) evaluations: (3161) 
#> Path [50] :Initial log joint density = -482767.869635 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.034e-02   1.957e-01    1.000e+00  1.000e+00      3025 -3.709e+03 -3.709e+03                   
#> Path [50] :Best Iter: [49] ELBO (-3708.691571) evaluations: (3025) 
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
#> Path [1] :Initial log joint density = -902.968750 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.547e+03      5.469e-03   4.581e+03    4.295e-03  4.295e-03      8579  2.282e+03 -1.504e+04                   
#> Path [1] :Best Iter: [41] ELBO (2281.755314) evaluations: (8579) 
#> Path [2] :Initial log joint density = -1161.963117 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.544e+03      1.790e-02   3.196e+03    1.307e-02  1.307e-02      8383  2.282e+03 -1.199e+04                   
#> Path [2] :Best Iter: [33] ELBO (2281.917961) evaluations: (8383) 
#> Path [3] :Initial log joint density = -482.168479 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.545e+03      2.372e-02   4.608e+03    2.556e-02  2.556e-02      8742  2.282e+03 -8.034e+04                   
#> Path [3] :Best Iter: [30] ELBO (2281.982669) evaluations: (8742) 
#> Path [4] :Initial log joint density = -459.558926 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.555e+03      1.444e-01   7.699e+03    4.006e-02  4.006e-02      8779  2.278e+03 -1.391e+04                   
#> Path [4] :Best Iter: [30] ELBO (2278.229108) evaluations: (8779) 
#> Path [5] :Initial log joint density = -383.381341 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.549e+03      6.764e-03   3.918e+03    2.849e-02  2.849e-02      9165  2.282e+03 -6.903e+03                   
#> Path [5] :Best Iter: [26] ELBO (2282.000385) evaluations: (9165) 
#> Path [6] :Initial log joint density = -365.395131 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.548e+03      1.030e-01   1.925e+03    3.104e-02  6.174e-02      8257  2.280e+03 -2.210e+04                   
#> Path [6] :Best Iter: [27] ELBO (2280.349779) evaluations: (8257) 
#> Path [7] :Initial log joint density = -428.442405 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.567e+03      1.008e-01   4.908e+03    3.536e-02  3.536e-02      8646  2.280e+03 -8.924e+03                   
#> Path [7] :Best Iter: [32] ELBO (2280.208431) evaluations: (8646) 
#> Path [8] :Initial log joint density = -220.070657 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.560e+03      4.852e-02   4.520e+03    1.177e-02  1.177e-02      8617  2.277e+03 -4.195e+10                   
#> Path [8] :Best Iter: [33] ELBO (2276.719359) evaluations: (8617) 
#> Path [9] :Initial log joint density = -557.215748 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.550e+03      4.432e-02   2.538e+03    3.787e-02  3.787e-02      8132  2.282e+03 -4.631e+03                   
#> Path [9] :Best Iter: [34] ELBO (2282.026798) evaluations: (8132) 
#> Path [10] :Initial log joint density = -524.013626 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.535e+03      3.022e-02   1.783e+03    3.389e-02  3.389e-02      8957  2.278e+03  5.884e+02                   
#> Path [10] :Best Iter: [35] ELBO (2277.544277) evaluations: (8957) 
#> Path [11] :Initial log joint density = -1079.727001 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.542e+03      3.335e-02   3.996e+03    1.907e-02  1.907e-02      8646  2.281e+03 -5.303e+04                   
#> Path [11] :Best Iter: [31] ELBO (2281.186277) evaluations: (8646) 
#> Path [12] :Initial log joint density = -465.690189 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.558e+03      4.113e-02   4.346e+03    1.451e-02  3.513e-02      8493  2.280e+03 -3.931e+04                   
#> Path [12] :Best Iter: [30] ELBO (2280.309476) evaluations: (8493) 
#> Path [13] :Initial log joint density = -178.211915 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      6.868e-02   1.368e+04    1.353e-02  2.792e-02      8674  2.282e+03 -4.164e+05                   
#> Path [13] :Best Iter: [30] ELBO (2281.639912) evaluations: (8674) 
#> Path [14] :Initial log joint density = -230.265684 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      8.442e-02   4.200e+03    4.441e-02  4.441e-02      8748  2.281e+03 -2.903e+04                   
#> Path [14] :Best Iter: [34] ELBO (2280.937740) evaluations: (8748) 
#> Path [15] :Initial log joint density = -793.051993 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.552e+03      1.684e-02   3.400e+03    2.031e-02  2.031e-02      8790  2.280e+03  6.368e+02                   
#> Path [15] :Best Iter: [25] ELBO (2280.073105) evaluations: (8790) 
#> Path [16] :Initial log joint density = -578.436960 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.564e+03      8.278e-02   5.328e+03    5.374e-02  5.374e-02      8377  2.280e+03 -6.433e+05                   
#> Path [16] :Best Iter: [33] ELBO (2280.497199) evaluations: (8377) 
#> Path [17] :Initial log joint density = -2028.358989 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.541e+03      3.679e-02   2.322e+03    1.005e-02  3.663e-02      8500  2.284e+03 -3.824e+06                   
#> Path [17] :Best Iter: [37] ELBO (2283.887134) evaluations: (8500) 
#> Path [18] :Initial log joint density = -1432.743671 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.535e+03      1.013e-01   2.024e+03    4.631e-02  8.259e-02      8283  2.281e+03 -7.973e+04                   
#> Path [18] :Best Iter: [32] ELBO (2281.183175) evaluations: (8283) 
#> Path [19] :Initial log joint density = -588.197326 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      8.220e-02   1.108e+04    1.056e-02  2.928e-02      8635  2.283e+03 -5.722e+06                   
#> Path [19] :Best Iter: [31] ELBO (2283.309332) evaluations: (8635) 
#> Path [20] :Initial log joint density = -623.767183 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.535e+03      7.503e-02   2.004e+03    5.798e-02  5.798e-02      8777  2.280e+03 -3.290e+04                   
#> Path [20] :Best Iter: [39] ELBO (2280.071008) evaluations: (8777) 
#> Path [21] :Initial log joint density = -634.509077 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.563e+03      1.285e-01   1.165e+04    3.788e-02  7.475e-02      8477  2.281e+03 -6.214e+05                   
#> Path [21] :Best Iter: [31] ELBO (2281.210964) evaluations: (8477) 
#> Path [22] :Initial log joint density = -2499.201417 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      7.821e-02   4.823e+03    3.161e-02  3.161e-02      8625  2.282e+03 -3.440e+04                   
#> Path [22] :Best Iter: [33] ELBO (2281.623174) evaluations: (8625) 
#> Path [23] :Initial log joint density = -411.836333 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      7.747e-02   4.734e+03    3.328e-02  6.266e-02      8636  2.277e+03 -1.151e+05                   
#> Path [23] :Best Iter: [28] ELBO (2277.209557) evaluations: (8636) 
#> Path [24] :Initial log joint density = -559.424689 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.546e+03      3.691e-02   3.027e+03    2.760e-02  2.760e-02      8434  2.275e+03 -1.843e+04                   
#> Path [24] :Best Iter: [33] ELBO (2275.382033) evaluations: (8434) 
#> Path [25] :Initial log joint density = -522.532230 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.533e+03      2.821e-02   1.835e+03    1.866e-02  1.866e-02      8243  2.279e+03 -1.575e+03                   
#> Path [25] :Best Iter: [40] ELBO (2278.849103) evaluations: (8243) 
#> Path [26] :Initial log joint density = -691.995282 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      2.518e-02   3.997e+03    1.741e-02  1.741e-02      8682  2.282e+03 -2.655e+04                   
#> Path [26] :Best Iter: [35] ELBO (2281.514117) evaluations: (8682) 
#> Path [27] :Initial log joint density = -512.403578 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.538e+03      2.831e-02   3.648e+03    2.495e-02  2.495e-02      8365  2.283e+03 -5.230e+04                   
#> Path [27] :Best Iter: [39] ELBO (2282.767566) evaluations: (8365) 
#> Path [28] :Initial log joint density = -1575.241488 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.528e+03      3.966e-02   2.725e+03    2.048e-02  2.048e-02      8323  2.282e+03 -1.524e+06                   
#> Path [28] :Best Iter: [40] ELBO (2282.192954) evaluations: (8323) 
#> Path [29] :Initial log joint density = -275.740487 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.558e+03      7.615e-02   5.330e+03    6.129e-02  6.129e-02      9127  2.283e+03 -3.109e+05                   
#> Path [29] :Best Iter: [31] ELBO (2282.530924) evaluations: (9127) 
#> Path [30] :Initial log joint density = -715.116811 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.544e+03      6.860e-02   2.617e+03    1.621e-02  3.373e-02      8722  2.279e+03 -1.400e+04                   
#> Path [30] :Best Iter: [38] ELBO (2279.198273) evaluations: (8722) 
#> Path [31] :Initial log joint density = -2684.787615 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      1.710e-02   8.064e+03    1.529e-02  1.529e-02      8611  2.282e+03 -4.490e+03                   
#> Path [31] :Best Iter: [29] ELBO (2281.915604) evaluations: (8611) 
#> Path [32] :Initial log joint density = -485.198740 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.565e+03      7.640e-02   1.268e+04    6.079e-02  6.079e-02      8925  2.283e+03  6.951e+02                   
#> Path [32] :Best Iter: [30] ELBO (2283.078447) evaluations: (8925) 
#> Path [33] :Initial log joint density = -6078.335094 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.541e+03      7.620e-02   3.937e+03    4.242e-02  4.242e-02      8665  2.279e+03 -3.366e+03                   
#> Path [33] :Best Iter: [36] ELBO (2279.453213) evaluations: (8665) 
#> Path [34] :Initial log joint density = -1173.869259 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.540e+03      4.209e-02   2.384e+03    1.676e-02  1.676e-02      8615  2.284e+03 -9.163e+03                   
#> Path [34] :Best Iter: [34] ELBO (2284.132407) evaluations: (8615) 
#> Path [35] :Initial log joint density = -332.811121 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      4.870e-02   5.141e+03    2.172e-02  2.172e-02      8255  2.279e+03 -3.792e+05                   
#> Path [35] :Best Iter: [30] ELBO (2278.778488) evaluations: (8255) 
#> Path [36] :Initial log joint density = -747.060435 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.559e+03      1.503e-02   3.095e+03    1.378e-02  1.378e-02      8201  2.278e+03 -2.279e+04                   
#> Path [36] :Best Iter: [32] ELBO (2277.511493) evaluations: (8201) 
#> Path [37] :Initial log joint density = -526.723191 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.561e+03      6.232e-02   5.949e+03    4.551e-02  4.551e-02      8650  2.281e+03 -2.269e+03                   
#> Path [37] :Best Iter: [30] ELBO (2280.558136) evaluations: (8650) 
#> Path [38] :Initial log joint density = -404.384347 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      1.375e-01   6.303e+03    5.706e-02  9.788e-02      8664  2.282e+03 -7.229e+04                   
#> Path [38] :Best Iter: [27] ELBO (2282.465134) evaluations: (8664) 
#> Path [39] :Initial log joint density = -632.941179 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.550e+03      4.273e-02   2.624e+03    1.527e-02  1.527e-02      9212  2.279e+03 -1.368e+04                   
#> Path [39] :Best Iter: [33] ELBO (2279.162542) evaluations: (9212) 
#> Path [40] :Initial log joint density = -775.346346 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.530e+03      6.272e-02   1.528e+03    8.953e-02  8.953e-02      8525  2.282e+03  5.353e+01                   
#> Path [40] :Best Iter: [36] ELBO (2282.472913) evaluations: (8525) 
#> Path [41] :Initial log joint density = -712.964437 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      1.158e-01   4.045e+03    2.051e-02  4.673e-02      8777  2.279e+03 -1.554e+06                   
#> Path [41] :Best Iter: [30] ELBO (2279.017463) evaluations: (8777) 
#> Path [42] :Initial log joint density = -742.481724 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.535e+03      4.342e-02   1.681e+03    2.675e-02  2.675e-02      8150  2.285e+03 -1.228e+04                   
#> Path [42] :Best Iter: [38] ELBO (2284.697530) evaluations: (8150) 
#> Path [43] :Initial log joint density = -896.167512 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.540e+03      5.332e-02   3.331e+03    1.454e-02  1.454e-02      8271  2.279e+03 -7.228e+06                   
#> Path [43] :Best Iter: [40] ELBO (2278.959754) evaluations: (8271) 
#> Path [44] :Initial log joint density = -685.456581 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.553e+03      1.849e-02   3.474e+03    2.683e-02  2.683e-02      8431  2.279e+03 -1.661e+04                   
#> Path [44] :Best Iter: [33] ELBO (2278.727162) evaluations: (8431) 
#> Path [45] :Initial log joint density = -643.163599 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.559e+03      7.774e-02   4.767e+03    5.292e-02  5.292e-02      8529  2.282e+03 -4.488e+05                   
#> Path [45] :Best Iter: [36] ELBO (2281.509429) evaluations: (8529) 
#> Path [46] :Initial log joint density = -1335.371318 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.537e+03      8.866e-02   2.322e+03    3.951e-02  3.951e-02      8568  2.283e+03 -3.178e+05                   
#> Path [46] :Best Iter: [38] ELBO (2282.753006) evaluations: (8568) 
#> Path [47] :Initial log joint density = -1070.599605 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.545e+03      7.592e-02   5.141e+03    2.477e-02  2.477e-02      8665  2.281e+03 -1.938e+05                   
#> Path [47] :Best Iter: [32] ELBO (2281.286466) evaluations: (8665) 
#> Path [48] :Initial log joint density = -616.101402 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.554e+03      3.762e-02   1.335e+04    1.550e-02  1.550e-02      8515  2.281e+03 -2.545e+04                   
#> Path [48] :Best Iter: [36] ELBO (2281.348816) evaluations: (8515) 
#> Path [49] :Initial log joint density = -515.544733 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.559e+03      4.951e-02   4.184e+03    6.724e-02  6.724e-02      8553  2.282e+03 -4.611e+03                   
#> Path [49] :Best Iter: [23] ELBO (2281.706488) evaluations: (8553) 
#> Path [50] :Initial log joint density = -461.375064 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.567e+03      5.172e-02   7.435e+03    2.628e-02  2.628e-02      8422  2.280e+03 -1.284e+06                   
#> Path [50] :Best Iter: [33] ELBO (2279.898143) evaluations: (8422) 
#> Finished in  15.7 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
