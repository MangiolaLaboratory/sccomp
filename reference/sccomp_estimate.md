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
  portable = TRUE,
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

- portable:

  Logical, whether to keep the result portable by caching required draws
  in memory and removing Stan draw CSV files after fitting. Default is
  TRUE to save disk space and move needed values into memory. Set to
  FALSE to keep draw CSV files on disk.

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

- c_rhat - R-hat convergence diagnostic for the composition (c)
  parameter; values close to 1.0 indicate convergence.

- c_ess_bulk - Bulk effective sample size for the composition (c)
  parameter; higher is better.

- c_ess_tail - Tail effective sample size for the composition (c)
  parameter; higher is better.

- v_lower - Lower (2.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_effect - Mean of the posterior distribution for a variability (v)
  parameter.

- v_upper - Upper (97.5%) quantile of the posterior distribution for a
  variability (v) parameter.

- v_rhat - R-hat convergence diagnostic for the variability (v)
  parameter.

- v_ess_bulk - Bulk effective sample size for the variability (v)
  parameter.

- v_ess_tail - Tail effective sample size for the variability (v)
  parameter.

Note: pH0 and FDR columns are not computed by `sccomp_estimate()`. Run
[`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md)
on the result to obtain hypothesis-test statistics.

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
#> Path [1] :Initial log joint density = -481751.550185 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.941e-03   2.010e-01    1.000e+00  1.000e+00      3162 -3.687e+03 -3.686e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3686.338709) evaluations: (3162) 
#> Path [2] :Initial log joint density = -481686.861673 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.172e-03   2.012e-01    1.000e+00  1.000e+00      3326 -3.682e+03 -3.685e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3682.240473) evaluations: (3326) 
#> Path [3] :Initial log joint density = -481464.681229 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.150e-02   1.693e-01    1.000e+00  1.000e+00      3530 -3.684e+03 -3.682e+03                   
#> Path [3] :Best Iter: [59] ELBO (-3682.157472) evaluations: (3530) 
#> Path [4] :Initial log joint density = -483810.991261 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.108e-03   2.560e-01    9.953e-01  9.953e-01      3208 -3.684e+03 -3.695e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3684.468762) evaluations: (3208) 
#> Path [5] :Initial log joint density = -482055.657247 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.056e-02   1.992e-01    1.000e+00  1.000e+00      3506 -3.684e+03 -3.684e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3684.085441) evaluations: (3506) 
#> Path [6] :Initial log joint density = -481656.335751 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.275e-02   2.731e-01    1.000e+00  1.000e+00      3085 -3.688e+03 -3.690e+03                   
#> Path [6] :Best Iter: [54] ELBO (-3688.391916) evaluations: (3085) 
#> Path [7] :Initial log joint density = -487949.565364 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.509e-03   2.064e-01    7.870e-01  7.870e-01      3384 -3.692e+03 -3.690e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3689.647185) evaluations: (3384) 
#> Path [8] :Initial log joint density = -481680.352033 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.156e-02   3.629e-01    1.000e+00  1.000e+00      3366 -3.682e+03 -3.695e+03                   
#> Path [8] :Best Iter: [58] ELBO (-3682.373152) evaluations: (3366) 
#> Path [9] :Initial log joint density = -482092.714677 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.871e-03   2.345e-01    7.868e-01  7.868e-01      3363 -3.684e+03 -3.695e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3684.216682) evaluations: (3363) 
#> Path [10] :Initial log joint density = -481748.115299 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.356e-03   2.276e-01    1.000e+00  1.000e+00      3088 -3.689e+03 -3.705e+03                   
#> Path [10] :Best Iter: [42] ELBO (-3689.474354) evaluations: (3088) 
#> Path [11] :Initial log joint density = -481920.931249 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.686e-03   2.480e-01    1.000e+00  1.000e+00      3262 -3.690e+03 -3.687e+03                   
#> Path [11] :Best Iter: [57] ELBO (-3687.435064) evaluations: (3262) 
#> Path [12] :Initial log joint density = -483814.169988 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.190e-03   2.260e-01    1.000e+00  1.000e+00      3580 -3.682e+03 -3.684e+03                   
#> Path [12] :Best Iter: [58] ELBO (-3682.211246) evaluations: (3580) 
#> Path [13] :Initial log joint density = -481993.074616 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.511e-03   2.312e-01    1.000e+00  1.000e+00      3363 -3.683e+03 -3.685e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3682.877475) evaluations: (3363) 
#> Path [14] :Initial log joint density = -481546.180125 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.160e-02   2.390e-01    1.000e+00  1.000e+00      3156 -3.683e+03 -3.686e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3683.490607) evaluations: (3156) 
#> Path [15] :Initial log joint density = -483586.445897 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.995e-03   3.065e-01    1.000e+00  1.000e+00      3210 -3.687e+03 -3.691e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3687.121067) evaluations: (3210) 
#> Path [16] :Initial log joint density = -481336.222752 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.425e-03   2.528e-01    7.433e-01  7.433e-01      2997 -3.690e+03 -3.706e+03                   
#> Path [16] :Best Iter: [52] ELBO (-3689.811993) evaluations: (2997) 
#> Path [17] :Initial log joint density = -484264.831140 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.165e-03   2.160e-01    1.000e+00  1.000e+00      3564 -3.684e+03 -3.685e+03                   
#> Path [17] :Best Iter: [57] ELBO (-3684.271151) evaluations: (3564) 
#> Path [18] :Initial log joint density = -485075.168993 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      5.873e-03   2.666e-01    3.994e-01  1.000e+00      3651 -3.684e+03 -3.694e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3684.187247) evaluations: (3651) 
#> Path [19] :Initial log joint density = -481721.085645 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.840e-02   2.843e-01    9.526e-01  9.526e-01      3650 -3.682e+03 -3.690e+03                   
#> Path [19] :Best Iter: [59] ELBO (-3681.688670) evaluations: (3650) 
#> Path [20] :Initial log joint density = -481704.409446 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.559e-03   2.974e-01    1.000e+00  1.000e+00      2830 -3.689e+03 -3.692e+03                   
#> Path [20] :Best Iter: [42] ELBO (-3689.277242) evaluations: (2830) 
#> Path [21] :Initial log joint density = -481567.357548 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.307e-03   2.571e-01    1.000e+00  1.000e+00      3157 -3.687e+03 -3.693e+03                   
#> Path [21] :Best Iter: [52] ELBO (-3686.561908) evaluations: (3157) 
#> Path [22] :Initial log joint density = -481628.056474 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.991e-03   2.111e-01    1.000e+00  1.000e+00      3245 -3.684e+03 -3.682e+03                   
#> Path [22] :Best Iter: [57] ELBO (-3682.140270) evaluations: (3245) 
#> Path [23] :Initial log joint density = -483418.271969 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.631e-03   1.970e-01    1.000e+00  1.000e+00      3355 -3.685e+03 -3.686e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3684.657076) evaluations: (3355) 
#> Path [24] :Initial log joint density = -481755.023643 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.511e-03   2.297e-01    1.000e+00  1.000e+00      3356 -3.685e+03 -3.685e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3684.687182) evaluations: (3356) 
#> Path [25] :Initial log joint density = -481808.875878 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.236e-02   3.425e-01    1.000e+00  1.000e+00      3210 -3.685e+03 -3.692e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3685.145823) evaluations: (3210) 
#> Path [26] :Initial log joint density = -481892.323209 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.492e-02   2.360e-01    1.000e+00  1.000e+00      3173 -3.682e+03 -3.684e+03                   
#> Path [26] :Best Iter: [56] ELBO (-3682.249949) evaluations: (3173) 
#> Path [27] :Initial log joint density = -481920.138444 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.856e-03   2.015e-01    1.000e+00  1.000e+00      2753 -3.694e+03 -3.710e+03                   
#> Path [27] :Best Iter: [37] ELBO (-3694.306622) evaluations: (2753) 
#> Path [28] :Initial log joint density = -481568.059327 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.119e-02   1.592e-01    1.000e+00  1.000e+00      3739 -3.681e+03 -3.684e+03                   
#> Path [28] :Best Iter: [59] ELBO (-3681.089195) evaluations: (3739) 
#> Path [29] :Initial log joint density = -481878.640423 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.958e-03   2.150e-01    9.226e-01  9.226e-01      3165 -3.684e+03 -3.691e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3683.797474) evaluations: (3165) 
#> Path [30] :Initial log joint density = -481470.723336 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.655e-03   2.551e-01    7.049e-01  7.049e-01      3071 -3.691e+03 -3.698e+03                   
#> Path [30] :Best Iter: [44] ELBO (-3690.778500) evaluations: (3071) 
#> Path [31] :Initial log joint density = -483936.217416 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.855e-03   1.500e-01    1.000e+00  1.000e+00      3269 -3.687e+03 -3.695e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3687.054685) evaluations: (3269) 
#> Path [32] :Initial log joint density = -481822.662240 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.643e-03   2.647e-01    6.388e-01  6.388e-01      3621 -3.683e+03 -3.698e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3682.833183) evaluations: (3621) 
#> Path [33] :Initial log joint density = -481527.117951 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.980e-03   2.516e-01    5.221e-01  5.221e-01      3321 -3.691e+03 -3.684e+03                   
#> Path [33] :Best Iter: [57] ELBO (-3684.207553) evaluations: (3321) 
#> Path [34] :Initial log joint density = -481612.925321 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.472e-02   3.064e-01    1.000e+00  1.000e+00      3233 -3.688e+03 -3.690e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3688.302135) evaluations: (3233) 
#> Path [35] :Initial log joint density = -482072.775006 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.867e-03   2.145e-01    8.794e-01  8.794e-01      3306 -3.685e+03 -3.690e+03                   
#> Path [35] :Best Iter: [56] ELBO (-3684.608110) evaluations: (3306) 
#> Path [36] :Initial log joint density = -482286.506819 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      2.956e-03   3.541e-01    5.430e-01  5.430e-01      3277 -3.684e+03 -3.700e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3683.861293) evaluations: (3277) 
#> Path [37] :Initial log joint density = -482510.712134 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      3.952e-03   1.072e-01    7.972e-01  7.972e-01      3394 -3.685e+03 -3.696e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3684.760754) evaluations: (3394) 
#> Path [38] :Initial log joint density = -481710.959788 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.602e-03   2.094e-01    1.000e+00  1.000e+00      2996 -3.693e+03 -3.702e+03                   
#> Path [38] :Best Iter: [50] ELBO (-3692.624025) evaluations: (2996) 
#> Path [39] :Initial log joint density = -482527.773636 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.075e-03   2.604e-01    6.720e-01  6.720e-01      3170 -3.689e+03 -3.694e+03                   
#> Path [39] :Best Iter: [52] ELBO (-3689.413764) evaluations: (3170) 
#> Path [40] :Initial log joint density = -481372.498938 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.721e-03   2.061e-01    7.387e-01  7.387e-01      3099 -3.687e+03 -3.699e+03                   
#> Path [40] :Best Iter: [55] ELBO (-3687.010606) evaluations: (3099) 
#> Path [41] :Initial log joint density = -481680.003892 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      2.975e-03   2.116e-01    6.696e-01  6.696e-01      3748 -3.685e+03 -3.696e+03                   
#> Path [41] :Best Iter: [61] ELBO (-3684.816194) evaluations: (3748) 
#> Path [42] :Initial log joint density = -481384.611341 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.805e-02   2.744e-01    1.000e+00  1.000e+00      2926 -3.694e+03 -3.695e+03                   
#> Path [42] :Best Iter: [52] ELBO (-3693.763719) evaluations: (2926) 
#> Path [43] :Initial log joint density = -481948.541685 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.397e-02   2.162e-01    1.000e+00  1.000e+00      3451 -3.686e+03 -3.685e+03                   
#> Path [43] :Best Iter: [60] ELBO (-3685.316768) evaluations: (3451) 
#> Path [44] :Initial log joint density = -481757.515673 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.508e-03   2.281e-01    1.000e+00  1.000e+00      3356 -3.684e+03 -3.697e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3684.060841) evaluations: (3356) 
#> Path [45] :Initial log joint density = -481662.285473 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.033e-03   3.160e-01    5.872e-01  5.872e-01      3087 -3.692e+03 -3.703e+03                   
#> Path [45] :Best Iter: [41] ELBO (-3691.521610) evaluations: (3087) 
#> Path [46] :Initial log joint density = -481755.909120 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.185e-03   2.042e-01    1.000e+00  1.000e+00      3304 -3.687e+03 -3.685e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3685.009095) evaluations: (3304) 
#> Path [47] :Initial log joint density = -484993.286446 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.335e-03   2.475e-01    9.128e-01  9.128e-01      3233 -3.684e+03 -3.693e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3684.289968) evaluations: (3233) 
#> Path [48] :Initial log joint density = -481638.652696 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.420e-03   1.874e-01    8.641e-01  8.641e-01      2915 -3.693e+03 -3.706e+03                   
#> Path [48] :Best Iter: [52] ELBO (-3692.615928) evaluations: (2915) 
#> Path [49] :Initial log joint density = -481591.010756 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.985e-03   2.751e-01    3.708e-01  1.000e+00      3414 -3.687e+03 -3.696e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3687.364409) evaluations: (3414) 
#> Path [50] :Initial log joint density = -482085.960104 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.048e-03   2.400e-01    7.387e-01  7.387e-01      3646 -3.685e+03 -3.697e+03                   
#> Path [50] :Best Iter: [59] ELBO (-3684.887737) evaluations: (3646) 
#> Finished in  13.7 seconds.
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
#> Path [1] :Initial log joint density = -279.788648 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.552e+03      3.533e-02   2.690e+03    2.446e-02  2.446e-02      8767  2.298e+03 -4.074e+04                   
#> Path [1] :Best Iter: [36] ELBO (2297.961563) evaluations: (8767) 
#> Path [2] :Initial log joint density = -1239.280811 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.553e+03      2.839e-02   1.825e+03    3.022e-02  3.022e-02      8637  2.297e+03 -5.758e+02                   
#> Path [2] :Best Iter: [34] ELBO (2297.141404) evaluations: (8637) 
#> Path [3] :Initial log joint density = -474.798212 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.570e+03      5.144e-02   4.588e+03    2.842e-02  2.842e-02      8763  2.289e+03 -1.118e+05                   
#> Path [3] :Best Iter: [35] ELBO (2289.302167) evaluations: (8763) 
#> Path [4] :Initial log joint density = -546.470418 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.571e+03      8.069e-02   6.416e+03    3.327e-02  3.327e-02      8665  2.295e+03 -1.994e+07                   
#> Path [4] :Best Iter: [33] ELBO (2295.452706) evaluations: (8665) 
#> Path [5] :Initial log joint density = -1144.623073 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      4.220e-02   2.386e+03    1.612e-02  1.612e-02      8456  2.291e+03 -2.364e+07                   
#> Path [5] :Best Iter: [34] ELBO (2291.351866) evaluations: (8456) 
#> Path [6] :Initial log joint density = -4658.191166 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      2.067e-02   6.002e+03    2.990e-02  2.990e-02      8284  2.298e+03 -1.936e+03                   
#> Path [6] :Best Iter: [31] ELBO (2297.579837) evaluations: (8284) 
#> Path [7] :Initial log joint density = -625.566915 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.553e+03      6.368e-02   2.837e+03    2.725e-02  5.260e-02      8285  2.298e+03 -1.765e+04                   
#> Path [7] :Best Iter: [41] ELBO (2298.112928) evaluations: (8285) 
#> Path [8] :Initial log joint density = -603.016745 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.552e+03      4.577e-02   2.367e+03    2.503e-02  6.251e-02      8507  2.298e+03 -3.405e+04                   
#> Path [8] :Best Iter: [39] ELBO (2297.706181) evaluations: (8507) 
#> Path [9] :Initial log joint density = -554.804910 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      6.788e-02   5.983e+03    5.876e-02  5.876e-02      8949  2.297e+03 -3.974e+04                   
#> Path [9] :Best Iter: [28] ELBO (2297.346423) evaluations: (8949) 
#> Path [10] :Initial log joint density = -1854.405476 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      5.259e-02   2.783e+03    2.460e-02  4.582e-02      8201  2.300e+03 -5.306e+05                   
#> Path [10] :Best Iter: [33] ELBO (2299.799545) evaluations: (8201) 
#> Path [11] :Initial log joint density = -859.568497 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      8.696e-02   2.299e+03    2.877e-02  6.060e-02      8390  2.296e+03 -3.574e+05                   
#> Path [11] :Best Iter: [33] ELBO (2296.417716) evaluations: (8390) 
#> Path [12] :Initial log joint density = -120.031735 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.576e+03      1.260e-01   4.823e+03    3.374e-02  3.374e-02      8625  2.296e+03 -2.130e+03                   
#> Path [12] :Best Iter: [29] ELBO (2296.261307) evaluations: (8625) 
#> Path [13] :Initial log joint density = -2610.885626 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.564e+03      1.134e-01   5.024e+03    6.069e-02  6.069e-02      8748  2.300e+03 -4.697e+03                   
#> Path [13] :Best Iter: [36] ELBO (2299.612792) evaluations: (8748) 
#> Path [14] :Initial log joint density = -1936.510076 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.565e+03      1.307e-02   2.239e+03    4.634e-02  4.634e-02      8945  2.296e+03  1.008e+03                   
#> Path [14] :Best Iter: [34] ELBO (2296.031651) evaluations: (8945) 
#> Path [15] :Initial log joint density = -587.622410 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.578e+03      4.295e-02   7.600e+03    2.041e-02  2.041e-02      8514  2.293e+03 -3.633e+04                   
#> Path [15] :Best Iter: [30] ELBO (2292.641518) evaluations: (8514) 
#> Path [16] :Initial log joint density = -630.946411 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.580e+03      3.075e-02   5.289e+03    1.834e-02  1.834e-02      8950  2.294e+03 -1.043e+05                   
#> Path [16] :Best Iter: [32] ELBO (2293.804440) evaluations: (8950) 
#> Path [17] :Initial log joint density = -486.296757 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.576e+03      1.411e-02   8.237e+03    1.832e-02  1.832e-02      8746  2.294e+03 -8.674e+04                   
#> Path [17] :Best Iter: [34] ELBO (2293.869015) evaluations: (8746) 
#> Path [18] :Initial log joint density = -1779.626518 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.573e+03      4.891e-02   6.357e+03    1.528e-02  1.528e-02      8525  2.299e+03 -3.819e+07                   
#> Path [18] :Best Iter: [36] ELBO (2298.797924) evaluations: (8525) 
#> Path [19] :Initial log joint density = -729.209888 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.577e+03      2.218e-02   5.096e+03    5.732e-02  5.732e-02      8587  2.296e+03 -1.217e+03                   
#> Path [19] :Best Iter: [29] ELBO (2296.498711) evaluations: (8587) 
#> Path [20] :Initial log joint density = -478.949746 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      1.909e-02   3.619e+03    3.461e-02  3.461e-02      8748  2.295e+03  2.505e+02                   
#> Path [20] :Best Iter: [35] ELBO (2295.156681) evaluations: (8748) 
#> Path [21] :Initial log joint density = -220.295149 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.572e+03      1.051e-01   3.820e+03    5.757e-02  5.757e-02      8759  2.296e+03 -4.267e+04                   
#> Path [21] :Best Iter: [30] ELBO (2296.126093) evaluations: (8759) 
#> Path [22] :Initial log joint density = -504.823892 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.562e+03      8.394e-02   7.953e+03    2.868e-02  7.715e-02      8592  2.296e+03 -1.410e+07                   
#> Path [22] :Best Iter: [35] ELBO (2296.432437) evaluations: (8592) 
#> Path [23] :Initial log joint density = -508.176318 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      2.574e-02   3.366e+03    1.186e-02  1.186e-02      8243  2.294e+03 -9.041e+04                   
#> Path [23] :Best Iter: [34] ELBO (2294.009583) evaluations: (8243) 
#> Path [24] :Initial log joint density = -513.940480 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.576e+03      1.249e-01   7.716e+03    4.327e-02  4.327e-02      9239  2.295e+03 -9.130e+04                   
#> Path [24] :Best Iter: [31] ELBO (2295.447870) evaluations: (9239) 
#> Path [25] :Initial log joint density = -531.711243 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.561e+03      7.276e-02   2.552e+03    5.909e-02  5.909e-02      8961  2.297e+03  1.840e+02                   
#> Path [25] :Best Iter: [33] ELBO (2297.459457) evaluations: (8961) 
#> Path [26] :Initial log joint density = -691.757357 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.582e+03      5.193e-02   8.918e+03    1.668e-02  1.668e-02      8696  2.298e+03 -2.451e+05                   
#> Path [26] :Best Iter: [28] ELBO (2298.163999) evaluations: (8696) 
#> Path [27] :Initial log joint density = -846.994103 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.557e+03      5.911e-02   4.087e+03    2.196e-02  2.196e-02      8335  2.298e+03 -3.052e+04                   
#> Path [27] :Best Iter: [41] ELBO (2297.595717) evaluations: (8335) 
#> Path [28] :Initial log joint density = -1564.235792 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      4.555e-02   3.252e+03    5.489e-02  5.489e-02      8549  2.294e+03 -7.234e+02                   
#> Path [28] :Best Iter: [35] ELBO (2293.553310) evaluations: (8549) 
#> Path [29] :Initial log joint density = -93.211367 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.584e+03      4.880e-02   1.099e+04    4.529e-02  4.529e-02      8290  2.299e+03 -2.828e+04                   
#> Path [29] :Best Iter: [35] ELBO (2299.398301) evaluations: (8290) 
#> Path [30] :Initial log joint density = -319.570028 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.574e+03      5.302e-02   4.803e+03    2.639e-02  2.639e-02      8690  2.299e+03 -3.639e+06                   
#> Path [30] :Best Iter: [34] ELBO (2299.386378) evaluations: (8690) 
#> Path [31] :Initial log joint density = -566.027925 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.577e+03      1.066e-01   6.162e+03    4.445e-02  4.445e-02      9043  2.296e+03 -1.001e+05                   
#> Path [31] :Best Iter: [29] ELBO (2295.738114) evaluations: (9043) 
#> Path [32] :Initial log joint density = -399.553190 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.565e+03      7.349e-02   5.880e+03    2.139e-02  2.139e-02      8839  2.296e+03 -1.790e+04                   
#> Path [32] :Best Iter: [33] ELBO (2295.721215) evaluations: (8839) 
#> Path [33] :Initial log joint density = -1365.210378 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.555e+03      1.808e-02   4.385e+03    2.442e-02  2.442e-02      8375  2.297e+03 -9.986e+03                   
#> Path [33] :Best Iter: [36] ELBO (2297.161132) evaluations: (8375) 
#> Path [34] :Initial log joint density = -654.434265 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.575e+03      5.353e-02   8.184e+03    1.546e-02  6.681e-02      8990  2.298e+03 -1.298e+05                   
#> Path [34] :Best Iter: [31] ELBO (2298.423265) evaluations: (8990) 
#> Path [35] :Initial log joint density = -277.621319 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.576e+03      9.111e-02   5.193e+03    5.207e-02  5.207e-02      9273  2.300e+03 -2.028e+04                   
#> Path [35] :Best Iter: [29] ELBO (2300.426698) evaluations: (9273) 
#> Path [36] :Initial log joint density = -421.305714 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.583e+03      6.989e-02   7.212e+03    1.479e-02  2.636e-02      8505  2.300e+03 -1.409e+06                   
#> Path [36] :Best Iter: [30] ELBO (2299.908184) evaluations: (8505) 
#> Path [37] :Initial log joint density = -773.887443 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.560e+03      9.011e-02   4.073e+03    7.255e-02  7.255e-02      8472  2.295e+03 -1.864e+04                   
#> Path [37] :Best Iter: [36] ELBO (2294.892453) evaluations: (8472) 
#> Path [38] :Initial log joint density = -419.464127 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.584e+03      3.222e-02   4.711e+03    3.065e-02  3.065e-02      8629  2.297e+03 -4.921e+04                   
#> Path [38] :Best Iter: [30] ELBO (2297.458387) evaluations: (8629) 
#> Path [39] :Initial log joint density = -583.091271 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      9.595e-02   2.474e+03    4.528e-02  8.627e-02      8340  2.293e+03 -6.689e+03                   
#> Path [39] :Best Iter: [32] ELBO (2293.102469) evaluations: (8340) 
#> Path [40] :Initial log joint density = -289.941090 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.565e+03      4.177e-02   4.647e+03    1.085e-02  1.085e-02      8825  2.295e+03 -2.509e+04                   
#> Path [40] :Best Iter: [33] ELBO (2295.055407) evaluations: (8825) 
#> Path [41] :Initial log joint density = -752.544326 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.559e+03      6.794e-02   2.384e+03    5.547e-02  5.547e-02      8688  2.297e+03 -5.084e+03                   
#> Path [41] :Best Iter: [32] ELBO (2297.479697) evaluations: (8688) 
#> Path [42] :Initial log joint density = -329.807205 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.577e+03      6.764e-02   7.579e+03    4.195e-02  4.195e-02      8902  2.301e+03 -4.858e+04                   
#> Path [42] :Best Iter: [31] ELBO (2300.579929) evaluations: (8902) 
#> Path [43] :Initial log joint density = -680.052224 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.568e+03      4.755e-02   5.233e+03    1.354e-02  1.354e-02      8980  2.298e+03 -4.471e+04                   
#> Path [43] :Best Iter: [33] ELBO (2297.505940) evaluations: (8980) 
#> Path [44] :Initial log joint density = -456.617749 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.569e+03      3.185e-02   7.020e+03    1.694e-02  1.694e-02      8497  2.302e+03 -3.447e+04                   
#> Path [44] :Best Iter: [35] ELBO (2301.924542) evaluations: (8497) 
#> Path [45] :Initial log joint density = -3153.224863 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.558e+03      2.312e-02   5.397e+03    5.860e-03  5.860e-03      8662  2.297e+03 -1.953e+06                   
#> Path [45] :Best Iter: [39] ELBO (2297.370215) evaluations: (8662) 
#> Path [46] :Initial log joint density = -7395.010781 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.555e+03      5.790e-02   1.740e+03    3.332e-02  3.332e-02      9036  2.301e+03 -2.245e+05                   
#> Path [46] :Best Iter: [35] ELBO (2301.075624) evaluations: (9036) 
#> Path [47] :Initial log joint density = -646.482249 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.566e+03      5.180e-02   3.418e+03    2.934e-02  2.934e-02      8702  2.294e+03 -3.341e+05                   
#> Path [47] :Best Iter: [30] ELBO (2294.073916) evaluations: (8702) 
#> Path [48] :Initial log joint density = -911.214366 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.553e+03      2.220e-02   3.683e+03    1.170e-02  1.170e-02      8321  2.297e+03 -2.566e+04                   
#> Path [48] :Best Iter: [41] ELBO (2297.401414) evaluations: (8321) 
#> Path [49] :Initial log joint density = -897.267417 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.556e+03      7.639e-02   3.570e+03    3.062e-02  3.062e-02      8193  2.296e+03 -8.567e+04                   
#> Path [49] :Best Iter: [33] ELBO (2296.167145) evaluations: (8193) 
#> Path [50] :Initial log joint density = -577.599259 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100       2.585e+03      7.832e-02   1.013e+04    1.406e-02  3.353e-02      8599  2.298e+03 -4.137e+08                   
#> Path [50] :Best Iter: [25] ELBO (2297.984471) evaluations: (8599) 
#> Finished in  16.0 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
