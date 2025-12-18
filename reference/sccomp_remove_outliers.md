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
#> Path [1] :Initial log joint density = -481896.467352 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.512e-03   2.101e-01    6.694e-01  6.694e-01      3450 -3.697e+03 -3.715e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3696.807202) evaluations: (3450) 
#> Path [2] :Initial log joint density = -484083.970253 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.071e-03   2.271e-01    8.427e-01  8.427e-01      3238 -3.704e+03 -3.711e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3703.553334) evaluations: (3238) 
#> Path [3] :Initial log joint density = -481222.846308 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.581e-02   2.352e-01    1.000e+00  1.000e+00      2964 -3.708e+03 -3.709e+03                   
#> Path [3] :Best Iter: [43] ELBO (-3708.289445) evaluations: (2964) 
#> Path [4] :Initial log joint density = -481773.234958 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.375e-02   2.651e-01    9.089e-01  9.089e-01      3000 -3.709e+03 -3.724e+03                   
#> Path [4] :Best Iter: [53] ELBO (-3709.305253) evaluations: (3000) 
#> Path [5] :Initial log joint density = -482892.999399 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.617e-03   2.319e-01    1.000e+00  1.000e+00      3236 -3.706e+03 -3.703e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3702.996901) evaluations: (3236) 
#> Path [6] :Initial log joint density = -483558.228932 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.555e-03   2.101e-01    8.697e-01  8.697e-01      3051 -3.706e+03 -3.721e+03                   
#> Path [6] :Best Iter: [52] ELBO (-3706.340078) evaluations: (3051) 
#> Path [7] :Initial log joint density = -481303.869532 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.165e-03   1.924e-01    1.000e+00  1.000e+00      3267 -3.706e+03 -3.710e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3705.578834) evaluations: (3267) 
#> Path [8] :Initial log joint density = -481686.263276 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.789e-03   2.467e-01    1.000e+00  1.000e+00      3568 -3.701e+03 -3.705e+03                   
#> Path [8] :Best Iter: [60] ELBO (-3701.034949) evaluations: (3568) 
#> Path [9] :Initial log joint density = -482383.940764 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.125e-03   2.594e-01    4.480e-01  1.000e+00      3382 -3.698e+03 -3.711e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3698.114163) evaluations: (3382) 
#> Path [10] :Initial log joint density = -481635.627523 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.280e-03   2.332e-01    1.000e+00  1.000e+00      3026 -3.708e+03 -3.702e+03                   
#> Path [10] :Best Iter: [55] ELBO (-3702.392418) evaluations: (3026) 
#> Path [11] :Initial log joint density = -481908.651841 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.559e-03   2.659e-01    6.441e-01  6.441e-01      3692 -3.699e+03 -3.711e+03                   
#> Path [11] :Best Iter: [60] ELBO (-3698.670819) evaluations: (3692) 
#> Path [12] :Initial log joint density = -481778.442824 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.360e-03   2.485e-01    1.000e+00  1.000e+00      3176 -3.705e+03 -3.705e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3704.517522) evaluations: (3176) 
#> Path [13] :Initial log joint density = -481786.031954 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.793e-03   2.057e-01    1.000e+00  1.000e+00      3212 -3.702e+03 -3.703e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3702.209203) evaluations: (3212) 
#> Path [14] :Initial log joint density = -481524.166240 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.523e-03   2.533e-01    1.000e+00  1.000e+00      2827 -3.709e+03 -3.715e+03                   
#> Path [14] :Best Iter: [48] ELBO (-3708.587090) evaluations: (2827) 
#> Path [15] :Initial log joint density = -484566.555895 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.244e-03   2.183e-01    1.000e+00  1.000e+00      3179 -3.701e+03 -3.700e+03                   
#> Path [15] :Best Iter: [57] ELBO (-3700.362695) evaluations: (3179) 
#> Path [16] :Initial log joint density = -481699.345883 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.834e-03   2.240e-01    6.887e-01  6.887e-01      2881 -3.708e+03 -3.720e+03                   
#> Path [16] :Best Iter: [51] ELBO (-3707.595965) evaluations: (2881) 
#> Path [17] :Initial log joint density = -481923.630681 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.858e-03   2.180e-01    7.100e-01  7.100e-01      3126 -3.709e+03 -3.717e+03                   
#> Path [17] :Best Iter: [48] ELBO (-3709.464405) evaluations: (3126) 
#> Path [18] :Initial log joint density = -481637.811270 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.806e-02   2.750e-01    1.000e+00  1.000e+00      3329 -3.701e+03 -3.700e+03                   
#> Path [18] :Best Iter: [58] ELBO (-3699.721646) evaluations: (3329) 
#> Path [19] :Initial log joint density = -481573.644062 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.350e-02   4.022e-01    1.000e+00  1.000e+00      3147 -3.712e+03 -3.705e+03                   
#> Path [19] :Best Iter: [55] ELBO (-3705.105296) evaluations: (3147) 
#> Path [20] :Initial log joint density = -481272.771647 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.241e-03   1.373e-01    1.000e+00  1.000e+00      3305 -3.699e+03 -3.712e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3699.089999) evaluations: (3305) 
#> Path [21] :Initial log joint density = -481920.990493 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.437e-02   2.811e-01    1.000e+00  1.000e+00      3461 -3.701e+03 -3.700e+03                   
#> Path [21] :Best Iter: [60] ELBO (-3700.482597) evaluations: (3461) 
#> Path [22] :Initial log joint density = -481317.850746 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.062e-02   2.637e-01    1.000e+00  1.000e+00      3157 -3.703e+03 -3.703e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3702.786997) evaluations: (3157) 
#> Path [23] :Initial log joint density = -481920.785305 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.656e-03   2.784e-01    1.000e+00  1.000e+00      3637 -3.701e+03 -3.704e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3701.470379) evaluations: (3637) 
#> Path [24] :Initial log joint density = -481664.566903 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.376e-03   2.411e-01    1.000e+00  1.000e+00      3079 -3.709e+03 -3.711e+03                   
#> Path [24] :Best Iter: [49] ELBO (-3709.097389) evaluations: (3079) 
#> Path [25] :Initial log joint density = -481715.228608 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.181e-02   2.195e-01    1.000e+00  1.000e+00      3396 -3.701e+03 -3.700e+03                   
#> Path [25] :Best Iter: [60] ELBO (-3700.467419) evaluations: (3396) 
#> Path [26] :Initial log joint density = -481617.243441 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.661e-03   1.846e-01    9.041e-01  9.041e-01      3561 -3.699e+03 -3.704e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3698.751167) evaluations: (3561) 
#> Path [27] :Initial log joint density = -481780.040886 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.780e-03   3.650e-01    3.670e-01  1.000e+00      3080 -3.708e+03 -3.715e+03                   
#> Path [27] :Best Iter: [49] ELBO (-3708.210777) evaluations: (3080) 
#> Path [28] :Initial log joint density = -481753.412444 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.416e-03   2.171e-01    1.000e+00  1.000e+00      3503 -3.697e+03 -3.704e+03                   
#> Path [28] :Best Iter: [57] ELBO (-3697.241979) evaluations: (3503) 
#> Path [29] :Initial log joint density = -485154.555876 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.396e-02   3.088e-01    1.000e+00  1.000e+00      3480 -3.696e+03 -3.706e+03                   
#> Path [29] :Best Iter: [58] ELBO (-3695.713075) evaluations: (3480) 
#> Path [30] :Initial log joint density = -481634.754658 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.883e-03   2.410e-01    1.000e+00  1.000e+00      3116 -3.707e+03 -3.703e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3703.316866) evaluations: (3116) 
#> Path [31] :Initial log joint density = -481534.485288 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.125e-02   3.355e-01    1.000e+00  1.000e+00      3169 -3.700e+03 -3.708e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3699.736811) evaluations: (3169) 
#> Path [32] :Initial log joint density = -485618.790138 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.000e-03   2.063e-01    1.000e+00  1.000e+00      3510 -3.702e+03 -3.704e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3701.928421) evaluations: (3510) 
#> Path [33] :Initial log joint density = -481557.121654 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.051e-03   2.690e-01    6.804e-01  6.804e-01      2914 -3.707e+03 -3.723e+03                   
#> Path [33] :Best Iter: [49] ELBO (-3706.863762) evaluations: (2914) 
#> Path [34] :Initial log joint density = -485051.908004 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.129e-02   3.728e-01    1.000e+00  1.000e+00      3392 -3.702e+03 -3.706e+03                   
#> Path [34] :Best Iter: [56] ELBO (-3701.716844) evaluations: (3392) 
#> Path [35] :Initial log joint density = -481666.395173 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.287e-02   1.850e-01    8.774e-01  8.774e-01      3980 -3.699e+03 -3.705e+03                   
#> Path [35] :Best Iter: [61] ELBO (-3699.279073) evaluations: (3980) 
#> Path [36] :Initial log joint density = -481370.012494 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.918e-02   2.577e-01    1.000e+00  1.000e+00      3489 -3.701e+03 -3.699e+03                   
#> Path [36] :Best Iter: [59] ELBO (-3698.869380) evaluations: (3489) 
#> Path [37] :Initial log joint density = -481472.497263 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.426e-02   2.420e-01    1.000e+00  1.000e+00      3158 -3.698e+03 -3.702e+03                   
#> Path [37] :Best Iter: [55] ELBO (-3697.732840) evaluations: (3158) 
#> Path [38] :Initial log joint density = -481531.521734 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.509e-02   3.466e-01    9.593e-01  9.593e-01      3503 -3.703e+03 -3.712e+03                   
#> Path [38] :Best Iter: [58] ELBO (-3702.999987) evaluations: (3503) 
#> Path [39] :Initial log joint density = -481593.608909 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.726e-03   1.832e-01    8.477e-01  8.477e-01      3080 -3.707e+03 -3.712e+03                   
#> Path [39] :Best Iter: [54] ELBO (-3706.922273) evaluations: (3080) 
#> Path [40] :Initial log joint density = -481466.586838 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.336e-03   2.282e-01    1.000e+00  1.000e+00      2957 -3.707e+03 -3.719e+03                   
#> Path [40] :Best Iter: [44] ELBO (-3706.614827) evaluations: (2957) 
#> Path [41] :Initial log joint density = -481447.711773 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.397e-02   2.186e-01    1.000e+00  1.000e+00      2866 -3.710e+03 -3.714e+03                   
#> Path [41] :Best Iter: [37] ELBO (-3709.965464) evaluations: (2866) 
#> Path [42] :Initial log joint density = -481879.509224 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      2.278e-03   1.989e-01    6.965e-01  6.965e-01      3745 -3.699e+03 -3.715e+03                   
#> Path [42] :Best Iter: [61] ELBO (-3698.616332) evaluations: (3745) 
#> Path [43] :Initial log joint density = -481622.394295 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.565e-03   1.590e-01    8.394e-01  8.394e-01      2961 -3.708e+03 -3.720e+03                   
#> Path [43] :Best Iter: [40] ELBO (-3708.435375) evaluations: (2961) 
#> Path [44] :Initial log joint density = -481488.688805 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.828e-03   2.527e-01    9.648e-01  9.648e-01      2669 -3.710e+03 -3.718e+03                   
#> Path [44] :Best Iter: [49] ELBO (-3709.822214) evaluations: (2669) 
#> Path [45] :Initial log joint density = -483315.963384 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.024e-02   2.202e-01    1.000e+00  1.000e+00      3607 -3.699e+03 -3.701e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3699.050640) evaluations: (3607) 
#> Path [46] :Initial log joint density = -481877.485840 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.522e-03   2.358e-01    1.000e+00  1.000e+00      2961 -3.708e+03 -3.724e+03                   
#> Path [46] :Best Iter: [50] ELBO (-3707.996840) evaluations: (2961) 
#> Path [47] :Initial log joint density = -481472.383508 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.238e-02   1.873e-01    1.000e+00  1.000e+00      3270 -3.700e+03 -3.703e+03                   
#> Path [47] :Best Iter: [56] ELBO (-3699.971115) evaluations: (3270) 
#> Path [48] :Initial log joint density = -481788.198136 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.152e-03   2.292e-01    1.000e+00  1.000e+00      2955 -3.707e+03 -3.713e+03                   
#> Path [48] :Best Iter: [46] ELBO (-3706.822298) evaluations: (2955) 
#> Path [49] :Initial log joint density = -481938.438955 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.132e-03   3.560e-01    9.130e-01  9.130e-01      3397 -3.702e+03 -3.712e+03                   
#> Path [49] :Best Iter: [56] ELBO (-3701.745714) evaluations: (3397) 
#> Path [50] :Initial log joint density = -481480.323530 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.104e-02   2.309e-01    1.000e+00  1.000e+00      3074 -3.706e+03 -3.705e+03                   
#> Path [50] :Best Iter: [55] ELBO (-3705.063488) evaluations: (3074) 
#> Finished in  13.6 seconds.
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
#> Path [1] :Initial log joint density = -433080.771678 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.634e-02   1.005e+03    2.491e-02  2.491e-02      8191 -3.324e+03 -6.933e+04                   
#> Path [1] :Best Iter: [55] ELBO (-3323.856544) evaluations: (8191) 
#> Path [2] :Initial log joint density = -431949.217341 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.040e-02   2.508e+03    2.158e-02  4.729e-02      8128 -3.327e+03 -4.726e+06                   
#> Path [2] :Best Iter: [41] ELBO (-3327.161313) evaluations: (8128) 
#> Path [3] :Initial log joint density = -432019.189686 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.315e-01   3.097e+03    1.949e-02  1.949e-02      8599 -3.325e+03 -1.630e+05                   
#> Path [3] :Best Iter: [41] ELBO (-3324.716260) evaluations: (8599) 
#> Path [4] :Initial log joint density = -431569.121141 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      7.140e-02   1.083e+03    7.080e-02  7.080e-02      7966 -3.324e+03 -1.449e+04                   
#> Path [4] :Best Iter: [52] ELBO (-3324.086912) evaluations: (7966) 
#> Path [5] :Initial log joint density = -432222.364311 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.282e-02   1.672e+03    2.115e-02  2.115e-02      8291 -3.325e+03 -7.465e+03                   
#> Path [5] :Best Iter: [49] ELBO (-3325.021823) evaluations: (8291) 
#> Path [6] :Initial log joint density = -431890.118415 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.878e-02   2.137e+03    3.086e-02  3.086e-02      8258 -3.323e+03 -2.249e+04                   
#> Path [6] :Best Iter: [40] ELBO (-3323.407594) evaluations: (8258) 
#> Path [7] :Initial log joint density = -431981.248349 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.610e-02   5.059e+03    1.371e-02  1.371e-02      8143 -3.327e+03 -1.698e+08                   
#> Path [7] :Best Iter: [43] ELBO (-3326.696557) evaluations: (8143) 
#> Path [8] :Initial log joint density = -431626.062948 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.964e-02   2.193e+03    1.686e-02  1.686e-02      8386 -3.325e+03 -8.135e+05                   
#> Path [8] :Best Iter: [47] ELBO (-3325.194569) evaluations: (8386) 
#> Path [9] :Initial log joint density = -431579.584206 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.369e-02   3.128e+03    1.744e-02  1.744e-02      8464 -3.325e+03 -1.025e+05                   
#> Path [9] :Best Iter: [43] ELBO (-3325.470401) evaluations: (8464) 
#> Path [10] :Initial log joint density = -432027.299682 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.763e-02   1.935e+03    3.841e-02  3.841e-02      7940 -3.316e+03 -6.983e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3316.080566) evaluations: (7940) 
#> Path [11] :Initial log joint density = -431774.053774 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.023e-01   5.848e+03    2.987e-02  2.987e-02      8308 -3.324e+03 -2.909e+06                   
#> Path [11] :Best Iter: [51] ELBO (-3323.559964) evaluations: (8308) 
#> Path [12] :Initial log joint density = -431665.278313 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.440e-02   2.014e+03    3.339e-02  3.339e-02      8217 -3.328e+03 -1.632e+04                   
#> Path [12] :Best Iter: [48] ELBO (-3328.035515) evaluations: (8217) 
#> Path [13] :Initial log joint density = -431943.842758 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.181e-02   9.060e+02    3.207e-02  8.216e-02      7976 -3.321e+03 -5.119e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3320.715768) evaluations: (7976) 
#> Path [14] :Initial log joint density = -432424.081777 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      8.557e-02   2.483e+03    4.174e-02  4.174e-02      8085 -3.319e+03 -1.025e+04                   
#> Path [14] :Best Iter: [58] ELBO (-3318.912969) evaluations: (8085) 
#> Path [15] :Initial log joint density = -431916.831548 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      9.646e-02   2.569e+03    2.652e-02  2.652e-02      8309 -3.322e+03 -9.905e+07                   
#> Path [15] :Best Iter: [45] ELBO (-3321.558289) evaluations: (8309) 
#> Path [16] :Initial log joint density = -431613.967463 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      7.636e-02   1.520e+03    5.066e-02  5.066e-02      8181 -3.320e+03 -3.889e+04                   
#> Path [16] :Best Iter: [49] ELBO (-3320.409466) evaluations: (8181) 
#> Path [17] :Initial log joint density = -435333.022189 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.323e-02   1.012e+03    3.920e-02  3.920e-02      8036 -3.325e+03 -1.417e+04                   
#> Path [17] :Best Iter: [53] ELBO (-3325.137235) evaluations: (8036) 
#> Path [18] :Initial log joint density = -433438.412158 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.074e-01   1.115e+03    2.409e-02  4.762e-02      8153 -3.326e+03 -9.080e+05                   
#> Path [18] :Best Iter: [48] ELBO (-3326.020790) evaluations: (8153) 
#> Path [19] :Initial log joint density = -432165.426871 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.548e-02   4.195e+03    2.888e-02  2.888e-02      8338 -3.326e+03 -4.106e+04                   
#> Path [19] :Best Iter: [47] ELBO (-3326.170924) evaluations: (8338) 
#> Path [20] :Initial log joint density = -431915.254627 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.125e-02   5.710e+03    2.569e-02  2.569e-02      8415 -3.330e+03 -3.149e+06                   
#> Path [20] :Best Iter: [44] ELBO (-3330.336218) evaluations: (8415) 
#> Path [21] :Initial log joint density = -431780.911061 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.007e-02   1.676e+03    4.320e-02  4.320e-02      8417 -3.330e+03 -5.555e+04                   
#> Path [21] :Best Iter: [44] ELBO (-3330.372516) evaluations: (8417) 
#> Path [22] :Initial log joint density = -434054.169591 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.356e-02   2.239e+03    1.205e-02  1.205e-02      8474 -3.323e+03 -2.516e+04                   
#> Path [22] :Best Iter: [49] ELBO (-3323.462665) evaluations: (8474) 
#> Path [23] :Initial log joint density = -431891.213400 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      5.913e-02   1.504e+03    2.470e-02  6.374e-02      7817 -3.319e+03 -1.006e+05                   
#> Path [23] :Best Iter: [57] ELBO (-3318.648636) evaluations: (7817) 
#> Path [24] :Initial log joint density = -434794.393636 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      4.104e-02   8.288e+02    4.014e-02  4.014e-02      8218 -3.319e+03 -1.002e+05                   
#> Path [24] :Best Iter: [57] ELBO (-3319.263082) evaluations: (8218) 
#> Path [25] :Initial log joint density = -431890.154889 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.400e-02   2.212e+03    1.871e-02  1.871e-02      8142 -3.327e+03 -1.003e+04                   
#> Path [25] :Best Iter: [48] ELBO (-3327.210831) evaluations: (8142) 
#> Path [26] :Initial log joint density = -433019.061482 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.014e-02   1.223e+03    2.025e-02  2.025e-02      8143 -3.321e+03 -7.395e+04                   
#> Path [26] :Best Iter: [56] ELBO (-3320.576506) evaluations: (8143) 
#> Path [27] :Initial log joint density = -431506.538666 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.023e-02   3.331e+03    2.239e-02  2.239e-02      8474 -3.326e+03 -1.743e+04                   
#> Path [27] :Best Iter: [43] ELBO (-3325.880243) evaluations: (8474) 
#> Path [28] :Initial log joint density = -432096.857394 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      8.816e-02   1.821e+03    5.257e-02  5.257e-02      8138 -3.320e+03 -2.152e+04                   
#> Path [28] :Best Iter: [57] ELBO (-3319.766831) evaluations: (8138) 
#> Path [29] :Initial log joint density = -433644.226219 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.936e-02   2.237e+03    2.941e-02  6.226e-02      8162 -3.322e+03 -6.358e+04                   
#> Path [29] :Best Iter: [55] ELBO (-3321.902901) evaluations: (8162) 
#> Path [30] :Initial log joint density = -431605.693849 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.048e-02   9.244e+03    1.223e-02  1.223e-02      8391 -3.329e+03 -4.184e+06                   
#> Path [30] :Best Iter: [41] ELBO (-3329.094064) evaluations: (8391) 
#> Path [31] :Initial log joint density = -432749.961475 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.335e-02   2.169e+03    3.714e-02  3.714e-02      8323 -3.326e+03 -5.081e+04                   
#> Path [31] :Best Iter: [49] ELBO (-3326.377036) evaluations: (8323) 
#> Path [32] :Initial log joint density = -431672.913398 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.208e-02   2.805e+03    2.809e-02  2.809e-02      8322 -3.329e+03 -4.868e+06                   
#> Path [32] :Best Iter: [48] ELBO (-3328.683704) evaluations: (8322) 
#> Path [33] :Initial log joint density = -431859.733699 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      6.967e-02   1.421e+03    3.466e-02  3.466e-02      8395 -3.327e+03 -1.469e+04                   
#> Path [33] :Best Iter: [56] ELBO (-3327.327368) evaluations: (8395) 
#> Path [34] :Initial log joint density = -431935.416593 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.787e-02   6.087e+03    2.306e-02  2.306e-02      8095 -3.329e+03 -1.447e+06                   
#> Path [34] :Best Iter: [41] ELBO (-3328.766023) evaluations: (8095) 
#> Path [35] :Initial log joint density = -431884.877672 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.265e-02   1.961e+03    1.472e-02  1.472e-02      7967 -3.322e+03 -1.160e+08                   
#> Path [35] :Best Iter: [56] ELBO (-3321.509849) evaluations: (7967) 
#> Path [36] :Initial log joint density = -432375.807700 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      2.325e-02   1.443e+03    2.656e-02  2.656e-02      8181 -3.318e+03 -4.224e+03                   
#> Path [36] :Best Iter: [57] ELBO (-3317.720129) evaluations: (8181) 
#> Path [37] :Initial log joint density = -431770.564516 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.903e-02   2.102e+03    1.419e-02  1.419e-02      8213 -3.326e+03 -1.356e+05                   
#> Path [37] :Best Iter: [44] ELBO (-3326.111522) evaluations: (8213) 
#> Path [38] :Initial log joint density = -432013.691962 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.016e-02   2.368e+03    1.806e-02  1.806e-02      8308 -3.327e+03 -1.276e+06                   
#> Path [38] :Best Iter: [51] ELBO (-3326.520923) evaluations: (8308) 
#> Path [39] :Initial log joint density = -431970.270966 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.262e-02   3.019e+03    1.899e-02  1.899e-02      8280 -3.328e+03 -4.563e+04                   
#> Path [39] :Best Iter: [47] ELBO (-3327.564592) evaluations: (8280) 
#> Path [40] :Initial log joint density = -433811.602858 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.816e-02   1.438e+03    4.514e-02  4.514e-02      8174 -3.318e+03 -1.947e+04                   
#> Path [40] :Best Iter: [55] ELBO (-3317.538379) evaluations: (8174) 
#> Path [41] :Initial log joint density = -432461.655685 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      4.283e-02   1.481e+03    1.923e-02  1.923e-02      8168 -3.323e+03 -3.223e+04                   
#> Path [41] :Best Iter: [53] ELBO (-3323.165222) evaluations: (8168) 
#> Path [42] :Initial log joint density = -433669.015119 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.662e-02   1.680e+03    2.845e-02  2.845e-02      8193 -3.323e+03 -1.639e+05                   
#> Path [42] :Best Iter: [56] ELBO (-3322.586945) evaluations: (8193) 
#> Path [43] :Initial log joint density = -432968.922252 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.415e-02   3.805e+03    2.685e-02  2.685e-02      8390 -3.325e+03 -3.081e+05                   
#> Path [43] :Best Iter: [46] ELBO (-3325.162424) evaluations: (8390) 
#> Path [44] :Initial log joint density = -433412.858946 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.931e-02   8.397e+02    3.209e-02  3.209e-02      8271 -3.324e+03 -6.685e+05                   
#> Path [44] :Best Iter: [55] ELBO (-3324.317125) evaluations: (8271) 
#> Path [45] :Initial log joint density = -431728.712645 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.506e-02   6.576e+03    1.349e-02  1.349e-02      8165 -3.325e+03 -5.641e+07                   
#> Path [45] :Best Iter: [44] ELBO (-3325.043214) evaluations: (8165) 
#> Path [46] :Initial log joint density = -431857.823816 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.763e-02   3.087e+03    2.643e-02  2.643e-02      8247 -3.329e+03 -1.495e+04                   
#> Path [46] :Best Iter: [53] ELBO (-3329.066212) evaluations: (8247) 
#> Path [47] :Initial log joint density = -434361.617286 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      9.612e-02   2.562e+03    1.819e-02  4.067e-02      8049 -3.321e+03 -1.070e+07                   
#> Path [47] :Best Iter: [57] ELBO (-3320.709397) evaluations: (8049) 
#> Path [48] :Initial log joint density = -431866.269737 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.167e-02   3.462e+03    1.943e-02  1.943e-02      8210 -3.326e+03 -2.014e+06                   
#> Path [48] :Best Iter: [42] ELBO (-3325.995014) evaluations: (8210) 
#> Path [49] :Initial log joint density = -436070.697148 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.161e-02   8.890e+02    3.416e-02  3.416e-02      8260 -3.321e+03 -7.345e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3320.816235) evaluations: (8260) 
#> Path [50] :Initial log joint density = -431968.407034 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.737e-02   2.390e+03    2.203e-02  2.203e-02      8436 -3.326e+03 -5.640e+04                   
#> Path [50] :Best Iter: [47] ELBO (-3326.208287) evaluations: (8436) 
#> Finished in  30.6 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -433080.771678 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.634e-02   1.005e+03    2.491e-02  2.491e-02      8191 -3.324e+03 -6.933e+04                   
#> Path [1] :Best Iter: [55] ELBO (-3323.856544) evaluations: (8191) 
#> Path [2] :Initial log joint density = -431949.217341 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.040e-02   2.508e+03    2.158e-02  4.729e-02      8128 -3.327e+03 -4.726e+06                   
#> Path [2] :Best Iter: [41] ELBO (-3327.161313) evaluations: (8128) 
#> Path [3] :Initial log joint density = -432019.189686 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.315e-01   3.097e+03    1.949e-02  1.949e-02      8599 -3.325e+03 -1.630e+05                   
#> Path [3] :Best Iter: [41] ELBO (-3324.716260) evaluations: (8599) 
#> Path [4] :Initial log joint density = -431569.121141 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      7.140e-02   1.083e+03    7.080e-02  7.080e-02      7966 -3.324e+03 -1.449e+04                   
#> Path [4] :Best Iter: [52] ELBO (-3324.086912) evaluations: (7966) 
#> Path [5] :Initial log joint density = -432222.364311 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.282e-02   1.672e+03    2.115e-02  2.115e-02      8291 -3.325e+03 -7.465e+03                   
#> Path [5] :Best Iter: [49] ELBO (-3325.021823) evaluations: (8291) 
#> Path [6] :Initial log joint density = -431890.118415 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.878e-02   2.137e+03    3.086e-02  3.086e-02      8258 -3.323e+03 -2.249e+04                   
#> Path [6] :Best Iter: [40] ELBO (-3323.407594) evaluations: (8258) 
#> Path [7] :Initial log joint density = -431981.248349 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.610e-02   5.059e+03    1.371e-02  1.371e-02      8143 -3.327e+03 -1.698e+08                   
#> Path [7] :Best Iter: [43] ELBO (-3326.696557) evaluations: (8143) 
#> Path [8] :Initial log joint density = -431626.062948 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.964e-02   2.193e+03    1.686e-02  1.686e-02      8386 -3.325e+03 -8.135e+05                   
#> Path [8] :Best Iter: [47] ELBO (-3325.194569) evaluations: (8386) 
#> Path [9] :Initial log joint density = -431579.584206 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.369e-02   3.128e+03    1.744e-02  1.744e-02      8464 -3.325e+03 -1.025e+05                   
#> Path [9] :Best Iter: [43] ELBO (-3325.470401) evaluations: (8464) 
#> Path [10] :Initial log joint density = -432027.299682 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.763e-02   1.935e+03    3.841e-02  3.841e-02      7940 -3.316e+03 -6.983e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3316.080566) evaluations: (7940) 
#> Path [11] :Initial log joint density = -431774.053774 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.023e-01   5.848e+03    2.987e-02  2.987e-02      8308 -3.324e+03 -2.909e+06                   
#> Path [11] :Best Iter: [51] ELBO (-3323.559964) evaluations: (8308) 
#> Path [12] :Initial log joint density = -431665.278313 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.440e-02   2.014e+03    3.339e-02  3.339e-02      8217 -3.328e+03 -1.632e+04                   
#> Path [12] :Best Iter: [48] ELBO (-3328.035515) evaluations: (8217) 
#> Path [13] :Initial log joint density = -431943.842758 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.181e-02   9.060e+02    3.207e-02  8.216e-02      7976 -3.321e+03 -5.119e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3320.715768) evaluations: (7976) 
#> Path [14] :Initial log joint density = -432424.081777 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      8.557e-02   2.483e+03    4.174e-02  4.174e-02      8085 -3.319e+03 -1.025e+04                   
#> Path [14] :Best Iter: [58] ELBO (-3318.912969) evaluations: (8085) 
#> Path [15] :Initial log joint density = -431916.831548 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      9.646e-02   2.569e+03    2.652e-02  2.652e-02      8309 -3.322e+03 -9.905e+07                   
#> Path [15] :Best Iter: [45] ELBO (-3321.558289) evaluations: (8309) 
#> Path [16] :Initial log joint density = -431613.967463 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      7.636e-02   1.520e+03    5.066e-02  5.066e-02      8181 -3.320e+03 -3.889e+04                   
#> Path [16] :Best Iter: [49] ELBO (-3320.409466) evaluations: (8181) 
#> Path [17] :Initial log joint density = -435333.022189 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.323e-02   1.012e+03    3.920e-02  3.920e-02      8036 -3.325e+03 -1.417e+04                   
#> Path [17] :Best Iter: [53] ELBO (-3325.137235) evaluations: (8036) 
#> Path [18] :Initial log joint density = -433438.412158 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.074e-01   1.115e+03    2.409e-02  4.762e-02      8153 -3.326e+03 -9.080e+05                   
#> Path [18] :Best Iter: [48] ELBO (-3326.020790) evaluations: (8153) 
#> Path [19] :Initial log joint density = -432165.426871 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.548e-02   4.195e+03    2.888e-02  2.888e-02      8338 -3.326e+03 -4.106e+04                   
#> Path [19] :Best Iter: [47] ELBO (-3326.170924) evaluations: (8338) 
#> Path [20] :Initial log joint density = -431915.254627 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.125e-02   5.710e+03    2.569e-02  2.569e-02      8415 -3.330e+03 -3.149e+06                   
#> Path [20] :Best Iter: [44] ELBO (-3330.336218) evaluations: (8415) 
#> Path [21] :Initial log joint density = -431780.911061 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.007e-02   1.676e+03    4.320e-02  4.320e-02      8417 -3.330e+03 -5.555e+04                   
#> Path [21] :Best Iter: [44] ELBO (-3330.372516) evaluations: (8417) 
#> Path [22] :Initial log joint density = -434054.169591 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.356e-02   2.239e+03    1.205e-02  1.205e-02      8474 -3.323e+03 -2.516e+04                   
#> Path [22] :Best Iter: [49] ELBO (-3323.462665) evaluations: (8474) 
#> Path [23] :Initial log joint density = -431891.213400 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      5.913e-02   1.504e+03    2.470e-02  6.374e-02      7817 -3.319e+03 -1.006e+05                   
#> Path [23] :Best Iter: [57] ELBO (-3318.648636) evaluations: (7817) 
#> Path [24] :Initial log joint density = -434794.393636 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      4.104e-02   8.288e+02    4.014e-02  4.014e-02      8218 -3.319e+03 -1.002e+05                   
#> Path [24] :Best Iter: [57] ELBO (-3319.263082) evaluations: (8218) 
#> Path [25] :Initial log joint density = -431890.154889 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.400e-02   2.212e+03    1.871e-02  1.871e-02      8142 -3.327e+03 -1.003e+04                   
#> Path [25] :Best Iter: [48] ELBO (-3327.210831) evaluations: (8142) 
#> Path [26] :Initial log joint density = -433019.061482 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.014e-02   1.223e+03    2.025e-02  2.025e-02      8143 -3.321e+03 -7.395e+04                   
#> Path [26] :Best Iter: [56] ELBO (-3320.576506) evaluations: (8143) 
#> Path [27] :Initial log joint density = -431506.538666 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.023e-02   3.331e+03    2.239e-02  2.239e-02      8474 -3.326e+03 -1.743e+04                   
#> Path [27] :Best Iter: [43] ELBO (-3325.880243) evaluations: (8474) 
#> Path [28] :Initial log joint density = -432096.857394 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      8.816e-02   1.821e+03    5.257e-02  5.257e-02      8138 -3.320e+03 -2.152e+04                   
#> Path [28] :Best Iter: [57] ELBO (-3319.766831) evaluations: (8138) 
#> Path [29] :Initial log joint density = -433644.226219 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.936e-02   2.237e+03    2.941e-02  6.226e-02      8162 -3.322e+03 -6.358e+04                   
#> Path [29] :Best Iter: [55] ELBO (-3321.902901) evaluations: (8162) 
#> Path [30] :Initial log joint density = -431605.693849 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.048e-02   9.244e+03    1.223e-02  1.223e-02      8391 -3.329e+03 -4.184e+06                   
#> Path [30] :Best Iter: [41] ELBO (-3329.094064) evaluations: (8391) 
#> Path [31] :Initial log joint density = -432749.961475 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.335e-02   2.169e+03    3.714e-02  3.714e-02      8323 -3.326e+03 -5.081e+04                   
#> Path [31] :Best Iter: [49] ELBO (-3326.377036) evaluations: (8323) 
#> Path [32] :Initial log joint density = -431672.913398 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.208e-02   2.805e+03    2.809e-02  2.809e-02      8322 -3.329e+03 -4.868e+06                   
#> Path [32] :Best Iter: [48] ELBO (-3328.683704) evaluations: (8322) 
#> Path [33] :Initial log joint density = -431859.733699 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      6.967e-02   1.421e+03    3.466e-02  3.466e-02      8395 -3.327e+03 -1.469e+04                   
#> Path [33] :Best Iter: [56] ELBO (-3327.327368) evaluations: (8395) 
#> Path [34] :Initial log joint density = -431935.416593 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.787e-02   6.087e+03    2.306e-02  2.306e-02      8095 -3.329e+03 -1.447e+06                   
#> Path [34] :Best Iter: [41] ELBO (-3328.766023) evaluations: (8095) 
#> Path [35] :Initial log joint density = -431884.877672 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.265e-02   1.961e+03    1.472e-02  1.472e-02      7967 -3.322e+03 -1.160e+08                   
#> Path [35] :Best Iter: [56] ELBO (-3321.509849) evaluations: (7967) 
#> Path [36] :Initial log joint density = -432375.807700 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      2.325e-02   1.443e+03    2.656e-02  2.656e-02      8181 -3.318e+03 -4.224e+03                   
#> Path [36] :Best Iter: [57] ELBO (-3317.720129) evaluations: (8181) 
#> Path [37] :Initial log joint density = -431770.564516 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      1.903e-02   2.102e+03    1.419e-02  1.419e-02      8213 -3.326e+03 -1.356e+05                   
#> Path [37] :Best Iter: [44] ELBO (-3326.111522) evaluations: (8213) 
#> Path [38] :Initial log joint density = -432013.691962 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.016e-02   2.368e+03    1.806e-02  1.806e-02      8308 -3.327e+03 -1.276e+06                   
#> Path [38] :Best Iter: [51] ELBO (-3326.520923) evaluations: (8308) 
#> Path [39] :Initial log joint density = -431970.270966 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.262e-02   3.019e+03    1.899e-02  1.899e-02      8280 -3.328e+03 -4.563e+04                   
#> Path [39] :Best Iter: [47] ELBO (-3327.564592) evaluations: (8280) 
#> Path [40] :Initial log joint density = -433811.602858 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.816e-02   1.438e+03    4.514e-02  4.514e-02      8174 -3.318e+03 -1.947e+04                   
#> Path [40] :Best Iter: [55] ELBO (-3317.538379) evaluations: (8174) 
#> Path [41] :Initial log joint density = -432461.655685 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      4.283e-02   1.481e+03    1.923e-02  1.923e-02      8168 -3.323e+03 -3.223e+04                   
#> Path [41] :Best Iter: [53] ELBO (-3323.165222) evaluations: (8168) 
#> Path [42] :Initial log joint density = -433669.015119 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.662e-02   1.680e+03    2.845e-02  2.845e-02      8193 -3.323e+03 -1.639e+05                   
#> Path [42] :Best Iter: [56] ELBO (-3322.586945) evaluations: (8193) 
#> Path [43] :Initial log joint density = -432968.922252 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      7.415e-02   3.805e+03    2.685e-02  2.685e-02      8390 -3.325e+03 -3.081e+05                   
#> Path [43] :Best Iter: [46] ELBO (-3325.162424) evaluations: (8390) 
#> Path [44] :Initial log joint density = -433412.858946 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      5.931e-02   8.397e+02    3.209e-02  3.209e-02      8271 -3.324e+03 -6.685e+05                   
#> Path [44] :Best Iter: [55] ELBO (-3324.317125) evaluations: (8271) 
#> Path [45] :Initial log joint density = -431728.712645 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      4.506e-02   6.576e+03    1.349e-02  1.349e-02      8165 -3.325e+03 -5.641e+07                   
#> Path [45] :Best Iter: [44] ELBO (-3325.043214) evaluations: (8165) 
#> Path [46] :Initial log joint density = -431857.823816 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      3.763e-02   3.087e+03    2.643e-02  2.643e-02      8247 -3.329e+03 -1.495e+04                   
#> Path [46] :Best Iter: [53] ELBO (-3329.066212) evaluations: (8247) 
#> Path [47] :Initial log joint density = -434361.617286 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      9.612e-02   2.562e+03    1.819e-02  4.067e-02      8049 -3.321e+03 -1.070e+07                   
#> Path [47] :Best Iter: [57] ELBO (-3320.709397) evaluations: (8049) 
#> Path [48] :Initial log joint density = -431866.269737 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      6.167e-02   3.462e+03    1.943e-02  1.943e-02      8210 -3.326e+03 -2.014e+06                   
#> Path [48] :Best Iter: [42] ELBO (-3325.995014) evaluations: (8210) 
#> Path [49] :Initial log joint density = -436070.697148 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.290e+05      3.161e-02   8.890e+02    3.416e-02  3.416e-02      8260 -3.321e+03 -7.345e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3320.816235) evaluations: (8260) 
#> Path [50] :Initial log joint density = -431968.407034 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.289e+05      2.737e-02   2.390e+03    2.203e-02  2.203e-02      8436 -3.326e+03 -5.640e+04                   
#> Path [50] :Best Iter: [47] ELBO (-3326.208287) evaluations: (8436) 
#> Finished in  27.6 seconds.
#> sccomp says: auto-cleanup removed 2 draw files from 'sccomp_draws_files'
# }
```
