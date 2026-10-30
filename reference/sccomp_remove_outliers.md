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
  portable = TRUE,
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

- portable:

  Logical, whether to keep the result portable by caching required draws
  in memory and removing Stan draw CSV files after fitting. Default is
  TRUE to save disk space and move needed values into memory. Set to
  FALSE to keep draw CSV files on disk.

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

Note: pH0 and FDR columns are not computed by
`sccomp_remove_outliers()`. Run
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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481339.020154 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      4.802e-02   2.267e-01    1.000e+00  1.000e+00      3800 -3.703e+03 -3.703e+03                   
#> Path [1] :Best Iter: [61] ELBO (-3702.888672) evaluations: (3800) 
#> Path [2] :Initial log joint density = -481499.039047 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.005e-02   2.711e-01    1.000e+00  1.000e+00      3772 -3.700e+03 -3.708e+03                   
#> Path [2] :Best Iter: [60] ELBO (-3700.473365) evaluations: (3772) 
#> Path [3] :Initial log joint density = -481395.907032 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.244e-02   3.486e-01    1.000e+00  1.000e+00      4059 -3.700e+03 -3.706e+03                   
#> Path [3] :Best Iter: [65] ELBO (-3700.340646) evaluations: (4059) 
#> Path [4] :Initial log joint density = -481303.271776 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.856e-02   3.504e-01    1.000e+00  1.000e+00      3549 -3.701e+03 -3.707e+03                   
#> Path [4] :Best Iter: [60] ELBO (-3700.973182) evaluations: (3549) 
#> Path [5] :Initial log joint density = -481694.592182 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.091e-02   2.136e-01    1.000e+00  1.000e+00      4040 -3.704e+03 -3.703e+03                   
#> Path [5] :Best Iter: [64] ELBO (-3702.964042) evaluations: (4040) 
#> Path [6] :Initial log joint density = -482169.194002 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      8.188e-03   2.115e-01    1.000e+00  1.000e+00      4859 -3.705e+03 -3.705e+03                   
#> Path [6] :Best Iter: [70] ELBO (-3704.938957) evaluations: (4859) 
#> Path [7] :Initial log joint density = -481509.557285 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.989e-02   4.209e-01    1.000e+00  1.000e+00      3924 -3.701e+03 -3.709e+03                   
#> Path [7] :Best Iter: [58] ELBO (-3700.703304) evaluations: (3924) 
#> Path [8] :Initial log joint density = -481691.580323 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.289e-03   1.577e-01    1.000e+00  1.000e+00      3825 -3.705e+03 -3.710e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3705.378191) evaluations: (3825) 
#> Path [9] :Initial log joint density = -482769.319344 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.031e-02   2.479e-01    1.000e+00  1.000e+00      4497 -3.700e+03 -3.705e+03                   
#> Path [9] :Best Iter: [67] ELBO (-3700.139282) evaluations: (4497) 
#> Path [10] :Initial log joint density = -481697.409105 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      7.585e-03   1.505e-01    8.224e-01  8.224e-01      3881 -3.700e+03 -3.716e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3700.259794) evaluations: (3881) 
#> Path [11] :Initial log joint density = -483556.607567 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      5.947e-03   2.339e-01    6.838e-01  6.838e-01      4340 -3.698e+03 -3.710e+03                   
#> Path [11] :Best Iter: [66] ELBO (-3698.056295) evaluations: (4340) 
#> Path [12] :Initial log joint density = -484568.892073 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.859e-03   1.662e-01    9.362e-01  9.362e-01      4300 -3.699e+03 -3.708e+03                   
#> Path [12] :Best Iter: [62] ELBO (-3698.781104) evaluations: (4300) 
#> Path [13] :Initial log joint density = -481712.030302 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.770e-03   1.647e-01    1.000e+00  1.000e+00      3852 -3.702e+03 -3.703e+03                   
#> Path [13] :Best Iter: [59] ELBO (-3702.330116) evaluations: (3852) 
#> Path [14] :Initial log joint density = -481348.755733 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.870e-02   1.972e-01    1.000e+00  1.000e+00      4163 -3.701e+03 -3.700e+03                   
#> Path [14] :Best Iter: [67] ELBO (-3699.663574) evaluations: (4163) 
#> Path [15] :Initial log joint density = -481616.467275 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.289e-02   2.148e-01    1.000e+00  1.000e+00      4442 -3.702e+03 -3.701e+03                   
#> Path [15] :Best Iter: [69] ELBO (-3700.559810) evaluations: (4442) 
#> Path [16] :Initial log joint density = -481743.320932 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.931e-02   1.751e-01    1.000e+00  1.000e+00      4310 -3.702e+03 -3.704e+03                   
#> Path [16] :Best Iter: [68] ELBO (-3702.381259) evaluations: (4310) 
#> Path [17] :Initial log joint density = -481529.629547 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.435e-02   3.584e-01    1.000e+00  1.000e+00      4519 -3.699e+03 -3.708e+03                   
#> Path [17] :Best Iter: [63] ELBO (-3699.115850) evaluations: (4519) 
#> Path [18] :Initial log joint density = -481570.379581 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.542e-02   3.398e-01    1.000e+00  1.000e+00      3778 -3.701e+03 -3.708e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3701.388010) evaluations: (3778) 
#> Path [19] :Initial log joint density = -481995.606694 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      9.716e-03   2.316e-01    1.000e+00  1.000e+00      5550 -3.697e+03 -3.708e+03                   
#> Path [19] :Best Iter: [76] ELBO (-3697.481672) evaluations: (5550) 
#> Path [20] :Initial log joint density = -481727.624433 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.463e-02   3.144e-01    4.581e-01  1.000e+00      4495 -3.699e+03 -3.708e+03                   
#> Path [20] :Best Iter: [66] ELBO (-3698.964554) evaluations: (4495) 
#> Path [21] :Initial log joint density = -481768.059657 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      4.102e-03   1.901e-01    8.228e-01  8.228e-01      5510 -3.699e+03 -3.709e+03                   
#> Path [21] :Best Iter: [78] ELBO (-3698.711855) evaluations: (5510) 
#> Path [22] :Initial log joint density = -481557.598542 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.277e-02   1.943e-01    1.000e+00  1.000e+00      4390 -3.701e+03 -3.705e+03                   
#> Path [22] :Best Iter: [61] ELBO (-3701.412357) evaluations: (4390) 
#> Path [23] :Initial log joint density = -481885.334251 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      8.478e-03   2.486e-01    1.000e+00  1.000e+00      3893 -3.703e+03 -3.712e+03                   
#> Path [23] :Best Iter: [61] ELBO (-3703.020768) evaluations: (3893) 
#> Path [24] :Initial log joint density = -481714.806795 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      3.024e-02   2.472e-01    1.000e+00  1.000e+00      3840 -3.706e+03 -3.706e+03                   
#> Path [24] :Best Iter: [63] ELBO (-3705.863757) evaluations: (3840) 
#> Path [25] :Initial log joint density = -481595.124033 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      9.785e-03   2.471e-01    1.000e+00  1.000e+00      4325 -3.704e+03 -3.711e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3704.433526) evaluations: (4325) 
#> Path [26] :Initial log joint density = -483071.418551 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.998e-03   2.621e-01    6.565e-01  6.565e-01      4375 -3.702e+03 -3.715e+03                   
#> Path [26] :Best Iter: [65] ELBO (-3702.054624) evaluations: (4375) 
#> Path [27] :Initial log joint density = -482539.858188 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.479e-03   2.767e-01    6.007e-01  6.007e-01      3867 -3.699e+03 -3.714e+03                   
#> Path [27] :Best Iter: [61] ELBO (-3698.800855) evaluations: (3867) 
#> Path [28] :Initial log joint density = -481433.820121 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.601e-02   1.704e-01    1.000e+00  1.000e+00      3594 -3.701e+03 -3.700e+03                   
#> Path [28] :Best Iter: [61] ELBO (-3699.589149) evaluations: (3594) 
#> Path [29] :Initial log joint density = -481551.629005 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.458e-02   2.461e-01    9.332e-01  9.332e-01      5081 -3.700e+03 -3.711e+03                   
#> Path [29] :Best Iter: [74] ELBO (-3700.235782) evaluations: (5081) 
#> Path [30] :Initial log joint density = -482096.729888 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      2.225e-02   9.856e-02    1.000e+00  1.000e+00      5398 -3.700e+03 -3.702e+03                   
#> Path [30] :Best Iter: [70] ELBO (-3699.620606) evaluations: (5398) 
#> Path [31] :Initial log joint density = -481422.064793 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      2.082e-02   3.121e-01    3.484e-01  1.000e+00      4700 -3.697e+03 -3.707e+03                   
#> Path [31] :Best Iter: [72] ELBO (-3696.794331) evaluations: (4700) 
#> Path [32] :Initial log joint density = -481681.834184 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.060e-03   2.119e-01    1.000e+00  1.000e+00      4799 -3.701e+03 -3.708e+03                   
#> Path [32] :Best Iter: [70] ELBO (-3700.560049) evaluations: (4799) 
#> Path [33] :Initial log joint density = -481518.742258 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      7.677e-03   1.754e-01    9.347e-01  9.347e-01      4768 -3.701e+03 -3.712e+03                   
#> Path [33] :Best Iter: [71] ELBO (-3700.733861) evaluations: (4768) 
#> Path [34] :Initial log joint density = -483675.631976 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      4.064e-03   1.877e-01    7.894e-01  7.894e-01      4599 -3.698e+03 -3.709e+03                   
#> Path [34] :Best Iter: [69] ELBO (-3697.680587) evaluations: (4599) 
#> Path [35] :Initial log joint density = -481502.609946 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      9.227e-03   1.852e-01    1.000e+00  1.000e+00      4303 -3.699e+03 -3.706e+03                   
#> Path [35] :Best Iter: [67] ELBO (-3698.770055) evaluations: (4303) 
#> Path [36] :Initial log joint density = -481458.597258 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.803e-02   2.864e-01    1.000e+00  1.000e+00      4283 -3.700e+03 -3.706e+03                   
#> Path [36] :Best Iter: [63] ELBO (-3699.761177) evaluations: (4283) 
#> Path [37] :Initial log joint density = -481835.103720 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      8.115e-03   2.195e-01    1.000e+00  1.000e+00      5255 -3.699e+03 -3.705e+03                   
#> Path [37] :Best Iter: [76] ELBO (-3698.990687) evaluations: (5255) 
#> Path [38] :Initial log joint density = -481829.248922 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.004e-02   1.810e-01    1.000e+00  1.000e+00      4909 -3.701e+03 -3.704e+03                   
#> Path [38] :Best Iter: [71] ELBO (-3700.788886) evaluations: (4909) 
#> Path [39] :Initial log joint density = -481507.533701 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.599e-02   1.763e-01    1.000e+00  1.000e+00      5019 -3.700e+03 -3.699e+03                   
#> Path [39] :Best Iter: [75] ELBO (-3699.257401) evaluations: (5019) 
#> Path [40] :Initial log joint density = -481902.666549 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.010e-02   1.632e-01    1.000e+00  1.000e+00      5795 -3.699e+03 -3.703e+03                   
#> Path [40] :Best Iter: [78] ELBO (-3698.832281) evaluations: (5795) 
#> Path [41] :Initial log joint density = -483652.694447 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.901e-03   2.299e-01    1.000e+00  1.000e+00      3239 -3.711e+03 -3.709e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3708.981766) evaluations: (3239) 
#> Path [42] :Initial log joint density = -481185.989486 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.125e-02   3.028e-01    1.000e+00  1.000e+00      3809 -3.703e+03 -3.709e+03                   
#> Path [42] :Best Iter: [63] ELBO (-3703.040311) evaluations: (3809) 
#> Path [43] :Initial log joint density = -485811.197912 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      8.468e-03   2.047e-01    1.000e+00  1.000e+00      4488 -3.699e+03 -3.702e+03                   
#> Path [43] :Best Iter: [63] ELBO (-3699.487360) evaluations: (4488) 
#> Path [44] :Initial log joint density = -481494.100348 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.527e-03   1.781e-01    1.000e+00  1.000e+00      3884 -3.699e+03 -3.704e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3699.217021) evaluations: (3884) 
#> Path [45] :Initial log joint density = -481375.245344 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.027e-02   3.304e-01    9.984e-01  9.984e-01      3835 -3.700e+03 -3.712e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3700.110921) evaluations: (3835) 
#> Path [46] :Initial log joint density = -481565.279755 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.298e-02   2.396e-01    4.986e-01  1.000e+00      4273 -3.702e+03 -3.709e+03                   
#> Path [46] :Best Iter: [62] ELBO (-3701.735854) evaluations: (4273) 
#> Path [47] :Initial log joint density = -481564.181514 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.673e-02   3.372e-01    1.000e+00  1.000e+00      3727 -3.700e+03 -3.707e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3700.425330) evaluations: (3727) 
#> Path [48] :Initial log joint density = -481450.598758 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      1.651e-02   2.597e-01    1.000e+00  1.000e+00      2915 -3.708e+03 -3.712e+03                   
#> Path [48] :Best Iter: [45] ELBO (-3708.344525) evaluations: (2915) 
#> Path [49] :Initial log joint density = -481580.771966 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      3.127e-03   2.930e-01    5.616e-01  5.616e-01      2948 -3.711e+03 -3.728e+03                   
#> Path [49] :Best Iter: [51] ELBO (-3711.039702) evaluations: (2948) 
#> Path [50] :Initial log joint density = -481328.804859 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      5.258e-03   2.011e-01    8.617e-01  8.617e-01      3628 -3.704e+03 -3.709e+03                   
#> Path [50] :Best Iter: [59] ELBO (-3703.728303) evaluations: (3628) 
#> Finished in  14.8 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 4.227 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier identification - step 1/2
#> Loading model from cache...
#> Path [1] :Initial log joint density = -428122.962980 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.254e+05      2.222e-02   2.149e-01    1.000e+00  1.000e+00      5993 -3.281e+03 -3.283e+03                   
#> Path [1] :Best Iter: [83] ELBO (-3281.199029) evaluations: (5993) 
#> Path [2] :Initial log joint density = -428107.531983 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.041e-02   2.116e-01    8.965e-01  8.965e-01      5770 -3.289e+03 -3.295e+03                   
#> Path [2] :Best Iter: [78] ELBO (-3289.494844) evaluations: (5770) 
#> Path [3] :Initial log joint density = -429547.273254 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.561e-03   2.447e-01    1.000e+00  1.000e+00      5226 -3.290e+03 -3.300e+03                   
#> Path [3] :Best Iter: [73] ELBO (-3290.007728) evaluations: (5226) 
#> Path [4] :Initial log joint density = -428544.905543 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.254e+05      1.427e-02   2.122e-01    1.000e+00  1.000e+00      7548 -3.283e+03 -3.289e+03                   
#> Path [4] :Best Iter: [95] ELBO (-3282.771753) evaluations: (7548) 
#> Path [5] :Initial log joint density = -430272.515279 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      9.975e-02   2.433e+02    8.262e-02  1.763e-01      8055 -3.291e+03 -5.262e+03                   
#> Path [5] :Best Iter: [60] ELBO (-3290.932205) evaluations: (8055) 
#> Path [6] :Initial log joint density = -427940.453802 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.254e+05      2.181e-03   2.386e-01    5.814e-01  5.814e-01      6396 -3.286e+03 -3.300e+03                   
#> Path [6] :Best Iter: [76] ELBO (-3286.426759) evaluations: (6396) 
#> Path [7] :Initial log joint density = -428164.690338 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.254e+05      1.367e-02   3.862e-01    1.000e+00  1.000e+00      5327 -3.284e+03 -3.294e+03                   
#> Path [7] :Best Iter: [77] ELBO (-3283.840384) evaluations: (5327) 
#> Path [8] :Initial log joint density = -428033.581197 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.311e-02   2.514e-01    1.000e+00  1.000e+00      5523 -3.286e+03 -3.289e+03                   
#> Path [8] :Best Iter: [78] ELBO (-3285.976091) evaluations: (5523) 
#> Path [9] :Initial log joint density = -428295.699394 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.254e+05      1.405e-02   3.346e-01    1.000e+00  1.000e+00      7189 -3.285e+03 -3.291e+03                   
#> Path [9] :Best Iter: [87] ELBO (-3285.144028) evaluations: (7189) 
#> Path [10] :Initial log joint density = -428764.864075 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      1.770e-02   2.110e-01    1.000e+00  1.000e+00      8117 -3.282e+03 -3.286e+03                   
#> Path [10] :Best Iter: [99] ELBO (-3282.143235) evaluations: (8117) 
#> Path [11] :Initial log joint density = -428665.195532 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.254e+05      8.459e-03   3.407e-01    4.491e-01  1.000e+00      5584 -3.287e+03 -3.301e+03                   
#> Path [11] :Best Iter: [80] ELBO (-3286.958834) evaluations: (5584) 
#> Path [12] :Initial log joint density = -428425.511832 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      2.160e-02   2.286e-01    1.000e+00  1.000e+00      5902 -3.288e+03 -3.292e+03                   
#> Path [12] :Best Iter: [74] ELBO (-3287.856787) evaluations: (5902) 
#> Path [13] :Initial log joint density = -428155.889473 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.254e+05      1.341e-02   2.243e-01    9.826e-01  9.826e-01      6212 -3.287e+03 -3.295e+03                   
#> Path [13] :Best Iter: [84] ELBO (-3286.754185) evaluations: (6212) 
#> Path [14] :Initial log joint density = -428251.587935 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.254e+05      1.539e-02   2.536e-01    1.000e+00  1.000e+00      6021 -3.286e+03 -3.294e+03                   
#> Path [14] :Best Iter: [83] ELBO (-3286.069204) evaluations: (6021) 
#> Path [15] :Initial log joint density = -427953.752843 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.254e+05      1.559e-02   2.783e-01    1.000e+00  1.000e+00      5110 -3.293e+03 -3.296e+03                   
#> Path [15] :Best Iter: [69] ELBO (-3292.584136) evaluations: (5110) 
#> Path [16] :Initial log joint density = -428877.057313 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.254e+05      1.186e-02   2.023e-01    1.000e+00  1.000e+00      7738 -3.286e+03 -3.290e+03                   
#> Path [16] :Best Iter: [91] ELBO (-3285.877195) evaluations: (7738) 
#> Path [17] :Initial log joint density = -427931.671337 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      1.647e-02   1.745e-01    1.000e+00  1.000e+00      5157 -3.289e+03 -3.291e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3289.424438) evaluations: (5157) 
#> Path [18] :Initial log joint density = -427866.504390 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      8.594e-03   3.036e-01    4.399e-01  1.000e+00      5697 -3.283e+03 -3.297e+03                   
#> Path [18] :Best Iter: [80] ELBO (-3283.331444) evaluations: (5697) 
#> Path [19] :Initial log joint density = -429400.271173 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      6.863e-03   2.648e-01    8.123e-01  8.123e-01      4738 -3.288e+03 -3.301e+03                   
#> Path [19] :Best Iter: [71] ELBO (-3288.309899) evaluations: (4738) 
#> Path [20] :Initial log joint density = -427775.008248 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.119e-02   1.991e-01    1.000e+00  1.000e+00      5465 -3.290e+03 -3.301e+03                   
#> Path [20] :Best Iter: [75] ELBO (-3290.145189) evaluations: (5465) 
#> Path [21] :Initial log joint density = -428204.576816 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.254e+05      1.252e-02   2.402e-01    1.000e+00  1.000e+00      5535 -3.286e+03 -3.292e+03                   
#> Path [21] :Best Iter: [78] ELBO (-3286.242970) evaluations: (5535) 
#> Path [22] :Initial log joint density = -428740.802435 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      1.169e-02   3.759e-01    1.000e+00  1.000e+00      7260 -3.288e+03 -3.293e+03                   
#> Path [22] :Best Iter: [93] ELBO (-3287.765594) evaluations: (7260) 
#> Path [23] :Initial log joint density = -428691.053039 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.397e-02   4.666e-01    1.000e+00  1.000e+00      7835 -3.283e+03 -3.293e+03                   
#> Path [23] :Best Iter: [97] ELBO (-3283.389705) evaluations: (7835) 
#> Path [24] :Initial log joint density = -428127.644664 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.254e+05      6.334e-03   2.130e-01    1.000e+00  1.000e+00      4664 -3.291e+03 -3.292e+03                   
#> Path [24] :Best Iter: [59] ELBO (-3290.870774) evaluations: (4664) 
#> Path [25] :Initial log joint density = -429941.733237 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      1.551e-01   2.043e+02    2.502e-01  2.502e-01      8229 -3.294e+03 -7.839e+03                   
#> Path [25] :Best Iter: [62] ELBO (-3293.539005) evaluations: (8229) 
#> Path [26] :Initial log joint density = -428792.567248 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.254e+05      8.833e-03   2.249e-01    8.877e-01  8.877e-01      7652 -3.287e+03 -3.293e+03                   
#> Path [26] :Best Iter: [86] ELBO (-3286.993292) evaluations: (7652) 
#> Path [27] :Initial log joint density = -428432.256001 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      1.164e-02   2.075e-01    1.000e+00  1.000e+00      7216 -3.284e+03 -3.287e+03                   
#> Path [27] :Best Iter: [91] ELBO (-3284.015574) evaluations: (7216) 
#> Path [28] :Initial log joint density = -428476.903441 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      7.620e-03   2.287e-01    7.459e-01  7.459e-01      7244 -3.283e+03 -3.299e+03                   
#> Path [28] :Best Iter: [93] ELBO (-3282.939604) evaluations: (7244) 
#> Path [29] :Initial log joint density = -428255.970281 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      7.428e-03   2.428e-01    1.000e+00  1.000e+00      8004 -3.287e+03 -3.288e+03                   
#> Path [29] :Best Iter: [96] ELBO (-3287.388562) evaluations: (8004) 
#> Path [30] :Initial log joint density = -428079.832984 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.443e-02   1.798e-01    1.000e+00  1.000e+00      5694 -3.286e+03 -3.291e+03                   
#> Path [30] :Best Iter: [78] ELBO (-3286.119216) evaluations: (5694) 
#> Path [31] :Initial log joint density = -428395.660176 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.254e+05      1.394e-02   2.788e-01    1.000e+00  1.000e+00      6551 -3.286e+03 -3.294e+03                   
#> Path [31] :Best Iter: [88] ELBO (-3286.314324) evaluations: (6551) 
#> Path [32] :Initial log joint density = -428100.907799 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.140e-02   2.415e-01    1.000e+00  1.000e+00      5384 -3.289e+03 -3.291e+03                   
#> Path [32] :Best Iter: [73] ELBO (-3289.222143) evaluations: (5384) 
#> Path [33] :Initial log joint density = -428142.647077 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      2.708e-03   2.691e-01    6.127e-01  6.127e-01      4769 -3.289e+03 -3.299e+03                   
#> Path [33] :Best Iter: [71] ELBO (-3289.012897) evaluations: (4769) 
#> Path [34] :Initial log joint density = -432300.808277 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.254e+05      1.220e-02   3.164e-01    1.000e+00  1.000e+00      5293 -3.288e+03 -3.295e+03                   
#> Path [34] :Best Iter: [74] ELBO (-3287.635650) evaluations: (5293) 
#> Path [35] :Initial log joint density = -429839.123076 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.254e+05      1.228e-02   2.825e-01    1.000e+00  1.000e+00      4486 -3.287e+03 -3.295e+03                   
#> Path [35] :Best Iter: [68] ELBO (-3287.073317) evaluations: (4486) 
#> Path [36] :Initial log joint density = -428055.951619 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      8.870e-03   2.386e-01    9.408e-01  9.408e-01      7942 -3.282e+03 -3.289e+03                   
#> Path [36] :Best Iter: [97] ELBO (-3282.045826) evaluations: (7942) 
#> Path [37] :Initial log joint density = -428057.672618 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.715e-02   2.466e-01    8.509e-01  8.509e-01      5742 -3.285e+03 -3.291e+03                   
#> Path [37] :Best Iter: [75] ELBO (-3285.140191) evaluations: (5742) 
#> Path [38] :Initial log joint density = -429494.978518 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      4.419e-03   2.093e-01    1.000e+00  1.000e+00      4806 -3.289e+03 -3.296e+03                   
#> Path [38] :Best Iter: [66] ELBO (-3288.832427) evaluations: (4806) 
#> Path [39] :Initial log joint density = -429257.052096 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.934e-03   2.357e-01    1.000e+00  1.000e+00      5282 -3.289e+03 -3.290e+03                   
#> Path [39] :Best Iter: [69] ELBO (-3289.079060) evaluations: (5282) 
#> Path [40] :Initial log joint density = -428724.116740 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.716e-02   5.297e-01    1.000e+00  1.000e+00      7973 -3.285e+03 -3.290e+03                   
#> Path [40] :Best Iter: [92] ELBO (-3284.756530) evaluations: (7973) 
#> Path [41] :Initial log joint density = -428330.290089 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.254e+05      1.600e-02   1.852e-01    1.000e+00  1.000e+00      6333 -3.287e+03 -3.287e+03                   
#> Path [41] :Best Iter: [84] ELBO (-3286.549833) evaluations: (6333) 
#> Path [42] :Initial log joint density = -428451.572529 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      8.143e-03   2.235e-01    1.000e+00  1.000e+00      7955 -3.284e+03 -3.287e+03                   
#> Path [42] :Best Iter: [99] ELBO (-3283.645953) evaluations: (7955) 
#> Path [43] :Initial log joint density = -428053.253264 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.254e+05      1.416e-02   1.892e-01    1.000e+00  1.000e+00      5808 -3.286e+03 -3.290e+03                   
#> Path [43] :Best Iter: [79] ELBO (-3286.114169) evaluations: (5808) 
#> Path [44] :Initial log joint density = -428205.112956 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.250e-03   2.241e-01    1.000e+00  1.000e+00      5157 -3.289e+03 -3.292e+03                   
#> Path [44] :Best Iter: [69] ELBO (-3289.439140) evaluations: (5157) 
#> Path [45] :Initial log joint density = -429524.680878 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      2.085e-02   2.492e-01    1.000e+00  1.000e+00      5371 -3.287e+03 -3.287e+03                   
#> Path [45] :Best Iter: [76] ELBO (-3287.325599) evaluations: (5371) 
#> Path [46] :Initial log joint density = -428006.503498 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.254e+05      1.013e-02   2.691e-01    1.000e+00  1.000e+00      4774 -3.292e+03 -3.299e+03                   
#> Path [46] :Best Iter: [66] ELBO (-3291.509284) evaluations: (4774) 
#> Path [47] :Initial log joint density = -428374.571150 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.154e-02   3.670e-01    1.000e+00  1.000e+00      8055 -3.284e+03 -3.293e+03                   
#> Path [47] :Best Iter: [98] ELBO (-3283.973141) evaluations: (8055) 
#> Path [48] :Initial log joint density = -428242.161206 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              90      -4.254e+05      1.168e-02   2.358e-01    1.000e+00  1.000e+00      6675 -3.285e+03 -3.291e+03                   
#> Path [48] :Best Iter: [83] ELBO (-3285.409605) evaluations: (6675) 
#> Path [49] :Initial log joint density = -428085.815051 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.254e+05      5.851e-03   2.080e-01    1.000e+00  1.000e+00      5648 -3.283e+03 -3.296e+03                   
#> Path [49] :Best Iter: [76] ELBO (-3283.392529) evaluations: (5648) 
#> Path [50] :Initial log joint density = -428407.901939 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.254e+05      9.377e-03   1.834e-01    1.000e+00  1.000e+00      5111 -3.289e+03 -3.290e+03                   
#> Path [50] :Best Iter: [69] ELBO (-3288.715162) evaluations: (5111) 
#> Finished in  23.7 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 21.19 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -428122.962980 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.254e+05      2.222e-02   2.149e-01    1.000e+00  1.000e+00      5993 -3.281e+03 -3.283e+03                   
#> Path [1] :Best Iter: [83] ELBO (-3281.199029) evaluations: (5993) 
#> Path [2] :Initial log joint density = -428107.531983 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.041e-02   2.116e-01    8.965e-01  8.965e-01      5770 -3.289e+03 -3.295e+03                   
#> Path [2] :Best Iter: [78] ELBO (-3289.494844) evaluations: (5770) 
#> Path [3] :Initial log joint density = -429547.273254 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.561e-03   2.447e-01    1.000e+00  1.000e+00      5226 -3.290e+03 -3.300e+03                   
#> Path [3] :Best Iter: [73] ELBO (-3290.007728) evaluations: (5226) 
#> Path [4] :Initial log joint density = -428544.905543 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.254e+05      1.427e-02   2.122e-01    1.000e+00  1.000e+00      7548 -3.283e+03 -3.289e+03                   
#> Path [4] :Best Iter: [95] ELBO (-3282.771753) evaluations: (7548) 
#> Path [5] :Initial log joint density = -430272.515279 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      9.975e-02   2.433e+02    8.262e-02  1.763e-01      8055 -3.291e+03 -5.262e+03                   
#> Path [5] :Best Iter: [60] ELBO (-3290.932205) evaluations: (8055) 
#> Path [6] :Initial log joint density = -427940.453802 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.254e+05      2.181e-03   2.386e-01    5.814e-01  5.814e-01      6396 -3.286e+03 -3.300e+03                   
#> Path [6] :Best Iter: [76] ELBO (-3286.426759) evaluations: (6396) 
#> Path [7] :Initial log joint density = -428164.690338 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.254e+05      1.367e-02   3.862e-01    1.000e+00  1.000e+00      5327 -3.284e+03 -3.294e+03                   
#> Path [7] :Best Iter: [77] ELBO (-3283.840384) evaluations: (5327) 
#> Path [8] :Initial log joint density = -428033.581197 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.311e-02   2.514e-01    1.000e+00  1.000e+00      5523 -3.286e+03 -3.289e+03                   
#> Path [8] :Best Iter: [78] ELBO (-3285.976091) evaluations: (5523) 
#> Path [9] :Initial log joint density = -428295.699394 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.254e+05      1.405e-02   3.346e-01    1.000e+00  1.000e+00      7189 -3.285e+03 -3.291e+03                   
#> Path [9] :Best Iter: [87] ELBO (-3285.144028) evaluations: (7189) 
#> Path [10] :Initial log joint density = -428764.864075 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      1.770e-02   2.110e-01    1.000e+00  1.000e+00      8117 -3.282e+03 -3.286e+03                   
#> Path [10] :Best Iter: [99] ELBO (-3282.143235) evaluations: (8117) 
#> Path [11] :Initial log joint density = -428665.195532 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.254e+05      8.459e-03   3.407e-01    4.491e-01  1.000e+00      5584 -3.287e+03 -3.301e+03                   
#> Path [11] :Best Iter: [80] ELBO (-3286.958834) evaluations: (5584) 
#> Path [12] :Initial log joint density = -428425.511832 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      2.160e-02   2.286e-01    1.000e+00  1.000e+00      5902 -3.288e+03 -3.292e+03                   
#> Path [12] :Best Iter: [74] ELBO (-3287.856787) evaluations: (5902) 
#> Path [13] :Initial log joint density = -428155.889473 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.254e+05      1.341e-02   2.243e-01    9.826e-01  9.826e-01      6212 -3.287e+03 -3.295e+03                   
#> Path [13] :Best Iter: [84] ELBO (-3286.754185) evaluations: (6212) 
#> Path [14] :Initial log joint density = -428251.587935 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.254e+05      1.539e-02   2.536e-01    1.000e+00  1.000e+00      6021 -3.286e+03 -3.294e+03                   
#> Path [14] :Best Iter: [83] ELBO (-3286.069204) evaluations: (6021) 
#> Path [15] :Initial log joint density = -427953.752843 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.254e+05      1.559e-02   2.783e-01    1.000e+00  1.000e+00      5110 -3.293e+03 -3.296e+03                   
#> Path [15] :Best Iter: [69] ELBO (-3292.584136) evaluations: (5110) 
#> Path [16] :Initial log joint density = -428877.057313 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.254e+05      1.186e-02   2.023e-01    1.000e+00  1.000e+00      7738 -3.286e+03 -3.290e+03                   
#> Path [16] :Best Iter: [91] ELBO (-3285.877195) evaluations: (7738) 
#> Path [17] :Initial log joint density = -427931.671337 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      1.647e-02   1.745e-01    1.000e+00  1.000e+00      5157 -3.289e+03 -3.291e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3289.424438) evaluations: (5157) 
#> Path [18] :Initial log joint density = -427866.504390 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      8.594e-03   3.036e-01    4.399e-01  1.000e+00      5697 -3.283e+03 -3.297e+03                   
#> Path [18] :Best Iter: [80] ELBO (-3283.331444) evaluations: (5697) 
#> Path [19] :Initial log joint density = -429400.271173 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      6.863e-03   2.648e-01    8.123e-01  8.123e-01      4738 -3.288e+03 -3.301e+03                   
#> Path [19] :Best Iter: [71] ELBO (-3288.309899) evaluations: (4738) 
#> Path [20] :Initial log joint density = -427775.008248 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.119e-02   1.991e-01    1.000e+00  1.000e+00      5465 -3.290e+03 -3.301e+03                   
#> Path [20] :Best Iter: [75] ELBO (-3290.145189) evaluations: (5465) 
#> Path [21] :Initial log joint density = -428204.576816 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.254e+05      1.252e-02   2.402e-01    1.000e+00  1.000e+00      5535 -3.286e+03 -3.292e+03                   
#> Path [21] :Best Iter: [78] ELBO (-3286.242970) evaluations: (5535) 
#> Path [22] :Initial log joint density = -428740.802435 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      1.169e-02   3.759e-01    1.000e+00  1.000e+00      7260 -3.288e+03 -3.293e+03                   
#> Path [22] :Best Iter: [93] ELBO (-3287.765594) evaluations: (7260) 
#> Path [23] :Initial log joint density = -428691.053039 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.397e-02   4.666e-01    1.000e+00  1.000e+00      7835 -3.283e+03 -3.293e+03                   
#> Path [23] :Best Iter: [97] ELBO (-3283.389705) evaluations: (7835) 
#> Path [24] :Initial log joint density = -428127.644664 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.254e+05      6.334e-03   2.130e-01    1.000e+00  1.000e+00      4664 -3.291e+03 -3.292e+03                   
#> Path [24] :Best Iter: [59] ELBO (-3290.870774) evaluations: (4664) 
#> Path [25] :Initial log joint density = -429941.733237 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      1.551e-01   2.043e+02    2.502e-01  2.502e-01      8229 -3.294e+03 -7.839e+03                   
#> Path [25] :Best Iter: [62] ELBO (-3293.539005) evaluations: (8229) 
#> Path [26] :Initial log joint density = -428792.567248 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.254e+05      8.833e-03   2.249e-01    8.877e-01  8.877e-01      7652 -3.287e+03 -3.293e+03                   
#> Path [26] :Best Iter: [86] ELBO (-3286.993292) evaluations: (7652) 
#> Path [27] :Initial log joint density = -428432.256001 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      1.164e-02   2.075e-01    1.000e+00  1.000e+00      7216 -3.284e+03 -3.287e+03                   
#> Path [27] :Best Iter: [91] ELBO (-3284.015574) evaluations: (7216) 
#> Path [28] :Initial log joint density = -428476.903441 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.254e+05      7.620e-03   2.287e-01    7.459e-01  7.459e-01      7244 -3.283e+03 -3.299e+03                   
#> Path [28] :Best Iter: [93] ELBO (-3282.939604) evaluations: (7244) 
#> Path [29] :Initial log joint density = -428255.970281 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      7.428e-03   2.428e-01    1.000e+00  1.000e+00      8004 -3.287e+03 -3.288e+03                   
#> Path [29] :Best Iter: [96] ELBO (-3287.388562) evaluations: (8004) 
#> Path [30] :Initial log joint density = -428079.832984 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.443e-02   1.798e-01    1.000e+00  1.000e+00      5694 -3.286e+03 -3.291e+03                   
#> Path [30] :Best Iter: [78] ELBO (-3286.119216) evaluations: (5694) 
#> Path [31] :Initial log joint density = -428395.660176 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.254e+05      1.394e-02   2.788e-01    1.000e+00  1.000e+00      6551 -3.286e+03 -3.294e+03                   
#> Path [31] :Best Iter: [88] ELBO (-3286.314324) evaluations: (6551) 
#> Path [32] :Initial log joint density = -428100.907799 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      1.140e-02   2.415e-01    1.000e+00  1.000e+00      5384 -3.289e+03 -3.291e+03                   
#> Path [32] :Best Iter: [73] ELBO (-3289.222143) evaluations: (5384) 
#> Path [33] :Initial log joint density = -428142.647077 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      2.708e-03   2.691e-01    6.127e-01  6.127e-01      4769 -3.289e+03 -3.299e+03                   
#> Path [33] :Best Iter: [71] ELBO (-3289.012897) evaluations: (4769) 
#> Path [34] :Initial log joint density = -432300.808277 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.254e+05      1.220e-02   3.164e-01    1.000e+00  1.000e+00      5293 -3.288e+03 -3.295e+03                   
#> Path [34] :Best Iter: [74] ELBO (-3287.635650) evaluations: (5293) 
#> Path [35] :Initial log joint density = -429839.123076 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.254e+05      1.228e-02   2.825e-01    1.000e+00  1.000e+00      4486 -3.287e+03 -3.295e+03                   
#> Path [35] :Best Iter: [68] ELBO (-3287.073317) evaluations: (4486) 
#> Path [36] :Initial log joint density = -428055.951619 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      8.870e-03   2.386e-01    9.408e-01  9.408e-01      7942 -3.282e+03 -3.289e+03                   
#> Path [36] :Best Iter: [97] ELBO (-3282.045826) evaluations: (7942) 
#> Path [37] :Initial log joint density = -428057.672618 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.254e+05      1.715e-02   2.466e-01    8.509e-01  8.509e-01      5742 -3.285e+03 -3.291e+03                   
#> Path [37] :Best Iter: [75] ELBO (-3285.140191) evaluations: (5742) 
#> Path [38] :Initial log joint density = -429494.978518 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.254e+05      4.419e-03   2.093e-01    1.000e+00  1.000e+00      4806 -3.289e+03 -3.296e+03                   
#> Path [38] :Best Iter: [66] ELBO (-3288.832427) evaluations: (4806) 
#> Path [39] :Initial log joint density = -429257.052096 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.934e-03   2.357e-01    1.000e+00  1.000e+00      5282 -3.289e+03 -3.290e+03                   
#> Path [39] :Best Iter: [69] ELBO (-3289.079060) evaluations: (5282) 
#> Path [40] :Initial log joint density = -428724.116740 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.716e-02   5.297e-01    1.000e+00  1.000e+00      7973 -3.285e+03 -3.290e+03                   
#> Path [40] :Best Iter: [92] ELBO (-3284.756530) evaluations: (7973) 
#> Path [41] :Initial log joint density = -428330.290089 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.254e+05      1.600e-02   1.852e-01    1.000e+00  1.000e+00      6333 -3.287e+03 -3.287e+03                   
#> Path [41] :Best Iter: [84] ELBO (-3286.549833) evaluations: (6333) 
#> Path [42] :Initial log joint density = -428451.572529 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      8.143e-03   2.235e-01    1.000e+00  1.000e+00      7955 -3.284e+03 -3.287e+03                   
#> Path [42] :Best Iter: [99] ELBO (-3283.645953) evaluations: (7955) 
#> Path [43] :Initial log joint density = -428053.253264 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.254e+05      1.416e-02   1.892e-01    1.000e+00  1.000e+00      5808 -3.286e+03 -3.290e+03                   
#> Path [43] :Best Iter: [79] ELBO (-3286.114169) evaluations: (5808) 
#> Path [44] :Initial log joint density = -428205.112956 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.254e+05      5.250e-03   2.241e-01    1.000e+00  1.000e+00      5157 -3.289e+03 -3.292e+03                   
#> Path [44] :Best Iter: [69] ELBO (-3289.439140) evaluations: (5157) 
#> Path [45] :Initial log joint density = -429524.680878 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.254e+05      2.085e-02   2.492e-01    1.000e+00  1.000e+00      5371 -3.287e+03 -3.287e+03                   
#> Path [45] :Best Iter: [76] ELBO (-3287.325599) evaluations: (5371) 
#> Path [46] :Initial log joint density = -428006.503498 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.254e+05      1.013e-02   2.691e-01    1.000e+00  1.000e+00      4774 -3.292e+03 -3.299e+03                   
#> Path [46] :Best Iter: [66] ELBO (-3291.509284) evaluations: (4774) 
#> Path [47] :Initial log joint density = -428374.571150 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.254e+05      2.154e-02   3.670e-01    1.000e+00  1.000e+00      8055 -3.284e+03 -3.293e+03                   
#> Path [47] :Best Iter: [98] ELBO (-3283.973141) evaluations: (8055) 
#> Path [48] :Initial log joint density = -428242.161206 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              90      -4.254e+05      1.168e-02   2.358e-01    1.000e+00  1.000e+00      6675 -3.285e+03 -3.291e+03                   
#> Path [48] :Best Iter: [83] ELBO (-3285.409605) evaluations: (6675) 
#> Path [49] :Initial log joint density = -428085.815051 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.254e+05      5.851e-03   2.080e-01    1.000e+00  1.000e+00      5648 -3.283e+03 -3.296e+03                   
#> Path [49] :Best Iter: [76] ELBO (-3283.392529) evaluations: (5648) 
#> Path [50] :Initial log joint density = -428407.901939 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.254e+05      9.377e-03   1.834e-01    1.000e+00  1.000e+00      5111 -3.289e+03 -3.290e+03                   
#> Path [50] :Best Iter: [69] ELBO (-3288.715162) evaluations: (5111) 
#> Finished in  21.3 seconds.
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
