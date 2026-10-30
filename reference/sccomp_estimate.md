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
  prior_mean = list(intercept = c(0, 2), coefficients = c(0, 1)),
  prior_overdispersion_mean_association = list(intercept = c(4, 2), slope = c(0, 2),
    standard_deviation = c(1, 0.5)),
  .sample_cell_group_pairs_to_exclude = NULL,
  output_directory = "sccomp_draws_files",
  verbose = TRUE,
  enable_loo = FALSE,
  noise_model = "multi_beta_binomial",
  exclude_mean_variability_association = FALSE,
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
  .abundance = NULL,
  exclude_priors = NULL
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

  A list with numeric length-2 vectors `intercept` and `coefficients`,
  giving the location and scale of the Student-t (3 degrees of freedom)
  priors on the CLR-scale intercept and on the remaining composition
  coefficients. The default intercept scale of 2 keeps rare cell groups
  (e.g. 0.1-1% abundance, CLR around -2 to -4) plausible a priori; the
  heavy tails let large effects escape shrinkage.

- prior_overdispersion_mean_association:

  A named list with numeric length-2 vectors `intercept`, `slope`, and
  `standard_deviation` passed to the Stan Student-t / Normal hyperpriors
  on `prec_intercept`, `prec_slope`, and `log_prec_sd`. Use `NULL` for
  package defaults. A scalar logical such as `FALSE` is not meaningful
  here and is treated as `NULL` after a message; to disable abundance
  dependence in the variability prior, use
  `exclude_mean_variability_association = TRUE` instead.

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

- exclude_mean_variability_association:

  Logical. When `TRUE`, the prior on the variability parameters does not
  depend on the abundance: the mean-variability regression is reduced to
  an intercept-only Normal (or a two-component mixture when
  `bimodal_mean_variability_association = TRUE`) while the rest of the
  hierarchical prior structure stays unchanged.

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
  FALSE to keep draw CSV files on disk. With `portable = FALSE`, CSVs
  remain for you to inspect or archive, but cmdstanr typically **still
  holds posterior draws in RAM** after fitting and summarisation
  (`fit$summary()`). The printed estimate table only calls
  `fit$summary()` on composition (`beta`, …) and variability (`alpha`,
  …), not on every saved parameter (e.g. `prec_sd`), yet cmdstanr still
  exposes all saved parameters from memory once output has been read, so
  `fit$draws(variables = "prec_sd")` can work after CSV deletion just
  like `beta`. Deleting CSVs does **not** reliably invalidate the fit in
  the same R session. Call
  [`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md)
  before deleting draw files, use `portable = TRUE` (draws cached then
  files removed), or run `incorporate_parameters_into_sccomp_object()`
  before deletion if you remove files manually.
  [`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md)
  stops with a clear error when recorded Stan output paths are missing
  unless draws were incorporated for portability as above.

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

- exclude_priors:

  **DEPRECATED**. Use `exclude_mean_variability_association` instead.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481665.172628 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.051e-02   2.098e-01    9.534e-01  9.534e-01      4765 -3.700e+03 -3.707e+03                   
#> Path [1] :Best Iter: [72] ELBO (-3700.456549) evaluations: (4765) 
#> Path [2] :Initial log joint density = -485532.045148 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      7.714e-03   2.491e-01    1.000e+00  1.000e+00      3685 -3.707e+03 -3.709e+03                   
#> Path [2] :Best Iter: [61] ELBO (-3707.345100) evaluations: (3685) 
#> Path [3] :Initial log joint density = -481337.512441 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      7.903e-03   2.476e-01    7.593e-01  7.593e-01      3243 -3.705e+03 -3.720e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3704.706504) evaluations: (3243) 
#> Path [4] :Initial log joint density = -481661.140316 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      2.598e-02   1.750e-01    1.000e+00  1.000e+00      5407 -3.699e+03 -3.702e+03                   
#> Path [4] :Best Iter: [72] ELBO (-3699.105739) evaluations: (5407) 
#> Path [5] :Initial log joint density = -481954.080376 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      9.725e-03   2.295e-01    1.000e+00  1.000e+00      3850 -3.705e+03 -3.706e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3704.523172) evaluations: (3850) 
#> Path [6] :Initial log joint density = -481155.960981 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.383e-02   2.650e-01    4.218e-01  1.000e+00      4175 -3.700e+03 -3.709e+03                   
#> Path [6] :Best Iter: [65] ELBO (-3699.694800) evaluations: (4175) 
#> Path [7] :Initial log joint density = -481479.909287 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.347e-03   2.109e-01    1.000e+00  1.000e+00      3900 -3.704e+03 -3.714e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3704.434816) evaluations: (3900) 
#> Path [8] :Initial log joint density = -481888.685393 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      8.785e-03   1.303e-01    1.000e+00  1.000e+00      5223 -3.699e+03 -3.699e+03                   
#> Path [8] :Best Iter: [78] ELBO (-3699.038875) evaluations: (5223) 
#> Path [9] :Initial log joint density = -484792.877981 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.485e-02   2.384e-01    1.000e+00  1.000e+00      4578 -3.703e+03 -3.702e+03                   
#> Path [9] :Best Iter: [71] ELBO (-3702.386034) evaluations: (4578) 
#> Path [10] :Initial log joint density = -482268.293579 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.574e-02   2.150e-01    1.000e+00  1.000e+00      5260 -3.702e+03 -3.705e+03                   
#> Path [10] :Best Iter: [77] ELBO (-3701.856216) evaluations: (5260) 
#> Path [11] :Initial log joint density = -482759.708603 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.444e-02   2.883e-01    1.000e+00  1.000e+00      3880 -3.701e+03 -3.702e+03                   
#> Path [11] :Best Iter: [58] ELBO (-3700.891131) evaluations: (3880) 
#> Path [12] :Initial log joint density = -482275.779468 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.826e-02   2.992e-01    1.000e+00  1.000e+00      3831 -3.705e+03 -3.709e+03                   
#> Path [12] :Best Iter: [61] ELBO (-3704.704452) evaluations: (3831) 
#> Path [13] :Initial log joint density = -481771.073283 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.752e-02   3.141e-01    1.000e+00  1.000e+00      4034 -3.701e+03 -3.708e+03                   
#> Path [13] :Best Iter: [65] ELBO (-3701.029886) evaluations: (4034) 
#> Path [14] :Initial log joint density = -481822.723574 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      9.057e-03   2.780e-01    9.370e-01  9.370e-01      3550 -3.706e+03 -3.719e+03                   
#> Path [14] :Best Iter: [57] ELBO (-3706.399865) evaluations: (3550) 
#> Path [15] :Initial log joint density = -481726.317253 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.237e-02   2.095e-01    1.000e+00  1.000e+00      4443 -3.700e+03 -3.701e+03                   
#> Path [15] :Best Iter: [68] ELBO (-3700.050716) evaluations: (4443) 
#> Path [16] :Initial log joint density = -482635.050559 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.603e-03   2.448e-01    6.223e-01  6.223e-01      3831 -3.704e+03 -3.717e+03                   
#> Path [16] :Best Iter: [61] ELBO (-3703.794733) evaluations: (3831) 
#> Path [17] :Initial log joint density = -481668.144499 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      4.451e-03   2.048e-01    7.295e-01  7.295e-01      4985 -3.702e+03 -3.714e+03                   
#> Path [17] :Best Iter: [73] ELBO (-3701.612098) evaluations: (4985) 
#> Path [18] :Initial log joint density = -481734.218449 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.809e-02   2.234e-01    1.000e+00  1.000e+00      4581 -3.700e+03 -3.707e+03                   
#> Path [18] :Best Iter: [66] ELBO (-3700.401038) evaluations: (4581) 
#> Path [19] :Initial log joint density = -481376.741624 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.718e-02   1.656e-01    1.000e+00  1.000e+00      3681 -3.702e+03 -3.712e+03                   
#> Path [19] :Best Iter: [57] ELBO (-3701.886471) evaluations: (3681) 
#> Path [20] :Initial log joint density = -481526.150943 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      4.654e-03   1.439e-01    7.505e-01  7.505e-01      5872 -3.697e+03 -3.708e+03                   
#> Path [20] :Best Iter: [78] ELBO (-3696.676600) evaluations: (5872) 
#> Path [21] :Initial log joint density = -481701.338872 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.194e-02   1.914e-01    1.000e+00  1.000e+00      5042 -3.699e+03 -3.702e+03                   
#> Path [21] :Best Iter: [69] ELBO (-3699.479747) evaluations: (5042) 
#> Path [22] :Initial log joint density = -481508.043598 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      7.129e-03   1.948e-01    9.842e-01  9.842e-01      3986 -3.702e+03 -3.707e+03                   
#> Path [22] :Best Iter: [63] ELBO (-3701.667981) evaluations: (3986) 
#> Path [23] :Initial log joint density = -481645.121071 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.358e-02   1.757e-01    1.000e+00  1.000e+00      3994 -3.701e+03 -3.700e+03                   
#> Path [23] :Best Iter: [66] ELBO (-3700.418650) evaluations: (3994) 
#> Path [24] :Initial log joint density = -482155.435873 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      1.663e-02   1.654e-01    1.000e+00  1.000e+00      6295 -3.701e+03 -3.700e+03                   
#> Path [24] :Best Iter: [86] ELBO (-3699.579170) evaluations: (6295) 
#> Path [25] :Initial log joint density = -481623.466923 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      3.329e-03   2.176e-01    7.477e-01  7.477e-01      4159 -3.704e+03 -3.715e+03                   
#> Path [25] :Best Iter: [58] ELBO (-3704.073153) evaluations: (4159) 
#> Path [26] :Initial log joint density = -481748.767205 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      8.693e-03   1.890e-01    1.000e+00  1.000e+00      5473 -3.699e+03 -3.711e+03                   
#> Path [26] :Best Iter: [77] ELBO (-3699.277231) evaluations: (5473) 
#> Path [27] :Initial log joint density = -481562.787043 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      8.984e-03   2.201e-01    1.000e+00  1.000e+00      5002 -3.701e+03 -3.702e+03                   
#> Path [27] :Best Iter: [71] ELBO (-3701.392273) evaluations: (5002) 
#> Path [28] :Initial log joint density = -481450.959569 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.118e-03   1.733e-01    1.000e+00  1.000e+00      4183 -3.700e+03 -3.709e+03                   
#> Path [28] :Best Iter: [65] ELBO (-3700.086361) evaluations: (4183) 
#> Path [29] :Initial log joint density = -481513.298986 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.702e-02   2.351e-01    1.000e+00  1.000e+00      3913 -3.702e+03 -3.707e+03                   
#> Path [29] :Best Iter: [62] ELBO (-3702.100861) evaluations: (3913) 
#> Path [30] :Initial log joint density = -481461.879262 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.221e-03   1.542e-01    1.000e+00  1.000e+00      3830 -3.699e+03 -3.707e+03                   
#> Path [30] :Best Iter: [60] ELBO (-3698.613308) evaluations: (3830) 
#> Path [31] :Initial log joint density = -481512.924447 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.391e-03   2.437e-01    6.447e-01  6.447e-01      4538 -3.699e+03 -3.713e+03                   
#> Path [31] :Best Iter: [70] ELBO (-3698.650232) evaluations: (4538) 
#> Path [32] :Initial log joint density = -481912.095216 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.187e-02   2.932e-01    9.185e-01  9.185e-01      4262 -3.702e+03 -3.712e+03                   
#> Path [32] :Best Iter: [68] ELBO (-3702.295602) evaluations: (4262) 
#> Path [33] :Initial log joint density = -481496.075251 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.629e-03   2.381e-01    8.705e-01  8.705e-01      3865 -3.698e+03 -3.713e+03                   
#> Path [33] :Best Iter: [62] ELBO (-3697.989977) evaluations: (3865) 
#> Path [34] :Initial log joint density = -481498.865790 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      4.148e-03   1.009e-01    1.000e+00  1.000e+00      3960 -3.703e+03 -3.716e+03                   
#> Path [34] :Best Iter: [61] ELBO (-3702.782610) evaluations: (3960) 
#> Path [35] :Initial log joint density = -481808.062395 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.103e-02   2.894e-01    9.167e-01  9.167e-01      3686 -3.699e+03 -3.712e+03                   
#> Path [35] :Best Iter: [60] ELBO (-3699.046177) evaluations: (3686) 
#> Path [36] :Initial log joint density = -481624.119640 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      7.240e-03   2.185e-01    1.000e+00  1.000e+00      3936 -3.705e+03 -3.713e+03                   
#> Path [36] :Best Iter: [62] ELBO (-3705.143007) evaluations: (3936) 
#> Path [37] :Initial log joint density = -481354.942954 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.420e-02   2.869e-01    1.000e+00  1.000e+00      4204 -3.698e+03 -3.706e+03                   
#> Path [37] :Best Iter: [65] ELBO (-3698.093428) evaluations: (4204) 
#> Path [38] :Initial log joint density = -483516.165371 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.603e-02   1.472e-01    1.000e+00  1.000e+00      4560 -3.703e+03 -3.708e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3703.406531) evaluations: (4560) 
#> Path [39] :Initial log joint density = -481291.029987 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.577e-02   3.298e-01    1.000e+00  1.000e+00      3968 -3.701e+03 -3.710e+03                   
#> Path [39] :Best Iter: [63] ELBO (-3701.171081) evaluations: (3968) 
#> Path [40] :Initial log joint density = -482063.344786 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      2.199e-02   1.350e-01    1.000e+00  1.000e+00      4834 -3.702e+03 -3.702e+03                   
#> Path [40] :Best Iter: [70] ELBO (-3701.573108) evaluations: (4834) 
#> Path [41] :Initial log joint density = -482391.145018 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.515e-02   2.543e-01    1.000e+00  1.000e+00      4727 -3.699e+03 -3.705e+03                   
#> Path [41] :Best Iter: [66] ELBO (-3698.624368) evaluations: (4727) 
#> Path [42] :Initial log joint density = -481334.514496 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.971e-03   2.026e-01    9.547e-01  9.547e-01      3866 -3.703e+03 -3.714e+03                   
#> Path [42] :Best Iter: [60] ELBO (-3703.302271) evaluations: (3866) 
#> Path [43] :Initial log joint density = -482253.828444 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.667e-02   2.423e-01    1.000e+00  1.000e+00      4340 -3.700e+03 -3.706e+03                   
#> Path [43] :Best Iter: [68] ELBO (-3699.724805) evaluations: (4340) 
#> Path [44] :Initial log joint density = -481462.198655 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.234e-02   1.984e-01    1.000e+00  1.000e+00      3972 -3.701e+03 -3.708e+03                   
#> Path [44] :Best Iter: [62] ELBO (-3700.540680) evaluations: (3972) 
#> Path [45] :Initial log joint density = -483016.398542 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.072e-02   2.710e-01    8.920e-01  8.920e-01      3629 -3.705e+03 -3.712e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3705.486977) evaluations: (3629) 
#> Path [46] :Initial log joint density = -481420.924668 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.875e-02   2.612e-01    1.000e+00  1.000e+00      5189 -3.699e+03 -3.704e+03                   
#> Path [46] :Best Iter: [77] ELBO (-3698.946528) evaluations: (5189) 
#> Path [47] :Initial log joint density = -481554.868411 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.122e-02   2.208e-01    1.000e+00  1.000e+00      3855 -3.706e+03 -3.703e+03                   
#> Path [47] :Best Iter: [63] ELBO (-3703.355362) evaluations: (3855) 
#> Path [48] :Initial log joint density = -481806.485870 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.058e-02   2.273e-01    1.000e+00  1.000e+00      4686 -3.700e+03 -3.703e+03                   
#> Path [48] :Best Iter: [68] ELBO (-3699.572805) evaluations: (4686) 
#> Path [49] :Initial log joint density = -483013.410916 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.431e-02   3.441e-01    1.000e+00  1.000e+00      3834 -3.706e+03 -3.708e+03                   
#> Path [49] :Best Iter: [61] ELBO (-3706.344479) evaluations: (3834) 
#> Path [50] :Initial log joint density = -482264.990797 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      4.038e-03   2.423e-01    6.598e-01  6.598e-01      5188 -3.700e+03 -3.713e+03                   
#> Path [50] :Best Iter: [74] ELBO (-3700.009993) evaluations: (5188) 
#> Finished in  15.2 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> sccomp says: proportion column is a proportion. The sum-constrained beta model will be used. When possible using counts is preferred as the binomial noise component is often dominating for rare groups (e.g. rare cell types).
#> sccomp says: estimation
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Warning: sccomp says: your proportion values include 0. Assuming that 0s derive from a precision threshold (e.g. deconvolution), 0s are converted to the smaller non 0 proportion value.
#> Path [1] :Initial log joint density = -2674.415490 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72       2.497e+03      1.187e-03   1.226e-02    1.000e+00  1.000e+00      4720  2.284e+03  2.279e+03                   
#> Path [1] :Best Iter: [69] ELBO (2284.460950) evaluations: (4720) 
#> Path [2] :Initial log joint density = -91.825862 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63       2.497e+03      6.397e-04   1.746e-02    4.472e-01  1.000e+00      3774  2.286e+03  2.274e+03                   
#> Path [2] :Best Iter: [56] ELBO (2286.177799) evaluations: (3774) 
#> Path [3] :Initial log joint density = -388.500784 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79       2.497e+03      4.240e-04   1.131e-02    1.000e+00  1.000e+00      5370  2.286e+03  2.283e+03                   
#> Path [3] :Best Iter: [77] ELBO (2286.197132) evaluations: (5370) 
#> Path [4] :Initial log joint density = -269.588501 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79       2.497e+03      5.821e-04   1.540e-02    1.000e+00  1.000e+00      5551  2.286e+03  2.280e+03                   
#> Path [4] :Best Iter: [77] ELBO (2285.627591) evaluations: (5551) 
#> Path [5] :Initial log joint density = -857.433598 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79       2.497e+03      1.505e-04   1.844e-02    6.023e-01  6.023e-01      5258  2.285e+03  2.273e+03                   
#> Path [5] :Best Iter: [63] ELBO (2284.716694) evaluations: (5258) 
#> Path [6] :Initial log joint density = -496.799575 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83       2.497e+03      4.812e-04   1.171e-02    8.193e-01  8.193e-01      5780  2.288e+03  2.278e+03                   
#> Path [6] :Best Iter: [77] ELBO (2287.540086) evaluations: (5780) 
#> Path [7] :Initial log joint density = -4706.090205 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73       2.497e+03      3.671e-04   1.747e-02    4.283e-01  1.000e+00      4808  2.287e+03  2.278e+03                   
#> Path [7] :Best Iter: [71] ELBO (2287.012188) evaluations: (4808) 
#> Path [8] :Initial log joint density = -731.484562 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      1.480e-03   1.859e-02    5.257e-01  1.000e+00      4273  2.287e+03  2.280e+03                   
#> Path [8] :Best Iter: [68] ELBO (2287.276418) evaluations: (4273) 
#> Path [9] :Initial log joint density = -522.406197 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      5.787e-04   2.476e-02    3.562e-01  1.000e+00      4403  2.285e+03  2.275e+03                   
#> Path [9] :Best Iter: [63] ELBO (2285.101200) evaluations: (4403) 
#> Path [10] :Initial log joint density = -525.771095 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      1.027e-03   1.811e-02    1.000e+00  1.000e+00      4148  2.284e+03  2.281e+03                   
#> Path [10] :Best Iter: [65] ELBO (2283.619124) evaluations: (4148) 
#> Path [11] :Initial log joint density = -755.279080 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70       2.497e+03      4.055e-04   1.418e-02    1.000e+00  1.000e+00      4516  2.284e+03  2.281e+03                   
#> Path [11] :Best Iter: [66] ELBO (2284.309300) evaluations: (4516) 
#> Path [12] :Initial log joint density = -400.702099 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      3.721e-04   1.503e-02    4.661e-01  1.000e+00      4334  2.286e+03  2.276e+03                   
#> Path [12] :Best Iter: [61] ELBO (2285.847427) evaluations: (4334) 
#> Path [13] :Initial log joint density = -395.740616 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78       2.497e+03      5.049e-04   1.636e-02    1.000e+00  1.000e+00      5192  2.288e+03  2.287e+03                   
#> Path [13] :Best Iter: [74] ELBO (2288.122226) evaluations: (5192) 
#> Path [14] :Initial log joint density = -1337.830125 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76       2.497e+03      1.003e-03   1.998e-02    1.000e+00  1.000e+00      5155  2.286e+03  2.282e+03                   
#> Path [14] :Best Iter: [67] ELBO (2285.742051) evaluations: (5155) 
#> Path [15] :Initial log joint density = -709.800846 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74       2.497e+03      6.311e-04   1.714e-02    1.000e+00  1.000e+00      4881  2.288e+03  2.285e+03                   
#> Path [15] :Best Iter: [72] ELBO (2287.806626) evaluations: (4881) 
#> Path [16] :Initial log joint density = -224.407642 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.497e+03      4.496e-04   1.259e-02    9.544e-01  9.544e-01      3958  2.285e+03  2.276e+03                   
#> Path [16] :Best Iter: [64] ELBO (2285.217954) evaluations: (3958) 
#> Path [17] :Initial log joint density = -322.315520 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      5.185e-04   8.748e-03    7.412e-01  7.412e-01      4464  2.288e+03  2.274e+03                   
#> Path [17] :Best Iter: [66] ELBO (2287.798911) evaluations: (4464) 
#> Path [18] :Initial log joint density = -350.195433 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      8.237e-04   1.859e-02    9.556e-01  9.556e-01      4274  2.286e+03  2.275e+03                   
#> Path [18] :Best Iter: [55] ELBO (2285.603219) evaluations: (4274) 
#> Path [19] :Initial log joint density = -1168.204862 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70       2.497e+03      4.330e-04   1.614e-02    4.534e-01  1.000e+00      4479  2.289e+03  2.277e+03                   
#> Path [19] :Best Iter: [69] ELBO (2288.723354) evaluations: (4479) 
#> Path [20] :Initial log joint density = -555.289658 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.497e+03      1.431e-03   1.669e-02    1.000e+00  1.000e+00      4106  2.285e+03  2.283e+03                   
#> Path [20] :Best Iter: [62] ELBO (2284.709341) evaluations: (4106) 
#> Path [21] :Initial log joint density = -1408.243971 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      5.421e-04   1.596e-02    1.000e+00  1.000e+00      4197  2.286e+03  2.284e+03                   
#> Path [21] :Best Iter: [60] ELBO (2286.358564) evaluations: (4197) 
#> Path [22] :Initial log joint density = -425.352357 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.497e+03      6.879e-04   2.557e-02    1.000e+00  1.000e+00      4282  2.285e+03  2.277e+03                   
#> Path [22] :Best Iter: [61] ELBO (2285.288362) evaluations: (4282) 
#> Path [23] :Initial log joint density = -473.352920 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.497e+03      1.966e-04   1.889e-02    5.372e-01  5.372e-01      4055  2.288e+03  2.277e+03                   
#> Path [23] :Best Iter: [64] ELBO (2288.109649) evaluations: (4055) 
#> Path [24] :Initial log joint density = -385.262279 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83       2.497e+03      2.341e-04   2.004e-02    5.761e-01  5.761e-01      5746  2.287e+03  2.278e+03                   
#> Path [24] :Best Iter: [71] ELBO (2287.160667) evaluations: (5746) 
#> Path [25] :Initial log joint density = -1810.246159 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      6.419e-04   1.711e-02    8.727e-01  8.727e-01      4457  2.287e+03  2.273e+03                   
#> Path [25] :Best Iter: [63] ELBO (2287.365521) evaluations: (4457) 
#> Path [26] :Initial log joint density = -483.360913 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      7.945e-04   1.798e-02    1.000e+00  1.000e+00      4465  2.283e+03  2.284e+03                   
#> Path [26] :Best Iter: [69] ELBO (2284.052812) evaluations: (4465) 
#> Path [27] :Initial log joint density = -188.905056 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64       2.497e+03      1.095e-03   1.678e-02    1.000e+00  1.000e+00      3992  2.283e+03  2.281e+03                   
#> Path [27] :Best Iter: [55] ELBO (2282.556142) evaluations: (3992) 
#> Path [28] :Initial log joint density = -3934.709072 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      3.843e-04   7.520e-03    1.000e+00  1.000e+00      4278  2.285e+03  2.277e+03                   
#> Path [28] :Best Iter: [56] ELBO (2285.487944) evaluations: (4278) 
#> Path [29] :Initial log joint density = -351.586070 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      5.400e-03   1.637e-02    1.000e+00  1.000e+00      4188  2.283e+03  2.282e+03                   
#> Path [29] :Best Iter: [57] ELBO (2282.593052) evaluations: (4188) 
#> Path [30] :Initial log joint density = -722.304370 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79       2.497e+03      5.531e-04   1.234e-02    1.000e+00  1.000e+00      5357  2.287e+03  2.285e+03                   
#> Path [30] :Best Iter: [75] ELBO (2287.055667) evaluations: (5357) 
#> Path [31] :Initial log joint density = -389.933935 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78       2.497e+03      3.214e-04   1.104e-02    8.549e-01  8.549e-01      5227  2.286e+03  2.281e+03                   
#> Path [31] :Best Iter: [76] ELBO (2286.392400) evaluations: (5227) 
#> Path [32] :Initial log joint density = -478.953762 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70       2.497e+03      7.627e-04   1.676e-02    1.000e+00  1.000e+00      4495  2.286e+03  2.285e+03                   
#> Path [32] :Best Iter: [69] ELBO (2286.270785) evaluations: (4495) 
#> Path [33] :Initial log joint density = -1061.330678 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82       2.497e+03      1.246e-03   1.631e-02    1.000e+00  1.000e+00      5754  2.288e+03  2.287e+03                   
#> Path [33] :Best Iter: [79] ELBO (2288.284462) evaluations: (5754) 
#> Path [34] :Initial log joint density = -541.397753 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73       2.497e+03      8.033e-04   1.463e-02    1.000e+00  1.000e+00      4803  2.288e+03  2.281e+03                   
#> Path [34] :Best Iter: [71] ELBO (2287.641474) evaluations: (4803) 
#> Path [35] :Initial log joint density = -1701.686326 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.497e+03      1.181e-03   2.413e-02    1.000e+00  1.000e+00      4306  2.287e+03  2.280e+03                   
#> Path [35] :Best Iter: [67] ELBO (2286.848106) evaluations: (4306) 
#> Path [36] :Initial log joint density = -1918.921900 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      1.971e-03   2.640e-02    1.000e+00  1.000e+00      4149  2.287e+03  2.281e+03                   
#> Path [36] :Best Iter: [65] ELBO (2286.987606) evaluations: (4149) 
#> Path [37] :Initial log joint density = -762.666444 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.497e+03      5.280e-04   1.347e-02    1.000e+00  1.000e+00      4365  2.284e+03  2.282e+03                   
#> Path [37] :Best Iter: [65] ELBO (2284.160463) evaluations: (4365) 
#> Path [38] :Initial log joint density = -5876.197446 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75       2.497e+03      9.547e-04   1.757e-02    1.000e+00  1.000e+00      5013  2.286e+03  2.283e+03                   
#> Path [38] :Best Iter: [73] ELBO (2285.919502) evaluations: (5013) 
#> Path [39] :Initial log joint density = -423.498393 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.497e+03      1.628e-04   2.094e-02    5.415e-01  5.415e-01      4190  2.283e+03  2.274e+03                   
#> Path [39] :Best Iter: [60] ELBO (2283.333820) evaluations: (4190) 
#> Path [40] :Initial log joint density = -3110.001971 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.497e+03      5.959e-04   1.483e-02    1.000e+00  1.000e+00      4025  2.285e+03  2.287e+03                   
#> Path [40] :Best Iter: [65] ELBO (2286.748233) evaluations: (4025) 
#> Path [41] :Initial log joint density = -378.208870 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.497e+03      8.785e-04   1.325e-02    1.000e+00  1.000e+00      4161  2.287e+03  2.285e+03                   
#> Path [41] :Best Iter: [61] ELBO (2286.900715) evaluations: (4161) 
#> Path [42] :Initial log joint density = -651.643745 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.497e+03      7.078e-04   2.320e-02    1.000e+00  1.000e+00      3965  2.288e+03  2.282e+03                   
#> Path [42] :Best Iter: [64] ELBO (2287.896303) evaluations: (3965) 
#> Path [43] :Initial log joint density = -321.818730 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63       2.497e+03      6.344e-04   1.504e-02    1.000e+00  1.000e+00      3775  2.283e+03  2.277e+03                   
#> Path [43] :Best Iter: [57] ELBO (2282.988701) evaluations: (3775) 
#> Path [44] :Initial log joint density = -879.488806 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81       2.497e+03      1.300e-03   1.705e-02    1.000e+00  1.000e+00      5962  2.288e+03  2.282e+03                   
#> Path [44] :Best Iter: [80] ELBO (2288.124345) evaluations: (5962) 
#> Path [45] :Initial log joint density = -606.423174 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.497e+03      7.278e-04   1.244e-02    1.000e+00  1.000e+00      4307  2.282e+03  2.284e+03                   
#> Path [45] :Best Iter: [68] ELBO (2283.583403) evaluations: (4307) 
#> Path [46] :Initial log joint density = -870.784118 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.497e+03      4.656e-04   1.114e-02    5.082e-01  1.000e+00      4414  2.285e+03  2.275e+03                   
#> Path [46] :Best Iter: [67] ELBO (2285.451660) evaluations: (4414) 
#> Path [47] :Initial log joint density = -286.340653 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72       2.497e+03      8.658e-04   2.723e-02    1.000e+00  1.000e+00      4623  2.287e+03  2.282e+03                   
#> Path [47] :Best Iter: [70] ELBO (2286.587290) evaluations: (4623) 
#> Path [48] :Initial log joint density = -304.018767 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55       2.497e+03      4.831e-04   2.389e-02    1.000e+00  1.000e+00      3166  2.275e+03  2.270e+03                   
#> Path [48] :Best Iter: [49] ELBO (2274.576912) evaluations: (3166) 
#> Path [49] :Initial log joint density = -860.592054 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78       2.497e+03      6.012e-04   1.064e-02    1.000e+00  1.000e+00      5458  2.286e+03  2.285e+03                   
#> Path [49] :Best Iter: [74] ELBO (2286.438870) evaluations: (5458) 
#> Path [50] :Initial log joint density = -323.338211 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76       2.497e+03      2.411e-04   2.042e-02    5.469e-01  5.469e-01      5006  2.287e+03  2.277e+03                   
#> Path [50] :Best Iter: [74] ELBO (2286.658199) evaluations: (5006) 
#> Finished in  9.8 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
