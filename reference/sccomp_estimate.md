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

  A list specifying prior knowledge about the mean distribution,
  including intercept and coefficients.

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
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481680.323498 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      7.115e-03   1.944e-01    1.000e+00  1.000e+00      4707 -3.688e+03 -3.697e+03                   
#> Path [1] :Best Iter: [64] ELBO (-3688.440187) evaluations: (4707) 
#> Path [2] :Initial log joint density = -485549.639229 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.577e-03   2.292e-01    8.636e-01  8.636e-01      3470 -3.690e+03 -3.697e+03                   
#> Path [2] :Best Iter: [57] ELBO (-3690.147890) evaluations: (3470) 
#> Path [3] :Initial log joint density = -481348.256047 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.072e-02   2.569e-01    1.000e+00  1.000e+00      3238 -3.689e+03 -3.695e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3689.267755) evaluations: (3238) 
#> Path [4] :Initial log joint density = -481676.347808 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      3.752e-03   2.049e-01    7.099e-01  7.099e-01      4727 -3.687e+03 -3.700e+03                   
#> Path [4] :Best Iter: [71] ELBO (-3686.703682) evaluations: (4727) 
#> Path [5] :Initial log joint density = -481967.121980 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.910e-03   3.544e-01    4.279e-01  1.000e+00      4113 -3.689e+03 -3.693e+03                   
#> Path [5] :Best Iter: [65] ELBO (-3688.833790) evaluations: (4113) 
#> Path [6] :Initial log joint density = -481170.130434 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.446e-03   1.734e-01    7.009e-01  7.009e-01      4373 -3.684e+03 -3.699e+03                   
#> Path [6] :Best Iter: [65] ELBO (-3683.540377) evaluations: (4373) 
#> Path [7] :Initial log joint density = -481498.016333 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.751e-02   3.233e-01    9.579e-01  9.579e-01      3775 -3.682e+03 -3.697e+03                   
#> Path [7] :Best Iter: [61] ELBO (-3682.129157) evaluations: (3775) 
#> Path [8] :Initial log joint density = -481907.104408 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      6.557e-03   2.009e-01    1.000e+00  1.000e+00      4559 -3.688e+03 -3.695e+03                   
#> Path [8] :Best Iter: [62] ELBO (-3687.837600) evaluations: (4559) 
#> Path [9] :Initial log joint density = -484806.760311 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.248e-02   3.369e-01    1.000e+00  1.000e+00      4153 -3.686e+03 -3.693e+03                   
#> Path [9] :Best Iter: [64] ELBO (-3686.209128) evaluations: (4153) 
#> Path [10] :Initial log joint density = -482283.410683 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.144e-02   2.372e-01    1.000e+00  1.000e+00      4825 -3.687e+03 -3.694e+03                   
#> Path [10] :Best Iter: [68] ELBO (-3687.159716) evaluations: (4825) 
#> Path [11] :Initial log joint density = -482779.782031 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      4.686e-03   1.542e-01    1.000e+00  1.000e+00      3327 -3.688e+03 -3.699e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3688.413698) evaluations: (3327) 
#> Path [12] :Initial log joint density = -482289.291925 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.724e-03   1.427e-01    1.000e+00  1.000e+00      3867 -3.690e+03 -3.699e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3690.042412) evaluations: (3867) 
#> Path [13] :Initial log joint density = -481784.758823 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      5.378e-03   2.203e-01    1.000e+00  1.000e+00      3582 -3.690e+03 -3.699e+03                   
#> Path [13] :Best Iter: [58] ELBO (-3690.181761) evaluations: (3582) 
#> Path [14] :Initial log joint density = -481837.913043 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      2.193e-02   2.722e-01    1.000e+00  1.000e+00      3900 -3.690e+03 -3.692e+03                   
#> Path [14] :Best Iter: [64] ELBO (-3689.969373) evaluations: (3900) 
#> Path [15] :Initial log joint density = -481745.695800 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      8.886e-03   1.572e-01    8.152e-01  8.152e-01      4440 -3.682e+03 -3.695e+03                   
#> Path [15] :Best Iter: [68] ELBO (-3682.026097) evaluations: (4440) 
#> Path [16] :Initial log joint density = -482648.892587 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.663e-03   1.793e-01    9.530e-01  9.530e-01      3554 -3.691e+03 -3.697e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3690.601581) evaluations: (3554) 
#> Path [17] :Initial log joint density = -481684.104535 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.088e-02   3.413e-01    1.000e+00  1.000e+00      3834 -3.688e+03 -3.693e+03                   
#> Path [17] :Best Iter: [62] ELBO (-3688.411483) evaluations: (3834) 
#> Path [18] :Initial log joint density = -481752.942477 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      7.685e-03   1.882e-01    1.000e+00  1.000e+00      3177 -3.693e+03 -3.699e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3693.499761) evaluations: (3177) 
#> Path [19] :Initial log joint density = -481391.232699 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.613e-03   2.549e-01    1.000e+00  1.000e+00      3349 -3.689e+03 -3.695e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3689.177490) evaluations: (3349) 
#> Path [20] :Initial log joint density = -481537.395125 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.273e-02   2.245e-01    1.000e+00  1.000e+00      3748 -3.687e+03 -3.688e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3687.225094) evaluations: (3748) 
#> Path [21] :Initial log joint density = -481718.530282 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      9.546e-03   2.744e-01    1.000e+00  1.000e+00      4627 -3.686e+03 -3.696e+03                   
#> Path [21] :Best Iter: [71] ELBO (-3685.821965) evaluations: (4627) 
#> Path [22] :Initial log joint density = -481520.673888 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.766e-02   2.376e-01    1.000e+00  1.000e+00      3630 -3.688e+03 -3.690e+03                   
#> Path [22] :Best Iter: [60] ELBO (-3687.652064) evaluations: (3630) 
#> Path [23] :Initial log joint density = -481659.460370 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.530e-02   1.988e-01    1.000e+00  1.000e+00      4088 -3.684e+03 -3.687e+03                   
#> Path [23] :Best Iter: [64] ELBO (-3684.450909) evaluations: (4088) 
#> Path [24] :Initial log joint density = -482175.906741 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.194e-02   1.920e-01    1.000e+00  1.000e+00      5638 -3.686e+03 -3.687e+03                   
#> Path [24] :Best Iter: [76] ELBO (-3685.852028) evaluations: (5638) 
#> Path [25] :Initial log joint density = -481638.216180 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      3.193e-03   2.873e-01    7.142e-01  7.142e-01      3209 -3.694e+03 -3.704e+03                   
#> Path [25] :Best Iter: [49] ELBO (-3693.605446) evaluations: (3209) 
#> Path [26] :Initial log joint density = -481769.739470 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.180e-02   2.293e-01    1.000e+00  1.000e+00      4603 -3.688e+03 -3.689e+03                   
#> Path [26] :Best Iter: [70] ELBO (-3687.788126) evaluations: (4603) 
#> Path [27] :Initial log joint density = -481580.217388 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.412e-02   2.067e-01    1.000e+00  1.000e+00      4332 -3.686e+03 -3.689e+03                   
#> Path [27] :Best Iter: [61] ELBO (-3685.864165) evaluations: (4332) 
#> Path [28] :Initial log joint density = -481464.654121 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.122e-02   2.840e-01    9.688e-01  9.688e-01      4297 -3.683e+03 -3.695e+03                   
#> Path [28] :Best Iter: [67] ELBO (-3682.767587) evaluations: (4297) 
#> Path [29] :Initial log joint density = -481526.639144 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.890e-03   1.610e-01    1.000e+00  1.000e+00      3201 -3.696e+03 -3.695e+03                   
#> Path [29] :Best Iter: [57] ELBO (-3695.324094) evaluations: (3201) 
#> Path [30] :Initial log joint density = -481475.101418 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      6.825e-03   1.431e-01    1.000e+00  1.000e+00      4208 -3.684e+03 -3.691e+03                   
#> Path [30] :Best Iter: [60] ELBO (-3684.127762) evaluations: (4208) 
#> Path [31] :Initial log joint density = -481524.422034 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.354e-02   3.466e-01    1.000e+00  1.000e+00      4115 -3.687e+03 -3.696e+03                   
#> Path [31] :Best Iter: [66] ELBO (-3686.666482) evaluations: (4115) 
#> Path [32] :Initial log joint density = -481928.770283 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.358e-02   3.037e-01    9.414e-01  9.414e-01      4216 -3.685e+03 -3.698e+03                   
#> Path [32] :Best Iter: [67] ELBO (-3684.819612) evaluations: (4216) 
#> Path [33] :Initial log joint density = -481511.602527 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      7.450e-03   2.229e-01    1.000e+00  1.000e+00      4111 -3.686e+03 -3.692e+03                   
#> Path [33] :Best Iter: [62] ELBO (-3685.650715) evaluations: (4111) 
#> Path [34] :Initial log joint density = -481514.811668 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.497e-02   2.041e-01    1.000e+00  1.000e+00      3940 -3.683e+03 -3.689e+03                   
#> Path [34] :Best Iter: [62] ELBO (-3682.856791) evaluations: (3940) 
#> Path [35] :Initial log joint density = -481822.452584 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.187e-02   2.104e-01    1.000e+00  1.000e+00      4493 -3.684e+03 -3.690e+03                   
#> Path [35] :Best Iter: [66] ELBO (-3684.409184) evaluations: (4493) 
#> Path [36] :Initial log joint density = -481640.109983 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.280e-03   1.497e-01    1.000e+00  1.000e+00      3587 -3.686e+03 -3.694e+03                   
#> Path [36] :Best Iter: [58] ELBO (-3686.134684) evaluations: (3587) 
#> Path [37] :Initial log joint density = -481369.066838 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.748e-02   3.083e-01    3.826e-01  1.000e+00      4419 -3.685e+03 -3.695e+03                   
#> Path [37] :Best Iter: [68] ELBO (-3685.244298) evaluations: (4419) 
#> Path [38] :Initial log joint density = -483532.928560 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      6.545e-03   2.217e-01    1.000e+00  1.000e+00      3551 -3.689e+03 -3.691e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3688.685930) evaluations: (3551) 
#> Path [39] :Initial log joint density = -481303.997263 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      1.586e-02   3.164e-01    8.844e-01  8.844e-01      2872 -3.696e+03 -3.711e+03                   
#> Path [39] :Best Iter: [51] ELBO (-3695.703595) evaluations: (2872) 
#> Path [40] :Initial log joint density = -482085.975130 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.076e-02   2.740e-01    1.000e+00  1.000e+00      4221 -3.689e+03 -3.691e+03                   
#> Path [40] :Best Iter: [58] ELBO (-3689.273039) evaluations: (4221) 
#> Path [41] :Initial log joint density = -482403.235504 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      7.757e-03   2.283e-01    1.000e+00  1.000e+00      4428 -3.684e+03 -3.686e+03                   
#> Path [41] :Best Iter: [66] ELBO (-3683.677061) evaluations: (4428) 
#> Path [42] :Initial log joint density = -481352.513840 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.585e-03   1.789e-01    1.000e+00  1.000e+00      4183 -3.685e+03 -3.688e+03                   
#> Path [42] :Best Iter: [60] ELBO (-3685.326011) evaluations: (4183) 
#> Path [43] :Initial log joint density = -482268.473169 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.718e-02   1.770e-01    1.000e+00  1.000e+00      4943 -3.682e+03 -3.690e+03                   
#> Path [43] :Best Iter: [69] ELBO (-3681.736432) evaluations: (4943) 
#> Path [44] :Initial log joint density = -481475.979482 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.410e-02   2.675e-01    9.881e-01  9.881e-01      3390 -3.688e+03 -3.699e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3688.481765) evaluations: (3390) 
#> Path [45] :Initial log joint density = -483030.540686 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      8.073e-03   2.751e-01    9.927e-01  9.927e-01      3093 -3.694e+03 -3.704e+03                   
#> Path [45] :Best Iter: [47] ELBO (-3694.219518) evaluations: (3093) 
#> Path [46] :Initial log joint density = -481435.807767 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      4.037e-03   2.144e-01    7.715e-01  7.715e-01      5408 -3.680e+03 -3.690e+03                   
#> Path [46] :Best Iter: [76] ELBO (-3680.126891) evaluations: (5408) 
#> Path [47] :Initial log joint density = -481573.142999 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.117e-02   4.261e-01    1.000e+00  1.000e+00      3191 -3.693e+03 -3.699e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3693.238153) evaluations: (3191) 
#> Path [48] :Initial log joint density = -481824.909464 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.730e-02   1.832e-01    1.000e+00  1.000e+00      4709 -3.684e+03 -3.690e+03                   
#> Path [48] :Best Iter: [65] ELBO (-3683.506725) evaluations: (4709) 
#> Path [49] :Initial log joint density = -483030.561114 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.683e-03   2.207e-01    1.000e+00  1.000e+00      3504 -3.688e+03 -3.695e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3688.282182) evaluations: (3504) 
#> Path [50] :Initial log joint density = -482285.755125 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.252e-02   2.210e-01    1.000e+00  1.000e+00      4674 -3.686e+03 -3.686e+03                   
#> Path [50] :Best Iter: [69] ELBO (-3686.204158) evaluations: (4674) 
#> Finished in  15.5 seconds.
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
#> Path [1] :Initial log joint density = -2692.351653 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.484e+03      4.035e-04   1.936e-02    2.988e-01  1.000e+00      4403  2.303e+03  2.293e+03                   
#> Path [1] :Best Iter: [63] ELBO (2303.372823) evaluations: (4403) 
#> Path [2] :Initial log joint density = -106.593791 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.484e+03      3.458e-04   1.684e-02    6.872e-01  6.872e-01      3958  2.302e+03  2.290e+03                   
#> Path [2] :Best Iter: [64] ELBO (2301.846403) evaluations: (3958) 
#> Path [3] :Initial log joint density = -407.463780 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80       2.484e+03      1.452e-03   1.124e-02    1.000e+00  1.000e+00      5495  2.302e+03  2.301e+03                   
#> Path [3] :Best Iter: [74] ELBO (2302.454576) evaluations: (5495) 
#> Path [4] :Initial log joint density = -290.548562 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      4.899e-04   1.598e-02    4.926e-01  1.000e+00      5376  2.305e+03  2.297e+03                   
#> Path [4] :Best Iter: [67] ELBO (2305.104411) evaluations: (5376) 
#> Path [5] :Initial log joint density = -871.093250 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      1.179e-03   1.161e-02    1.000e+00  1.000e+00      5180  2.305e+03  2.304e+03                   
#> Path [5] :Best Iter: [63] ELBO (2305.398511) evaluations: (5180) 
#> Path [6] :Initial log joint density = -510.965226 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81       2.484e+03      6.616e-04   1.142e-02    1.000e+00  1.000e+00      5515  2.306e+03  2.304e+03                   
#> Path [6] :Best Iter: [80] ELBO (2306.045580) evaluations: (5515) 
#> Path [7] :Initial log joint density = -4722.584330 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75       2.484e+03      8.513e-04   1.141e-02    1.000e+00  1.000e+00      5147  2.303e+03  2.303e+03                   
#> Path [7] :Best Iter: [71] ELBO (2303.126889) evaluations: (5147) 
#> Path [8] :Initial log joint density = -754.387532 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70       2.484e+03      3.930e-04   1.393e-02    3.377e-01  1.000e+00      4384  2.305e+03  2.297e+03                   
#> Path [8] :Best Iter: [68] ELBO (2305.359244) evaluations: (4384) 
#> Path [9] :Initial log joint density = -537.419811 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      7.800e-04   1.554e-02    1.000e+00  1.000e+00      4368  2.304e+03  2.303e+03                   
#> Path [9] :Best Iter: [65] ELBO (2304.101460) evaluations: (4368) 
#> Path [10] :Initial log joint density = -548.084941 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.484e+03      1.097e-03   2.331e-02    1.000e+00  1.000e+00      3994  2.302e+03  2.297e+03                   
#> Path [10] :Best Iter: [65] ELBO (2301.897403) evaluations: (3994) 
#> Path [11] :Initial log joint density = -773.647807 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71       2.484e+03      5.682e-04   1.566e-02    1.000e+00  1.000e+00      4685  2.301e+03  2.297e+03                   
#> Path [11] :Best Iter: [70] ELBO (2301.497685) evaluations: (4685) 
#> Path [12] :Initial log joint density = -417.660624 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      5.072e-04   2.656e-02    1.000e+00  1.000e+00      4369  2.302e+03  2.294e+03                   
#> Path [12] :Best Iter: [65] ELBO (2302.157403) evaluations: (4369) 
#> Path [13] :Initial log joint density = -409.107791 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      9.189e-04   2.270e-02    1.000e+00  1.000e+00      5119  2.305e+03  2.301e+03                   
#> Path [13] :Best Iter: [71] ELBO (2305.395353) evaluations: (5119) 
#> Path [14] :Initial log joint density = -1352.089967 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78       2.484e+03      8.133e-04   2.938e-02    1.000e+00  1.000e+00      5416  2.307e+03  2.298e+03                   
#> Path [14] :Best Iter: [77] ELBO (2306.688956) evaluations: (5416) 
#> Path [15] :Initial log joint density = -722.583235 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74       2.484e+03      3.998e-04   2.635e-02    6.451e-01  6.451e-01      4849  2.305e+03  2.291e+03                   
#> Path [15] :Best Iter: [72] ELBO (2304.580600) evaluations: (4849) 
#> Path [16] :Initial log joint density = -241.591763 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.484e+03      1.994e-04   6.564e-03    4.296e-01  9.340e-01      4116  2.303e+03  2.291e+03                   
#> Path [16] :Best Iter: [64] ELBO (2302.771627) evaluations: (4116) 
#> Path [17] :Initial log joint density = -337.978689 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      3.680e-04   1.168e-02    1.000e+00  1.000e+00      4389  2.304e+03  2.293e+03                   
#> Path [17] :Best Iter: [60] ELBO (2304.056250) evaluations: (4389) 
#> Path [18] :Initial log joint density = -365.208534 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.484e+03      4.740e-04   1.224e-02    1.000e+00  1.000e+00      4096  2.304e+03  2.303e+03                   
#> Path [18] :Best Iter: [63] ELBO (2303.586315) evaluations: (4096) 
#> Path [19] :Initial log joint density = -1178.105879 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.484e+03      1.133e-03   1.265e-02    1.000e+00  1.000e+00      4401  2.304e+03  2.305e+03                   
#> Path [19] :Best Iter: [69] ELBO (2305.269012) evaluations: (4401) 
#> Path [20] :Initial log joint density = -575.129397 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.484e+03      1.017e-03   1.765e-02    1.000e+00  1.000e+00      4018  2.303e+03  2.303e+03                   
#> Path [20] :Best Iter: [62] ELBO (2303.056698) evaluations: (4018) 
#> Path [21] :Initial log joint density = -1424.892115 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69       2.484e+03      2.471e-04   1.204e-02    3.242e-01  1.000e+00      4342  2.304e+03  2.297e+03                   
#> Path [21] :Best Iter: [62] ELBO (2304.198698) evaluations: (4342) 
#> Path [22] :Initial log joint density = -442.139090 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66       2.484e+03      4.851e-04   2.039e-02    9.546e-01  9.546e-01      4051  2.303e+03  2.291e+03                   
#> Path [22] :Best Iter: [61] ELBO (2302.902724) evaluations: (4051) 
#> Path [23] :Initial log joint density = -489.087476 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.484e+03      2.002e-04   1.366e-02    6.332e-01  6.332e-01      4018  2.305e+03  2.294e+03                   
#> Path [23] :Best Iter: [63] ELBO (2304.525682) evaluations: (4018) 
#> Path [24] :Initial log joint density = -400.001010 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83       2.484e+03      4.393e-04   8.624e-03    1.000e+00  1.000e+00      5760  2.305e+03  2.297e+03                   
#> Path [24] :Best Iter: [80] ELBO (2305.098598) evaluations: (5760) 
#> Path [25] :Initial log joint density = -1827.168706 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      1.179e-03   1.824e-02    1.000e+00  1.000e+00      4346  2.306e+03  2.305e+03                   
#> Path [25] :Best Iter: [63] ELBO (2306.064689) evaluations: (4346) 
#> Path [26] :Initial log joint density = -501.389258 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      4.160e-04   1.633e-02    1.000e+00  1.000e+00      4367  2.301e+03  2.295e+03                   
#> Path [26] :Best Iter: [56] ELBO (2301.286502) evaluations: (4367) 
#> Path [27] :Initial log joint density = -198.874091 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63       2.484e+03      4.463e-04   2.114e-02    4.092e-01  1.000e+00      3896  2.300e+03  2.291e+03                   
#> Path [27] :Best Iter: [57] ELBO (2300.440220) evaluations: (3896) 
#> Path [28] :Initial log joint density = -3952.104421 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      5.490e-04   1.493e-02    1.000e+00  1.000e+00      4243  2.304e+03  2.301e+03                   
#> Path [28] :Best Iter: [67] ELBO (2304.162495) evaluations: (4243) 
#> Path [29] :Initial log joint density = -368.250359 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65       2.484e+03      1.823e-03   1.116e-02    1.000e+00  1.000e+00      4047  2.302e+03  2.302e+03                   
#> Path [29] :Best Iter: [56] ELBO (2302.132467) evaluations: (4047) 
#> Path [30] :Initial log joint density = -733.880284 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      8.214e-04   1.946e-02    1.000e+00  1.000e+00      5087  2.304e+03  2.303e+03                   
#> Path [30] :Best Iter: [75] ELBO (2303.968565) evaluations: (5087) 
#> Path [31] :Initial log joint density = -408.509628 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      3.785e-04   1.299e-02    7.775e-01  7.775e-01      5157  2.306e+03  2.296e+03                   
#> Path [31] :Best Iter: [74] ELBO (2305.676237) evaluations: (5157) 
#> Path [32] :Initial log joint density = -495.290416 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70       2.484e+03      2.533e-04   2.363e-02    5.589e-01  5.589e-01      4495  2.303e+03  2.291e+03                   
#> Path [32] :Best Iter: [63] ELBO (2302.523563) evaluations: (4495) 
#> Path [33] :Initial log joint density = -1078.837384 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83       2.484e+03      1.645e-03   2.098e-02    1.000e+00  1.000e+00      5947  2.307e+03  2.300e+03                   
#> Path [33] :Best Iter: [81] ELBO (2306.744441) evaluations: (5947) 
#> Path [34] :Initial log joint density = -555.785261 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      6.211e-04   1.511e-02    1.000e+00  1.000e+00      5182  2.303e+03  2.298e+03                   
#> Path [34] :Best Iter: [76] ELBO (2303.471428) evaluations: (5182) 
#> Path [35] :Initial log joint density = -1720.893990 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71       2.484e+03      5.113e-04   1.531e-02    8.962e-01  8.962e-01      4537  2.306e+03  2.294e+03                   
#> Path [35] :Best Iter: [70] ELBO (2306.434503) evaluations: (4537) 
#> Path [36] :Initial log joint density = -1939.804211 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      8.138e-04   1.125e-02    1.000e+00  1.000e+00      4304  2.304e+03  2.302e+03                   
#> Path [36] :Best Iter: [65] ELBO (2304.229480) evaluations: (4304) 
#> Path [37] :Initial log joint density = -777.826662 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71       2.484e+03      8.217e-04   1.287e-02    1.000e+00  1.000e+00      4598  2.302e+03  2.299e+03                   
#> Path [37] :Best Iter: [63] ELBO (2302.199904) evaluations: (4598) 
#> Path [38] :Initial log joint density = -5892.802493 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74       2.484e+03      8.793e-04   1.430e-02    1.000e+00  1.000e+00      4963  2.305e+03  2.303e+03                   
#> Path [38] :Best Iter: [67] ELBO (2305.032421) evaluations: (4963) 
#> Path [39] :Initial log joint density = -439.164595 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64       2.484e+03      7.700e-04   1.343e-02    1.000e+00  1.000e+00      3836  2.305e+03  2.301e+03                   
#> Path [39] :Best Iter: [60] ELBO (2305.491547) evaluations: (3836) 
#> Path [40] :Initial log joint density = -3125.359820 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63       2.484e+03      3.735e-04   9.913e-03    1.000e+00  1.000e+00      3896  2.303e+03  2.299e+03                   
#> Path [40] :Best Iter: [60] ELBO (2302.980233) evaluations: (3896) 
#> Path [41] :Initial log joint density = -399.830650 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.484e+03      9.443e-04   7.172e-03    1.000e+00  1.000e+00      4335  2.302e+03  2.298e+03                   
#> Path [41] :Best Iter: [57] ELBO (2302.208949) evaluations: (4335) 
#> Path [42] :Initial log joint density = -672.943573 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.484e+03      6.275e-04   1.932e-02    1.000e+00  1.000e+00      4219  2.302e+03  2.298e+03                   
#> Path [42] :Best Iter: [65] ELBO (2302.192908) evaluations: (4219) 
#> Path [43] :Initial log joint density = -341.514703 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64       2.484e+03      1.122e-03   1.319e-02    1.000e+00  1.000e+00      3839  2.301e+03  2.300e+03                   
#> Path [43] :Best Iter: [60] ELBO (2300.879781) evaluations: (3839) 
#> Path [44] :Initial log joint density = -894.139852 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81       2.484e+03      4.536e-04   1.403e-02    8.550e-01  8.550e-01      5913  2.305e+03  2.299e+03                   
#> Path [44] :Best Iter: [79] ELBO (2305.136023) evaluations: (5913) 
#> Path [45] :Initial log joint density = -622.175635 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67       2.484e+03      2.898e-04   8.320e-03    4.516e-01  1.000e+00      4178  2.303e+03  2.295e+03                   
#> Path [45] :Best Iter: [66] ELBO (2302.878275) evaluations: (4178) 
#> Path [46] :Initial log joint density = -887.295946 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68       2.484e+03      1.045e-03   1.836e-02    9.639e-01  9.639e-01      4355  2.304e+03  2.297e+03                   
#> Path [46] :Best Iter: [67] ELBO (2303.880886) evaluations: (4355) 
#> Path [47] :Initial log joint density = -300.626687 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73       2.484e+03      1.331e-03   1.794e-02    1.000e+00  1.000e+00      4680  2.306e+03  2.300e+03                   
#> Path [47] :Best Iter: [69] ELBO (2305.702461) evaluations: (4680) 
#> Path [48] :Initial log joint density = -319.287628 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54       2.484e+03      3.616e-04   2.129e-02    9.369e-01  9.369e-01      3047  2.295e+03  2.261e+03                   
#> Path [48] :Best Iter: [49] ELBO (2294.930407) evaluations: (3047) 
#> Path [49] :Initial log joint density = -876.773197 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75       2.484e+03      1.068e-03   1.261e-02    1.000e+00  1.000e+00      5140  2.304e+03  2.303e+03                   
#> Path [49] :Best Iter: [71] ELBO (2304.036677) evaluations: (5140) 
#> Path [50] :Initial log joint density = -343.019780 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77       2.484e+03      6.555e-04   1.142e-02    1.000e+00  1.000e+00      5092  2.304e+03  2.299e+03                   
#> Path [50] :Best Iter: [76] ELBO (2303.997225) evaluations: (5092) 
#> Finished in  11.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
