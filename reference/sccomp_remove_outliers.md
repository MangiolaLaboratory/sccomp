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
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481353.719698 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.566e-02   2.321e-01    7.502e-01  7.502e-01      4120 -3.685e+03 -3.696e+03                   
#> Path [1] :Best Iter: [58] ELBO (-3685.085371) evaluations: (4120) 
#> Path [2] :Initial log joint density = -481513.420270 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.849e-02   1.775e-01    1.000e+00  1.000e+00      4037 -3.688e+03 -3.692e+03                   
#> Path [2] :Best Iter: [62] ELBO (-3687.908869) evaluations: (4037) 
#> Path [3] :Initial log joint density = -481408.821776 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.076e-03   2.380e-01    1.000e+00  1.000e+00      4081 -3.684e+03 -3.695e+03                   
#> Path [3] :Best Iter: [63] ELBO (-3684.433329) evaluations: (4081) 
#> Path [4] :Initial log joint density = -481318.404052 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      3.314e-03   1.879e-01    8.388e-01  8.388e-01      3539 -3.686e+03 -3.699e+03                   
#> Path [4] :Best Iter: [59] ELBO (-3686.256939) evaluations: (3539) 
#> Path [5] :Initial log joint density = -481706.124609 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.594e-03   2.378e-01    5.562e-01  5.562e-01      4322 -3.683e+03 -3.698e+03                   
#> Path [5] :Best Iter: [58] ELBO (-3683.067748) evaluations: (4322) 
#> Path [6] :Initial log joint density = -482185.180770 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.137e-02   2.088e-01    1.000e+00  1.000e+00      5268 -3.685e+03 -3.687e+03                   
#> Path [6] :Best Iter: [75] ELBO (-3685.028490) evaluations: (5268) 
#> Path [7] :Initial log joint density = -481522.595269 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      5.956e-03   2.357e-01    4.852e-01  1.000e+00      3200 -3.696e+03 -3.699e+03                   
#> Path [7] :Best Iter: [53] ELBO (-3696.295806) evaluations: (3200) 
#> Path [8] :Initial log joint density = -481707.823652 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      3.493e-03   3.798e-01    5.427e-01  5.427e-01      3080 -3.694e+03 -3.706e+03                   
#> Path [8] :Best Iter: [53] ELBO (-3693.664774) evaluations: (3080) 
#> Path [9] :Initial log joint density = -482785.458933 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.563e-02   3.553e-01    1.000e+00  1.000e+00      5149 -3.682e+03 -3.691e+03                   
#> Path [9] :Best Iter: [64] ELBO (-3681.851577) evaluations: (5149) 
#> Path [10] :Initial log joint density = -481712.816292 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      6.257e-03   2.152e-01    1.000e+00  1.000e+00      2887 -3.697e+03 -3.700e+03                   
#> Path [10] :Best Iter: [50] ELBO (-3696.662658) evaluations: (2887) 
#> Path [11] :Initial log joint density = -483570.389284 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      7.844e-03   1.298e-01    8.991e-01  8.991e-01      4438 -3.682e+03 -3.694e+03                   
#> Path [11] :Best Iter: [64] ELBO (-3681.974270) evaluations: (4438) 
#> Path [12] :Initial log joint density = -484585.371215 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.083e-02   1.589e-01    1.000e+00  1.000e+00      4463 -3.683e+03 -3.684e+03                   
#> Path [12] :Best Iter: [64] ELBO (-3682.923868) evaluations: (4463) 
#> Path [13] :Initial log joint density = -481732.145239 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.080e-02   2.167e-01    1.000e+00  1.000e+00      3125 -3.691e+03 -3.686e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3685.524595) evaluations: (3125) 
#> Path [14] :Initial log joint density = -481365.525034 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.978e-02   2.361e-01    1.000e+00  1.000e+00      4201 -3.685e+03 -3.691e+03                   
#> Path [14] :Best Iter: [62] ELBO (-3685.396723) evaluations: (4201) 
#> Path [15] :Initial log joint density = -481635.174219 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.747e-02   1.341e-01    1.000e+00  1.000e+00      4312 -3.686e+03 -3.686e+03                   
#> Path [15] :Best Iter: [68] ELBO (-3685.622004) evaluations: (4312) 
#> Path [16] :Initial log joint density = -481759.906029 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      9.332e-03   1.952e-01    8.180e-01  8.180e-01      4419 -3.683e+03 -3.695e+03                   
#> Path [16] :Best Iter: [68] ELBO (-3683.407445) evaluations: (4419) 
#> Path [17] :Initial log joint density = -481541.292122 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      9.138e-03   2.339e-01    1.000e+00  1.000e+00      4538 -3.685e+03 -3.686e+03                   
#> Path [17] :Best Iter: [66] ELBO (-3684.535332) evaluations: (4538) 
#> Path [18] :Initial log joint density = -481591.527381 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      9.240e-03   2.555e-01    1.000e+00  1.000e+00      2944 -3.695e+03 -3.698e+03                   
#> Path [18] :Best Iter: [47] ELBO (-3694.614052) evaluations: (2944) 
#> Path [19] :Initial log joint density = -482011.355551 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.168e-02   3.377e-01    1.000e+00  1.000e+00      4823 -3.685e+03 -3.697e+03                   
#> Path [19] :Best Iter: [71] ELBO (-3685.249917) evaluations: (4823) 
#> Path [20] :Initial log joint density = -481743.721396 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      8.596e-03   2.519e-01    1.000e+00  1.000e+00      3513 -3.688e+03 -3.696e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3688.406414) evaluations: (3513) 
#> Path [21] :Initial log joint density = -481779.053237 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.196e-02   2.141e-01    1.000e+00  1.000e+00      4599 -3.685e+03 -3.686e+03                   
#> Path [21] :Best Iter: [71] ELBO (-3685.258243) evaluations: (4599) 
#> Path [22] :Initial log joint density = -481573.863095 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.936e-03   2.500e-01    3.086e-01  3.086e-01      4334 -3.687e+03 -3.690e+03                   
#> Path [22] :Best Iter: [67] ELBO (-3686.814225) evaluations: (4334) 
#> Path [23] :Initial log joint density = -481902.001993 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      2.272e-02   2.852e-01    1.000e+00  1.000e+00      4688 -3.683e+03 -3.687e+03                   
#> Path [23] :Best Iter: [72] ELBO (-3683.060827) evaluations: (4688) 
#> Path [24] :Initial log joint density = -481732.272823 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      6.043e-03   2.305e-01    1.000e+00  1.000e+00      3214 -3.693e+03 -3.690e+03                   
#> Path [24] :Best Iter: [56] ELBO (-3690.175239) evaluations: (3214) 
#> Path [25] :Initial log joint density = -481609.790473 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.869e-03   1.930e-01    1.000e+00  1.000e+00      3509 -3.685e+03 -3.694e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3684.625873) evaluations: (3509) 
#> Path [26] :Initial log joint density = -483088.960508 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.140e-02   4.295e-01    1.000e+00  1.000e+00      3456 -3.689e+03 -3.699e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3689.424736) evaluations: (3456) 
#> Path [27] :Initial log joint density = -482556.740832 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      5.742e-03   2.073e-01    1.000e+00  1.000e+00      3330 -3.689e+03 -3.693e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3689.395883) evaluations: (3330) 
#> Path [28] :Initial log joint density = -481452.236604 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.987e-03   2.045e-01    1.000e+00  1.000e+00      3109 -3.692e+03 -3.702e+03                   
#> Path [28] :Best Iter: [52] ELBO (-3691.592158) evaluations: (3109) 
#> Path [29] :Initial log joint density = -481565.379052 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.025e-02   2.214e-01    1.000e+00  1.000e+00      4332 -3.685e+03 -3.690e+03                   
#> Path [29] :Best Iter: [65] ELBO (-3685.072403) evaluations: (4332) 
#> Path [30] :Initial log joint density = -482111.151556 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      6.654e-03   2.577e-01    1.000e+00  1.000e+00      4179 -3.687e+03 -3.698e+03                   
#> Path [30] :Best Iter: [64] ELBO (-3687.005771) evaluations: (4179) 
#> Path [31] :Initial log joint density = -481439.505562 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.146e-02   1.994e-01    1.000e+00  1.000e+00      4316 -3.688e+03 -3.691e+03                   
#> Path [31] :Best Iter: [64] ELBO (-3688.121401) evaluations: (4316) 
#> Path [32] :Initial log joint density = -481695.751006 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      8.283e-03   2.529e-01    1.000e+00  1.000e+00      4195 -3.688e+03 -3.695e+03                   
#> Path [32] :Best Iter: [64] ELBO (-3688.394435) evaluations: (4195) 
#> Path [33] :Initial log joint density = -481534.139293 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      4.501e-03   2.319e-01    6.858e-01  6.858e-01      4197 -3.685e+03 -3.705e+03                   
#> Path [33] :Best Iter: [65] ELBO (-3684.707349) evaluations: (4197) 
#> Path [34] :Initial log joint density = -483691.724293 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      3.332e-02   2.701e-01    1.000e+00  1.000e+00      4473 -3.682e+03 -3.688e+03                   
#> Path [34] :Best Iter: [69] ELBO (-3681.962295) evaluations: (4473) 
#> Path [35] :Initial log joint density = -481517.038180 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      4.219e-02   2.068e-01    1.000e+00  1.000e+00      4079 -3.685e+03 -3.688e+03                   
#> Path [35] :Best Iter: [59] ELBO (-3684.913230) evaluations: (4079) 
#> Path [36] :Initial log joint density = -481473.829272 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      2.951e-02   2.419e-01    1.000e+00  1.000e+00      4313 -3.687e+03 -3.688e+03                   
#> Path [36] :Best Iter: [63] ELBO (-3686.652351) evaluations: (4313) 
#> Path [37] :Initial log joint density = -481850.507464 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      8.469e-03   1.673e-01    1.000e+00  1.000e+00      4814 -3.687e+03 -3.699e+03                   
#> Path [37] :Best Iter: [71] ELBO (-3687.207731) evaluations: (4814) 
#> Path [38] :Initial log joint density = -481843.920829 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.200e-02   2.535e-01    1.000e+00  1.000e+00      4602 -3.684e+03 -3.686e+03                   
#> Path [38] :Best Iter: [70] ELBO (-3684.146122) evaluations: (4602) 
#> Path [39] :Initial log joint density = -481520.593820 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.165e-02   2.663e-01    1.000e+00  1.000e+00      4200 -3.687e+03 -3.686e+03                   
#> Path [39] :Best Iter: [67] ELBO (-3686.081509) evaluations: (4200) 
#> Path [40] :Initial log joint density = -481919.236636 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.260e-02   3.422e-01    1.000e+00  1.000e+00      5278 -3.684e+03 -3.699e+03                   
#> Path [40] :Best Iter: [76] ELBO (-3683.614660) evaluations: (5278) 
#> Path [41] :Initial log joint density = -483674.364479 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.250e-02   4.543e-01    1.000e+00  1.000e+00      3363 -3.691e+03 -3.696e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3690.814427) evaluations: (3363) 
#> Path [42] :Initial log joint density = -481201.764970 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.292e-02   3.033e-01    1.000e+00  1.000e+00      3735 -3.685e+03 -3.692e+03                   
#> Path [42] :Best Iter: [61] ELBO (-3685.233207) evaluations: (3735) 
#> Path [43] :Initial log joint density = -485826.136157 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.970e-02   2.254e-01    1.000e+00  1.000e+00      5050 -3.680e+03 -3.684e+03                   
#> Path [43] :Best Iter: [75] ELBO (-3680.138378) evaluations: (5050) 
#> Path [44] :Initial log joint density = -481508.117302 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      3.975e-03   1.813e-01    7.658e-01  7.658e-01      4304 -3.686e+03 -3.692e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3685.587052) evaluations: (4304) 
#> Path [45] :Initial log joint density = -481391.337336 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      2.283e-02   2.862e-01    1.000e+00  1.000e+00      4435 -3.684e+03 -3.689e+03                   
#> Path [45] :Best Iter: [59] ELBO (-3683.847581) evaluations: (4435) 
#> Path [46] :Initial log joint density = -481583.206700 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      6.200e-03   1.795e-01    1.000e+00  1.000e+00      3328 -3.687e+03 -3.687e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3686.650531) evaluations: (3328) 
#> Path [47] :Initial log joint density = -481582.297862 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      8.840e-03   3.886e-01    4.265e-01  1.000e+00      3020 -3.693e+03 -3.698e+03                   
#> Path [47] :Best Iter: [53] ELBO (-3693.142842) evaluations: (3020) 
#> Path [48] :Initial log joint density = -481464.903770 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.227e-02   2.815e-01    8.857e-01  8.857e-01      3000 -3.693e+03 -3.703e+03                   
#> Path [48] :Best Iter: [51] ELBO (-3692.522942) evaluations: (3000) 
#> Path [49] :Initial log joint density = -481598.814745 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      7.512e-03   1.682e-01    9.082e-01  9.082e-01      3131 -3.695e+03 -3.698e+03                   
#> Path [49] :Best Iter: [53] ELBO (-3695.499233) evaluations: (3131) 
#> Path [50] :Initial log joint density = -481344.977960 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.441e-02   2.392e-01    1.000e+00  1.000e+00      2988 -3.694e+03 -3.694e+03                   
#> Path [50] :Best Iter: [54] ELBO (-3693.549752) evaluations: (2988) 
#> Finished in  15.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 4.182 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier identification - step 1/2
#> Loading model from cache...
#> Path [1] :Initial log joint density = -428930.251263 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      1.907e-02   1.662e-01    1.000e+00  1.000e+00      5369 -3.287e+03 -3.289e+03                   
#> Path [1] :Best Iter: [76] ELBO (-3287.497684) evaluations: (5369) 
#> Path [2] :Initial log joint density = -428914.583986 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      9.396e-03   2.630e-01    1.000e+00  1.000e+00      4183 -3.293e+03 -3.305e+03                   
#> Path [2] :Best Iter: [64] ELBO (-3293.402220) evaluations: (4183) 
#> Path [3] :Initial log joint density = -430354.859811 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.262e+05      2.013e-02   2.640e-01    4.676e-01  1.000e+00      5599 -3.291e+03 -3.300e+03                   
#> Path [3] :Best Iter: [78] ELBO (-3290.607673) evaluations: (5599) 
#> Path [4] :Initial log joint density = -429355.554519 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      3.895e-03   2.079e-01    6.796e-01  6.796e-01      6146 -3.290e+03 -3.309e+03                   
#> Path [4] :Best Iter: [82] ELBO (-3289.921421) evaluations: (6146) 
#> Path [5] :Initial log joint density = -431078.660627 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      8.088e-03   3.151e-01    4.431e-01  1.000e+00      4340 -3.292e+03 -3.299e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3292.144341) evaluations: (4340) 
#> Path [6] :Initial log joint density = -428752.521719 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      8.297e-03   2.492e-01    1.000e+00  1.000e+00      4342 -3.294e+03 -3.302e+03                   
#> Path [6] :Best Iter: [66] ELBO (-3293.805300) evaluations: (4342) 
#> Path [7] :Initial log joint density = -428973.299849 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.262e+05      2.210e-02   2.648e-01    1.000e+00  1.000e+00      5220 -3.287e+03 -3.297e+03                   
#> Path [7] :Best Iter: [71] ELBO (-3286.759207) evaluations: (5220) 
#> Path [8] :Initial log joint density = -428843.125882 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      1.274e-02   3.439e-01    1.000e+00  1.000e+00      4227 -3.296e+03 -3.299e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3295.951597) evaluations: (4227) 
#> Path [9] :Initial log joint density = -429110.490191 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.262e+05      9.184e-03   2.163e-01    1.000e+00  1.000e+00      5755 -3.293e+03 -3.305e+03                   
#> Path [9] :Best Iter: [77] ELBO (-3292.974876) evaluations: (5755) 
#> Path [10] :Initial log joint density = -429582.055966 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      3.138e-03   2.127e-01    7.241e-01  7.241e-01      6145 -3.290e+03 -3.301e+03                   
#> Path [10] :Best Iter: [83] ELBO (-3290.477389) evaluations: (6145) 
#> Path [11] :Initial log joint density = -429473.757820 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.262e+05      7.246e-03   2.444e-01    1.000e+00  1.000e+00      4972 -3.293e+03 -3.298e+03                   
#> Path [11] :Best Iter: [72] ELBO (-3293.246618) evaluations: (4972) 
#> Path [12] :Initial log joint density = -429235.427887 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.262e+05      7.164e-03   3.159e-01    4.324e-01  1.000e+00      5646 -3.292e+03 -3.304e+03                   
#> Path [12] :Best Iter: [80] ELBO (-3291.946072) evaluations: (5646) 
#> Path [13] :Initial log joint density = -428960.019234 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      7.287e-03   2.447e-01    1.000e+00  1.000e+00      4183 -3.295e+03 -3.304e+03                   
#> Path [13] :Best Iter: [65] ELBO (-3294.747470) evaluations: (4183) 
#> Path [14] :Initial log joint density = -429065.666965 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.262e+05      2.408e-03   1.606e-01    7.308e-01  7.308e-01      5122 -3.292e+03 -3.309e+03                   
#> Path [14] :Best Iter: [75] ELBO (-3291.591358) evaluations: (5122) 
#> Path [15] :Initial log joint density = -428765.253947 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      8.342e-03   2.209e-01    1.000e+00  1.000e+00      5242 -3.292e+03 -3.296e+03                   
#> Path [15] :Best Iter: [75] ELBO (-3291.560001) evaluations: (5242) 
#> Path [16] :Initial log joint density = -429686.434714 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.262e+05      6.633e-03   2.000e-01    1.000e+00  1.000e+00      6408 -3.291e+03 -3.291e+03                   
#> Path [16] :Best Iter: [84] ELBO (-3290.905251) evaluations: (6408) 
#> Path [17] :Initial log joint density = -428740.439809 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      6.625e-03   2.051e-01    1.000e+00  1.000e+00      4204 -3.296e+03 -3.301e+03                   
#> Path [17] :Best Iter: [63] ELBO (-3296.142987) evaluations: (4204) 
#> Path [18] :Initial log joint density = -428670.254897 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.262e+05      1.056e-02   1.362e-01    1.000e+00  1.000e+00      5696 -3.288e+03 -3.294e+03                   
#> Path [18] :Best Iter: [77] ELBO (-3288.113995) evaluations: (5696) 
#> Path [19] :Initial log joint density = -430205.214053 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      1.028e-02   2.951e-01    1.000e+00  1.000e+00      4346 -3.292e+03 -3.295e+03                   
#> Path [19] :Best Iter: [68] ELBO (-3292.331434) evaluations: (4346) 
#> Path [20] :Initial log joint density = -428579.292940 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.262e+05      4.089e-03   2.700e-01    6.048e-01  6.048e-01      4907 -3.292e+03 -3.302e+03                   
#> Path [20] :Best Iter: [67] ELBO (-3292.210719) evaluations: (4907) 
#> Path [21] :Initial log joint density = -429015.705603 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.262e+05      7.948e-03   1.676e-01    1.000e+00  1.000e+00      4077 -3.295e+03 -3.299e+03                   
#> Path [21] :Best Iter: [61] ELBO (-3295.346505) evaluations: (4077) 
#> Path [22] :Initial log joint density = -429548.061087 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      3.364e-03   1.584e-01    6.835e-01  6.835e-01      7091 -3.290e+03 -3.304e+03                   
#> Path [22] :Best Iter: [90] ELBO (-3290.060592) evaluations: (7091) 
#> Path [23] :Initial log joint density = -429506.987691 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      1.646e-02   3.729e-01    1.000e+00  1.000e+00      7028 -3.288e+03 -3.293e+03                   
#> Path [23] :Best Iter: [92] ELBO (-3288.236588) evaluations: (7028) 
#> Path [24] :Initial log joint density = -428935.909162 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      1.361e-02   2.861e-01    1.000e+00  1.000e+00      4347 -3.295e+03 -3.299e+03                   
#> Path [24] :Best Iter: [66] ELBO (-3295.322630) evaluations: (4347) 
#> Path [25] :Initial log joint density = -430754.500999 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.262e+05      6.656e-03   2.605e-01    7.030e-01  7.030e-01      4540 -3.293e+03 -3.309e+03                   
#> Path [25] :Best Iter: [69] ELBO (-3293.273167) evaluations: (4540) 
#> Path [26] :Initial log joint density = -429597.751329 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      6.932e-03   2.320e-01    1.000e+00  1.000e+00      6066 -3.289e+03 -3.296e+03                   
#> Path [26] :Best Iter: [81] ELBO (-3289.282312) evaluations: (6066) 
#> Path [27] :Initial log joint density = -429236.745400 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.262e+05      1.107e-02   3.626e-01    1.000e+00  1.000e+00      6438 -3.293e+03 -3.300e+03                   
#> Path [27] :Best Iter: [83] ELBO (-3293.009607) evaluations: (6438) 
#> Path [28] :Initial log joint density = -429287.119198 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.262e+05      6.132e-03   2.310e-01    1.000e+00  1.000e+00      5739 -3.292e+03 -3.306e+03                   
#> Path [28] :Best Iter: [74] ELBO (-3292.463058) evaluations: (5739) 
#> Path [29] :Initial log joint density = -429067.071925 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.262e+05      3.363e-03   2.821e-01    7.480e-01  7.480e-01      5875 -3.291e+03 -3.306e+03                   
#> Path [29] :Best Iter: [81] ELBO (-3291.485481) evaluations: (5875) 
#> Path [30] :Initial log joint density = -428891.401145 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      1.406e-02   2.327e-01    1.000e+00  1.000e+00      4436 -3.292e+03 -3.296e+03                   
#> Path [30] :Best Iter: [66] ELBO (-3291.601205) evaluations: (4436) 
#> Path [31] :Initial log joint density = -429208.178933 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.262e+05      1.245e-02   2.920e-01    1.000e+00  1.000e+00      5719 -3.292e+03 -3.302e+03                   
#> Path [31] :Best Iter: [81] ELBO (-3292.211556) evaluations: (5719) 
#> Path [32] :Initial log joint density = -428901.903684 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      4.200e-03   2.958e-01    6.693e-01  6.693e-01      4238 -3.296e+03 -3.308e+03                   
#> Path [32] :Best Iter: [66] ELBO (-3296.109845) evaluations: (4238) 
#> Path [33] :Initial log joint density = -428942.826799 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      1.301e-02   4.417e-01    1.000e+00  1.000e+00      4261 -3.291e+03 -3.301e+03                   
#> Path [33] :Best Iter: [67] ELBO (-3290.510257) evaluations: (4261) 
#> Path [34] :Initial log joint density = -433111.476340 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      8.873e-03   2.739e-01    4.989e-01  1.000e+00      4505 -3.292e+03 -3.305e+03                   
#> Path [34] :Best Iter: [67] ELBO (-3292.159768) evaluations: (4505) 
#> Path [35] :Initial log joint density = -430649.980518 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      6.447e-03   2.416e-01    6.663e-01  6.663e-01      5382 -3.286e+03 -3.305e+03                   
#> Path [35] :Best Iter: [76] ELBO (-3286.496477) evaluations: (5382) 
#> Path [36] :Initial log joint density = -428861.654734 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.262e+05      1.057e-02   3.172e-01    1.000e+00  1.000e+00      6143 -3.291e+03 -3.294e+03                   
#> Path [36] :Best Iter: [83] ELBO (-3291.429717) evaluations: (6143) 
#> Path [37] :Initial log joint density = -428870.847631 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      5.562e-03   2.064e-01    1.000e+00  1.000e+00      4146 -3.296e+03 -3.300e+03                   
#> Path [37] :Best Iter: [64] ELBO (-3296.393854) evaluations: (4146) 
#> Path [38] :Initial log joint density = -430297.582363 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.262e+05      1.069e-02   2.947e-01    1.000e+00  1.000e+00      4035 -3.298e+03 -3.305e+03                   
#> Path [38] :Best Iter: [62] ELBO (-3297.539701) evaluations: (4035) 
#> Path [39] :Initial log joint density = -430071.919788 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      6.740e-03   1.999e-01    1.000e+00  1.000e+00      4150 -3.295e+03 -3.301e+03                   
#> Path [39] :Best Iter: [63] ELBO (-3295.439090) evaluations: (4150) 
#> Path [40] :Initial log joint density = -429536.784005 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      3.211e-03   1.678e-01    7.145e-01  7.145e-01      7189 -3.290e+03 -3.303e+03                   
#> Path [40] :Best Iter: [91] ELBO (-3289.578334) evaluations: (7189) 
#> Path [41] :Initial log joint density = -429150.903194 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.262e+05      8.236e-03   2.094e-01    1.000e+00  1.000e+00      4941 -3.293e+03 -3.303e+03                   
#> Path [41] :Best Iter: [72] ELBO (-3292.642946) evaluations: (4941) 
#> Path [42] :Initial log joint density = -429254.287623 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.262e+05      1.411e-02   3.420e-01    1.000e+00  1.000e+00      6455 -3.289e+03 -3.296e+03                   
#> Path [42] :Best Iter: [86] ELBO (-3289.449731) evaluations: (6455) 
#> Path [43] :Initial log joint density = -428858.867993 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.262e+05      1.204e-02   4.500e-01    1.000e+00  1.000e+00      3838 -3.296e+03 -3.308e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3295.972923) evaluations: (3838) 
#> Path [44] :Initial log joint density = -429019.829078 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      7.013e-03   1.776e-01    1.000e+00  1.000e+00      4244 -3.292e+03 -3.301e+03                   
#> Path [44] :Best Iter: [65] ELBO (-3291.998088) evaluations: (4244) 
#> Path [45] :Initial log joint density = -430332.557992 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.262e+05      2.964e-02   2.447e-01    1.000e+00  1.000e+00      4971 -3.292e+03 -3.293e+03                   
#> Path [45] :Best Iter: [71] ELBO (-3291.507757) evaluations: (4971) 
#> Path [46] :Initial log joint density = -428809.592495 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      5.287e-03   2.202e-01    1.000e+00  1.000e+00      4449 -3.296e+03 -3.308e+03                   
#> Path [46] :Best Iter: [67] ELBO (-3296.347109) evaluations: (4449) 
#> Path [47] :Initial log joint density = -429190.080380 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              88      -4.262e+05      7.348e-03   2.180e-01    1.000e+00  1.000e+00      6446 -3.291e+03 -3.296e+03                   
#> Path [47] :Best Iter: [85] ELBO (-3290.752147) evaluations: (6446) 
#> Path [48] :Initial log joint density = -429049.541562 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.262e+05      7.733e-03   2.219e-01    1.000e+00  1.000e+00      4858 -3.294e+03 -3.296e+03                   
#> Path [48] :Best Iter: [69] ELBO (-3294.129989) evaluations: (4858) 
#> Path [49] :Initial log joint density = -428899.064737 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.262e+05      1.300e-02   1.834e-01    1.000e+00  1.000e+00      5122 -3.287e+03 -3.291e+03                   
#> Path [49] :Best Iter: [74] ELBO (-3287.482130) evaluations: (5122) 
#> Path [50] :Initial log joint density = -429220.221086 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      1.385e-02   2.135e-01    1.000e+00  1.000e+00      4493 -3.292e+03 -3.296e+03                   
#> Path [50] :Best Iter: [68] ELBO (-3291.640161) evaluations: (4493) 
#> Finished in  24.1 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 20.819 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -428930.251263 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      1.907e-02   1.662e-01    1.000e+00  1.000e+00      5369 -3.287e+03 -3.289e+03                   
#> Path [1] :Best Iter: [76] ELBO (-3287.497684) evaluations: (5369) 
#> Path [2] :Initial log joint density = -428914.583986 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      9.396e-03   2.630e-01    1.000e+00  1.000e+00      4183 -3.293e+03 -3.305e+03                   
#> Path [2] :Best Iter: [64] ELBO (-3293.402220) evaluations: (4183) 
#> Path [3] :Initial log joint density = -430354.859811 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.262e+05      2.013e-02   2.640e-01    4.676e-01  1.000e+00      5599 -3.291e+03 -3.300e+03                   
#> Path [3] :Best Iter: [78] ELBO (-3290.607673) evaluations: (5599) 
#> Path [4] :Initial log joint density = -429355.554519 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      3.895e-03   2.079e-01    6.796e-01  6.796e-01      6146 -3.290e+03 -3.309e+03                   
#> Path [4] :Best Iter: [82] ELBO (-3289.921421) evaluations: (6146) 
#> Path [5] :Initial log joint density = -431078.660627 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      8.088e-03   3.151e-01    4.431e-01  1.000e+00      4340 -3.292e+03 -3.299e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3292.144341) evaluations: (4340) 
#> Path [6] :Initial log joint density = -428752.521719 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      8.297e-03   2.492e-01    1.000e+00  1.000e+00      4342 -3.294e+03 -3.302e+03                   
#> Path [6] :Best Iter: [66] ELBO (-3293.805300) evaluations: (4342) 
#> Path [7] :Initial log joint density = -428973.299849 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.262e+05      2.210e-02   2.648e-01    1.000e+00  1.000e+00      5220 -3.287e+03 -3.297e+03                   
#> Path [7] :Best Iter: [71] ELBO (-3286.759207) evaluations: (5220) 
#> Path [8] :Initial log joint density = -428843.125882 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      1.274e-02   3.439e-01    1.000e+00  1.000e+00      4227 -3.296e+03 -3.299e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3295.951597) evaluations: (4227) 
#> Path [9] :Initial log joint density = -429110.490191 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.262e+05      9.184e-03   2.163e-01    1.000e+00  1.000e+00      5755 -3.293e+03 -3.305e+03                   
#> Path [9] :Best Iter: [77] ELBO (-3292.974876) evaluations: (5755) 
#> Path [10] :Initial log joint density = -429582.055966 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      3.138e-03   2.127e-01    7.241e-01  7.241e-01      6145 -3.290e+03 -3.301e+03                   
#> Path [10] :Best Iter: [83] ELBO (-3290.477389) evaluations: (6145) 
#> Path [11] :Initial log joint density = -429473.757820 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.262e+05      7.246e-03   2.444e-01    1.000e+00  1.000e+00      4972 -3.293e+03 -3.298e+03                   
#> Path [11] :Best Iter: [72] ELBO (-3293.246618) evaluations: (4972) 
#> Path [12] :Initial log joint density = -429235.427887 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.262e+05      7.164e-03   3.159e-01    4.324e-01  1.000e+00      5646 -3.292e+03 -3.304e+03                   
#> Path [12] :Best Iter: [80] ELBO (-3291.946072) evaluations: (5646) 
#> Path [13] :Initial log joint density = -428960.019234 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      7.287e-03   2.447e-01    1.000e+00  1.000e+00      4183 -3.295e+03 -3.304e+03                   
#> Path [13] :Best Iter: [65] ELBO (-3294.747470) evaluations: (4183) 
#> Path [14] :Initial log joint density = -429065.666965 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.262e+05      2.408e-03   1.606e-01    7.308e-01  7.308e-01      5122 -3.292e+03 -3.309e+03                   
#> Path [14] :Best Iter: [75] ELBO (-3291.591358) evaluations: (5122) 
#> Path [15] :Initial log joint density = -428765.253947 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      8.342e-03   2.209e-01    1.000e+00  1.000e+00      5242 -3.292e+03 -3.296e+03                   
#> Path [15] :Best Iter: [75] ELBO (-3291.560001) evaluations: (5242) 
#> Path [16] :Initial log joint density = -429686.434714 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.262e+05      6.633e-03   2.000e-01    1.000e+00  1.000e+00      6408 -3.291e+03 -3.291e+03                   
#> Path [16] :Best Iter: [84] ELBO (-3290.905251) evaluations: (6408) 
#> Path [17] :Initial log joint density = -428740.439809 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      6.625e-03   2.051e-01    1.000e+00  1.000e+00      4204 -3.296e+03 -3.301e+03                   
#> Path [17] :Best Iter: [63] ELBO (-3296.142987) evaluations: (4204) 
#> Path [18] :Initial log joint density = -428670.254897 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.262e+05      1.056e-02   1.362e-01    1.000e+00  1.000e+00      5696 -3.288e+03 -3.294e+03                   
#> Path [18] :Best Iter: [77] ELBO (-3288.113995) evaluations: (5696) 
#> Path [19] :Initial log joint density = -430205.214053 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      1.028e-02   2.951e-01    1.000e+00  1.000e+00      4346 -3.292e+03 -3.295e+03                   
#> Path [19] :Best Iter: [68] ELBO (-3292.331434) evaluations: (4346) 
#> Path [20] :Initial log joint density = -428579.292940 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.262e+05      4.089e-03   2.700e-01    6.048e-01  6.048e-01      4907 -3.292e+03 -3.302e+03                   
#> Path [20] :Best Iter: [67] ELBO (-3292.210719) evaluations: (4907) 
#> Path [21] :Initial log joint density = -429015.705603 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.262e+05      7.948e-03   1.676e-01    1.000e+00  1.000e+00      4077 -3.295e+03 -3.299e+03                   
#> Path [21] :Best Iter: [61] ELBO (-3295.346505) evaluations: (4077) 
#> Path [22] :Initial log joint density = -429548.061087 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      3.364e-03   1.584e-01    6.835e-01  6.835e-01      7091 -3.290e+03 -3.304e+03                   
#> Path [22] :Best Iter: [90] ELBO (-3290.060592) evaluations: (7091) 
#> Path [23] :Initial log joint density = -429506.987691 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      1.646e-02   3.729e-01    1.000e+00  1.000e+00      7028 -3.288e+03 -3.293e+03                   
#> Path [23] :Best Iter: [92] ELBO (-3288.236588) evaluations: (7028) 
#> Path [24] :Initial log joint density = -428935.909162 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.262e+05      1.361e-02   2.861e-01    1.000e+00  1.000e+00      4347 -3.295e+03 -3.299e+03                   
#> Path [24] :Best Iter: [66] ELBO (-3295.322630) evaluations: (4347) 
#> Path [25] :Initial log joint density = -430754.500999 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.262e+05      6.656e-03   2.605e-01    7.030e-01  7.030e-01      4540 -3.293e+03 -3.309e+03                   
#> Path [25] :Best Iter: [69] ELBO (-3293.273167) evaluations: (4540) 
#> Path [26] :Initial log joint density = -429597.751329 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.262e+05      6.932e-03   2.320e-01    1.000e+00  1.000e+00      6066 -3.289e+03 -3.296e+03                   
#> Path [26] :Best Iter: [81] ELBO (-3289.282312) evaluations: (6066) 
#> Path [27] :Initial log joint density = -429236.745400 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.262e+05      1.107e-02   3.626e-01    1.000e+00  1.000e+00      6438 -3.293e+03 -3.300e+03                   
#> Path [27] :Best Iter: [83] ELBO (-3293.009607) evaluations: (6438) 
#> Path [28] :Initial log joint density = -429287.119198 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.262e+05      6.132e-03   2.310e-01    1.000e+00  1.000e+00      5739 -3.292e+03 -3.306e+03                   
#> Path [28] :Best Iter: [74] ELBO (-3292.463058) evaluations: (5739) 
#> Path [29] :Initial log joint density = -429067.071925 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.262e+05      3.363e-03   2.821e-01    7.480e-01  7.480e-01      5875 -3.291e+03 -3.306e+03                   
#> Path [29] :Best Iter: [81] ELBO (-3291.485481) evaluations: (5875) 
#> Path [30] :Initial log joint density = -428891.401145 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      1.406e-02   2.327e-01    1.000e+00  1.000e+00      4436 -3.292e+03 -3.296e+03                   
#> Path [30] :Best Iter: [66] ELBO (-3291.601205) evaluations: (4436) 
#> Path [31] :Initial log joint density = -429208.178933 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.262e+05      1.245e-02   2.920e-01    1.000e+00  1.000e+00      5719 -3.292e+03 -3.302e+03                   
#> Path [31] :Best Iter: [81] ELBO (-3292.211556) evaluations: (5719) 
#> Path [32] :Initial log joint density = -428901.903684 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      4.200e-03   2.958e-01    6.693e-01  6.693e-01      4238 -3.296e+03 -3.308e+03                   
#> Path [32] :Best Iter: [66] ELBO (-3296.109845) evaluations: (4238) 
#> Path [33] :Initial log joint density = -428942.826799 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      1.301e-02   4.417e-01    1.000e+00  1.000e+00      4261 -3.291e+03 -3.301e+03                   
#> Path [33] :Best Iter: [67] ELBO (-3290.510257) evaluations: (4261) 
#> Path [34] :Initial log joint density = -433111.476340 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      8.873e-03   2.739e-01    4.989e-01  1.000e+00      4505 -3.292e+03 -3.305e+03                   
#> Path [34] :Best Iter: [67] ELBO (-3292.159768) evaluations: (4505) 
#> Path [35] :Initial log joint density = -430649.980518 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.262e+05      6.447e-03   2.416e-01    6.663e-01  6.663e-01      5382 -3.286e+03 -3.305e+03                   
#> Path [35] :Best Iter: [76] ELBO (-3286.496477) evaluations: (5382) 
#> Path [36] :Initial log joint density = -428861.654734 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.262e+05      1.057e-02   3.172e-01    1.000e+00  1.000e+00      6143 -3.291e+03 -3.294e+03                   
#> Path [36] :Best Iter: [83] ELBO (-3291.429717) evaluations: (6143) 
#> Path [37] :Initial log joint density = -428870.847631 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      5.562e-03   2.064e-01    1.000e+00  1.000e+00      4146 -3.296e+03 -3.300e+03                   
#> Path [37] :Best Iter: [64] ELBO (-3296.393854) evaluations: (4146) 
#> Path [38] :Initial log joint density = -430297.582363 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.262e+05      1.069e-02   2.947e-01    1.000e+00  1.000e+00      4035 -3.298e+03 -3.305e+03                   
#> Path [38] :Best Iter: [62] ELBO (-3297.539701) evaluations: (4035) 
#> Path [39] :Initial log joint density = -430071.919788 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.262e+05      6.740e-03   1.999e-01    1.000e+00  1.000e+00      4150 -3.295e+03 -3.301e+03                   
#> Path [39] :Best Iter: [63] ELBO (-3295.439090) evaluations: (4150) 
#> Path [40] :Initial log joint density = -429536.784005 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.262e+05      3.211e-03   1.678e-01    7.145e-01  7.145e-01      7189 -3.290e+03 -3.303e+03                   
#> Path [40] :Best Iter: [91] ELBO (-3289.578334) evaluations: (7189) 
#> Path [41] :Initial log joint density = -429150.903194 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.262e+05      8.236e-03   2.094e-01    1.000e+00  1.000e+00      4941 -3.293e+03 -3.303e+03                   
#> Path [41] :Best Iter: [72] ELBO (-3292.642946) evaluations: (4941) 
#> Path [42] :Initial log joint density = -429254.287623 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.262e+05      1.411e-02   3.420e-01    1.000e+00  1.000e+00      6455 -3.289e+03 -3.296e+03                   
#> Path [42] :Best Iter: [86] ELBO (-3289.449731) evaluations: (6455) 
#> Path [43] :Initial log joint density = -428858.867993 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.262e+05      1.204e-02   4.500e-01    1.000e+00  1.000e+00      3838 -3.296e+03 -3.308e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3295.972923) evaluations: (3838) 
#> Path [44] :Initial log joint density = -429019.829078 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.262e+05      7.013e-03   1.776e-01    1.000e+00  1.000e+00      4244 -3.292e+03 -3.301e+03                   
#> Path [44] :Best Iter: [65] ELBO (-3291.998088) evaluations: (4244) 
#> Path [45] :Initial log joint density = -430332.557992 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.262e+05      2.964e-02   2.447e-01    1.000e+00  1.000e+00      4971 -3.292e+03 -3.293e+03                   
#> Path [45] :Best Iter: [71] ELBO (-3291.507757) evaluations: (4971) 
#> Path [46] :Initial log joint density = -428809.592495 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      5.287e-03   2.202e-01    1.000e+00  1.000e+00      4449 -3.296e+03 -3.308e+03                   
#> Path [46] :Best Iter: [67] ELBO (-3296.347109) evaluations: (4449) 
#> Path [47] :Initial log joint density = -429190.080380 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              88      -4.262e+05      7.348e-03   2.180e-01    1.000e+00  1.000e+00      6446 -3.291e+03 -3.296e+03                   
#> Path [47] :Best Iter: [85] ELBO (-3290.752147) evaluations: (6446) 
#> Path [48] :Initial log joint density = -429049.541562 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.262e+05      7.733e-03   2.219e-01    1.000e+00  1.000e+00      4858 -3.294e+03 -3.296e+03                   
#> Path [48] :Best Iter: [69] ELBO (-3294.129989) evaluations: (4858) 
#> Path [49] :Initial log joint density = -428899.064737 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.262e+05      1.300e-02   1.834e-01    1.000e+00  1.000e+00      5122 -3.287e+03 -3.291e+03                   
#> Path [49] :Best Iter: [74] ELBO (-3287.482130) evaluations: (5122) 
#> Path [50] :Initial log joint density = -429220.221086 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.262e+05      1.385e-02   2.135e-01    1.000e+00  1.000e+00      4493 -3.292e+03 -3.296e+03                   
#> Path [50] :Best Iter: [68] ELBO (-3291.640161) evaluations: (4493) 
#> Finished in  21.1 seconds.
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
