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
#> Path [1] :Initial log joint density = -481678.231158 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.901e-03   2.167e-01    1.000e+00  1.000e+00      2988 -3.689e+03 -3.701e+03                   
#> Path [1] :Best Iter: [47] ELBO (-3689.370971) evaluations: (2988) 
#> Path [2] :Initial log joint density = -481542.163837 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.331e-02   2.516e-01    1.000e+00  1.000e+00      3701 -3.679e+03 -3.685e+03                   
#> Path [2] :Best Iter: [60] ELBO (-3679.400461) evaluations: (3701) 
#> Path [3] :Initial log joint density = -481563.290368 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.323e-03   2.537e-01    1.000e+00  1.000e+00      3201 -3.684e+03 -3.685e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3684.010662) evaluations: (3201) 
#> Path [4] :Initial log joint density = -481102.925071 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      2.905e-03   3.520e-01    5.625e-01  5.625e-01      2910 -3.691e+03 -3.704e+03                   
#> Path [4] :Best Iter: [41] ELBO (-3690.706549) evaluations: (2910) 
#> Path [5] :Initial log joint density = -481772.334889 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.164e-03   2.265e-01    1.000e+00  1.000e+00      3258 -3.684e+03 -3.685e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3684.304842) evaluations: (3258) 
#> Path [6] :Initial log joint density = -481923.142402 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      7.937e-03   2.736e-01    1.000e+00  1.000e+00      3888 -3.681e+03 -3.697e+03                   
#> Path [6] :Best Iter: [61] ELBO (-3681.390202) evaluations: (3888) 
#> Path [7] :Initial log joint density = -481876.508555 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.028e-02   3.026e-01    9.983e-01  9.983e-01      3649 -3.682e+03 -3.686e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3682.126325) evaluations: (3649) 
#> Path [8] :Initial log joint density = -481689.022118 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.964e-03   2.672e-01    1.000e+00  1.000e+00      3193 -3.683e+03 -3.693e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3682.864517) evaluations: (3193) 
#> Path [9] :Initial log joint density = -481844.554164 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.521e-03   2.608e-01    8.338e-01  8.338e-01      2825 -3.691e+03 -3.704e+03                   
#> Path [9] :Best Iter: [51] ELBO (-3691.153013) evaluations: (2825) 
#> Path [10] :Initial log joint density = -482217.438572 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.616e-02   3.331e-01    1.000e+00  1.000e+00      3241 -3.684e+03 -3.684e+03                   
#> Path [10] :Best Iter: [56] ELBO (-3683.740005) evaluations: (3241) 
#> Path [11] :Initial log joint density = -481494.046713 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.288e-02   1.883e-01    1.000e+00  1.000e+00      3339 -3.690e+03 -3.688e+03                   
#> Path [11] :Best Iter: [56] ELBO (-3687.988574) evaluations: (3339) 
#> Path [12] :Initial log joint density = -481678.904165 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.833e-03   2.055e-01    9.460e-01  9.460e-01      3356 -3.683e+03 -3.693e+03                   
#> Path [12] :Best Iter: [56] ELBO (-3683.408574) evaluations: (3356) 
#> Path [13] :Initial log joint density = -481540.093775 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.042e-02   2.524e-01    1.000e+00  1.000e+00      2917 -3.688e+03 -3.697e+03                   
#> Path [13] :Best Iter: [41] ELBO (-3688.442889) evaluations: (2917) 
#> Path [14] :Initial log joint density = -481631.800459 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.205e-02   2.117e-01    1.000e+00  1.000e+00      3589 -3.682e+03 -3.683e+03                   
#> Path [14] :Best Iter: [60] ELBO (-3681.731222) evaluations: (3589) 
#> Path [15] :Initial log joint density = -482486.174106 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      9.549e-03   2.565e-01    1.000e+00  1.000e+00      3296 -3.688e+03 -3.686e+03                   
#> Path [15] :Best Iter: [57] ELBO (-3685.660094) evaluations: (3296) 
#> Path [16] :Initial log joint density = -481498.297375 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.172e-02   3.263e-01    1.000e+00  1.000e+00      2829 -3.693e+03 -3.705e+03                   
#> Path [16] :Best Iter: [45] ELBO (-3693.103115) evaluations: (2829) 
#> Path [17] :Initial log joint density = -484936.480722 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.201e-03   1.857e-01    1.000e+00  1.000e+00      3406 -3.687e+03 -3.684e+03                   
#> Path [17] :Best Iter: [59] ELBO (-3683.664174) evaluations: (3406) 
#> Path [18] :Initial log joint density = -481537.473035 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.096e-03   2.773e-01    5.992e-01  5.992e-01      3399 -3.684e+03 -3.691e+03                   
#> Path [18] :Best Iter: [58] ELBO (-3684.446198) evaluations: (3399) 
#> Path [19] :Initial log joint density = -481827.093152 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.551e-03   2.029e-01    8.636e-01  8.636e-01      3320 -3.683e+03 -3.694e+03                   
#> Path [19] :Best Iter: [55] ELBO (-3683.074552) evaluations: (3320) 
#> Path [20] :Initial log joint density = -481960.805451 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.040e-02   1.968e-01    1.000e+00  1.000e+00      2953 -3.687e+03 -3.700e+03                   
#> Path [20] :Best Iter: [51] ELBO (-3686.760042) evaluations: (2953) 
#> Path [21] :Initial log joint density = -481796.400852 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.662e-03   2.048e-01    8.644e-01  8.644e-01      3263 -3.686e+03 -3.691e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3685.998576) evaluations: (3263) 
#> Path [22] :Initial log joint density = -482421.750492 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.393e-02   3.407e-01    1.000e+00  1.000e+00      3241 -3.684e+03 -3.694e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3683.674269) evaluations: (3241) 
#> Path [23] :Initial log joint density = -483300.412757 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.988e-03   2.461e-01    1.000e+00  1.000e+00      3077 -3.691e+03 -3.683e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3683.276285) evaluations: (3077) 
#> Path [24] :Initial log joint density = -481476.904789 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.563e-03   2.112e-01    1.000e+00  1.000e+00      2879 -3.690e+03 -3.695e+03                   
#> Path [24] :Best Iter: [44] ELBO (-3690.271102) evaluations: (2879) 
#> Path [25] :Initial log joint density = -481926.667281 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.521e-03   2.839e-01    8.981e-01  8.981e-01      3237 -3.684e+03 -3.694e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3684.225750) evaluations: (3237) 
#> Path [26] :Initial log joint density = -481634.869795 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.781e-03   3.937e-01    5.423e-01  5.423e-01      3425 -3.685e+03 -3.698e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3685.244871) evaluations: (3425) 
#> Path [27] :Initial log joint density = -481479.923878 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      6.910e-03   2.247e-01    6.021e-01  6.021e-01      3902 -3.680e+03 -3.691e+03                   
#> Path [27] :Best Iter: [61] ELBO (-3680.344727) evaluations: (3902) 
#> Path [28] :Initial log joint density = -481537.164150 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.007e-02   1.960e-01    9.650e-01  9.650e-01      3675 -3.682e+03 -3.687e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3681.979964) evaluations: (3675) 
#> Path [29] :Initial log joint density = -481486.073937 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.448e-02   3.196e-01    1.000e+00  1.000e+00      3185 -3.687e+03 -3.692e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3686.584473) evaluations: (3185) 
#> Path [30] :Initial log joint density = -481767.323486 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.077e-02   2.611e-01    1.000e+00  1.000e+00      3086 -3.691e+03 -3.701e+03                   
#> Path [30] :Best Iter: [52] ELBO (-3691.011522) evaluations: (3086) 
#> Path [31] :Initial log joint density = -481849.041897 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.344e-03   2.155e-01    1.000e+00  1.000e+00      3290 -3.687e+03 -3.696e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3686.813227) evaluations: (3290) 
#> Path [32] :Initial log joint density = -481657.841355 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      8.691e-03   2.565e-01    9.133e-01  9.133e-01      3897 -3.683e+03 -3.693e+03                   
#> Path [32] :Best Iter: [61] ELBO (-3683.470967) evaluations: (3897) 
#> Path [33] :Initial log joint density = -482147.556691 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.501e-02   2.801e-01    1.000e+00  1.000e+00      3333 -3.683e+03 -3.684e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3682.830359) evaluations: (3333) 
#> Path [34] :Initial log joint density = -481845.645030 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.333e-03   2.937e-01    5.890e-01  5.890e-01      3356 -3.681e+03 -3.693e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3681.349294) evaluations: (3356) 
#> Path [35] :Initial log joint density = -483904.874570 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.286e-03   3.309e-01    5.662e-01  5.662e-01      3193 -3.685e+03 -3.695e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3685.434855) evaluations: (3193) 
#> Path [36] :Initial log joint density = -483398.369735 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.395e-02   3.479e-01    1.000e+00  1.000e+00      3208 -3.685e+03 -3.694e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3684.891822) evaluations: (3208) 
#> Path [37] :Initial log joint density = -481472.491622 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.589e-02   2.284e-01    1.000e+00  1.000e+00      3446 -3.681e+03 -3.688e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3681.149632) evaluations: (3446) 
#> Path [38] :Initial log joint density = -481558.638655 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.126e-02   2.883e-01    1.000e+00  1.000e+00      3187 -3.683e+03 -3.688e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3682.704578) evaluations: (3187) 
#> Path [39] :Initial log joint density = -481557.723373 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.598e-03   2.277e-01    1.000e+00  1.000e+00      3220 -3.682e+03 -3.687e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3681.727682) evaluations: (3220) 
#> Path [40] :Initial log joint density = -483190.127441 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.386e-03   2.708e-01    1.000e+00  1.000e+00      3192 -3.692e+03 -3.694e+03                   
#> Path [40] :Best Iter: [49] ELBO (-3692.356493) evaluations: (3192) 
#> Path [41] :Initial log joint density = -481729.965220 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.024e-03   2.597e-01    1.000e+00  1.000e+00      2828 -3.691e+03 -3.703e+03                   
#> Path [41] :Best Iter: [43] ELBO (-3690.754698) evaluations: (2828) 
#> Path [42] :Initial log joint density = -481747.909159 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.195e-03   2.997e-01    6.613e-01  6.613e-01      3073 -3.692e+03 -3.701e+03                   
#> Path [42] :Best Iter: [48] ELBO (-3692.150652) evaluations: (3073) 
#> Path [43] :Initial log joint density = -483006.968488 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.122e-03   3.433e-01    7.423e-01  7.423e-01      3464 -3.682e+03 -3.697e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3681.881056) evaluations: (3464) 
#> Path [44] :Initial log joint density = -481682.000121 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.681e-03   2.677e-01    4.977e-01  1.000e+00      2752 -3.690e+03 -3.703e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3690.492797) evaluations: (2752) 
#> Path [45] :Initial log joint density = -482063.959020 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.129e-03   1.761e-01    1.000e+00  1.000e+00      3710 -3.685e+03 -3.693e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3684.754154) evaluations: (3710) 
#> Path [46] :Initial log joint density = -481364.351217 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.684e-03   1.780e-01    9.628e-01  9.628e-01      3171 -3.692e+03 -3.691e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3690.775840) evaluations: (3171) 
#> Path [47] :Initial log joint density = -482175.520169 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.103e-02   3.971e-01    1.000e+00  1.000e+00      3738 -3.681e+03 -3.689e+03                   
#> Path [47] :Best Iter: [60] ELBO (-3681.369873) evaluations: (3738) 
#> Path [48] :Initial log joint density = -482079.072666 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.445e-02   2.980e-01    1.000e+00  1.000e+00      3760 -3.681e+03 -3.684e+03                   
#> Path [48] :Best Iter: [61] ELBO (-3680.858659) evaluations: (3760) 
#> Path [49] :Initial log joint density = -481559.357510 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.194e-03   1.535e-01    8.886e-01  8.886e-01      3117 -3.694e+03 -3.691e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3690.661483) evaluations: (3117) 
#> Path [50] :Initial log joint density = -481686.759965 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.629e-03   2.541e-01    1.000e+00  1.000e+00      3105 -3.686e+03 -3.686e+03                   
#> Path [50] :Best Iter: [55] ELBO (-3686.080241) evaluations: (3105) 
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
#> Path [1] :Initial log joint density = -427704.525907 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.191e-01   1.285e+03    2.761e-02  5.108e-02      8069 -3.287e+03 -5.923e+07                   
#> Path [1] :Best Iter: [45] ELBO (-3286.992056) evaluations: (8069) 
#> Path [2] :Initial log joint density = -427473.558206 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.038e-01   1.914e+03    2.883e-02  2.883e-02      8239 -3.283e+03 -2.615e+04                   
#> Path [2] :Best Iter: [50] ELBO (-3283.301404) evaluations: (8239) 
#> Path [3] :Initial log joint density = -427482.726710 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.982e-02   5.136e+03    2.093e-02  2.093e-02      8308 -3.288e+03 -6.194e+05                   
#> Path [3] :Best Iter: [44] ELBO (-3288.205943) evaluations: (8308) 
#> Path [4] :Initial log joint density = -427496.809011 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.221e-02   1.470e+03    4.245e-02  4.245e-02      8369 -3.281e+03 -8.926e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3281.473696) evaluations: (8369) 
#> Path [5] :Initial log joint density = -427640.519208 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.772e-02   5.396e+03    1.596e-02  1.596e-02      8115 -3.285e+03 -2.777e+06                   
#> Path [5] :Best Iter: [41] ELBO (-3284.979444) evaluations: (8115) 
#> Path [6] :Initial log joint density = -427859.913385 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.496e-02   2.698e+03    3.185e-02  3.185e-02      8153 -3.285e+03 -3.020e+05                   
#> Path [6] :Best Iter: [47] ELBO (-3284.869420) evaluations: (8153) 
#> Path [7] :Initial log joint density = -428640.639482 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.075e-02   1.465e+03    2.865e-02  5.048e-02      8461 -3.287e+03 -3.647e+04                   
#> Path [7] :Best Iter: [46] ELBO (-3286.598523) evaluations: (8461) 
#> Path [8] :Initial log joint density = -427518.247966 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.583e-02   3.171e+03    2.080e-02  2.080e-02      8351 -3.287e+03 -2.570e+08                   
#> Path [8] :Best Iter: [43] ELBO (-3286.800789) evaluations: (8351) 
#> Path [9] :Initial log joint density = -427520.947823 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.216e-02   2.304e+03    1.298e-02  2.636e-02      8204 -3.292e+03 -1.735e+05                   
#> Path [9] :Best Iter: [44] ELBO (-3291.978646) evaluations: (8204) 
#> Path [10] :Initial log joint density = -427363.978792 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.313e-02   3.972e+03    2.133e-02  2.133e-02      8299 -3.288e+03 -2.155e+05                   
#> Path [10] :Best Iter: [46] ELBO (-3287.963936) evaluations: (8299) 
#> Path [11] :Initial log joint density = -427284.904811 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.118e-02   1.903e+03    3.626e-02  3.626e-02      8308 -3.287e+03 -1.900e+04                   
#> Path [11] :Best Iter: [42] ELBO (-3287.103018) evaluations: (8308) 
#> Path [12] :Initial log joint density = -427540.342981 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.056e-02   2.024e+03    2.623e-02  2.623e-02      8140 -3.287e+03 -1.037e+05                   
#> Path [12] :Best Iter: [45] ELBO (-3287.206490) evaluations: (8140) 
#> Path [13] :Initial log joint density = -427838.374160 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.514e-02   1.942e+03    1.748e-02  6.219e-02      8332 -3.287e+03 -1.191e+06                   
#> Path [13] :Best Iter: [38] ELBO (-3287.293923) evaluations: (8332) 
#> Path [14] :Initial log joint density = -430292.123340 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.786e-02   3.080e+03    1.402e-02  1.402e-02      8191 -3.287e+03 -5.592e+05                   
#> Path [14] :Best Iter: [52] ELBO (-3286.943107) evaluations: (8191) 
#> Path [15] :Initial log joint density = -428442.843832 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.565e-02   1.607e+03    5.802e-02  5.802e-02      8123 -3.281e+03 -1.937e+04                   
#> Path [15] :Best Iter: [57] ELBO (-3281.350992) evaluations: (8123) 
#> Path [16] :Initial log joint density = -427819.013065 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.659e-02   3.019e+03    1.274e-02  1.274e-02      8424 -3.280e+03 -2.232e+05                   
#> Path [16] :Best Iter: [56] ELBO (-3279.717250) evaluations: (8424) 
#> Path [17] :Initial log joint density = -428020.332443 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.863e-02   1.968e+03    3.569e-02  3.569e-02      8156 -3.277e+03 -2.053e+04                   
#> Path [17] :Best Iter: [56] ELBO (-3276.690225) evaluations: (8156) 
#> Path [18] :Initial log joint density = -427797.173710 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.071e-02   3.716e+03    2.024e-02  2.024e-02      8186 -3.288e+03 -1.876e+04                   
#> Path [18] :Best Iter: [46] ELBO (-3287.805950) evaluations: (8186) 
#> Path [19] :Initial log joint density = -432390.726917 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.520e-02   1.528e+03    3.012e-02  5.134e-02      8383 -3.284e+03 -2.112e+04                   
#> Path [19] :Best Iter: [55] ELBO (-3283.525998) evaluations: (8383) 
#> Path [20] :Initial log joint density = -427939.878905 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.435e-02   1.749e+03    1.340e-02  3.394e-02      7960 -3.279e+03 -4.268e+06                   
#> Path [20] :Best Iter: [56] ELBO (-3279.039940) evaluations: (7960) 
#> Path [21] :Initial log joint density = -427956.586031 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.001e-01   7.573e+03    1.314e-02  2.586e-02      8388 -3.285e+03 -3.125e+08                   
#> Path [21] :Best Iter: [49] ELBO (-3285.050425) evaluations: (8388) 
#> Path [22] :Initial log joint density = -427779.453204 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.564e-02   1.507e+03    2.532e-02  4.443e-02      8284 -3.291e+03 -3.103e+05                   
#> Path [22] :Best Iter: [47] ELBO (-3290.940000) evaluations: (8284) 
#> Path [23] :Initial log joint density = -428806.829805 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.279e-02   1.578e+03    4.388e-02  4.388e-02      8109 -3.283e+03 -3.820e+04                   
#> Path [23] :Best Iter: [55] ELBO (-3283.085338) evaluations: (8109) 
#> Path [24] :Initial log joint density = -427413.526396 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.062e-02   1.312e+03    2.748e-02  2.748e-02      8083 -3.288e+03 -4.561e+04                   
#> Path [24] :Best Iter: [48] ELBO (-3287.584360) evaluations: (8083) 
#> Path [25] :Initial log joint density = -427471.877431 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.382e-02   2.053e+03    2.240e-02  4.548e-02      8060 -3.285e+03 -6.910e+04                   
#> Path [25] :Best Iter: [49] ELBO (-3285.255839) evaluations: (8060) 
#> Path [26] :Initial log joint density = -427854.924419 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.657e-02   2.782e+03    2.543e-02  2.543e-02      8298 -3.287e+03 -3.402e+04                   
#> Path [26] :Best Iter: [44] ELBO (-3286.688246) evaluations: (8298) 
#> Path [27] :Initial log joint density = -427800.207539 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.930e-02   2.756e+03    1.741e-02  1.741e-02      7931 -3.284e+03 -1.183e+04                   
#> Path [27] :Best Iter: [46] ELBO (-3283.905547) evaluations: (7931) 
#> Path [28] :Initial log joint density = -428256.530822 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.366e-02   1.911e+03    2.770e-02  2.770e-02      8111 -3.284e+03 -5.940e+03                   
#> Path [28] :Best Iter: [47] ELBO (-3283.809776) evaluations: (8111) 
#> Path [29] :Initial log joint density = -427833.949622 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.297e-02   6.511e+03    1.204e-02  1.204e-02      8452 -3.284e+03 -3.110e+06                   
#> Path [29] :Best Iter: [44] ELBO (-3283.732331) evaluations: (8452) 
#> Path [30] :Initial log joint density = -429050.247954 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.262e-02   2.857e+03    3.619e-02  3.619e-02      8156 -3.284e+03 -6.999e+04                   
#> Path [30] :Best Iter: [50] ELBO (-3284.498966) evaluations: (8156) 
#> Path [31] :Initial log joint density = -427589.478775 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.165e-02   1.879e+03    6.352e-02  6.352e-02      8169 -3.285e+03 -1.338e+04                   
#> Path [31] :Best Iter: [55] ELBO (-3285.193843) evaluations: (8169) 
#> Path [32] :Initial log joint density = -427505.314937 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.267e-01   6.550e+03    5.228e-02  5.228e-02      8310 -3.289e+03 -1.718e+04                   
#> Path [32] :Best Iter: [40] ELBO (-3288.746408) evaluations: (8310) 
#> Path [33] :Initial log joint density = -427823.157241 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.535e-02   1.201e+03    4.656e-02  4.656e-02      7961 -3.281e+03 -8.945e+04                   
#> Path [33] :Best Iter: [55] ELBO (-3281.194104) evaluations: (7961) 
#> Path [34] :Initial log joint density = -428787.526344 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.904e-02   1.907e+03    3.984e-02  3.984e-02      8247 -3.287e+03 -4.771e+04                   
#> Path [34] :Best Iter: [49] ELBO (-3287.123475) evaluations: (8247) 
#> Path [35] :Initial log joint density = -428150.915481 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.485e-02   8.907e+02    2.503e-02  2.503e-02      8180 -3.276e+03 -3.273e+05                   
#> Path [35] :Best Iter: [58] ELBO (-3275.653398) evaluations: (8180) 
#> Path [36] :Initial log joint density = -427521.121801 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.362e-02   1.778e+03    4.187e-02  1.031e-01      8189 -3.284e+03 -1.033e+05                   
#> Path [36] :Best Iter: [47] ELBO (-3284.370986) evaluations: (8189) 
#> Path [37] :Initial log joint density = -427502.037677 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.272e-02   2.092e+03    3.262e-02  3.262e-02      8385 -3.283e+03 -1.177e+05                   
#> Path [37] :Best Iter: [44] ELBO (-3282.721884) evaluations: (8385) 
#> Path [38] :Initial log joint density = -427578.034741 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.742e-02   3.798e+03    2.287e-02  2.287e-02      8501 -3.286e+03 -1.011e+05                   
#> Path [38] :Best Iter: [43] ELBO (-3286.394075) evaluations: (8501) 
#> Path [39] :Initial log joint density = -427576.992114 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.894e-02   2.937e+03    2.416e-02  2.416e-02      8349 -3.286e+03 -7.952e+07                   
#> Path [39] :Best Iter: [48] ELBO (-3285.667022) evaluations: (8349) 
#> Path [40] :Initial log joint density = -428595.241928 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.338e-02   4.444e+03    2.620e-02  2.620e-02      8501 -3.279e+03 -7.421e+05                   
#> Path [40] :Best Iter: [46] ELBO (-3279.323019) evaluations: (8501) 
#> Path [41] :Initial log joint density = -427649.859707 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.753e-02   9.369e+02    5.285e-02  1.123e-01      8104 -3.281e+03 -1.165e+04                   
#> Path [41] :Best Iter: [59] ELBO (-3280.767890) evaluations: (8104) 
#> Path [42] :Initial log joint density = -427789.441832 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.751e-02   3.550e+03    1.439e-02  3.340e-02      8287 -3.287e+03 -3.354e+05                   
#> Path [42] :Best Iter: [39] ELBO (-3287.342680) evaluations: (8287) 
#> Path [43] :Initial log joint density = -427762.127651 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.796e-02   1.465e+03    3.173e-02  3.173e-02      8035 -3.281e+03 -4.773e+06                   
#> Path [43] :Best Iter: [55] ELBO (-3281.383206) evaluations: (8035) 
#> Path [44] :Initial log joint density = -431070.741479 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      2.476e-02   3.286e+03    1.723e-02  1.723e-02      8431 -3.283e+03 -9.882e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3282.704473) evaluations: (8431) 
#> Path [45] :Initial log joint density = -427718.695429 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.114e-02   3.895e+03    1.208e-02  2.248e-02      8218 -3.287e+03 -1.675e+04                   
#> Path [45] :Best Iter: [46] ELBO (-3286.614427) evaluations: (8218) 
#> Path [46] :Initial log joint density = -427668.810895 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.093e-02   2.728e+03    1.733e-02  1.733e-02      8465 -3.284e+03 -6.609e+06                   
#> Path [46] :Best Iter: [44] ELBO (-3283.699578) evaluations: (8465) 
#> Path [47] :Initial log joint density = -427603.507783 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.162e-02   1.148e+03    5.875e-02  5.875e-02      8122 -3.283e+03 -4.136e+04                   
#> Path [47] :Best Iter: [43] ELBO (-3282.829316) evaluations: (8122) 
#> Path [48] :Initial log joint density = -429943.554657 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      2.188e-02   2.014e+03    3.552e-02  3.552e-02      8214 -3.283e+03 -6.380e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3283.295969) evaluations: (8214) 
#> Path [49] :Initial log joint density = -427618.576695 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.104e-01   3.494e+03    2.495e-02  2.495e-02      8228 -3.288e+03 -1.191e+06                   
#> Path [49] :Best Iter: [43] ELBO (-3288.128183) evaluations: (8228) 
#> Path [50] :Initial log joint density = -428258.780652 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.763e-02   1.345e+03    3.071e-02  5.876e-02      8338 -3.277e+03 -1.992e+04                   
#> Path [50] :Best Iter: [56] ELBO (-3276.600617) evaluations: (8338) 
#> Finished in  30.2 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -427704.525907 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.191e-01   1.285e+03    2.761e-02  5.108e-02      8069 -3.287e+03 -5.923e+07                   
#> Path [1] :Best Iter: [45] ELBO (-3286.992056) evaluations: (8069) 
#> Path [2] :Initial log joint density = -427473.558206 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.038e-01   1.914e+03    2.883e-02  2.883e-02      8239 -3.283e+03 -2.615e+04                   
#> Path [2] :Best Iter: [50] ELBO (-3283.301404) evaluations: (8239) 
#> Path [3] :Initial log joint density = -427482.726710 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.982e-02   5.136e+03    2.093e-02  2.093e-02      8308 -3.288e+03 -6.194e+05                   
#> Path [3] :Best Iter: [44] ELBO (-3288.205943) evaluations: (8308) 
#> Path [4] :Initial log joint density = -427496.809011 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.221e-02   1.470e+03    4.245e-02  4.245e-02      8369 -3.281e+03 -8.926e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3281.473696) evaluations: (8369) 
#> Path [5] :Initial log joint density = -427640.519208 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.772e-02   5.396e+03    1.596e-02  1.596e-02      8115 -3.285e+03 -2.777e+06                   
#> Path [5] :Best Iter: [41] ELBO (-3284.979444) evaluations: (8115) 
#> Path [6] :Initial log joint density = -427859.913385 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.496e-02   2.698e+03    3.185e-02  3.185e-02      8153 -3.285e+03 -3.020e+05                   
#> Path [6] :Best Iter: [47] ELBO (-3284.869420) evaluations: (8153) 
#> Path [7] :Initial log joint density = -428640.639482 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.075e-02   1.465e+03    2.865e-02  5.048e-02      8461 -3.287e+03 -3.647e+04                   
#> Path [7] :Best Iter: [46] ELBO (-3286.598523) evaluations: (8461) 
#> Path [8] :Initial log joint density = -427518.247966 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.583e-02   3.171e+03    2.080e-02  2.080e-02      8351 -3.287e+03 -2.570e+08                   
#> Path [8] :Best Iter: [43] ELBO (-3286.800789) evaluations: (8351) 
#> Path [9] :Initial log joint density = -427520.947823 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.216e-02   2.304e+03    1.298e-02  2.636e-02      8204 -3.292e+03 -1.735e+05                   
#> Path [9] :Best Iter: [44] ELBO (-3291.978646) evaluations: (8204) 
#> Path [10] :Initial log joint density = -427363.978792 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.313e-02   3.972e+03    2.133e-02  2.133e-02      8299 -3.288e+03 -2.155e+05                   
#> Path [10] :Best Iter: [46] ELBO (-3287.963936) evaluations: (8299) 
#> Path [11] :Initial log joint density = -427284.904811 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.118e-02   1.903e+03    3.626e-02  3.626e-02      8308 -3.287e+03 -1.900e+04                   
#> Path [11] :Best Iter: [42] ELBO (-3287.103018) evaluations: (8308) 
#> Path [12] :Initial log joint density = -427540.342981 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.056e-02   2.024e+03    2.623e-02  2.623e-02      8140 -3.287e+03 -1.037e+05                   
#> Path [12] :Best Iter: [45] ELBO (-3287.206490) evaluations: (8140) 
#> Path [13] :Initial log joint density = -427838.374160 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.514e-02   1.942e+03    1.748e-02  6.219e-02      8332 -3.287e+03 -1.191e+06                   
#> Path [13] :Best Iter: [38] ELBO (-3287.293923) evaluations: (8332) 
#> Path [14] :Initial log joint density = -430292.123340 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.786e-02   3.080e+03    1.402e-02  1.402e-02      8191 -3.287e+03 -5.592e+05                   
#> Path [14] :Best Iter: [52] ELBO (-3286.943107) evaluations: (8191) 
#> Path [15] :Initial log joint density = -428442.843832 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.565e-02   1.607e+03    5.802e-02  5.802e-02      8123 -3.281e+03 -1.937e+04                   
#> Path [15] :Best Iter: [57] ELBO (-3281.350992) evaluations: (8123) 
#> Path [16] :Initial log joint density = -427819.013065 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.659e-02   3.019e+03    1.274e-02  1.274e-02      8424 -3.280e+03 -2.232e+05                   
#> Path [16] :Best Iter: [56] ELBO (-3279.717250) evaluations: (8424) 
#> Path [17] :Initial log joint density = -428020.332443 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.863e-02   1.968e+03    3.569e-02  3.569e-02      8156 -3.277e+03 -2.053e+04                   
#> Path [17] :Best Iter: [56] ELBO (-3276.690225) evaluations: (8156) 
#> Path [18] :Initial log joint density = -427797.173710 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.071e-02   3.716e+03    2.024e-02  2.024e-02      8186 -3.288e+03 -1.876e+04                   
#> Path [18] :Best Iter: [46] ELBO (-3287.805950) evaluations: (8186) 
#> Path [19] :Initial log joint density = -432390.726917 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.520e-02   1.528e+03    3.012e-02  5.134e-02      8383 -3.284e+03 -2.112e+04                   
#> Path [19] :Best Iter: [55] ELBO (-3283.525998) evaluations: (8383) 
#> Path [20] :Initial log joint density = -427939.878905 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.435e-02   1.749e+03    1.340e-02  3.394e-02      7960 -3.279e+03 -4.268e+06                   
#> Path [20] :Best Iter: [56] ELBO (-3279.039940) evaluations: (7960) 
#> Path [21] :Initial log joint density = -427956.586031 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.001e-01   7.573e+03    1.314e-02  2.586e-02      8388 -3.285e+03 -3.125e+08                   
#> Path [21] :Best Iter: [49] ELBO (-3285.050425) evaluations: (8388) 
#> Path [22] :Initial log joint density = -427779.453204 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.564e-02   1.507e+03    2.532e-02  4.443e-02      8284 -3.291e+03 -3.103e+05                   
#> Path [22] :Best Iter: [47] ELBO (-3290.940000) evaluations: (8284) 
#> Path [23] :Initial log joint density = -428806.829805 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.279e-02   1.578e+03    4.388e-02  4.388e-02      8109 -3.283e+03 -3.820e+04                   
#> Path [23] :Best Iter: [55] ELBO (-3283.085338) evaluations: (8109) 
#> Path [24] :Initial log joint density = -427413.526396 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.062e-02   1.312e+03    2.748e-02  2.748e-02      8083 -3.288e+03 -4.561e+04                   
#> Path [24] :Best Iter: [48] ELBO (-3287.584360) evaluations: (8083) 
#> Path [25] :Initial log joint density = -427471.877431 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.382e-02   2.053e+03    2.240e-02  4.548e-02      8060 -3.285e+03 -6.910e+04                   
#> Path [25] :Best Iter: [49] ELBO (-3285.255839) evaluations: (8060) 
#> Path [26] :Initial log joint density = -427854.924419 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.657e-02   2.782e+03    2.543e-02  2.543e-02      8298 -3.287e+03 -3.402e+04                   
#> Path [26] :Best Iter: [44] ELBO (-3286.688246) evaluations: (8298) 
#> Path [27] :Initial log joint density = -427800.207539 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.930e-02   2.756e+03    1.741e-02  1.741e-02      7931 -3.284e+03 -1.183e+04                   
#> Path [27] :Best Iter: [46] ELBO (-3283.905547) evaluations: (7931) 
#> Path [28] :Initial log joint density = -428256.530822 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.366e-02   1.911e+03    2.770e-02  2.770e-02      8111 -3.284e+03 -5.940e+03                   
#> Path [28] :Best Iter: [47] ELBO (-3283.809776) evaluations: (8111) 
#> Path [29] :Initial log joint density = -427833.949622 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.297e-02   6.511e+03    1.204e-02  1.204e-02      8452 -3.284e+03 -3.110e+06                   
#> Path [29] :Best Iter: [44] ELBO (-3283.732331) evaluations: (8452) 
#> Path [30] :Initial log joint density = -429050.247954 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.262e-02   2.857e+03    3.619e-02  3.619e-02      8156 -3.284e+03 -6.999e+04                   
#> Path [30] :Best Iter: [50] ELBO (-3284.498966) evaluations: (8156) 
#> Path [31] :Initial log joint density = -427589.478775 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      7.165e-02   1.879e+03    6.352e-02  6.352e-02      8169 -3.285e+03 -1.338e+04                   
#> Path [31] :Best Iter: [55] ELBO (-3285.193843) evaluations: (8169) 
#> Path [32] :Initial log joint density = -427505.314937 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.267e-01   6.550e+03    5.228e-02  5.228e-02      8310 -3.289e+03 -1.718e+04                   
#> Path [32] :Best Iter: [40] ELBO (-3288.746408) evaluations: (8310) 
#> Path [33] :Initial log joint density = -427823.157241 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.535e-02   1.201e+03    4.656e-02  4.656e-02      7961 -3.281e+03 -8.945e+04                   
#> Path [33] :Best Iter: [55] ELBO (-3281.194104) evaluations: (7961) 
#> Path [34] :Initial log joint density = -428787.526344 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.904e-02   1.907e+03    3.984e-02  3.984e-02      8247 -3.287e+03 -4.771e+04                   
#> Path [34] :Best Iter: [49] ELBO (-3287.123475) evaluations: (8247) 
#> Path [35] :Initial log joint density = -428150.915481 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.485e-02   8.907e+02    2.503e-02  2.503e-02      8180 -3.276e+03 -3.273e+05                   
#> Path [35] :Best Iter: [58] ELBO (-3275.653398) evaluations: (8180) 
#> Path [36] :Initial log joint density = -427521.121801 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.362e-02   1.778e+03    4.187e-02  1.031e-01      8189 -3.284e+03 -1.033e+05                   
#> Path [36] :Best Iter: [47] ELBO (-3284.370986) evaluations: (8189) 
#> Path [37] :Initial log joint density = -427502.037677 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.272e-02   2.092e+03    3.262e-02  3.262e-02      8385 -3.283e+03 -1.177e+05                   
#> Path [37] :Best Iter: [44] ELBO (-3282.721884) evaluations: (8385) 
#> Path [38] :Initial log joint density = -427578.034741 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      4.742e-02   3.798e+03    2.287e-02  2.287e-02      8501 -3.286e+03 -1.011e+05                   
#> Path [38] :Best Iter: [43] ELBO (-3286.394075) evaluations: (8501) 
#> Path [39] :Initial log joint density = -427576.992114 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.894e-02   2.937e+03    2.416e-02  2.416e-02      8349 -3.286e+03 -7.952e+07                   
#> Path [39] :Best Iter: [48] ELBO (-3285.667022) evaluations: (8349) 
#> Path [40] :Initial log joint density = -428595.241928 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      9.338e-02   4.444e+03    2.620e-02  2.620e-02      8501 -3.279e+03 -7.421e+05                   
#> Path [40] :Best Iter: [46] ELBO (-3279.323019) evaluations: (8501) 
#> Path [41] :Initial log joint density = -427649.859707 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.753e-02   9.369e+02    5.285e-02  1.123e-01      8104 -3.281e+03 -1.165e+04                   
#> Path [41] :Best Iter: [59] ELBO (-3280.767890) evaluations: (8104) 
#> Path [42] :Initial log joint density = -427789.441832 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.751e-02   3.550e+03    1.439e-02  3.340e-02      8287 -3.287e+03 -3.354e+05                   
#> Path [42] :Best Iter: [39] ELBO (-3287.342680) evaluations: (8287) 
#> Path [43] :Initial log joint density = -427762.127651 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      6.796e-02   1.465e+03    3.173e-02  3.173e-02      8035 -3.281e+03 -4.773e+06                   
#> Path [43] :Best Iter: [55] ELBO (-3281.383206) evaluations: (8035) 
#> Path [44] :Initial log joint density = -431070.741479 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      2.476e-02   3.286e+03    1.723e-02  1.723e-02      8431 -3.283e+03 -9.882e+03                   
#> Path [44] :Best Iter: [45] ELBO (-3282.704473) evaluations: (8431) 
#> Path [45] :Initial log joint density = -427718.695429 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      3.114e-02   3.895e+03    1.208e-02  2.248e-02      8218 -3.287e+03 -1.675e+04                   
#> Path [45] :Best Iter: [46] ELBO (-3286.614427) evaluations: (8218) 
#> Path [46] :Initial log joint density = -427668.810895 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.093e-02   2.728e+03    1.733e-02  1.733e-02      8465 -3.284e+03 -6.609e+06                   
#> Path [46] :Best Iter: [44] ELBO (-3283.699578) evaluations: (8465) 
#> Path [47] :Initial log joint density = -427603.507783 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      8.162e-02   1.148e+03    5.875e-02  5.875e-02      8122 -3.283e+03 -4.136e+04                   
#> Path [47] :Best Iter: [43] ELBO (-3282.829316) evaluations: (8122) 
#> Path [48] :Initial log joint density = -429943.554657 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      2.188e-02   2.014e+03    3.552e-02  3.552e-02      8214 -3.283e+03 -6.380e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3283.295969) evaluations: (8214) 
#> Path [49] :Initial log joint density = -427618.576695 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      1.104e-01   3.494e+03    2.495e-02  2.495e-02      8228 -3.288e+03 -1.191e+06                   
#> Path [49] :Best Iter: [43] ELBO (-3288.128183) evaluations: (8228) 
#> Path [50] :Initial log joint density = -428258.780652 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.247e+05      5.763e-02   1.345e+03    3.071e-02  5.876e-02      8338 -3.277e+03 -1.992e+04                   
#> Path [50] :Best Iter: [56] ELBO (-3276.600617) evaluations: (8338) 
#> Finished in  27.4 seconds.
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
