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
#> Chain 1  Elapsed Time: 4.605 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier identification - step 1/2
#> Loading model from cache...
#> Path [1] :Initial log joint density = -425167.600455 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.078e-01   8.233e+02    3.160e-02  8.215e-02      8219 -3.270e+03 -1.798e+05                   
#> Path [1] :Best Iter: [57] ELBO (-3269.930244) evaluations: (8219) 
#> Path [2] :Initial log joint density = -424929.009268 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.018e-02   2.257e+03    4.463e-02  4.463e-02      8259 -3.279e+03 -5.540e+05                   
#> Path [2] :Best Iter: [45] ELBO (-3278.905239) evaluations: (8259) 
#> Path [3] :Initial log joint density = -425091.927347 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.443e-02   3.595e+03    1.282e-02  1.282e-02      8144 -3.276e+03 -6.425e+05                   
#> Path [3] :Best Iter: [44] ELBO (-3276.385489) evaluations: (8144) 
#> Path [4] :Initial log joint density = -425480.485526 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.701e-02   2.100e+03    3.129e-02  3.129e-02      8363 -3.277e+03 -1.027e+04                   
#> Path [4] :Best Iter: [41] ELBO (-3277.314365) evaluations: (8363) 
#> Path [5] :Initial log joint density = -425399.494654 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.061e-01   3.871e+03    3.094e-02  3.094e-02      8257 -3.275e+03 -6.898e+05                   
#> Path [5] :Best Iter: [48] ELBO (-3275.244909) evaluations: (8257) 
#> Path [6] :Initial log joint density = -425625.722353 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.759e-02   8.921e+02    3.837e-02  3.837e-02      8076 -3.266e+03 -5.318e+04                   
#> Path [6] :Best Iter: [58] ELBO (-3266.095730) evaluations: (8076) 
#> Path [7] :Initial log joint density = -425076.027190 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.623e-01   1.983e+03    5.422e-02  5.422e-02      8431 -3.276e+03 -5.987e+04                   
#> Path [7] :Best Iter: [50] ELBO (-3276.399202) evaluations: (8431) 
#> Path [8] :Initial log joint density = -429647.890004 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.514e-02   1.107e+03    1.985e-02  1.985e-02      8173 -3.273e+03 -3.276e+04                   
#> Path [8] :Best Iter: [55] ELBO (-3272.525399) evaluations: (8173) 
#> Path [9] :Initial log joint density = -425083.075988 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.969e-02   1.686e+03    2.641e-02  6.656e-02      8081 -3.268e+03 -1.030e+04                   
#> Path [9] :Best Iter: [57] ELBO (-3267.979248) evaluations: (8081) 
#> Path [10] :Initial log joint density = -425146.073911 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.546e-02   1.701e+03    2.180e-02  2.180e-02      8321 -3.282e+03 -3.941e+04                   
#> Path [10] :Best Iter: [43] ELBO (-3281.558167) evaluations: (8321) 
#> Path [11] :Initial log joint density = -425196.739827 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.704e-02   2.059e+03    3.279e-02  3.279e-02      8346 -3.278e+03 -6.944e+04                   
#> Path [11] :Best Iter: [45] ELBO (-3278.393583) evaluations: (8346) 
#> Path [12] :Initial log joint density = -425544.809591 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.859e-02   8.921e+02    1.789e-02  1.789e-02      8067 -3.267e+03 -4.508e+04                   
#> Path [12] :Best Iter: [55] ELBO (-3267.310464) evaluations: (8067) 
#> Path [13] :Initial log joint density = -425111.987875 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.204e-02   2.910e+03    3.684e-02  3.684e-02      8196 -3.278e+03 -1.841e+04                   
#> Path [13] :Best Iter: [42] ELBO (-3278.132009) evaluations: (8196) 
#> Path [14] :Initial log joint density = -425427.622423 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.252e-01   1.688e+03    4.269e-02  4.269e-02      8096 -3.271e+03 -2.405e+05                   
#> Path [14] :Best Iter: [61] ELBO (-3271.247131) evaluations: (8096) 
#> Path [15] :Initial log joint density = -425217.083745 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.640e-02   1.979e+03    4.289e-02  4.289e-02      8084 -3.280e+03 -1.653e+04                   
#> Path [15] :Best Iter: [45] ELBO (-3280.483533) evaluations: (8084) 
#> Path [16] :Initial log joint density = -425070.207399 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.685e-02   2.046e+03    1.155e-02  2.903e-02      8416 -3.273e+03 -1.749e+10                   
#> Path [16] :Best Iter: [55] ELBO (-3272.765746) evaluations: (8416) 
#> Path [17] :Initial log joint density = -425109.414170 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.280e-02   2.471e+03    1.710e-02  1.710e-02      8484 -3.281e+03 -1.661e+06                   
#> Path [17] :Best Iter: [47] ELBO (-3280.626813) evaluations: (8484) 
#> Path [18] :Initial log joint density = -426712.294113 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.982e-02   9.034e+02    3.236e-02  3.236e-02      8106 -3.269e+03 -6.156e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3268.830040) evaluations: (8106) 
#> Path [19] :Initial log joint density = -426356.465513 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.457e-01   1.440e+03    2.537e-02  5.748e-02      8264 -3.275e+03 -1.010e+08                   
#> Path [19] :Best Iter: [49] ELBO (-3274.616654) evaluations: (8264) 
#> Path [20] :Initial log joint density = -425064.741198 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      4.098e-02   2.107e+03    1.685e-02  1.685e-02      8438 -3.276e+03 -4.380e+05                   
#> Path [20] :Best Iter: [43] ELBO (-3276.493440) evaluations: (8438) 
#> Path [21] :Initial log joint density = -425129.230143 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.743e-02   2.416e+03    2.165e-02  2.165e-02      8456 -3.277e+03 -4.578e+06                   
#> Path [21] :Best Iter: [47] ELBO (-3276.794124) evaluations: (8456) 
#> Path [22] :Initial log joint density = -426779.460887 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.217e-02   1.867e+03    3.248e-02  3.248e-02      8090 -3.276e+03 -7.320e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3276.315152) evaluations: (8090) 
#> Path [23] :Initial log joint density = -428751.524216 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.136e-01   8.141e+02    3.577e-02  3.577e-02      8242 -3.267e+03 -1.022e+05                   
#> Path [23] :Best Iter: [55] ELBO (-3267.463148) evaluations: (8242) 
#> Path [24] :Initial log joint density = -425035.826219 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.127e-02   2.844e+03    2.492e-02  4.616e-02      8305 -3.279e+03 -6.370e+05                   
#> Path [24] :Best Iter: [40] ELBO (-3279.313885) evaluations: (8305) 
#> Path [25] :Initial log joint density = -425049.711298 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.019e-02   2.094e+03    2.093e-02  2.093e-02      7951 -3.276e+03 -1.526e+04                   
#> Path [25] :Best Iter: [55] ELBO (-3275.684482) evaluations: (7951) 
#> Path [26] :Initial log joint density = -425187.673867 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.415e-02   2.083e+03    2.862e-02  2.862e-02      7992 -3.277e+03 -9.165e+03                   
#> Path [26] :Best Iter: [48] ELBO (-3276.500864) evaluations: (7992) 
#> Path [27] :Initial log joint density = -425653.274204 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.488e-02   2.978e+03    4.269e-02  4.269e-02      8265 -3.270e+03 -1.536e+04                   
#> Path [27] :Best Iter: [61] ELBO (-3270.448662) evaluations: (8265) 
#> Path [28] :Initial log joint density = -427906.376757 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.200e-01   1.514e+03    2.332e-02  2.332e-02      8044 -3.271e+03 -3.340e+05                   
#> Path [28] :Best Iter: [56] ELBO (-3271.035147) evaluations: (8044) 
#> Path [29] :Initial log joint density = -425197.486246 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.879e-02   2.940e+03    2.613e-02  2.613e-02      8076 -3.278e+03 -3.090e+05                   
#> Path [29] :Best Iter: [45] ELBO (-3278.067981) evaluations: (8076) 
#> Path [30] :Initial log joint density = -425506.357236 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.743e-02   7.099e+02    3.947e-02  3.947e-02      8158 -3.271e+03 -1.876e+04                   
#> Path [30] :Best Iter: [57] ELBO (-3270.657641) evaluations: (8158) 
#> Path [31] :Initial log joint density = -425313.260086 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.065e-02   2.684e+03    4.257e-02  4.257e-02      8315 -3.275e+03 -1.149e+06                   
#> Path [31] :Best Iter: [39] ELBO (-3275.272492) evaluations: (8315) 
#> Path [32] :Initial log joint density = -425397.663863 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.520e-02   4.399e+03    1.569e-02  2.870e-02      8246 -3.272e+03 -4.501e+06                   
#> Path [32] :Best Iter: [55] ELBO (-3272.209232) evaluations: (8246) 
#> Path [33] :Initial log joint density = -428404.377487 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.455e-02   1.647e+03    2.522e-02  2.522e-02      8049 -3.270e+03 -6.139e+04                   
#> Path [33] :Best Iter: [57] ELBO (-3269.702908) evaluations: (8049) 
#> Path [34] :Initial log joint density = -425317.204225 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.122e-01   3.151e+03    2.909e-02  2.909e-02      8340 -3.276e+03 -3.705e+07                   
#> Path [34] :Best Iter: [48] ELBO (-3275.629808) evaluations: (8340) 
#> Path [35] :Initial log joint density = -425487.563038 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.036e-02   1.299e+03    3.962e-02  3.962e-02      7994 -3.273e+03 -9.131e+04                   
#> Path [35] :Best Iter: [59] ELBO (-3272.906230) evaluations: (7994) 
#> Path [36] :Initial log joint density = -425192.252393 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.649e-02   4.114e+03    3.212e-02  6.211e-02      8199 -3.274e+03 -3.910e+06                   
#> Path [36] :Best Iter: [46] ELBO (-3273.635565) evaluations: (8199) 
#> Path [37] :Initial log joint density = -425068.833182 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.511e-02   1.704e+03    2.129e-02  4.544e-02      8027 -3.269e+03 -1.593e+05                   
#> Path [37] :Best Iter: [58] ELBO (-3268.702968) evaluations: (8027) 
#> Path [38] :Initial log joint density = -425219.739918 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.168e-01   9.357e+02    5.580e-02  9.294e-02      8183 -3.267e+03 -1.465e+06                   
#> Path [38] :Best Iter: [60] ELBO (-3267.394392) evaluations: (8183) 
#> Path [39] :Initial log joint density = -425158.132710 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      9.410e-02   4.523e+03    1.666e-02  1.666e-02      8233 -3.276e+03 -3.804e+11                   
#> Path [39] :Best Iter: [48] ELBO (-3276.026330) evaluations: (8233) 
#> Path [40] :Initial log joint density = -425220.022429 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.175e-02   2.825e+03    3.069e-02  3.069e-02      8394 -3.278e+03 -6.119e+06                   
#> Path [40] :Best Iter: [44] ELBO (-3277.782265) evaluations: (8394) 
#> Path [41] :Initial log joint density = -425183.518875 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.478e-02   1.113e+03    2.116e-02  4.095e-02      8206 -3.271e+03 -9.680e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3270.625599) evaluations: (8206) 
#> Path [42] :Initial log joint density = -425194.696697 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.962e-02   2.197e+03    2.780e-02  2.780e-02      8221 -3.279e+03 -8.988e+04                   
#> Path [42] :Best Iter: [45] ELBO (-3278.848280) evaluations: (8221) 
#> Path [43] :Initial log joint density = -425576.682640 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.103e-01   8.533e+02    4.296e-02  4.296e-02      8026 -3.272e+03 -1.254e+04                   
#> Path [43] :Best Iter: [56] ELBO (-3271.575119) evaluations: (8026) 
#> Path [44] :Initial log joint density = -425351.085126 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.621e-02   2.265e+03    2.070e-02  4.704e-02      8184 -3.274e+03 -1.069e+07                   
#> Path [44] :Best Iter: [54] ELBO (-3273.640880) evaluations: (8184) 
#> Path [45] :Initial log joint density = -425130.389110 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.193e-02   7.574e+03    1.630e-02  1.630e-02      8492 -3.278e+03 -1.378e+06                   
#> Path [45] :Best Iter: [47] ELBO (-3277.593672) evaluations: (8492) 
#> Path [46] :Initial log joint density = -426990.109503 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      4.323e-02   8.759e+02    2.822e-02  2.822e-02      8192 -3.269e+03 -2.322e+04                   
#> Path [46] :Best Iter: [58] ELBO (-3269.142015) evaluations: (8192) 
#> Path [47] :Initial log joint density = -425083.636177 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.155e-02   1.965e+03    2.204e-02  2.204e-02      8219 -3.278e+03 -2.092e+05                   
#> Path [47] :Best Iter: [49] ELBO (-3277.758546) evaluations: (8219) 
#> Path [48] :Initial log joint density = -424986.953034 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.458e-02   2.932e+03    2.889e-02  2.889e-02      8213 -3.277e+03 -3.016e+04                   
#> Path [48] :Best Iter: [55] ELBO (-3277.155453) evaluations: (8213) 
#> Path [49] :Initial log joint density = -425166.660795 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.974e-02   3.061e+03    2.461e-02  2.461e-02      8218 -3.279e+03 -2.593e+05                   
#> Path [49] :Best Iter: [42] ELBO (-3278.946856) evaluations: (8218) 
#> Path [50] :Initial log joint density = -425197.603897 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      9.307e-02   2.747e+03    1.869e-02  1.869e-02      8126 -3.278e+03 -1.232e+07                   
#> Path [50] :Best Iter: [45] ELBO (-3278.094882) evaluations: (8126) 
#> Finished in  30.7 seconds.
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 23.249 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> sccomp says: outlier-free model fitting - step 2/2
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -425167.600455 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.078e-01   8.233e+02    3.160e-02  8.215e-02      8219 -3.270e+03 -1.798e+05                   
#> Path [1] :Best Iter: [57] ELBO (-3269.930244) evaluations: (8219) 
#> Path [2] :Initial log joint density = -424929.009268 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.018e-02   2.257e+03    4.463e-02  4.463e-02      8259 -3.279e+03 -5.540e+05                   
#> Path [2] :Best Iter: [45] ELBO (-3278.905239) evaluations: (8259) 
#> Path [3] :Initial log joint density = -425091.927347 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.443e-02   3.595e+03    1.282e-02  1.282e-02      8144 -3.276e+03 -6.425e+05                   
#> Path [3] :Best Iter: [44] ELBO (-3276.385489) evaluations: (8144) 
#> Path [4] :Initial log joint density = -425480.485526 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.701e-02   2.100e+03    3.129e-02  3.129e-02      8363 -3.277e+03 -1.027e+04                   
#> Path [4] :Best Iter: [41] ELBO (-3277.314365) evaluations: (8363) 
#> Path [5] :Initial log joint density = -425399.494654 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.061e-01   3.871e+03    3.094e-02  3.094e-02      8257 -3.275e+03 -6.898e+05                   
#> Path [5] :Best Iter: [48] ELBO (-3275.244909) evaluations: (8257) 
#> Path [6] :Initial log joint density = -425625.722353 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.759e-02   8.921e+02    3.837e-02  3.837e-02      8076 -3.266e+03 -5.318e+04                   
#> Path [6] :Best Iter: [58] ELBO (-3266.095730) evaluations: (8076) 
#> Path [7] :Initial log joint density = -425076.027190 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.623e-01   1.983e+03    5.422e-02  5.422e-02      8431 -3.276e+03 -5.987e+04                   
#> Path [7] :Best Iter: [50] ELBO (-3276.399202) evaluations: (8431) 
#> Path [8] :Initial log joint density = -429647.890004 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.514e-02   1.107e+03    1.985e-02  1.985e-02      8173 -3.273e+03 -3.276e+04                   
#> Path [8] :Best Iter: [55] ELBO (-3272.525399) evaluations: (8173) 
#> Path [9] :Initial log joint density = -425083.075988 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.969e-02   1.686e+03    2.641e-02  6.656e-02      8081 -3.268e+03 -1.030e+04                   
#> Path [9] :Best Iter: [57] ELBO (-3267.979248) evaluations: (8081) 
#> Path [10] :Initial log joint density = -425146.073911 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.546e-02   1.701e+03    2.180e-02  2.180e-02      8321 -3.282e+03 -3.941e+04                   
#> Path [10] :Best Iter: [43] ELBO (-3281.558167) evaluations: (8321) 
#> Path [11] :Initial log joint density = -425196.739827 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.704e-02   2.059e+03    3.279e-02  3.279e-02      8346 -3.278e+03 -6.944e+04                   
#> Path [11] :Best Iter: [45] ELBO (-3278.393583) evaluations: (8346) 
#> Path [12] :Initial log joint density = -425544.809591 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.859e-02   8.921e+02    1.789e-02  1.789e-02      8067 -3.267e+03 -4.508e+04                   
#> Path [12] :Best Iter: [55] ELBO (-3267.310464) evaluations: (8067) 
#> Path [13] :Initial log joint density = -425111.987875 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.204e-02   2.910e+03    3.684e-02  3.684e-02      8196 -3.278e+03 -1.841e+04                   
#> Path [13] :Best Iter: [42] ELBO (-3278.132009) evaluations: (8196) 
#> Path [14] :Initial log joint density = -425427.622423 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.252e-01   1.688e+03    4.269e-02  4.269e-02      8096 -3.271e+03 -2.405e+05                   
#> Path [14] :Best Iter: [61] ELBO (-3271.247131) evaluations: (8096) 
#> Path [15] :Initial log joint density = -425217.083745 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.640e-02   1.979e+03    4.289e-02  4.289e-02      8084 -3.280e+03 -1.653e+04                   
#> Path [15] :Best Iter: [45] ELBO (-3280.483533) evaluations: (8084) 
#> Path [16] :Initial log joint density = -425070.207399 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.685e-02   2.046e+03    1.155e-02  2.903e-02      8416 -3.273e+03 -1.749e+10                   
#> Path [16] :Best Iter: [55] ELBO (-3272.765746) evaluations: (8416) 
#> Path [17] :Initial log joint density = -425109.414170 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.280e-02   2.471e+03    1.710e-02  1.710e-02      8484 -3.281e+03 -1.661e+06                   
#> Path [17] :Best Iter: [47] ELBO (-3280.626813) evaluations: (8484) 
#> Path [18] :Initial log joint density = -426712.294113 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.982e-02   9.034e+02    3.236e-02  3.236e-02      8106 -3.269e+03 -6.156e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3268.830040) evaluations: (8106) 
#> Path [19] :Initial log joint density = -426356.465513 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.457e-01   1.440e+03    2.537e-02  5.748e-02      8264 -3.275e+03 -1.010e+08                   
#> Path [19] :Best Iter: [49] ELBO (-3274.616654) evaluations: (8264) 
#> Path [20] :Initial log joint density = -425064.741198 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      4.098e-02   2.107e+03    1.685e-02  1.685e-02      8438 -3.276e+03 -4.380e+05                   
#> Path [20] :Best Iter: [43] ELBO (-3276.493440) evaluations: (8438) 
#> Path [21] :Initial log joint density = -425129.230143 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.743e-02   2.416e+03    2.165e-02  2.165e-02      8456 -3.277e+03 -4.578e+06                   
#> Path [21] :Best Iter: [47] ELBO (-3276.794124) evaluations: (8456) 
#> Path [22] :Initial log joint density = -426779.460887 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.217e-02   1.867e+03    3.248e-02  3.248e-02      8090 -3.276e+03 -7.320e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3276.315152) evaluations: (8090) 
#> Path [23] :Initial log joint density = -428751.524216 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.136e-01   8.141e+02    3.577e-02  3.577e-02      8242 -3.267e+03 -1.022e+05                   
#> Path [23] :Best Iter: [55] ELBO (-3267.463148) evaluations: (8242) 
#> Path [24] :Initial log joint density = -425035.826219 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.127e-02   2.844e+03    2.492e-02  4.616e-02      8305 -3.279e+03 -6.370e+05                   
#> Path [24] :Best Iter: [40] ELBO (-3279.313885) evaluations: (8305) 
#> Path [25] :Initial log joint density = -425049.711298 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      2.019e-02   2.094e+03    2.093e-02  2.093e-02      7951 -3.276e+03 -1.526e+04                   
#> Path [25] :Best Iter: [55] ELBO (-3275.684482) evaluations: (7951) 
#> Path [26] :Initial log joint density = -425187.673867 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.415e-02   2.083e+03    2.862e-02  2.862e-02      7992 -3.277e+03 -9.165e+03                   
#> Path [26] :Best Iter: [48] ELBO (-3276.500864) evaluations: (7992) 
#> Path [27] :Initial log joint density = -425653.274204 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.488e-02   2.978e+03    4.269e-02  4.269e-02      8265 -3.270e+03 -1.536e+04                   
#> Path [27] :Best Iter: [61] ELBO (-3270.448662) evaluations: (8265) 
#> Path [28] :Initial log joint density = -427906.376757 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.200e-01   1.514e+03    2.332e-02  2.332e-02      8044 -3.271e+03 -3.340e+05                   
#> Path [28] :Best Iter: [56] ELBO (-3271.035147) evaluations: (8044) 
#> Path [29] :Initial log joint density = -425197.486246 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.879e-02   2.940e+03    2.613e-02  2.613e-02      8076 -3.278e+03 -3.090e+05                   
#> Path [29] :Best Iter: [45] ELBO (-3278.067981) evaluations: (8076) 
#> Path [30] :Initial log joint density = -425506.357236 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.743e-02   7.099e+02    3.947e-02  3.947e-02      8158 -3.271e+03 -1.876e+04                   
#> Path [30] :Best Iter: [57] ELBO (-3270.657641) evaluations: (8158) 
#> Path [31] :Initial log joint density = -425313.260086 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.065e-02   2.684e+03    4.257e-02  4.257e-02      8315 -3.275e+03 -1.149e+06                   
#> Path [31] :Best Iter: [39] ELBO (-3275.272492) evaluations: (8315) 
#> Path [32] :Initial log joint density = -425397.663863 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.520e-02   4.399e+03    1.569e-02  2.870e-02      8246 -3.272e+03 -4.501e+06                   
#> Path [32] :Best Iter: [55] ELBO (-3272.209232) evaluations: (8246) 
#> Path [33] :Initial log joint density = -428404.377487 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.455e-02   1.647e+03    2.522e-02  2.522e-02      8049 -3.270e+03 -6.139e+04                   
#> Path [33] :Best Iter: [57] ELBO (-3269.702908) evaluations: (8049) 
#> Path [34] :Initial log joint density = -425317.204225 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.122e-01   3.151e+03    2.909e-02  2.909e-02      8340 -3.276e+03 -3.705e+07                   
#> Path [34] :Best Iter: [48] ELBO (-3275.629808) evaluations: (8340) 
#> Path [35] :Initial log joint density = -425487.563038 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.036e-02   1.299e+03    3.962e-02  3.962e-02      7994 -3.273e+03 -9.131e+04                   
#> Path [35] :Best Iter: [59] ELBO (-3272.906230) evaluations: (7994) 
#> Path [36] :Initial log joint density = -425192.252393 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.649e-02   4.114e+03    3.212e-02  6.211e-02      8199 -3.274e+03 -3.910e+06                   
#> Path [36] :Best Iter: [46] ELBO (-3273.635565) evaluations: (8199) 
#> Path [37] :Initial log joint density = -425068.833182 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      5.511e-02   1.704e+03    2.129e-02  4.544e-02      8027 -3.269e+03 -1.593e+05                   
#> Path [37] :Best Iter: [58] ELBO (-3268.702968) evaluations: (8027) 
#> Path [38] :Initial log joint density = -425219.739918 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.168e-01   9.357e+02    5.580e-02  9.294e-02      8183 -3.267e+03 -1.465e+06                   
#> Path [38] :Best Iter: [60] ELBO (-3267.394392) evaluations: (8183) 
#> Path [39] :Initial log joint density = -425158.132710 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      9.410e-02   4.523e+03    1.666e-02  1.666e-02      8233 -3.276e+03 -3.804e+11                   
#> Path [39] :Best Iter: [48] ELBO (-3276.026330) evaluations: (8233) 
#> Path [40] :Initial log joint density = -425220.022429 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.175e-02   2.825e+03    3.069e-02  3.069e-02      8394 -3.278e+03 -6.119e+06                   
#> Path [40] :Best Iter: [44] ELBO (-3277.782265) evaluations: (8394) 
#> Path [41] :Initial log joint density = -425183.518875 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      3.478e-02   1.113e+03    2.116e-02  4.095e-02      8206 -3.271e+03 -9.680e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3270.625599) evaluations: (8206) 
#> Path [42] :Initial log joint density = -425194.696697 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.962e-02   2.197e+03    2.780e-02  2.780e-02      8221 -3.279e+03 -8.988e+04                   
#> Path [42] :Best Iter: [45] ELBO (-3278.848280) evaluations: (8221) 
#> Path [43] :Initial log joint density = -425576.682640 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      1.103e-01   8.533e+02    4.296e-02  4.296e-02      8026 -3.272e+03 -1.254e+04                   
#> Path [43] :Best Iter: [56] ELBO (-3271.575119) evaluations: (8026) 
#> Path [44] :Initial log joint density = -425351.085126 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.621e-02   2.265e+03    2.070e-02  4.704e-02      8184 -3.274e+03 -1.069e+07                   
#> Path [44] :Best Iter: [54] ELBO (-3273.640880) evaluations: (8184) 
#> Path [45] :Initial log joint density = -425130.389110 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      8.193e-02   7.574e+03    1.630e-02  1.630e-02      8492 -3.278e+03 -1.378e+06                   
#> Path [45] :Best Iter: [47] ELBO (-3277.593672) evaluations: (8492) 
#> Path [46] :Initial log joint density = -426990.109503 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      4.323e-02   8.759e+02    2.822e-02  2.822e-02      8192 -3.269e+03 -2.322e+04                   
#> Path [46] :Best Iter: [58] ELBO (-3269.142015) evaluations: (8192) 
#> Path [47] :Initial log joint density = -425083.636177 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      7.155e-02   1.965e+03    2.204e-02  2.204e-02      8219 -3.278e+03 -2.092e+05                   
#> Path [47] :Best Iter: [49] ELBO (-3277.758546) evaluations: (8219) 
#> Path [48] :Initial log joint density = -424986.953034 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.458e-02   2.932e+03    2.889e-02  2.889e-02      8213 -3.277e+03 -3.016e+04                   
#> Path [48] :Best Iter: [55] ELBO (-3277.155453) evaluations: (8213) 
#> Path [49] :Initial log joint density = -425166.660795 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      6.974e-02   3.061e+03    2.461e-02  2.461e-02      8218 -3.279e+03 -2.593e+05                   
#> Path [49] :Best Iter: [42] ELBO (-3278.946856) evaluations: (8218) 
#> Path [50] :Initial log joint density = -425197.603897 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.223e+05      9.307e-02   2.747e+03    1.869e-02  1.869e-02      8126 -3.278e+03 -1.232e+07                   
#> Path [50] :Best Iter: [45] ELBO (-3278.094882) evaluations: (8126) 
#> Finished in  27.5 seconds.
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
