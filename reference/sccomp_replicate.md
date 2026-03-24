# sccomp_replicate

This function replicates counts from a real-world dataset.

## Usage

``` r
sccomp_replicate(
  fit,
  formula_composition = NULL,
  formula_variability = NULL,
  number_of_draws = 1,
  mcmc_seed = sample_seed(),
  cache_stan_model = sccomp_stan_models_cache_dir
)
```

## Arguments

- fit:

  The result of sccomp_estimate.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment. This formula can be a sub-formula
  of your estimated model; in this case all other factor will be
  factored out.

- formula_variability:

  A formula. The formula describing the model for differential
  variability, for example ~treatment. In most cases, if differentially
  variability is of interest, the formula should only include the factor
  of interest as a large anount of data is needed to define variability
  depending to each factors. This formula can be a sub-formula of your
  estimated model; in this case all other factor will be factored out.

- number_of_draws:

  An integer. How may copies of the data you want to draw from the model
  joint posterior distribution.

- mcmc_seed:

  An integer. Used for Markov-chain Monte Carlo reproducibility. By
  default a random number is sampled from 1 to 999999. This itself can
  be controlled by set.seed()

- cache_stan_model:

  A character string specifying the cache directory for compiled Stan
  models. The sccomp version will be automatically appended to ensure
  version isolation. Default is `sccomp_stan_models_cache_dir` which
  points to `~/.sccomp_models`.

## Value

A tibble `tbl` with cell_group-wise statistics

A tibble (`tbl`), with the following columns:

- **cell_group** - A character column representing the cell group being
  tested.

- **sample** - A factor column representing the sample name from which
  data was generated.

- **generated_proportions** - A numeric column representing the
  proportions generated from the model.

- **generated_counts** - An integer column representing the counts
  generated from the model.

- **replicate** - An integer column representing the replicate number,
  where each row corresponds to a different replicate of the data.

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
  if (instantiate::stan_cmdstan_exists() && .Platform$OS.type == "unix") {
    data("counts_obj")

    sccomp_estimate(
      counts_obj,
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    ) |>
    sccomp_replicate()
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481599.229303 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      8.154e-03   2.861e-01    1.000e+00  1.000e+00      2675 -3.710e+03 -3.718e+03                   
#> Path [1] :Best Iter: [43] ELBO (-3710.113106) evaluations: (2675) 
#> Path [2] :Initial log joint density = -481832.110628 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.202e-02   3.838e-01    1.000e+00  1.000e+00      3220 -3.701e+03 -3.708e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3700.845165) evaluations: (3220) 
#> Path [3] :Initial log joint density = -482290.569816 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.061e-02   2.746e-01    1.000e+00  1.000e+00      3405 -3.699e+03 -3.704e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3699.499487) evaluations: (3405) 
#> Path [4] :Initial log joint density = -481529.126314 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.060e-03   2.883e-01    7.007e-01  7.007e-01      3149 -3.707e+03 -3.714e+03                   
#> Path [4] :Best Iter: [52] ELBO (-3706.934954) evaluations: (3149) 
#> Path [5] :Initial log joint density = -481483.613787 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.103e-02   3.625e-01    1.000e+00  1.000e+00      3270 -3.701e+03 -3.707e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3701.285203) evaluations: (3270) 
#> Path [6] :Initial log joint density = -481636.772429 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.381e-03   1.961e-01    1.000e+00  1.000e+00      2958 -3.705e+03 -3.716e+03                   
#> Path [6] :Best Iter: [42] ELBO (-3705.250403) evaluations: (2958) 
#> Path [7] :Initial log joint density = -481772.750542 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.247e-02   3.056e-01    1.000e+00  1.000e+00      3120 -3.706e+03 -3.705e+03                   
#> Path [7] :Best Iter: [56] ELBO (-3705.346805) evaluations: (3120) 
#> Path [8] :Initial log joint density = -481693.386754 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.946e-02   2.442e-01    1.000e+00  1.000e+00      3789 -3.698e+03 -3.704e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3697.905196) evaluations: (3789) 
#> Path [9] :Initial log joint density = -481212.652211 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      8.694e-03   2.678e-01    1.000e+00  1.000e+00      2711 -3.706e+03 -3.713e+03                   
#> Path [9] :Best Iter: [43] ELBO (-3706.301998) evaluations: (2711) 
#> Path [10] :Initial log joint density = -481218.168884 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              48      -4.788e+05      9.323e-03   1.783e-01    1.000e+00  1.000e+00      2473 -3.704e+03 -3.707e+03                   
#> Path [10] :Best Iter: [45] ELBO (-3704.038969) evaluations: (2473) 
#> Path [11] :Initial log joint density = -481631.859786 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.011e-03   1.458e-01    1.000e+00  1.000e+00      2827 -3.707e+03 -3.721e+03                   
#> Path [11] :Best Iter: [48] ELBO (-3707.063259) evaluations: (2827) 
#> Path [12] :Initial log joint density = -481318.467001 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.654e-03   2.275e-01    1.000e+00  1.000e+00      3391 -3.700e+03 -3.704e+03                   
#> Path [12] :Best Iter: [58] ELBO (-3699.714811) evaluations: (3391) 
#> Path [13] :Initial log joint density = -481715.856332 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.527e-03   1.658e-01    1.000e+00  1.000e+00      3391 -3.702e+03 -3.710e+03                   
#> Path [13] :Best Iter: [57] ELBO (-3701.608754) evaluations: (3391) 
#> Path [14] :Initial log joint density = -481607.523310 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.318e-02   2.877e-01    9.001e-01  9.001e-01      3180 -3.709e+03 -3.709e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3709.264348) evaluations: (3180) 
#> Path [15] :Initial log joint density = -481800.049111 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.317e-03   2.273e-01    9.832e-01  9.832e-01      3612 -3.696e+03 -3.711e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3696.330941) evaluations: (3612) 
#> Path [16] :Initial log joint density = -481425.411338 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.857e-03   2.038e-01    1.000e+00  1.000e+00      3628 -3.697e+03 -3.697e+03                   
#> Path [16] :Best Iter: [62] ELBO (-3696.747201) evaluations: (3628) 
#> Path [17] :Initial log joint density = -482350.006490 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.554e-03   2.839e-01    1.000e+00  1.000e+00      3157 -3.710e+03 -3.713e+03                   
#> Path [17] :Best Iter: [46] ELBO (-3709.933401) evaluations: (3157) 
#> Path [18] :Initial log joint density = -481415.871719 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.102e-02   2.649e-01    8.005e-01  8.005e-01      3136 -3.701e+03 -3.711e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3701.388414) evaluations: (3136) 
#> Path [19] :Initial log joint density = -481468.172253 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.005e-02   2.355e-01    1.000e+00  1.000e+00      3599 -3.698e+03 -3.697e+03                   
#> Path [19] :Best Iter: [60] ELBO (-3696.896064) evaluations: (3599) 
#> Path [20] :Initial log joint density = -481361.930747 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.805e-03   3.146e-01    5.754e-01  5.754e-01      2944 -3.708e+03 -3.718e+03                   
#> Path [20] :Best Iter: [46] ELBO (-3708.193119) evaluations: (2944) 
#> Path [21] :Initial log joint density = -481588.138827 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.385e-02   2.790e-01    1.000e+00  1.000e+00      3119 -3.705e+03 -3.701e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3700.733216) evaluations: (3119) 
#> Path [22] :Initial log joint density = -482619.202039 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.621e-03   2.377e-01    1.000e+00  1.000e+00      3272 -3.701e+03 -3.704e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3701.274945) evaluations: (3272) 
#> Path [23] :Initial log joint density = -481665.393475 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.981e-03   1.533e-01    1.000e+00  1.000e+00      3386 -3.700e+03 -3.713e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3700.202987) evaluations: (3386) 
#> Path [24] :Initial log joint density = -481628.413750 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.718e-03   2.567e-01    1.000e+00  1.000e+00      3026 -3.711e+03 -3.705e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3705.122699) evaluations: (3026) 
#> Path [25] :Initial log joint density = -481655.985432 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      7.407e-03   1.971e-01    5.034e-01  1.000e+00      3444 -3.700e+03 -3.712e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3699.613329) evaluations: (3444) 
#> Path [26] :Initial log joint density = -481548.348375 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.500e-02   2.558e-01    1.000e+00  1.000e+00      3364 -3.700e+03 -3.697e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3697.448489) evaluations: (3364) 
#> Path [27] :Initial log joint density = -481475.660344 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.888e-03   3.249e-01    3.915e-01  1.000e+00      2919 -3.702e+03 -3.720e+03                   
#> Path [27] :Best Iter: [50] ELBO (-3702.040267) evaluations: (2919) 
#> Path [28] :Initial log joint density = -481659.592578 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.197e-02   2.207e-01    1.000e+00  1.000e+00      3479 -3.701e+03 -3.706e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3701.366352) evaluations: (3479) 
#> Path [29] :Initial log joint density = -482511.811599 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.603e-03   1.827e-01    1.000e+00  1.000e+00      3313 -3.700e+03 -3.709e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3700.052106) evaluations: (3313) 
#> Path [30] :Initial log joint density = -481660.932527 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.863e-03   2.023e-01    8.898e-01  8.898e-01      2938 -3.707e+03 -3.719e+03                   
#> Path [30] :Best Iter: [52] ELBO (-3707.430574) evaluations: (2938) 
#> Path [31] :Initial log joint density = -482102.274212 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.090e-02   2.969e-01    1.000e+00  1.000e+00      3270 -3.699e+03 -3.708e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3699.373754) evaluations: (3270) 
#> Path [32] :Initial log joint density = -481518.601412 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.252e-03   1.913e-01    1.000e+00  1.000e+00      2969 -3.707e+03 -3.706e+03                   
#> Path [32] :Best Iter: [53] ELBO (-3706.059933) evaluations: (2969) 
#> Path [33] :Initial log joint density = -482154.097310 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.941e-03   1.918e-01    8.922e-01  8.922e-01      3424 -3.701e+03 -3.708e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3700.614346) evaluations: (3424) 
#> Path [34] :Initial log joint density = -482005.164631 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      1.048e-02   2.156e-01    1.000e+00  1.000e+00      3874 -3.700e+03 -3.706e+03                   
#> Path [34] :Best Iter: [63] ELBO (-3699.971078) evaluations: (3874) 
#> Path [35] :Initial log joint density = -483022.346752 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.074e-02   3.154e-01    9.152e-01  9.152e-01      3239 -3.700e+03 -3.713e+03                   
#> Path [35] :Best Iter: [56] ELBO (-3700.002074) evaluations: (3239) 
#> Path [36] :Initial log joint density = -483661.900050 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.094e-02   3.451e-01    1.000e+00  1.000e+00      3595 -3.702e+03 -3.706e+03                   
#> Path [36] :Best Iter: [58] ELBO (-3702.320959) evaluations: (3595) 
#> Path [37] :Initial log joint density = -481672.601664 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.822e-03   2.311e-01    8.379e-01  8.379e-01      2944 -3.706e+03 -3.721e+03                   
#> Path [37] :Best Iter: [42] ELBO (-3706.013462) evaluations: (2944) 
#> Path [38] :Initial log joint density = -481477.589195 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.266e-03   1.435e-01    8.297e-01  8.297e-01      3128 -3.707e+03 -3.715e+03                   
#> Path [38] :Best Iter: [47] ELBO (-3706.712319) evaluations: (3128) 
#> Path [39] :Initial log joint density = -481645.636482 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.309e-02   4.106e-01    1.000e+00  1.000e+00      2919 -3.707e+03 -3.718e+03                   
#> Path [39] :Best Iter: [52] ELBO (-3706.814575) evaluations: (2919) 
#> Path [40] :Initial log joint density = -481785.864334 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.643e-02   3.545e-01    1.000e+00  1.000e+00      3262 -3.700e+03 -3.702e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3700.044105) evaluations: (3262) 
#> Path [41] :Initial log joint density = -481667.000338 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      8.294e-03   1.985e-01    1.000e+00  1.000e+00      3745 -3.700e+03 -3.707e+03                   
#> Path [41] :Best Iter: [60] ELBO (-3700.390520) evaluations: (3745) 
#> Path [42] :Initial log joint density = -481778.244727 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.726e-02   3.135e-01    1.000e+00  1.000e+00      3451 -3.696e+03 -3.703e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3696.204380) evaluations: (3451) 
#> Path [43] :Initial log joint density = -482315.315808 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.169e-02   2.926e-01    1.000e+00  1.000e+00      3561 -3.700e+03 -3.708e+03                   
#> Path [43] :Best Iter: [58] ELBO (-3700.281210) evaluations: (3561) 
#> Path [44] :Initial log joint density = -481553.635347 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.058e-03   1.676e-01    1.000e+00  1.000e+00      2827 -3.708e+03 -3.717e+03                   
#> Path [44] :Best Iter: [49] ELBO (-3708.149894) evaluations: (2827) 
#> Path [45] :Initial log joint density = -482856.720464 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.323e-02   2.887e-01    1.000e+00  1.000e+00      3363 -3.700e+03 -3.703e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3699.812348) evaluations: (3363) 
#> Path [46] :Initial log joint density = -481597.634721 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.089e-02   3.193e-01    1.000e+00  1.000e+00      3556 -3.699e+03 -3.700e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3698.959525) evaluations: (3556) 
#> Path [47] :Initial log joint density = -481772.098000 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.273e-03   2.787e-01    5.530e-01  5.530e-01      3642 -3.706e+03 -3.717e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3706.266594) evaluations: (3642) 
#> Path [48] :Initial log joint density = -481654.920218 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.612e-03   2.389e-01    1.000e+00  1.000e+00      3305 -3.703e+03 -3.713e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3703.053320) evaluations: (3305) 
#> Path [49] :Initial log joint density = -482018.401436 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.961e-03   1.912e-01    1.000e+00  1.000e+00      3411 -3.700e+03 -3.704e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3699.692226) evaluations: (3411) 
#> Path [50] :Initial log joint density = -481609.690605 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.595e-02   3.050e-01    1.000e+00  1.000e+00      3347 -3.697e+03 -3.699e+03                   
#> Path [50] :Best Iter: [58] ELBO (-3696.676846) evaluations: (3347) 
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
#> # A tibble: 720 × 5
#>    cell_group sample       generated_proportions generated_counts replicate
#>    <chr>      <fct>                        <dbl>            <int>     <int>
#>  1 B1         10x_6K                      0.0712              356         1
#>  2 B1         10x_8K                      0.0625              312         1
#>  3 B1         GSE115189                   0.0366              183         1
#>  4 B1         SCP345_580                  0.0621              310         1
#>  5 B1         SCP345_860                  0.0538              269         1
#>  6 B1         SCP424_pbmc1                0.0739              369         1
#>  7 B1         SCP424_pbmc2                0.0580              290         1
#>  8 B1         SCP591                      0.0657              328         1
#>  9 B1         SI-GA-E5                    0.0244              122         1
#> 10 B1         SI-GA-E7                    0.0486              242         1
#> # ℹ 710 more rows
# }
```
