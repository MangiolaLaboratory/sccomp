# simulate_data

This function simulates data from a fitted model.

## Usage

``` r
simulate_data(
  .data,
  .estimate_object,
  formula_composition,
  formula_variability = NULL,
  .sample = NULL,
  .cell_group = NULL,
  .coefficients = NULL,
  variability_multiplier = 5,
  number_of_draws = 1,
  mcmc_seed = sample_seed(),
  cores = detectCores(),
  sig_figs = 9,
  cache_stan_model = sccomp_stan_models_cache_dir
)
```

## Arguments

- .data:

  A tibble including a cell_group name column \| sample name column \|
  read counts column \| factor columns \| Pvalue column \| a
  significance column

- .estimate_object:

  The result of sccomp_estimate execution. This is used for sampling
  from real-data properties.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment

- formula_variability:

  A formula. The formula describing the model for differential
  variability, for example ~treatment

- .sample:

  A column name as symbol. The sample identifier

- .cell_group:

  A column name as symbol. The cell_group identifier

- .coefficients:

  The column names for coefficients, for example, c(b_0, b_1)

- variability_multiplier:

  A real scalar. This can be used for artificially increasing the
  variability of the simulation for benchmarking purposes.

- number_of_draws:

  An integer. How may copies of the data you want to draw from the model
  joint posterior distribution.

- mcmc_seed:

  An integer. Used for Markov-chain Monte Carlo reproducibility. By
  default a random number is sampled from 1 to 999999. This itself can
  be controlled by set.seed()#' @param cores Integer, the number of
  cores to be used for parallel calculations.

- cores:

  Integer, the number of cores to be used for parallel calculations.

- sig_figs:

  Number of significant figures to use for Stan model output. Default is
  9.

- cache_stan_model:

  A character string specifying the cache directory for compiled Stan
  models. The sccomp version will be automatically appended to ensure
  version isolation. Default is `sccomp_stan_models_cache_dir` which
  points to `~/.sccomp_models`.

## Value

A tibble (`tbl`) with the following columns:

- **sample** - A character column representing the sample name.

- **type** - A factor column representing the type of the sample.

- **phenotype** - A factor column representing the phenotype in the
  data.

- **count** - An integer column representing the original cell counts.

- **cell_group** - A character column representing the cell group
  identifier.

- **b_0** - A numeric column representing the first coefficient used for
  simulation.

- **b_1** - A numeric column representing the second coefficient used
  for simulation.

- **generated_proportions** - A numeric column representing the
  generated proportions from the simulation.

- **generated_counts** - An integer column representing the generated
  cell counts from the simulation.

- **replicate** - An integer column representing the replicate number
  for each draw from the posterior distribution.

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
    library(dplyr)

    estimate = sccomp_estimate(
      counts_obj,
      ~ type, ~1, "sample", "cell_group", "count",
      cores = 1
    )

    # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
    counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)

    # Simulate data
    simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
  }
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481605.483530 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      3.923e-03   1.249e-01    7.334e-01  7.334e-01      3538 -3.701e+03 -3.716e+03                   
#> Path [1] :Best Iter: [58] ELBO (-3701.022059) evaluations: (3538) 
#> Path [2] :Initial log joint density = -482018.910838 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.426e-03   2.835e-01    6.691e-01  6.691e-01      3599 -3.698e+03 -3.713e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3697.863230) evaluations: (3599) 
#> Path [3] :Initial log joint density = -481672.736906 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.423e-02   2.700e-01    7.221e-01  7.221e-01      2673 -3.710e+03 -3.726e+03                   
#> Path [3] :Best Iter: [39] ELBO (-3709.712537) evaluations: (2673) 
#> Path [4] :Initial log joint density = -482518.275451 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      6.811e-03   2.325e-01    1.000e+00  1.000e+00      3330 -3.703e+03 -3.706e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3702.626938) evaluations: (3330) 
#> Path [5] :Initial log joint density = -481528.167528 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.815e-03   2.282e-01    9.315e-01  9.315e-01      2854 -3.708e+03 -3.723e+03                   
#> Path [5] :Best Iter: [46] ELBO (-3707.660487) evaluations: (2854) 
#> Path [6] :Initial log joint density = -481181.060738 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.745e-03   2.696e-01    1.000e+00  1.000e+00      3074 -3.708e+03 -3.704e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3704.424896) evaluations: (3074) 
#> Path [7] :Initial log joint density = -481499.586766 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.688e-02   2.476e-01    1.000e+00  1.000e+00      3688 -3.700e+03 -3.707e+03                   
#> Path [7] :Best Iter: [59] ELBO (-3699.633278) evaluations: (3688) 
#> Path [8] :Initial log joint density = -481857.889907 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      3.237e-03   2.516e-01    7.544e-01  7.544e-01      3566 -3.697e+03 -3.710e+03                   
#> Path [8] :Best Iter: [58] ELBO (-3697.498399) evaluations: (3566) 
#> Path [9] :Initial log joint density = -482390.148268 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.691e-03   1.589e-01    7.658e-01  7.658e-01      3465 -3.701e+03 -3.714e+03                   
#> Path [9] :Best Iter: [56] ELBO (-3700.722291) evaluations: (3465) 
#> Path [10] :Initial log joint density = -481441.830994 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.616e-03   2.445e-01    1.000e+00  1.000e+00      2890 -3.709e+03 -3.716e+03                   
#> Path [10] :Best Iter: [45] ELBO (-3708.742167) evaluations: (2890) 
#> Path [11] :Initial log joint density = -481413.527660 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.643e-03   2.284e-01    1.000e+00  1.000e+00      3382 -3.699e+03 -3.704e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3699.176865) evaluations: (3382) 
#> Path [12] :Initial log joint density = -481763.841204 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.999e-03   2.207e-01    8.240e-01  8.240e-01      3155 -3.702e+03 -3.708e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3701.880708) evaluations: (3155) 
#> Path [13] :Initial log joint density = -481613.037874 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.305e-02   3.228e-01    9.113e-01  9.113e-01      3060 -3.707e+03 -3.710e+03                   
#> Path [13] :Best Iter: [49] ELBO (-3707.236007) evaluations: (3060) 
#> Path [14] :Initial log joint density = -481968.766991 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.086e-03   3.139e-01    6.990e-01  6.990e-01      3139 -3.707e+03 -3.716e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3707.211051) evaluations: (3139) 
#> Path [15] :Initial log joint density = -481538.015455 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.587e-03   2.085e-01    8.136e-01  8.136e-01      2994 -3.707e+03 -3.719e+03                   
#> Path [15] :Best Iter: [52] ELBO (-3707.467470) evaluations: (2994) 
#> Path [16] :Initial log joint density = -481768.007787 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.865e-03   2.250e-01    1.000e+00  1.000e+00      3393 -3.705e+03 -3.698e+03                   
#> Path [16] :Best Iter: [57] ELBO (-3698.113305) evaluations: (3393) 
#> Path [17] :Initial log joint density = -482426.919861 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.787e-03   2.153e-01    9.040e-01  9.040e-01      3209 -3.706e+03 -3.711e+03                   
#> Path [17] :Best Iter: [44] ELBO (-3705.951979) evaluations: (3209) 
#> Path [18] :Initial log joint density = -481826.410735 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.267e-03   2.711e-01    1.000e+00  1.000e+00      2788 -3.710e+03 -3.708e+03                   
#> Path [18] :Best Iter: [51] ELBO (-3708.069959) evaluations: (2788) 
#> Path [19] :Initial log joint density = -481852.973594 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.826e-03   2.078e-01    1.000e+00  1.000e+00      3037 -3.711e+03 -3.726e+03                   
#> Path [19] :Best Iter: [46] ELBO (-3711.365952) evaluations: (3037) 
#> Path [20] :Initial log joint density = -481740.712137 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.091e-03   2.001e-01    1.000e+00  1.000e+00      3380 -3.703e+03 -3.699e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3699.156727) evaluations: (3380) 
#> Path [21] :Initial log joint density = -481612.642379 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.402e-02   4.342e-01    1.000e+00  1.000e+00      3404 -3.702e+03 -3.708e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3702.369602) evaluations: (3404) 
#> Path [22] :Initial log joint density = -481626.376609 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.766e-03   2.098e-01    8.585e-01  8.585e-01      3109 -3.703e+03 -3.708e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3702.676476) evaluations: (3109) 
#> Path [23] :Initial log joint density = -481405.486642 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.835e-02   2.141e-01    1.000e+00  1.000e+00      3006 -3.705e+03 -3.707e+03                   
#> Path [23] :Best Iter: [45] ELBO (-3704.626733) evaluations: (3006) 
#> Path [24] :Initial log joint density = -481821.789938 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.426e-03   1.891e-01    1.000e+00  1.000e+00      3327 -3.700e+03 -3.703e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3700.324346) evaluations: (3327) 
#> Path [25] :Initial log joint density = -481873.779955 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.122e-03   2.534e-01    6.647e-01  6.647e-01      3409 -3.700e+03 -3.716e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3699.865192) evaluations: (3409) 
#> Path [26] :Initial log joint density = -481837.967478 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.210e-02   2.427e-01    8.880e-01  8.880e-01      3425 -3.700e+03 -3.709e+03                   
#> Path [26] :Best Iter: [57] ELBO (-3700.294008) evaluations: (3425) 
#> Path [27] :Initial log joint density = -485624.133260 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.598e-03   1.561e-01    8.164e-01  8.164e-01      3599 -3.702e+03 -3.712e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3702.107641) evaluations: (3599) 
#> Path [28] :Initial log joint density = -481702.441039 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.072e-02   3.130e-01    1.000e+00  1.000e+00      3217 -3.702e+03 -3.703e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3702.146161) evaluations: (3217) 
#> Path [29] :Initial log joint density = -481413.057650 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.531e-03   2.660e-01    7.143e-01  7.143e-01      3338 -3.700e+03 -3.710e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3700.120611) evaluations: (3338) 
#> Path [30] :Initial log joint density = -481472.136102 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.841e-03   2.537e-01    7.768e-01  7.768e-01      2751 -3.708e+03 -3.725e+03                   
#> Path [30] :Best Iter: [45] ELBO (-3707.755177) evaluations: (2751) 
#> Path [31] :Initial log joint density = -481828.416258 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.013e-02   2.908e-01    1.000e+00  1.000e+00      3109 -3.708e+03 -3.703e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3703.314159) evaluations: (3109) 
#> Path [32] :Initial log joint density = -481603.067698 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.399e-03   2.548e-01    1.000e+00  1.000e+00      3075 -3.708e+03 -3.703e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3702.576548) evaluations: (3075) 
#> Path [33] :Initial log joint density = -481839.725261 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.222e-03   1.686e-01    1.000e+00  1.000e+00      3359 -3.700e+03 -3.710e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3699.997407) evaluations: (3359) 
#> Path [34] :Initial log joint density = -481749.598996 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.340e-03   2.648e-01    1.000e+00  1.000e+00      3008 -3.707e+03 -3.715e+03                   
#> Path [34] :Best Iter: [42] ELBO (-3707.123557) evaluations: (3008) 
#> Path [35] :Initial log joint density = -481576.777300 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      2.482e-03   2.861e-01    6.111e-01  6.111e-01      3222 -3.706e+03 -3.712e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3706.254917) evaluations: (3222) 
#> Path [36] :Initial log joint density = -481547.691475 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.301e-03   2.749e-01    5.632e-01  5.632e-01      3003 -3.708e+03 -3.718e+03                   
#> Path [36] :Best Iter: [46] ELBO (-3708.295589) evaluations: (3003) 
#> Path [37] :Initial log joint density = -481601.740641 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.334e-02   4.070e-01    1.000e+00  1.000e+00      3038 -3.706e+03 -3.716e+03                   
#> Path [37] :Best Iter: [44] ELBO (-3705.956297) evaluations: (3038) 
#> Path [38] :Initial log joint density = -482807.688497 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.903e-03   2.409e-01    6.803e-01  6.803e-01      3289 -3.700e+03 -3.708e+03                   
#> Path [38] :Best Iter: [56] ELBO (-3700.465755) evaluations: (3289) 
#> Path [39] :Initial log joint density = -481541.149401 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.156e-02   2.590e-01    1.000e+00  1.000e+00      3185 -3.697e+03 -3.710e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3697.120529) evaluations: (3185) 
#> Path [40] :Initial log joint density = -481983.661422 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.504e-02   3.491e-01    1.000e+00  1.000e+00      3710 -3.698e+03 -3.701e+03                   
#> Path [40] :Best Iter: [61] ELBO (-3697.592714) evaluations: (3710) 
#> Path [41] :Initial log joint density = -482532.043722 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.230e-03   1.604e-01    8.375e-01  8.375e-01      3651 -3.700e+03 -3.713e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3700.036674) evaluations: (3651) 
#> Path [42] :Initial log joint density = -482036.619021 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.448e-02   2.507e-01    1.000e+00  1.000e+00      3478 -3.703e+03 -3.704e+03                   
#> Path [42] :Best Iter: [59] ELBO (-3702.677467) evaluations: (3478) 
#> Path [43] :Initial log joint density = -481581.289773 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.046e-03   2.179e-01    1.000e+00  1.000e+00      3074 -3.711e+03 -3.702e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3702.420780) evaluations: (3074) 
#> Path [44] :Initial log joint density = -483022.967263 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.748e-03   1.635e-01    1.000e+00  1.000e+00      3326 -3.701e+03 -3.707e+03                   
#> Path [44] :Best Iter: [55] ELBO (-3700.829425) evaluations: (3326) 
#> Path [45] :Initial log joint density = -481583.757311 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.243e-03   2.022e-01    1.000e+00  1.000e+00      3213 -3.707e+03 -3.704e+03                   
#> Path [45] :Best Iter: [56] ELBO (-3704.198172) evaluations: (3213) 
#> Path [46] :Initial log joint density = -481490.575329 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.994e-03   1.941e-01    1.000e+00  1.000e+00      2971 -3.707e+03 -3.702e+03                   
#> Path [46] :Best Iter: [55] ELBO (-3701.709424) evaluations: (2971) 
#> Path [47] :Initial log joint density = -481424.929248 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.703e-03   2.153e-01    9.819e-01  9.819e-01      2914 -3.712e+03 -3.726e+03                   
#> Path [47] :Best Iter: [40] ELBO (-3711.833841) evaluations: (2914) 
#> Path [48] :Initial log joint density = -481902.945041 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.267e-03   2.338e-01    7.562e-01  7.562e-01      2995 -3.705e+03 -3.719e+03                   
#> Path [48] :Best Iter: [43] ELBO (-3705.367488) evaluations: (2995) 
#> Path [49] :Initial log joint density = -481635.334636 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.550e-03   2.100e-01    8.291e-01  8.291e-01      3124 -3.706e+03 -3.716e+03                   
#> Path [49] :Best Iter: [52] ELBO (-3706.446333) evaluations: (3124) 
#> Path [50] :Initial log joint density = -481683.614549 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.914e-03   3.256e-01    6.147e-01  6.147e-01      2995 -3.707e+03 -3.723e+03                   
#> Path [50] :Best Iter: [51] ELBO (-3706.592818) evaluations: (2995) 
#> Finished in  13.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> Error in column_to_rownames(., quo_name(.cell_type)): Can't find column `cell_group` in `.data`.
# }
```
