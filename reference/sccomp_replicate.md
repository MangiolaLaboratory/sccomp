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
#> Path [1] :Initial log joint density = -481857.352512 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.072e-02   2.772e-01    5.264e-01  1.000e+00      3071 -3.708e+03 -3.708e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3707.571304) evaluations: (3071) 
#> Path [2] :Initial log joint density = -481263.536047 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.827e-03   2.019e-01    9.344e-01  9.344e-01      2832 -3.706e+03 -3.721e+03                   
#> Path [2] :Best Iter: [51] ELBO (-3706.222412) evaluations: (2832) 
#> Path [3] :Initial log joint density = -481682.009281 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.229e-02   2.397e-01    1.000e+00  1.000e+00      2995 -3.706e+03 -3.709e+03                   
#> Path [3] :Best Iter: [53] ELBO (-3705.821053) evaluations: (2995) 
#> Path [4] :Initial log joint density = -482140.254034 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      7.059e-03   2.271e-01    1.000e+00  1.000e+00      3674 -3.699e+03 -3.703e+03                   
#> Path [4] :Best Iter: [58] ELBO (-3699.287060) evaluations: (3674) 
#> Path [5] :Initial log joint density = -481199.527570 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.464e-03   2.567e-01    1.000e+00  1.000e+00      3262 -3.706e+03 -3.701e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3701.142608) evaluations: (3262) 
#> Path [6] :Initial log joint density = -481725.049984 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.072e-02   2.604e-01    1.000e+00  1.000e+00      3535 -3.704e+03 -3.702e+03                   
#> Path [6] :Best Iter: [60] ELBO (-3702.251819) evaluations: (3535) 
#> Path [7] :Initial log joint density = -482410.654722 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.476e-02   2.443e-01    1.000e+00  1.000e+00      3028 -3.708e+03 -3.702e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3701.535344) evaluations: (3028) 
#> Path [8] :Initial log joint density = -481492.383749 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.169e-03   1.952e-01    8.410e-01  8.410e-01      2797 -3.705e+03 -3.730e+03                   
#> Path [8] :Best Iter: [49] ELBO (-3705.492961) evaluations: (2797) 
#> Path [9] :Initial log joint density = -481764.651131 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.950e-03   2.214e-01    8.903e-01  8.903e-01      3195 -3.708e+03 -3.709e+03                   
#> Path [9] :Best Iter: [55] ELBO (-3707.559354) evaluations: (3195) 
#> Path [10] :Initial log joint density = -481473.116263 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.949e-03   2.365e-01    9.878e-01  9.878e-01      3106 -3.708e+03 -3.708e+03                   
#> Path [10] :Best Iter: [55] ELBO (-3708.196652) evaluations: (3106) 
#> Path [11] :Initial log joint density = -481891.986901 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.274e-03   3.115e-01    6.483e-01  6.483e-01      3434 -3.700e+03 -3.712e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3700.273471) evaluations: (3434) 
#> Path [12] :Initial log joint density = -485313.233904 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.410e-02   2.247e-01    1.000e+00  1.000e+00      3533 -3.701e+03 -3.703e+03                   
#> Path [12] :Best Iter: [58] ELBO (-3701.117858) evaluations: (3533) 
#> Path [13] :Initial log joint density = -481941.573464 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.342e-03   2.283e-01    1.000e+00  1.000e+00      3620 -3.701e+03 -3.699e+03                   
#> Path [13] :Best Iter: [62] ELBO (-3698.763517) evaluations: (3620) 
#> Path [14] :Initial log joint density = -483622.857925 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.939e-03   2.252e-01    1.000e+00  1.000e+00      2957 -3.707e+03 -3.718e+03                   
#> Path [14] :Best Iter: [48] ELBO (-3706.710122) evaluations: (2957) 
#> Path [15] :Initial log joint density = -481675.727968 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.896e-03   3.260e-01    1.000e+00  1.000e+00      3122 -3.710e+03 -3.706e+03                   
#> Path [15] :Best Iter: [55] ELBO (-3706.202511) evaluations: (3122) 
#> Path [16] :Initial log joint density = -482002.369050 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.629e-03   1.253e-01    8.056e-01  8.056e-01      3591 -3.702e+03 -3.713e+03                   
#> Path [16] :Best Iter: [59] ELBO (-3702.227555) evaluations: (3591) 
#> Path [17] :Initial log joint density = -481920.783513 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.711e-03   2.751e-01    7.667e-01  7.667e-01      3265 -3.702e+03 -3.712e+03                   
#> Path [17] :Best Iter: [55] ELBO (-3702.310919) evaluations: (3265) 
#> Path [18] :Initial log joint density = -482555.134240 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.475e-02   2.071e-01    1.000e+00  1.000e+00      3423 -3.703e+03 -3.702e+03                   
#> Path [18] :Best Iter: [59] ELBO (-3702.209607) evaluations: (3423) 
#> Path [19] :Initial log joint density = -483990.222971 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.523e-02   2.670e-01    1.000e+00  1.000e+00      2986 -3.709e+03 -3.713e+03                   
#> Path [19] :Best Iter: [45] ELBO (-3708.950406) evaluations: (2986) 
#> Path [20] :Initial log joint density = -482014.813191 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.657e-02   3.120e-01    5.113e-01  1.000e+00      4007 -3.697e+03 -3.704e+03                   
#> Path [20] :Best Iter: [63] ELBO (-3697.371682) evaluations: (4007) 
#> Path [21] :Initial log joint density = -481670.996419 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.855e-02   2.302e-01    1.000e+00  1.000e+00      3257 -3.704e+03 -3.701e+03                   
#> Path [21] :Best Iter: [56] ELBO (-3701.343670) evaluations: (3257) 
#> Path [22] :Initial log joint density = -481388.051225 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      2.757e-03   2.124e-01    7.642e-01  7.642e-01      2961 -3.702e+03 -3.722e+03                   
#> Path [22] :Best Iter: [51] ELBO (-3701.973754) evaluations: (2961) 
#> Path [23] :Initial log joint density = -481692.347937 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.635e-02   4.033e-01    1.000e+00  1.000e+00      2854 -3.705e+03 -3.710e+03                   
#> Path [23] :Best Iter: [51] ELBO (-3704.532218) evaluations: (2854) 
#> Path [24] :Initial log joint density = -485362.943854 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.433e-02   3.109e-01    1.000e+00  1.000e+00      3298 -3.700e+03 -3.699e+03                   
#> Path [24] :Best Iter: [57] ELBO (-3699.150335) evaluations: (3298) 
#> Path [25] :Initial log joint density = -482234.963368 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      7.803e-03   2.127e-01    1.000e+00  1.000e+00      3772 -3.700e+03 -3.707e+03                   
#> Path [25] :Best Iter: [61] ELBO (-3699.817680) evaluations: (3772) 
#> Path [26] :Initial log joint density = -481655.402253 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.741e-03   2.851e-01    1.000e+00  1.000e+00      3008 -3.706e+03 -3.715e+03                   
#> Path [26] :Best Iter: [48] ELBO (-3705.940808) evaluations: (3008) 
#> Path [27] :Initial log joint density = -481701.413051 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.501e-03   2.471e-01    7.866e-01  7.866e-01      3160 -3.702e+03 -3.717e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3702.091432) evaluations: (3160) 
#> Path [28] :Initial log joint density = -481836.262825 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.455e-02   3.698e-01    1.000e+00  1.000e+00      3229 -3.704e+03 -3.711e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3704.314872) evaluations: (3229) 
#> Path [29] :Initial log joint density = -481614.553185 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.436e-03   2.493e-01    1.000e+00  1.000e+00      3036 -3.706e+03 -3.712e+03                   
#> Path [29] :Best Iter: [42] ELBO (-3705.607002) evaluations: (3036) 
#> Path [30] :Initial log joint density = -481713.050002 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.191e-02   3.243e-01    1.000e+00  1.000e+00      2991 -3.707e+03 -3.711e+03                   
#> Path [30] :Best Iter: [43] ELBO (-3707.342788) evaluations: (2991) 
#> Path [31] :Initial log joint density = -481273.168557 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.180e-02   3.834e-01    1.000e+00  1.000e+00      2991 -3.708e+03 -3.724e+03                   
#> Path [31] :Best Iter: [46] ELBO (-3708.080784) evaluations: (2991) 
#> Path [32] :Initial log joint density = -482887.927356 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.487e-03   2.747e-01    8.353e-01  8.353e-01      3240 -3.701e+03 -3.712e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3701.287790) evaluations: (3240) 
#> Path [33] :Initial log joint density = -481609.144559 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.730e-03   1.959e-01    1.000e+00  1.000e+00      3109 -3.709e+03 -3.710e+03                   
#> Path [33] :Best Iter: [52] ELBO (-3708.953872) evaluations: (3109) 
#> Path [34] :Initial log joint density = -481418.022135 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.530e-03   3.964e-01    5.809e-01  5.809e-01      3105 -3.705e+03 -3.716e+03                   
#> Path [34] :Best Iter: [44] ELBO (-3705.113596) evaluations: (3105) 
#> Path [35] :Initial log joint density = -482064.760528 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      9.423e-03   2.904e-01    1.000e+00  1.000e+00      3769 -3.701e+03 -3.703e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3700.742872) evaluations: (3769) 
#> Path [36] :Initial log joint density = -481767.303765 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.773e-02   2.734e-01    9.345e-01  9.345e-01      3385 -3.698e+03 -3.706e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3697.731451) evaluations: (3385) 
#> Path [37] :Initial log joint density = -481623.902326 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.502e-03   2.003e-01    1.000e+00  1.000e+00      3175 -3.705e+03 -3.705e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3704.765860) evaluations: (3175) 
#> Path [38] :Initial log joint density = -481752.940232 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.280e-03   3.445e-01    4.339e-01  1.000e+00      2995 -3.708e+03 -3.724e+03                   
#> Path [38] :Best Iter: [43] ELBO (-3707.820251) evaluations: (2995) 
#> Path [39] :Initial log joint density = -484948.503389 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.485e-02   2.450e-01    1.000e+00  1.000e+00      3421 -3.699e+03 -3.700e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3698.527984) evaluations: (3421) 
#> Path [40] :Initial log joint density = -481848.889677 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.201e-02   2.555e-01    1.000e+00  1.000e+00      2833 -3.707e+03 -3.715e+03                   
#> Path [40] :Best Iter: [50] ELBO (-3707.268683) evaluations: (2833) 
#> Path [41] :Initial log joint density = -481662.350201 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.771e-03   2.258e-01    8.547e-01  8.547e-01      3057 -3.706e+03 -3.729e+03                   
#> Path [41] :Best Iter: [46] ELBO (-3705.677210) evaluations: (3057) 
#> Path [42] :Initial log joint density = -481495.192395 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.654e-03   2.112e-01    1.000e+00  1.000e+00      3136 -3.706e+03 -3.704e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3703.954043) evaluations: (3136) 
#> Path [43] :Initial log joint density = -481617.743788 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.126e-02   2.771e-01    1.000e+00  1.000e+00      3345 -3.702e+03 -3.703e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3701.908004) evaluations: (3345) 
#> Path [44] :Initial log joint density = -481347.701139 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.230e-02   3.035e-01    1.000e+00  1.000e+00      2874 -3.707e+03 -3.711e+03                   
#> Path [44] :Best Iter: [41] ELBO (-3706.825189) evaluations: (2874) 
#> Path [45] :Initial log joint density = -481616.840901 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.657e-03   2.598e-01    1.000e+00  1.000e+00      3093 -3.708e+03 -3.711e+03                   
#> Path [45] :Best Iter: [52] ELBO (-3707.982397) evaluations: (3093) 
#> Path [46] :Initial log joint density = -481830.077846 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.780e-03   1.912e-01    1.000e+00  1.000e+00      3158 -3.708e+03 -3.708e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3707.677671) evaluations: (3158) 
#> Path [47] :Initial log joint density = -481725.120970 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.041e-03   1.511e-01    1.000e+00  1.000e+00      2925 -3.707e+03 -3.716e+03                   
#> Path [47] :Best Iter: [50] ELBO (-3706.618925) evaluations: (2925) 
#> Path [48] :Initial log joint density = -481673.716321 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.061e-02   2.269e-01    1.000e+00  1.000e+00      3343 -3.700e+03 -3.699e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3698.712217) evaluations: (3343) 
#> Path [49] :Initial log joint density = -481548.724803 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.468e-03   3.233e-01    6.704e-01  6.704e-01      3121 -3.712e+03 -3.717e+03                   
#> Path [49] :Best Iter: [40] ELBO (-3712.300190) evaluations: (3121) 
#> Path [50] :Initial log joint density = -482713.963241 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.371e-02   2.222e-01    1.000e+00  1.000e+00      3437 -3.700e+03 -3.698e+03                   
#> Path [50] :Best Iter: [58] ELBO (-3697.583389) evaluations: (3437) 
#> Finished in  13.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 5
#>    cell_group sample       generated_proportions generated_counts replicate
#>    <chr>      <fct>                        <dbl>            <int>     <int>
#>  1 B1         10x_6K                      0.0387              193         1
#>  2 B1         10x_8K                      0.0222              110         1
#>  3 B1         GSE115189                   0.0721              360         1
#>  4 B1         SCP345_580                  0.0738              368         1
#>  5 B1         SCP345_860                  0.0384              191         1
#>  6 B1         SCP424_pbmc1                0.0525              262         1
#>  7 B1         SCP424_pbmc2                0.0701              350         1
#>  8 B1         SCP591                      0.0676              337         1
#>  9 B1         SI-GA-E5                    0.0422              210         1
#> 10 B1         SI-GA-E7                    0.0185               92         1
#> # ℹ 710 more rows
# }
```
