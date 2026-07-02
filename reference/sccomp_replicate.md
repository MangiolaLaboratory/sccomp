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
#> Path [1] :Initial log joint density = -482182.306870 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.118e-02   2.747e-01    1.000e+00  1.000e+00      3277 -3.682e+03 -3.688e+03                   
#> Path [1] :Best Iter: [57] ELBO (-3681.997940) evaluations: (3277) 
#> Path [2] :Initial log joint density = -481682.144145 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      5.123e-03   1.672e-01    1.000e+00  1.000e+00      2828 -3.689e+03 -3.698e+03                   
#> Path [2] :Best Iter: [49] ELBO (-3688.687132) evaluations: (2828) 
#> Path [3] :Initial log joint density = -482061.279901 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.282e-03   2.644e-01    9.880e-01  9.880e-01      3503 -3.681e+03 -3.694e+03                   
#> Path [3] :Best Iter: [59] ELBO (-3681.108570) evaluations: (3503) 
#> Path [4] :Initial log joint density = -481223.991759 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.099e-02   3.281e-01    9.420e-01  9.420e-01      3081 -3.691e+03 -3.693e+03                   
#> Path [4] :Best Iter: [48] ELBO (-3691.107722) evaluations: (3081) 
#> Path [5] :Initial log joint density = -481948.597140 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.071e-02   2.073e-01    1.000e+00  1.000e+00      3206 -3.684e+03 -3.682e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3681.702168) evaluations: (3206) 
#> Path [6] :Initial log joint density = -481704.934335 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.166e-02   2.442e-01    1.000e+00  1.000e+00      3555 -3.683e+03 -3.687e+03                   
#> Path [6] :Best Iter: [57] ELBO (-3682.635915) evaluations: (3555) 
#> Path [7] :Initial log joint density = -481468.232257 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.101e-02   2.583e-01    1.000e+00  1.000e+00      3278 -3.684e+03 -3.689e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3683.805929) evaluations: (3278) 
#> Path [8] :Initial log joint density = -481652.271741 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.260e-03   2.099e-01    7.835e-01  7.835e-01      3109 -3.685e+03 -3.694e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3684.969936) evaluations: (3109) 
#> Path [9] :Initial log joint density = -481438.582667 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.234e-02   3.331e-01    1.000e+00  1.000e+00      3488 -3.680e+03 -3.688e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3680.194101) evaluations: (3488) 
#> Path [10] :Initial log joint density = -481515.543933 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      5.711e-03   2.120e-01    6.464e-01  6.464e-01      3719 -3.683e+03 -3.693e+03                   
#> Path [10] :Best Iter: [56] ELBO (-3682.558013) evaluations: (3719) 
#> Path [11] :Initial log joint density = -481986.215757 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.682e-03   2.487e-01    1.000e+00  1.000e+00      3123 -3.691e+03 -3.694e+03                   
#> Path [11] :Best Iter: [47] ELBO (-3690.789022) evaluations: (3123) 
#> Path [12] :Initial log joint density = -481483.428290 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.774e-03   2.661e-01    9.216e-01  9.216e-01      3149 -3.685e+03 -3.700e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3685.054496) evaluations: (3149) 
#> Path [13] :Initial log joint density = -481880.656771 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.486e-03   1.856e-01    9.881e-01  9.881e-01      3519 -3.684e+03 -3.689e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3683.868374) evaluations: (3519) 
#> Path [14] :Initial log joint density = -481387.392234 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.136e-03   2.325e-01    7.456e-01  7.456e-01      2895 -3.690e+03 -3.708e+03                   
#> Path [14] :Best Iter: [46] ELBO (-3689.994684) evaluations: (2895) 
#> Path [15] :Initial log joint density = -481597.705598 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.823e-02   1.667e-01    1.000e+00  1.000e+00      3772 -3.683e+03 -3.683e+03                   
#> Path [15] :Best Iter: [60] ELBO (-3682.608840) evaluations: (3772) 
#> Path [16] :Initial log joint density = -481773.069430 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.339e-02   3.778e-01    1.000e+00  1.000e+00      2957 -3.692e+03 -3.703e+03                   
#> Path [16] :Best Iter: [52] ELBO (-3692.011897) evaluations: (2957) 
#> Path [17] :Initial log joint density = -481441.886518 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.862e-03   2.044e-01    1.000e+00  1.000e+00      3150 -3.691e+03 -3.683e+03                   
#> Path [17] :Best Iter: [55] ELBO (-3682.958444) evaluations: (3150) 
#> Path [18] :Initial log joint density = -481655.396459 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.698e-03   2.555e-01    9.221e-01  9.221e-01      3566 -3.684e+03 -3.699e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3683.776035) evaluations: (3566) 
#> Path [19] :Initial log joint density = -481568.174986 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.173e-02   2.209e-01    1.000e+00  1.000e+00      3041 -3.690e+03 -3.699e+03                   
#> Path [19] :Best Iter: [48] ELBO (-3690.370360) evaluations: (3041) 
#> Path [20] :Initial log joint density = -481560.280775 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.060e-02   2.135e-01    8.853e-01  8.853e-01      2958 -3.693e+03 -3.705e+03                   
#> Path [20] :Best Iter: [45] ELBO (-3693.368085) evaluations: (2958) 
#> Path [21] :Initial log joint density = -482364.378631 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      9.306e-03   2.684e-01    9.898e-01  9.898e-01      3108 -3.684e+03 -3.694e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3684.333774) evaluations: (3108) 
#> Path [22] :Initial log joint density = -481644.902919 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.051e-02   4.145e-01    1.000e+00  1.000e+00      2869 -3.691e+03 -3.704e+03                   
#> Path [22] :Best Iter: [46] ELBO (-3690.654481) evaluations: (2869) 
#> Path [23] :Initial log joint density = -482173.688641 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.187e-03   2.306e-01    6.914e-01  6.914e-01      3326 -3.683e+03 -3.697e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3683.264611) evaluations: (3326) 
#> Path [24] :Initial log joint density = -481679.493829 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.346e-03   2.004e-01    7.716e-01  7.716e-01      3237 -3.682e+03 -3.694e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3682.488869) evaluations: (3237) 
#> Path [25] :Initial log joint density = -482157.341335 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.240e-02   2.830e-01    1.000e+00  1.000e+00      3108 -3.683e+03 -3.689e+03                   
#> Path [25] :Best Iter: [55] ELBO (-3683.264286) evaluations: (3108) 
#> Path [26] :Initial log joint density = -481391.788903 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.616e-03   1.915e-01    8.293e-01  8.293e-01      2909 -3.688e+03 -3.701e+03                   
#> Path [26] :Best Iter: [46] ELBO (-3688.458627) evaluations: (2909) 
#> Path [27] :Initial log joint density = -481683.229877 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      5.678e-03   1.969e-01    1.000e+00  1.000e+00      2672 -3.688e+03 -3.696e+03                   
#> Path [27] :Best Iter: [42] ELBO (-3687.844647) evaluations: (2672) 
#> Path [28] :Initial log joint density = -482904.587636 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.017e-02   2.636e-01    1.000e+00  1.000e+00      3450 -3.682e+03 -3.686e+03                   
#> Path [28] :Best Iter: [58] ELBO (-3681.733645) evaluations: (3450) 
#> Path [29] :Initial log joint density = -481582.382139 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      2.041e-02   1.305e-01    1.000e+00  1.000e+00      3911 -3.682e+03 -3.682e+03                   
#> Path [29] :Best Iter: [61] ELBO (-3681.852311) evaluations: (3911) 
#> Path [30] :Initial log joint density = -481635.652722 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.103e-03   1.863e-01    1.000e+00  1.000e+00      3356 -3.685e+03 -3.685e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3684.696429) evaluations: (3356) 
#> Path [31] :Initial log joint density = -481886.741682 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.793e-03   1.681e-01    1.000e+00  1.000e+00      3192 -3.689e+03 -3.695e+03                   
#> Path [31] :Best Iter: [55] ELBO (-3689.307082) evaluations: (3192) 
#> Path [32] :Initial log joint density = -481456.061955 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.203e-02   1.620e-01    1.000e+00  1.000e+00      3240 -3.685e+03 -3.687e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3685.206605) evaluations: (3240) 
#> Path [33] :Initial log joint density = -482984.691089 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.503e-03   2.315e-01    1.000e+00  1.000e+00      3118 -3.690e+03 -3.685e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3685.195467) evaluations: (3118) 
#> Path [34] :Initial log joint density = -481375.421133 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.967e-03   2.366e-01    9.343e-01  9.343e-01      2807 -3.692e+03 -3.697e+03                   
#> Path [34] :Best Iter: [45] ELBO (-3692.183080) evaluations: (2807) 
#> Path [35] :Initial log joint density = -481533.192196 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.380e-02   1.573e-01    1.000e+00  1.000e+00      3591 -3.685e+03 -3.687e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3685.111517) evaluations: (3591) 
#> Path [36] :Initial log joint density = -481712.263048 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      4.959e-03   2.059e-01    8.564e-01  8.564e-01      3852 -3.681e+03 -3.695e+03                   
#> Path [36] :Best Iter: [62] ELBO (-3680.800245) evaluations: (3852) 
#> Path [37] :Initial log joint density = -481606.360096 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.262e-03   2.358e-01    1.000e+00  1.000e+00      3561 -3.681e+03 -3.687e+03                   
#> Path [37] :Best Iter: [56] ELBO (-3680.713703) evaluations: (3561) 
#> Path [38] :Initial log joint density = -481379.776898 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.147e-02   3.514e-01    9.506e-01  9.506e-01      2936 -3.691e+03 -3.709e+03                   
#> Path [38] :Best Iter: [52] ELBO (-3691.441761) evaluations: (2936) 
#> Path [39] :Initial log joint density = -481667.578076 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.707e-02   2.458e-01    1.000e+00  1.000e+00      2751 -3.690e+03 -3.696e+03                   
#> Path [39] :Best Iter: [41] ELBO (-3690.094260) evaluations: (2751) 
#> Path [40] :Initial log joint density = -481664.292969 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.077e-02   2.987e-01    8.999e-01  8.999e-01      3119 -3.690e+03 -3.695e+03                   
#> Path [40] :Best Iter: [48] ELBO (-3690.459686) evaluations: (3119) 
#> Path [41] :Initial log joint density = -486208.297554 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.005e-02   3.616e-01    1.000e+00  1.000e+00      3278 -3.686e+03 -3.687e+03                   
#> Path [41] :Best Iter: [57] ELBO (-3686.069972) evaluations: (3278) 
#> Path [42] :Initial log joint density = -481564.390156 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.294e-03   2.250e-01    1.000e+00  1.000e+00      3125 -3.689e+03 -3.683e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3683.153862) evaluations: (3125) 
#> Path [43] :Initial log joint density = -483115.612111 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.999e-03   3.420e-01    3.238e-01  1.000e+00      2990 -3.694e+03 -3.712e+03                   
#> Path [43] :Best Iter: [36] ELBO (-3694.084284) evaluations: (2990) 
#> Path [44] :Initial log joint density = -481372.209754 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.390e-02   2.747e-01    9.545e-01  9.545e-01      3446 -3.683e+03 -3.692e+03                   
#> Path [44] :Best Iter: [58] ELBO (-3683.462845) evaluations: (3446) 
#> Path [45] :Initial log joint density = -481556.815572 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      7.535e-03   2.000e-01    1.000e+00  1.000e+00      3675 -3.683e+03 -3.685e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3683.391548) evaluations: (3675) 
#> Path [46] :Initial log joint density = -481675.688415 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      8.492e-03   2.064e-01    5.180e-01  1.000e+00      2992 -3.689e+03 -3.703e+03                   
#> Path [46] :Best Iter: [47] ELBO (-3688.568813) evaluations: (2992) 
#> Path [47] :Initial log joint density = -481804.539648 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      9.303e-03   2.636e-01    5.225e-01  1.000e+00      3600 -3.680e+03 -3.691e+03                   
#> Path [47] :Best Iter: [59] ELBO (-3679.960256) evaluations: (3600) 
#> Path [48] :Initial log joint density = -481416.659963 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.915e-03   2.255e-01    3.040e-01  1.000e+00      3209 -3.686e+03 -3.695e+03                   
#> Path [48] :Best Iter: [55] ELBO (-3685.944421) evaluations: (3209) 
#> Path [49] :Initial log joint density = -481962.667815 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.109e-03   1.829e-01    1.000e+00  1.000e+00      3503 -3.685e+03 -3.697e+03                   
#> Path [49] :Best Iter: [57] ELBO (-3684.919768) evaluations: (3503) 
#> Path [50] :Initial log joint density = -481707.831191 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.877e-03   1.650e-01    4.696e-01  1.000e+00      3121 -3.689e+03 -3.694e+03                   
#> Path [50] :Best Iter: [43] ELBO (-3688.622091) evaluations: (3121) 
#> Finished in  13.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.001 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 720 × 5
#>    cell_group sample       generated_proportions generated_counts replicate
#>    <chr>      <fct>                        <dbl>            <int>     <int>
#>  1 B1         10x_6K                      0.0767              383         1
#>  2 B1         10x_8K                      0.0708              353         1
#>  3 B1         GSE115189                   0.0710              355         1
#>  4 B1         SCP345_580                  0.0786              393         1
#>  5 B1         SCP345_860                  0.0422              211         1
#>  6 B1         SCP424_pbmc1                0.0847              423         1
#>  7 B1         SCP424_pbmc2                0.0450              225         1
#>  8 B1         SCP591                      0.0373              186         1
#>  9 B1         SI-GA-E5                    0.0178               89         1
#> 10 B1         SI-GA-E7                    0.0445              222         1
#> # ℹ 710 more rows
# }
```
