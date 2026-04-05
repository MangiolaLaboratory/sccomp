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
#> Path [1] :Initial log joint density = -481597.254942 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      1.152e-02   2.101e-01    1.000e+00  1.000e+00      2626 -3.692e+03 -3.698e+03                   
#> Path [1] :Best Iter: [35] ELBO (-3691.998724) evaluations: (2626) 
#> Path [2] :Initial log joint density = -481830.911459 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.558e-02   3.196e-01    1.000e+00  1.000e+00      3220 -3.684e+03 -3.688e+03                   
#> Path [2] :Best Iter: [56] ELBO (-3683.731628) evaluations: (3220) 
#> Path [3] :Initial log joint density = -482289.220340 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.008e-02   2.094e-01    1.000e+00  1.000e+00      3569 -3.680e+03 -3.685e+03                   
#> Path [3] :Best Iter: [55] ELBO (-3680.268798) evaluations: (3569) 
#> Path [4] :Initial log joint density = -481527.949786 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.777e-02   2.037e-01    1.000e+00  1.000e+00      2928 -3.691e+03 -3.687e+03                   
#> Path [4] :Best Iter: [52] ELBO (-3686.502369) evaluations: (2928) 
#> Path [5] :Initial log joint density = -481482.913442 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.132e-02   2.560e-01    1.000e+00  1.000e+00      3270 -3.684e+03 -3.686e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3684.330827) evaluations: (3270) 
#> Path [6] :Initial log joint density = -481634.362330 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.286e-03   2.329e-01    1.000e+00  1.000e+00      3125 -3.687e+03 -3.686e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3686.361010) evaluations: (3125) 
#> Path [7] :Initial log joint density = -481770.749116 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.586e-03   1.558e-01    1.000e+00  1.000e+00      3109 -3.691e+03 -3.695e+03                   
#> Path [7] :Best Iter: [43] ELBO (-3691.084659) evaluations: (3109) 
#> Path [8] :Initial log joint density = -481690.819879 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      9.349e-03   1.827e-01    1.000e+00  1.000e+00      3789 -3.682e+03 -3.686e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3681.921153) evaluations: (3789) 
#> Path [9] :Initial log joint density = -481210.930886 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.544e-02   2.923e-01    1.000e+00  1.000e+00      2950 -3.692e+03 -3.694e+03                   
#> Path [9] :Best Iter: [44] ELBO (-3691.613591) evaluations: (2950) 
#> Path [10] :Initial log joint density = -481217.866231 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              49      -4.788e+05      1.120e-02   2.513e-01    1.000e+00  1.000e+00      2593 -3.687e+03 -3.690e+03                   
#> Path [10] :Best Iter: [45] ELBO (-3686.590734) evaluations: (2593) 
#> Path [11] :Initial log joint density = -481630.639221 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.394e-03   2.665e-01    1.000e+00  1.000e+00      2936 -3.689e+03 -3.694e+03                   
#> Path [11] :Best Iter: [48] ELBO (-3689.296708) evaluations: (2936) 
#> Path [12] :Initial log joint density = -481317.580695 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.107e-02   2.973e-01    1.000e+00  1.000e+00      3391 -3.684e+03 -3.686e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3683.633153) evaluations: (3391) 
#> Path [13] :Initial log joint density = -481715.068358 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      7.605e-03   2.075e-01    1.000e+00  1.000e+00      3394 -3.686e+03 -3.684e+03                   
#> Path [13] :Best Iter: [59] ELBO (-3683.852656) evaluations: (3394) 
#> Path [14] :Initial log joint density = -481606.182992 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.712e-03   1.828e-01    1.000e+00  1.000e+00      2880 -3.690e+03 -3.705e+03                   
#> Path [14] :Best Iter: [48] ELBO (-3689.758772) evaluations: (2880) 
#> Path [15] :Initial log joint density = -481799.550882 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.978e-03   2.396e-01    1.000e+00  1.000e+00      3503 -3.680e+03 -3.685e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3680.379895) evaluations: (3503) 
#> Path [16] :Initial log joint density = -481424.443096 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.925e-03   2.465e-01    1.000e+00  1.000e+00      3243 -3.690e+03 -3.684e+03                   
#> Path [16] :Best Iter: [57] ELBO (-3684.149135) evaluations: (3243) 
#> Path [17] :Initial log joint density = -482348.876413 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      3.637e-03   2.168e-01    7.134e-01  7.134e-01      3393 -3.685e+03 -3.694e+03                   
#> Path [17] :Best Iter: [56] ELBO (-3685.172762) evaluations: (3393) 
#> Path [18] :Initial log joint density = -481414.738360 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      4.089e-03   1.503e-01    7.885e-01  7.885e-01      3519 -3.684e+03 -3.696e+03                   
#> Path [18] :Best Iter: [59] ELBO (-3683.759551) evaluations: (3519) 
#> Path [19] :Initial log joint density = -481466.357223 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.211e-02   3.288e-01    1.000e+00  1.000e+00      3202 -3.689e+03 -3.690e+03                   
#> Path [19] :Best Iter: [53] ELBO (-3689.414810) evaluations: (3202) 
#> Path [20] :Initial log joint density = -481361.734937 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.519e-02   4.122e-01    1.000e+00  1.000e+00      2910 -3.692e+03 -3.699e+03                   
#> Path [20] :Best Iter: [41] ELBO (-3692.033254) evaluations: (2910) 
#> Path [21] :Initial log joint density = -481586.613733 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      3.568e-03   2.895e-01    6.427e-01  6.427e-01      3132 -3.689e+03 -3.699e+03                   
#> Path [21] :Best Iter: [49] ELBO (-3689.124571) evaluations: (3132) 
#> Path [22] :Initial log joint density = -482618.442015 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.556e-03   3.291e-01    5.691e-01  5.691e-01      3417 -3.682e+03 -3.695e+03                   
#> Path [22] :Best Iter: [56] ELBO (-3682.229301) evaluations: (3417) 
#> Path [23] :Initial log joint density = -481664.238446 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.958e-03   2.241e-01    9.953e-01  9.953e-01      3143 -3.687e+03 -3.691e+03                   
#> Path [23] :Best Iter: [53] ELBO (-3687.242615) evaluations: (3143) 
#> Path [24] :Initial log joint density = -481626.863480 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.338e-03   3.200e-01    1.000e+00  1.000e+00      3037 -3.691e+03 -3.701e+03                   
#> Path [24] :Best Iter: [49] ELBO (-3691.056352) evaluations: (3037) 
#> Path [25] :Initial log joint density = -481653.219018 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.979e-03   2.142e-01    1.000e+00  1.000e+00      3581 -3.683e+03 -3.693e+03                   
#> Path [25] :Best Iter: [60] ELBO (-3682.603626) evaluations: (3581) 
#> Path [26] :Initial log joint density = -481546.647424 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      8.164e-03   1.689e-01    8.981e-01  8.981e-01      3683 -3.681e+03 -3.689e+03                   
#> Path [26] :Best Iter: [60] ELBO (-3680.842582) evaluations: (3683) 
#> Path [27] :Initial log joint density = -481474.200998 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.648e-03   3.673e-01    4.965e-01  1.000e+00      3046 -3.689e+03 -3.705e+03                   
#> Path [27] :Best Iter: [47] ELBO (-3689.169303) evaluations: (3046) 
#> Path [28] :Initial log joint density = -481659.435879 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.724e-02   2.566e-01    9.647e-01  9.647e-01      3207 -3.684e+03 -3.687e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3683.798448) evaluations: (3207) 
#> Path [29] :Initial log joint density = -482509.897564 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.214e-02   2.886e-01    1.000e+00  1.000e+00      3414 -3.684e+03 -3.683e+03                   
#> Path [29] :Best Iter: [59] ELBO (-3683.245984) evaluations: (3414) 
#> Path [30] :Initial log joint density = -481658.957825 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      2.498e-03   3.189e-01    6.248e-01  6.248e-01      2853 -3.689e+03 -3.709e+03                   
#> Path [30] :Best Iter: [47] ELBO (-3689.082932) evaluations: (2853) 
#> Path [31] :Initial log joint density = -482102.216436 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.556e-03   2.134e-01    8.121e-01  8.121e-01      3415 -3.684e+03 -3.693e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3683.650118) evaluations: (3415) 
#> Path [32] :Initial log joint density = -481516.347477 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      6.692e-03   2.771e-01    1.000e+00  1.000e+00      2715 -3.690e+03 -3.703e+03                   
#> Path [32] :Best Iter: [46] ELBO (-3689.833338) evaluations: (2715) 
#> Path [33] :Initial log joint density = -482153.339426 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      4.315e-03   1.730e-01    1.000e+00  1.000e+00      3467 -3.682e+03 -3.692e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3682.109380) evaluations: (3467) 
#> Path [34] :Initial log joint density = -482003.687435 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.399e-02   2.348e-01    1.000e+00  1.000e+00      3710 -3.682e+03 -3.684e+03                   
#> Path [34] :Best Iter: [58] ELBO (-3682.287823) evaluations: (3710) 
#> Path [35] :Initial log joint density = -483019.400711 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.135e-02   2.624e-01    1.000e+00  1.000e+00      3239 -3.682e+03 -3.686e+03                   
#> Path [35] :Best Iter: [56] ELBO (-3682.452873) evaluations: (3239) 
#> Path [36] :Initial log joint density = -483658.873725 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.115e-03   2.011e-01    1.000e+00  1.000e+00      3461 -3.685e+03 -3.687e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3684.839510) evaluations: (3461) 
#> Path [37] :Initial log joint density = -481671.301859 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.198e-03   3.679e-01    1.000e+00  1.000e+00      2950 -3.690e+03 -3.704e+03                   
#> Path [37] :Best Iter: [42] ELBO (-3690.171216) evaluations: (2950) 
#> Path [38] :Initial log joint density = -481477.336006 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.382e-02   4.039e-01    1.000e+00  1.000e+00      3299 -3.686e+03 -3.690e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3686.244334) evaluations: (3299) 
#> Path [39] :Initial log joint density = -481644.005787 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.959e-03   3.038e-01    9.569e-01  9.569e-01      2909 -3.689e+03 -3.699e+03                   
#> Path [39] :Best Iter: [50] ELBO (-3689.053853) evaluations: (2909) 
#> Path [40] :Initial log joint density = -481783.712309 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.042e-02   2.229e-01    1.000e+00  1.000e+00      3308 -3.685e+03 -3.688e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3684.523336) evaluations: (3308) 
#> Path [41] :Initial log joint density = -481666.816700 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.788e+05      9.521e-03   3.246e-01    1.000e+00  1.000e+00      3955 -3.683e+03 -3.688e+03                   
#> Path [41] :Best Iter: [64] ELBO (-3683.154372) evaluations: (3955) 
#> Path [42] :Initial log joint density = -481777.982027 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.950e-02   2.235e-01    1.000e+00  1.000e+00      3397 -3.680e+03 -3.687e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3680.150734) evaluations: (3397) 
#> Path [43] :Initial log joint density = -482312.335514 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.976e-03   1.867e-01    1.000e+00  1.000e+00      3537 -3.683e+03 -3.685e+03                   
#> Path [43] :Best Iter: [58] ELBO (-3682.663089) evaluations: (3537) 
#> Path [44] :Initial log joint density = -481553.523877 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      2.748e-03   2.377e-01    6.736e-01  6.736e-01      3072 -3.690e+03 -3.699e+03                   
#> Path [44] :Best Iter: [53] ELBO (-3690.397036) evaluations: (3072) 
#> Path [45] :Initial log joint density = -482855.608201 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.287e-02   3.283e-01    1.000e+00  1.000e+00      3363 -3.683e+03 -3.694e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3682.534745) evaluations: (3363) 
#> Path [46] :Initial log joint density = -481596.512626 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.365e-02   3.582e-01    9.905e-01  9.905e-01      3554 -3.682e+03 -3.687e+03                   
#> Path [46] :Best Iter: [58] ELBO (-3681.500654) evaluations: (3554) 
#> Path [47] :Initial log joint density = -481770.873208 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      4.757e-03   1.534e-01    8.530e-01  8.530e-01      3739 -3.683e+03 -3.693e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3682.862198) evaluations: (3739) 
#> Path [48] :Initial log joint density = -481652.837057 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.072e-03   2.422e-01    8.010e-01  8.010e-01      3391 -3.683e+03 -3.695e+03                   
#> Path [48] :Best Iter: [58] ELBO (-3682.934723) evaluations: (3391) 
#> Path [49] :Initial log joint density = -482017.729929 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.743e-03   3.449e-01    6.381e-01  6.381e-01      3330 -3.683e+03 -3.699e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3683.438318) evaluations: (3330) 
#> Path [50] :Initial log joint density = -481609.166997 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.788e+05      2.000e-02   1.139e-01    1.000e+00  1.000e+00      4178 -3.681e+03 -3.685e+03                   
#> Path [50] :Best Iter: [60] ELBO (-3681.400939) evaluations: (4178) 
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
#>  1 B1         10x_6K                      0.0383              191         1
#>  2 B1         10x_8K                      0.0433              216         1
#>  3 B1         GSE115189                   0.0483              241         1
#>  4 B1         SCP345_580                  0.0373              186         1
#>  5 B1         SCP345_860                  0.0458              229         1
#>  6 B1         SCP424_pbmc1                0.0973              486         1
#>  7 B1         SCP424_pbmc2                0.0474              236         1
#>  8 B1         SCP591                      0.0996              498         1
#>  9 B1         SI-GA-E5                    0.0138               68         1
#> 10 B1         SI-GA-E7                    0.0132               65         1
#> # ℹ 710 more rows
# }
```
