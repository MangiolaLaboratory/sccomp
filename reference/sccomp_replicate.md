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
#> Path [1] :Initial log joint density = -481509.193778 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.369e-03   2.865e-01    1.000e+00  1.000e+00      3212 -3.688e+03 -3.696e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3688.291667) evaluations: (3212) 
#> Path [2] :Initial log joint density = -481650.295269 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      7.132e-03   2.389e-01    8.617e-01  8.617e-01      3793 -3.688e+03 -3.698e+03                   
#> Path [2] :Best Iter: [62] ELBO (-3687.828259) evaluations: (3793) 
#> Path [3] :Initial log joint density = -482211.014011 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.919e-03   1.435e-01    1.000e+00  1.000e+00      4087 -3.685e+03 -3.689e+03                   
#> Path [3] :Best Iter: [57] ELBO (-3684.501836) evaluations: (4087) 
#> Path [4] :Initial log joint density = -484179.852763 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.818e-03   2.480e-01    1.000e+00  1.000e+00      3275 -3.689e+03 -3.699e+03                   
#> Path [4] :Best Iter: [55] ELBO (-3689.043632) evaluations: (3275) 
#> Path [5] :Initial log joint density = -481691.124897 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.261e-02   2.281e-01    4.904e-01  1.000e+00      4377 -3.687e+03 -3.695e+03                   
#> Path [5] :Best Iter: [60] ELBO (-3686.749927) evaluations: (4377) 
#> Path [6] :Initial log joint density = -482211.338837 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.255e-02   3.484e-01    1.000e+00  1.000e+00      5151 -3.683e+03 -3.692e+03                   
#> Path [6] :Best Iter: [74] ELBO (-3683.131296) evaluations: (5151) 
#> Path [7] :Initial log joint density = -482117.003664 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.291e-02   2.257e-01    6.919e-01  6.919e-01      4052 -3.687e+03 -3.697e+03                   
#> Path [7] :Best Iter: [64] ELBO (-3687.028598) evaluations: (4052) 
#> Path [8] :Initial log joint density = -482765.434415 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.155e-02   1.942e-01    8.380e-01  8.380e-01      4839 -3.681e+03 -3.695e+03                   
#> Path [8] :Best Iter: [72] ELBO (-3680.924775) evaluations: (4839) 
#> Path [9] :Initial log joint density = -482025.469614 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.086e-02   1.979e-01    9.001e-01  9.001e-01      4491 -3.686e+03 -3.697e+03                   
#> Path [9] :Best Iter: [67] ELBO (-3686.140945) evaluations: (4491) 
#> Path [10] :Initial log joint density = -481574.375501 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      6.689e-03   1.549e-01    1.000e+00  1.000e+00      4021 -3.685e+03 -3.696e+03                   
#> Path [10] :Best Iter: [63] ELBO (-3685.287791) evaluations: (4021) 
#> Path [11] :Initial log joint density = -481749.569949 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.158e-02   2.175e-01    1.000e+00  1.000e+00      3504 -3.687e+03 -3.685e+03                   
#> Path [11] :Best Iter: [60] ELBO (-3684.574855) evaluations: (3504) 
#> Path [12] :Initial log joint density = -482069.422027 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      6.302e-03   2.167e-01    7.412e-01  7.412e-01      4842 -3.684e+03 -3.702e+03                   
#> Path [12] :Best Iter: [72] ELBO (-3684.412046) evaluations: (4842) 
#> Path [13] :Initial log joint density = -481428.670112 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      5.909e-03   1.961e-01    1.000e+00  1.000e+00      3294 -3.696e+03 -3.697e+03                   
#> Path [13] :Best Iter: [53] ELBO (-3696.403918) evaluations: (3294) 
#> Path [14] :Initial log joint density = -481871.853925 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      8.467e-03   2.600e-01    9.562e-01  9.562e-01      4736 -3.688e+03 -3.699e+03                   
#> Path [14] :Best Iter: [59] ELBO (-3688.278126) evaluations: (4736) 
#> Path [15] :Initial log joint density = -481448.165613 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      6.394e-03   2.154e-01    1.000e+00  1.000e+00      3203 -3.697e+03 -3.696e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3696.192733) evaluations: (3203) 
#> Path [16] :Initial log joint density = -481785.919013 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.172e-02   1.904e-01    1.000e+00  1.000e+00      4305 -3.687e+03 -3.688e+03                   
#> Path [16] :Best Iter: [66] ELBO (-3686.698789) evaluations: (4305) 
#> Path [17] :Initial log joint density = -481256.857059 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.367e-02   1.789e-01    1.000e+00  1.000e+00      4003 -3.685e+03 -3.690e+03                   
#> Path [17] :Best Iter: [61] ELBO (-3684.750495) evaluations: (4003) 
#> Path [18] :Initial log joint density = -481617.029242 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.471e-02   1.957e-01    1.000e+00  1.000e+00      4536 -3.685e+03 -3.687e+03                   
#> Path [18] :Best Iter: [62] ELBO (-3685.453902) evaluations: (4536) 
#> Path [19] :Initial log joint density = -481569.884508 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.797e-03   1.961e-01    1.000e+00  1.000e+00      3221 -3.697e+03 -3.694e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3693.887939) evaluations: (3221) 
#> Path [20] :Initial log joint density = -482114.949798 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      7.915e-03   2.186e-01    1.000e+00  1.000e+00      3382 -3.689e+03 -3.690e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3689.309628) evaluations: (3382) 
#> Path [21] :Initial log joint density = -481884.661343 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      8.793e-03   2.727e-01    7.578e-01  7.578e-01      3217 -3.689e+03 -3.702e+03                   
#> Path [21] :Best Iter: [55] ELBO (-3689.360592) evaluations: (3217) 
#> Path [22] :Initial log joint density = -481606.427423 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.158e-02   2.631e-01    9.675e-01  9.675e-01      4385 -3.683e+03 -3.692e+03                   
#> Path [22] :Best Iter: [68] ELBO (-3682.827812) evaluations: (4385) 
#> Path [23] :Initial log joint density = -481608.362142 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.474e-02   3.461e-01    9.834e-01  9.834e-01      5242 -3.684e+03 -3.693e+03                   
#> Path [23] :Best Iter: [74] ELBO (-3683.687317) evaluations: (5242) 
#> Path [24] :Initial log joint density = -482286.523954 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.078e-02   1.825e-01    1.000e+00  1.000e+00      4986 -3.686e+03 -3.687e+03                   
#> Path [24] :Best Iter: [71] ELBO (-3685.815438) evaluations: (4986) 
#> Path [25] :Initial log joint density = -481650.740155 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      4.055e-02   2.044e-01    1.000e+00  1.000e+00      4250 -3.685e+03 -3.685e+03                   
#> Path [25] :Best Iter: [59] ELBO (-3685.092691) evaluations: (4250) 
#> Path [26] :Initial log joint density = -481765.116634 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      7.215e-03   2.430e-01    1.000e+00  1.000e+00      4500 -3.689e+03 -3.691e+03                   
#> Path [26] :Best Iter: [69] ELBO (-3689.033718) evaluations: (4500) 
#> Path [27] :Initial log joint density = -481730.676175 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      7.452e-03   2.845e-01    4.864e-01  1.000e+00      3235 -3.692e+03 -3.697e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3692.426977) evaluations: (3235) 
#> Path [28] :Initial log joint density = -482526.014471 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.815e-03   1.562e-01    1.000e+00  1.000e+00      3417 -3.687e+03 -3.691e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3687.108624) evaluations: (3417) 
#> Path [29] :Initial log joint density = -481428.764794 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      9.280e-03   3.011e-01    4.114e-01  1.000e+00      2917 -3.692e+03 -3.710e+03                   
#> Path [29] :Best Iter: [47] ELBO (-3692.402644) evaluations: (2917) 
#> Path [30] :Initial log joint density = -481878.264340 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      7.923e-03   2.011e-01    1.000e+00  1.000e+00      3877 -3.690e+03 -3.695e+03                   
#> Path [30] :Best Iter: [56] ELBO (-3690.261691) evaluations: (3877) 
#> Path [31] :Initial log joint density = -481375.733166 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.046e-02   2.041e-01    1.000e+00  1.000e+00      3646 -3.688e+03 -3.689e+03                   
#> Path [31] :Best Iter: [58] ELBO (-3688.124337) evaluations: (3646) 
#> Path [32] :Initial log joint density = -482693.492943 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.850e-03   2.875e-01    7.824e-01  7.824e-01      3160 -3.691e+03 -3.702e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3691.283641) evaluations: (3160) 
#> Path [33] :Initial log joint density = -481806.879752 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.301e-02   2.428e-01    1.000e+00  1.000e+00      4438 -3.684e+03 -3.687e+03                   
#> Path [33] :Best Iter: [68] ELBO (-3683.687827) evaluations: (4438) 
#> Path [34] :Initial log joint density = -481519.877944 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.320e-02   1.894e-01    1.000e+00  1.000e+00      4411 -3.686e+03 -3.687e+03                   
#> Path [34] :Best Iter: [60] ELBO (-3685.694622) evaluations: (4411) 
#> Path [35] :Initial log joint density = -481454.285805 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.202e-02   2.057e-01    1.000e+00  1.000e+00      4354 -3.683e+03 -3.691e+03                   
#> Path [35] :Best Iter: [64] ELBO (-3683.344511) evaluations: (4354) 
#> Path [36] :Initial log joint density = -481492.473830 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.389e-02   1.828e-01    1.000e+00  1.000e+00      4496 -3.685e+03 -3.684e+03                   
#> Path [36] :Best Iter: [70] ELBO (-3684.021427) evaluations: (4496) 
#> Path [37] :Initial log joint density = -481890.652354 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      5.006e-03   2.220e-01    7.160e-01  7.160e-01      4407 -3.684e+03 -3.694e+03                   
#> Path [37] :Best Iter: [66] ELBO (-3684.470912) evaluations: (4407) 
#> Path [38] :Initial log joint density = -483882.461296 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.063e-02   2.324e-01    8.057e-01  8.057e-01      4054 -3.682e+03 -3.698e+03                   
#> Path [38] :Best Iter: [63] ELBO (-3682.386243) evaluations: (4054) 
#> Path [39] :Initial log joint density = -481525.760723 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.826e-02   2.776e-01    1.000e+00  1.000e+00      4538 -3.683e+03 -3.689e+03                   
#> Path [39] :Best Iter: [69] ELBO (-3683.028774) evaluations: (4538) 
#> Path [40] :Initial log joint density = -481494.936187 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      9.809e-03   2.163e-01    1.000e+00  1.000e+00      3842 -3.685e+03 -3.688e+03                   
#> Path [40] :Best Iter: [61] ELBO (-3684.828798) evaluations: (3842) 
#> Path [41] :Initial log joint density = -482052.229448 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      6.765e-03   2.086e-01    1.000e+00  1.000e+00      4191 -3.688e+03 -3.699e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3687.576931) evaluations: (4191) 
#> Path [42] :Initial log joint density = -481405.825353 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      7.560e-03   2.066e-01    8.528e-01  8.528e-01      4684 -3.682e+03 -3.690e+03                   
#> Path [42] :Best Iter: [70] ELBO (-3681.708579) evaluations: (4684) 
#> Path [43] :Initial log joint density = -481353.272603 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.112e-02   3.905e-01    1.000e+00  1.000e+00      4073 -3.686e+03 -3.693e+03                   
#> Path [43] :Best Iter: [64] ELBO (-3685.595818) evaluations: (4073) 
#> Path [44] :Initial log joint density = -481627.677556 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.327e-02   2.221e-01    1.000e+00  1.000e+00      4116 -3.685e+03 -3.695e+03                   
#> Path [44] :Best Iter: [65] ELBO (-3684.866540) evaluations: (4116) 
#> Path [45] :Initial log joint density = -482514.655888 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      8.391e-03   2.260e-01    1.000e+00  1.000e+00      3649 -3.688e+03 -3.694e+03                   
#> Path [45] :Best Iter: [56] ELBO (-3687.838168) evaluations: (3649) 
#> Path [46] :Initial log joint density = -482267.896191 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.294e-02   3.024e-01    1.000e+00  1.000e+00      5146 -3.686e+03 -3.690e+03                   
#> Path [46] :Best Iter: [76] ELBO (-3685.745741) evaluations: (5146) 
#> Path [47] :Initial log joint density = -481693.943593 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      9.511e-03   1.729e-01    1.000e+00  1.000e+00      3828 -3.689e+03 -3.686e+03                   
#> Path [47] :Best Iter: [63] ELBO (-3685.988695) evaluations: (3828) 
#> Path [48] :Initial log joint density = -481726.663113 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.701e-02   1.794e-01    1.000e+00  1.000e+00      4151 -3.684e+03 -3.687e+03                   
#> Path [48] :Best Iter: [60] ELBO (-3683.693974) evaluations: (4151) 
#> Path [49] :Initial log joint density = -481155.819549 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      8.155e-03   2.533e-01    7.038e-01  7.038e-01      4020 -3.683e+03 -3.697e+03                   
#> Path [49] :Best Iter: [58] ELBO (-3682.885674) evaluations: (4020) 
#> Path [50] :Initial log joint density = -481845.832140 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      9.050e-03   2.577e-01    1.000e+00  1.000e+00      4980 -3.688e+03 -3.690e+03                   
#> Path [50] :Best Iter: [73] ELBO (-3688.443319) evaluations: (4980) 
#> Finished in  15.6 seconds.
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
#>  1 B1         10x_6K                      0.0563              281         1
#>  2 B1         10x_8K                      0.0359              179         1
#>  3 B1         GSE115189                   0.0338              169         1
#>  4 B1         SCP345_580                  0.0493              246         1
#>  5 B1         SCP345_860                  0.0609              304         1
#>  6 B1         SCP424_pbmc1                0.0429              214         1
#>  7 B1         SCP424_pbmc2                0.0123               61         1
#>  8 B1         SCP591                      0.0566              283         1
#>  9 B1         SI-GA-E5                    0.0154               76         1
#> 10 B1         SI-GA-E7                    0.0202              101         1
#> # ℹ 710 more rows
# }
```
