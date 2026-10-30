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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481490.726881 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.760e-02   2.578e-01    1.000e+00  1.000e+00      3747 -3.702e+03 -3.706e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3702.367911) evaluations: (3747) 
#> Path [2] :Initial log joint density = -481633.079460 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.206e-02   2.450e-01    1.000e+00  1.000e+00      4562 -3.699e+03 -3.707e+03                   
#> Path [2] :Best Iter: [69] ELBO (-3698.520929) evaluations: (4562) 
#> Path [3] :Initial log joint density = -482196.786309 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.018e-02   2.659e-01    1.000e+00  1.000e+00      4148 -3.700e+03 -3.699e+03                   
#> Path [3] :Best Iter: [67] ELBO (-3698.758243) evaluations: (4148) 
#> Path [4] :Initial log joint density = -484166.697174 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      3.115e-03   1.984e-01    7.428e-01  7.428e-01      4009 -3.701e+03 -3.712e+03                   
#> Path [4] :Best Iter: [60] ELBO (-3701.360365) evaluations: (4009) 
#> Path [5] :Initial log joint density = -481673.579370 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.963e-02   2.721e-01    1.000e+00  1.000e+00      3998 -3.701e+03 -3.706e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3700.514785) evaluations: (3998) 
#> Path [6] :Initial log joint density = -482195.164524 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.268e-02   2.412e-01    9.550e-01  9.550e-01      4659 -3.700e+03 -3.711e+03                   
#> Path [6] :Best Iter: [68] ELBO (-3699.832539) evaluations: (4659) 
#> Path [7] :Initial log joint density = -482100.691308 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      9.774e-03   1.892e-01    9.211e-01  9.211e-01      4132 -3.706e+03 -3.712e+03                   
#> Path [7] :Best Iter: [63] ELBO (-3706.183579) evaluations: (4132) 
#> Path [8] :Initial log joint density = -482746.041182 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.395e-03   1.942e-01    1.000e+00  1.000e+00      4654 -3.705e+03 -3.713e+03                   
#> Path [8] :Best Iter: [68] ELBO (-3705.176605) evaluations: (4654) 
#> Path [9] :Initial log joint density = -482010.973115 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      8.164e-03   1.915e-01    1.000e+00  1.000e+00      5581 -3.698e+03 -3.708e+03                   
#> Path [9] :Best Iter: [74] ELBO (-3698.319157) evaluations: (5581) 
#> Path [10] :Initial log joint density = -481556.885910 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      8.454e-03   2.151e-01    1.000e+00  1.000e+00      4723 -3.703e+03 -3.704e+03                   
#> Path [10] :Best Iter: [72] ELBO (-3702.602810) evaluations: (4723) 
#> Path [11] :Initial log joint density = -481736.394616 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      2.546e-03   2.145e-01    6.673e-01  6.673e-01      3416 -3.705e+03 -3.718e+03                   
#> Path [11] :Best Iter: [56] ELBO (-3704.886169) evaluations: (3416) 
#> Path [12] :Initial log joint density = -482056.783431 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      2.022e-03   2.394e-01    5.882e-01  5.882e-01      5154 -3.699e+03 -3.716e+03                   
#> Path [12] :Best Iter: [75] ELBO (-3698.660636) evaluations: (5154) 
#> Path [13] :Initial log joint density = -481413.966385 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      8.429e-03   1.972e-01    1.000e+00  1.000e+00      3539 -3.704e+03 -3.705e+03                   
#> Path [13] :Best Iter: [58] ELBO (-3704.139128) evaluations: (3539) 
#> Path [14] :Initial log joint density = -481852.735217 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      2.261e-02   3.507e-01    1.000e+00  1.000e+00      5286 -3.700e+03 -3.710e+03                   
#> Path [14] :Best Iter: [76] ELBO (-3699.525132) evaluations: (5286) 
#> Path [15] :Initial log joint density = -481434.973319 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.244e-03   1.953e-01    1.000e+00  1.000e+00      3872 -3.706e+03 -3.707e+03                   
#> Path [15] :Best Iter: [60] ELBO (-3705.621043) evaluations: (3872) 
#> Path [16] :Initial log joint density = -481769.133125 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      7.497e-03   1.277e-01    1.000e+00  1.000e+00      5307 -3.700e+03 -3.708e+03                   
#> Path [16] :Best Iter: [75] ELBO (-3699.712821) evaluations: (5307) 
#> Path [17] :Initial log joint density = -481240.399967 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      9.293e-03   1.917e-01    1.000e+00  1.000e+00      3686 -3.700e+03 -3.700e+03                   
#> Path [17] :Best Iter: [62] ELBO (-3700.301028) evaluations: (3686) 
#> Path [18] :Initial log joint density = -481595.015351 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.664e-02   2.321e-01    1.000e+00  1.000e+00      5302 -3.700e+03 -3.703e+03                   
#> Path [18] :Best Iter: [75] ELBO (-3700.171083) evaluations: (5302) 
#> Path [19] :Initial log joint density = -481554.066685 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      8.544e-03   2.196e-01    1.000e+00  1.000e+00      3711 -3.702e+03 -3.711e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3701.665170) evaluations: (3711) 
#> Path [20] :Initial log joint density = -482103.001917 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.965e-02   3.701e-01    1.000e+00  1.000e+00      4083 -3.702e+03 -3.708e+03                   
#> Path [20] :Best Iter: [65] ELBO (-3702.394397) evaluations: (4083) 
#> Path [21] :Initial log joint density = -481866.700527 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.373e-02   1.397e-01    1.000e+00  1.000e+00      4081 -3.698e+03 -3.701e+03                   
#> Path [21] :Best Iter: [61] ELBO (-3698.479828) evaluations: (4081) 
#> Path [22] :Initial log joint density = -481585.482280 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.030e-02   2.262e-01    1.000e+00  1.000e+00      4507 -3.700e+03 -3.709e+03                   
#> Path [22] :Best Iter: [65] ELBO (-3700.158059) evaluations: (4507) 
#> Path [23] :Initial log joint density = -481595.396711 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      4.441e-03   2.105e-01    7.276e-01  7.276e-01      6027 -3.697e+03 -3.711e+03                   
#> Path [23] :Best Iter: [81] ELBO (-3696.853301) evaluations: (6027) 
#> Path [24] :Initial log joint density = -482270.768959 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      3.074e-03   3.037e-01    5.602e-01  5.602e-01      5439 -3.697e+03 -3.715e+03                   
#> Path [24] :Best Iter: [76] ELBO (-3697.489295) evaluations: (5439) 
#> Path [25] :Initial log joint density = -481629.946715 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.305e-02   2.151e-01    9.834e-01  9.834e-01      4308 -3.701e+03 -3.706e+03                   
#> Path [25] :Best Iter: [66] ELBO (-3701.187702) evaluations: (4308) 
#> Path [26] :Initial log joint density = -481750.950650 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      8.981e-03   1.856e-01    1.000e+00  1.000e+00      5220 -3.701e+03 -3.700e+03                   
#> Path [26] :Best Iter: [78] ELBO (-3699.914635) evaluations: (5220) 
#> Path [27] :Initial log joint density = -481718.168671 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.769e-02   2.594e-01    1.000e+00  1.000e+00      3641 -3.702e+03 -3.712e+03                   
#> Path [27] :Best Iter: [59] ELBO (-3701.631255) evaluations: (3641) 
#> Path [28] :Initial log joint density = -482513.876068 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      3.084e-03   2.235e-01    6.992e-01  6.992e-01      3277 -3.704e+03 -3.713e+03                   
#> Path [28] :Best Iter: [55] ELBO (-3704.038916) evaluations: (3277) 
#> Path [29] :Initial log joint density = -481411.683769 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      3.892e-03   1.942e-01    6.526e-01  6.526e-01      3882 -3.702e+03 -3.716e+03                   
#> Path [29] :Best Iter: [61] ELBO (-3702.278230) evaluations: (3882) 
#> Path [30] :Initial log joint density = -481858.395500 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      4.893e-03   1.919e-01    4.100e-01  1.000e+00      6199 -3.698e+03 -3.711e+03                   
#> Path [30] :Best Iter: [85] ELBO (-3698.005392) evaluations: (6199) 
#> Path [31] :Initial log joint density = -481360.044438 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      7.426e-03   2.103e-01    9.318e-01  9.318e-01      3556 -3.702e+03 -3.712e+03                   
#> Path [31] :Best Iter: [58] ELBO (-3701.603324) evaluations: (3556) 
#> Path [32] :Initial log joint density = -482676.917300 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      7.960e-03   2.013e-01    1.000e+00  1.000e+00      3642 -3.702e+03 -3.711e+03                   
#> Path [32] :Best Iter: [57] ELBO (-3701.891729) evaluations: (3642) 
#> Path [33] :Initial log joint density = -481791.689124 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.897e-02   2.709e-01    1.000e+00  1.000e+00      5328 -3.699e+03 -3.701e+03                   
#> Path [33] :Best Iter: [74] ELBO (-3698.850234) evaluations: (5328) 
#> Path [34] :Initial log joint density = -481501.682066 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.516e-02   2.307e-01    1.000e+00  1.000e+00      4799 -3.699e+03 -3.704e+03                   
#> Path [34] :Best Iter: [66] ELBO (-3698.944227) evaluations: (4799) 
#> Path [35] :Initial log joint density = -481440.262124 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.398e-02   2.473e-01    1.000e+00  1.000e+00      4261 -3.700e+03 -3.704e+03                   
#> Path [35] :Best Iter: [67] ELBO (-3700.451219) evaluations: (4261) 
#> Path [36] :Initial log joint density = -481476.752258 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      4.335e-03   2.163e-01    8.654e-01  8.654e-01      3952 -3.703e+03 -3.715e+03                   
#> Path [36] :Best Iter: [62] ELBO (-3703.251859) evaluations: (3952) 
#> Path [37] :Initial log joint density = -481871.332788 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.097e-02   2.928e-01    9.172e-01  9.172e-01      4007 -3.698e+03 -3.711e+03                   
#> Path [37] :Best Iter: [63] ELBO (-3697.969794) evaluations: (4007) 
#> Path [38] :Initial log joint density = -483870.131498 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.301e-02   2.997e-01    7.584e-01  7.584e-01      4256 -3.700e+03 -3.714e+03                   
#> Path [38] :Best Iter: [62] ELBO (-3699.648050) evaluations: (4256) 
#> Path [39] :Initial log joint density = -481514.172893 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      7.523e-03   2.175e-01    1.000e+00  1.000e+00      4088 -3.703e+03 -3.706e+03                   
#> Path [39] :Best Iter: [64] ELBO (-3702.701283) evaluations: (4088) 
#> Path [40] :Initial log joint density = -481481.463540 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.665e-03   1.319e-01    1.000e+00  1.000e+00      4782 -3.699e+03 -3.709e+03                   
#> Path [40] :Best Iter: [70] ELBO (-3698.812197) evaluations: (4782) 
#> Path [41] :Initial log joint density = -482032.105563 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      6.103e-03   1.786e-01    1.000e+00  1.000e+00      4501 -3.704e+03 -3.704e+03                   
#> Path [41] :Best Iter: [70] ELBO (-3703.833459) evaluations: (4501) 
#> Path [42] :Initial log joint density = -481391.120786 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.314e-03   2.695e-01    1.000e+00  1.000e+00      3160 -3.709e+03 -3.716e+03                   
#> Path [42] :Best Iter: [49] ELBO (-3709.276562) evaluations: (3160) 
#> Path [43] :Initial log joint density = -481336.101711 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.337e-02   2.088e-01    1.000e+00  1.000e+00      4958 -3.700e+03 -3.703e+03                   
#> Path [43] :Best Iter: [72] ELBO (-3700.428206) evaluations: (4958) 
#> Path [44] :Initial log joint density = -481612.726892 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.626e-02   2.668e-01    1.000e+00  1.000e+00      4803 -3.697e+03 -3.701e+03                   
#> Path [44] :Best Iter: [73] ELBO (-3697.118721) evaluations: (4803) 
#> Path [45] :Initial log joint density = -482500.122862 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.043e-02   2.502e-01    1.000e+00  1.000e+00      3730 -3.702e+03 -3.711e+03                   
#> Path [45] :Best Iter: [60] ELBO (-3702.261715) evaluations: (3730) 
#> Path [46] :Initial log joint density = -482251.808969 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      8.999e-03   2.218e-01    1.000e+00  1.000e+00      5862 -3.699e+03 -3.708e+03                   
#> Path [46] :Best Iter: [80] ELBO (-3698.571007) evaluations: (5862) 
#> Path [47] :Initial log joint density = -481676.682581 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      7.142e-03   2.176e-01    1.000e+00  1.000e+00      3746 -3.704e+03 -3.711e+03                   
#> Path [47] :Best Iter: [58] ELBO (-3704.227300) evaluations: (3746) 
#> Path [48] :Initial log joint density = -481709.822419 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      5.924e-03   2.073e-01    1.000e+00  1.000e+00      3832 -3.699e+03 -3.701e+03                   
#> Path [48] :Best Iter: [60] ELBO (-3699.472709) evaluations: (3832) 
#> Path [49] :Initial log joint density = -481143.820822 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      5.911e-03   2.328e-01    1.000e+00  1.000e+00      3926 -3.700e+03 -3.712e+03                   
#> Path [49] :Best Iter: [61] ELBO (-3700.325393) evaluations: (3926) 
#> Path [50] :Initial log joint density = -481828.758357 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      5.747e-03   1.950e-01    7.310e-01  7.310e-01      5651 -3.698e+03 -3.712e+03                   
#> Path [50] :Best Iter: [79] ELBO (-3698.201855) evaluations: (5651) 
#> Finished in  14.8 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
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
#>  1 B1         10x_6K                      0.0741              370         1
#>  2 B1         10x_8K                      0.0616              308         1
#>  3 B1         GSE115189                   0.0953              476         1
#>  4 B1         SCP345_580                  0.0511              255         1
#>  5 B1         SCP345_860                  0.0644              322         1
#>  6 B1         SCP424_pbmc1                0.0307              153         1
#>  7 B1         SCP424_pbmc2                0.0432              216         1
#>  8 B1         SCP591                      0.0351              175         1
#>  9 B1         SI-GA-E5                    0.0238              118         1
#> 10 B1         SI-GA-E7                    0.0327              163         1
#> # ℹ 710 more rows
# }
```
