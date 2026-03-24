# Plot 1D Intervals for Cell-group Effects

This function creates a series of 1D interval plots for cell-group
effects, highlighting significant differences based on a given
significance threshold.

## Usage

``` r
plot_1D_intervals(
  .data,
  significance_threshold = 0.05,
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
  show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR")
)
```

## Arguments

- .data:

  Data frame containing the main data.

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences.

- test_composition_above_logit_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test. A value of 0.2 correspond to a change in cell proportion of 10%
  for a cell type with baseline proportion of 50%. That is, a cell type
  goes from 45% to 50%. When the baseline proportion is closer to 0 or 1
  this effect thrshold has consistent value in the logit uncontrained
  scale.

- show_fdr_message:

  Logical. Whether to show the Bayesian FDR interpretation message on
  the plot. Default is TRUE.

- significance_statistic:

  Character vector indicating which statistic to highlight. Default is
  "pH0".

## Value

A combined plot of 1D interval plots.

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."

# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
    ) |> 
    sccomp_test()
    
  # Example usage:
  my_plot = plot_1D_intervals(estimate)
    
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481571.451993 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.473e-02   3.238e-01    1.000e+00  1.000e+00      2783 -3.707e+03 -3.720e+03                   
#> Path [1] :Best Iter: [43] ELBO (-3707.368926) evaluations: (2783) 
#> Path [2] :Initial log joint density = -481457.890373 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.345e-03   1.480e-01    1.000e+00  1.000e+00      3756 -3.699e+03 -3.705e+03                   
#> Path [2] :Best Iter: [57] ELBO (-3698.518774) evaluations: (3756) 
#> Path [3] :Initial log joint density = -481784.308961 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.312e-02   4.003e-01    1.000e+00  1.000e+00      3487 -3.701e+03 -3.708e+03                   
#> Path [3] :Best Iter: [58] ELBO (-3700.776256) evaluations: (3487) 
#> Path [4] :Initial log joint density = -481850.228131 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.284e-02   3.416e-01    1.000e+00  1.000e+00      3073 -3.705e+03 -3.708e+03                   
#> Path [4] :Best Iter: [53] ELBO (-3705.218442) evaluations: (3073) 
#> Path [5] :Initial log joint density = -481879.321797 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.196e-02   2.238e-01    1.000e+00  1.000e+00      3344 -3.706e+03 -3.701e+03                   
#> Path [5] :Best Iter: [57] ELBO (-3701.363232) evaluations: (3344) 
#> Path [6] :Initial log joint density = -482218.240844 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.257e-02   3.054e-01    1.000e+00  1.000e+00      3420 -3.702e+03 -3.703e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3702.316986) evaluations: (3420) 
#> Path [7] :Initial log joint density = -481858.201534 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.105e-02   3.665e-01    1.000e+00  1.000e+00      3355 -3.700e+03 -3.706e+03                   
#> Path [7] :Best Iter: [58] ELBO (-3699.695502) evaluations: (3355) 
#> Path [8] :Initial log joint density = -485047.333381 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.167e-03   2.068e-01    1.000e+00  1.000e+00      3475 -3.700e+03 -3.700e+03                   
#> Path [8] :Best Iter: [55] ELBO (-3700.001727) evaluations: (3475) 
#> Path [9] :Initial log joint density = -481445.158826 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.234e-02   3.360e-01    1.000e+00  1.000e+00      3242 -3.699e+03 -3.703e+03                   
#> Path [9] :Best Iter: [56] ELBO (-3699.362517) evaluations: (3242) 
#> Path [10] :Initial log joint density = -481814.339317 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.681e-03   2.488e-01    7.718e-01  7.718e-01      2989 -3.708e+03 -3.724e+03                   
#> Path [10] :Best Iter: [52] ELBO (-3708.305504) evaluations: (2989) 
#> Path [11] :Initial log joint density = -481797.978056 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.276e-02   2.478e-01    1.000e+00  1.000e+00      3387 -3.704e+03 -3.702e+03                   
#> Path [11] :Best Iter: [59] ELBO (-3702.334824) evaluations: (3387) 
#> Path [12] :Initial log joint density = -481586.241410 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.420e-03   1.806e-01    6.342e-01  6.342e-01      3445 -3.702e+03 -3.715e+03                   
#> Path [12] :Best Iter: [58] ELBO (-3701.936028) evaluations: (3445) 
#> Path [13] :Initial log joint density = -482052.280732 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.285e-03   2.309e-01    1.000e+00  1.000e+00      3393 -3.701e+03 -3.706e+03                   
#> Path [13] :Best Iter: [56] ELBO (-3701.121320) evaluations: (3393) 
#> Path [14] :Initial log joint density = -481617.859144 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.902e-03   2.166e-01    8.742e-01  8.742e-01      2948 -3.706e+03 -3.717e+03                   
#> Path [14] :Best Iter: [50] ELBO (-3705.558884) evaluations: (2948) 
#> Path [15] :Initial log joint density = -481708.419726 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.004e-03   2.226e-01    1.000e+00  1.000e+00      3449 -3.701e+03 -3.707e+03                   
#> Path [15] :Best Iter: [56] ELBO (-3700.865012) evaluations: (3449) 
#> Path [16] :Initial log joint density = -481711.854545 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.264e-03   1.999e-01    7.465e-01  7.465e-01      3136 -3.698e+03 -3.713e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3698.298539) evaluations: (3136) 
#> Path [17] :Initial log joint density = -481722.610456 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.275e-03   2.035e-01    8.615e-01  8.615e-01      3123 -3.709e+03 -3.711e+03                   
#> Path [17] :Best Iter: [52] ELBO (-3708.513909) evaluations: (3123) 
#> Path [18] :Initial log joint density = -482013.691356 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      3.524e-03   2.729e-01    6.296e-01  6.296e-01      3624 -3.701e+03 -3.717e+03                   
#> Path [18] :Best Iter: [56] ELBO (-3700.514082) evaluations: (3624) 
#> Path [19] :Initial log joint density = -481551.132532 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.531e-03   2.998e-01    6.093e-01  6.093e-01      2995 -3.707e+03 -3.723e+03                   
#> Path [19] :Best Iter: [52] ELBO (-3706.881834) evaluations: (2995) 
#> Path [20] :Initial log joint density = -482335.329990 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      5.554e-03   2.374e-01    1.000e+00  1.000e+00      3591 -3.700e+03 -3.702e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3699.945307) evaluations: (3591) 
#> Path [21] :Initial log joint density = -481487.266382 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.772e-03   2.102e-01    1.000e+00  1.000e+00      3062 -3.708e+03 -3.720e+03                   
#> Path [21] :Best Iter: [47] ELBO (-3708.498871) evaluations: (3062) 
#> Path [22] :Initial log joint density = -481304.386302 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      9.025e-03   2.801e-01    1.000e+00  1.000e+00      2684 -3.710e+03 -3.722e+03                   
#> Path [22] :Best Iter: [42] ELBO (-3710.408082) evaluations: (2684) 
#> Path [23] :Initial log joint density = -481739.406122 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      6.555e-03   1.603e-01    9.098e-01  9.098e-01      3561 -3.699e+03 -3.712e+03                   
#> Path [23] :Best Iter: [59] ELBO (-3699.351238) evaluations: (3561) 
#> Path [24] :Initial log joint density = -481598.051446 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.487e-03   2.112e-01    1.000e+00  1.000e+00      2944 -3.705e+03 -3.711e+03                   
#> Path [24] :Best Iter: [52] ELBO (-3705.212862) evaluations: (2944) 
#> Path [25] :Initial log joint density = -481925.678637 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.403e-02   3.546e-01    1.000e+00  1.000e+00      2984 -3.709e+03 -3.713e+03                   
#> Path [25] :Best Iter: [45] ELBO (-3708.885211) evaluations: (2984) 
#> Path [26] :Initial log joint density = -481543.185107 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.631e-03   2.898e-01    5.570e-01  5.570e-01      3072 -3.706e+03 -3.717e+03                   
#> Path [26] :Best Iter: [52] ELBO (-3706.334296) evaluations: (3072) 
#> Path [27] :Initial log joint density = -482144.909337 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.389e-03   1.839e-01    8.168e-01  8.168e-01      3192 -3.702e+03 -3.713e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3702.124925) evaluations: (3192) 
#> Path [28] :Initial log joint density = -482685.795208 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.198e-03   2.232e-01    8.896e-01  8.896e-01      3514 -3.697e+03 -3.705e+03                   
#> Path [28] :Best Iter: [57] ELBO (-3697.429491) evaluations: (3514) 
#> Path [29] :Initial log joint density = -481665.030408 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.435e-02   2.328e-01    1.000e+00  1.000e+00      3183 -3.710e+03 -3.703e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3703.491376) evaluations: (3183) 
#> Path [30] :Initial log joint density = -481509.970499 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.594e-03   1.605e-01    1.000e+00  1.000e+00      2828 -3.707e+03 -3.713e+03                   
#> Path [30] :Best Iter: [48] ELBO (-3707.125721) evaluations: (2828) 
#> Path [31] :Initial log joint density = -481413.528867 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.971e-03   3.071e-01    1.000e+00  1.000e+00      3064 -3.706e+03 -3.707e+03                   
#> Path [31] :Best Iter: [53] ELBO (-3705.593938) evaluations: (3064) 
#> Path [32] :Initial log joint density = -482676.631307 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.261e-03   2.339e-01    4.904e-01  1.000e+00      3120 -3.702e+03 -3.709e+03                   
#> Path [32] :Best Iter: [55] ELBO (-3702.073204) evaluations: (3120) 
#> Path [33] :Initial log joint density = -481607.180632 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.010e-03   2.850e-01    1.000e+00  1.000e+00      3220 -3.701e+03 -3.711e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3700.825855) evaluations: (3220) 
#> Path [34] :Initial log joint density = -485001.027596 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.176e-02   2.596e-01    1.000e+00  1.000e+00      3569 -3.700e+03 -3.700e+03                   
#> Path [34] :Best Iter: [59] ELBO (-3699.914855) evaluations: (3569) 
#> Path [35] :Initial log joint density = -481619.837115 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.459e-03   2.320e-01    8.387e-01  8.387e-01      3355 -3.702e+03 -3.717e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3701.655848) evaluations: (3355) 
#> Path [36] :Initial log joint density = -481747.776837 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              50      -4.788e+05      7.378e-03   2.686e-01    8.618e-01  8.618e-01      2764 -3.708e+03 -3.722e+03                   
#> Path [36] :Best Iter: [43] ELBO (-3707.768276) evaluations: (2764) 
#> Path [37] :Initial log joint density = -481394.618150 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.877e-02   3.680e-01    1.000e+00  1.000e+00      3451 -3.698e+03 -3.707e+03                   
#> Path [37] :Best Iter: [59] ELBO (-3698.357623) evaluations: (3451) 
#> Path [38] :Initial log joint density = -481642.662554 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.033e-03   2.490e-01    1.000e+00  1.000e+00      2867 -3.707e+03 -3.717e+03                   
#> Path [38] :Best Iter: [47] ELBO (-3707.302841) evaluations: (2867) 
#> Path [39] :Initial log joint density = -481662.640328 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.865e-03   2.122e-01    8.711e-01  8.711e-01      3101 -3.705e+03 -3.711e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3704.948326) evaluations: (3101) 
#> Path [40] :Initial log joint density = -481797.608644 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.105e-02   3.030e-01    5.146e-01  1.000e+00      3161 -3.701e+03 -3.710e+03                   
#> Path [40] :Best Iter: [55] ELBO (-3700.610684) evaluations: (3161) 
#> Path [41] :Initial log joint density = -483018.924959 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.217e-03   2.369e-01    1.000e+00  1.000e+00      3072 -3.703e+03 -3.700e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3700.375535) evaluations: (3072) 
#> Path [42] :Initial log joint density = -481809.528801 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.260e-02   3.903e-01    1.000e+00  1.000e+00      3165 -3.705e+03 -3.715e+03                   
#> Path [42] :Best Iter: [49] ELBO (-3704.690359) evaluations: (3165) 
#> Path [43] :Initial log joint density = -481441.377757 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.295e-02   1.499e-01    1.000e+00  1.000e+00      3227 -3.701e+03 -3.696e+03                   
#> Path [43] :Best Iter: [56] ELBO (-3695.731676) evaluations: (3227) 
#> Path [44] :Initial log joint density = -482574.742180 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.115e-03   2.749e-01    1.000e+00  1.000e+00      3088 -3.708e+03 -3.709e+03                   
#> Path [44] :Best Iter: [51] ELBO (-3707.910271) evaluations: (3088) 
#> Path [45] :Initial log joint density = -481736.173632 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.732e-03   2.278e-01    8.189e-01  8.189e-01      3246 -3.701e+03 -3.708e+03                   
#> Path [45] :Best Iter: [55] ELBO (-3701.198785) evaluations: (3246) 
#> Path [46] :Initial log joint density = -481798.831014 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.228e-02   3.129e-01    9.802e-01  9.802e-01      2992 -3.708e+03 -3.719e+03                   
#> Path [46] :Best Iter: [43] ELBO (-3707.841922) evaluations: (2992) 
#> Path [47] :Initial log joint density = -481819.814260 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      2.822e-03   2.377e-01    7.351e-01  7.351e-01      2993 -3.707e+03 -3.724e+03                   
#> Path [47] :Best Iter: [43] ELBO (-3707.121653) evaluations: (2993) 
#> Path [48] :Initial log joint density = -481554.132195 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      4.590e-03   2.646e-01    7.070e-01  7.070e-01      3044 -3.706e+03 -3.719e+03                   
#> Path [48] :Best Iter: [52] ELBO (-3705.961135) evaluations: (3044) 
#> Path [49] :Initial log joint density = -481960.677038 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.182e-03   2.279e-01    1.000e+00  1.000e+00      3332 -3.708e+03 -3.710e+03                   
#> Path [49] :Best Iter: [51] ELBO (-3707.566825) evaluations: (3332) 
#> Path [50] :Initial log joint density = -481713.532048 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.446e-02   3.803e-01    1.000e+00  1.000e+00      3042 -3.704e+03 -3.716e+03                   
#> Path [50] :Best Iter: [52] ELBO (-3704.003838) evaluations: (3042) 
#> Finished in  13.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }

```
