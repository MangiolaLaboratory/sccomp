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
  significance_statistic = c("FDR", "pH0")
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
  "FDR".

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
#> Path [1] :Initial log joint density = -481538.497909 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.106e-03   2.266e-01    7.183e-01  7.183e-01      3327 -3.703e+03 -3.710e+03                   
#> Path [1] :Best Iter: [56] ELBO (-3702.618174) evaluations: (3327) 
#> Path [2] :Initial log joint density = -481504.392291 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.444e-02   2.581e-01    1.000e+00  1.000e+00      3109 -3.699e+03 -3.704e+03                   
#> Path [2] :Best Iter: [55] ELBO (-3699.193124) evaluations: (3109) 
#> Path [3] :Initial log joint density = -481738.716172 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.197e-03   2.091e-01    1.000e+00  1.000e+00      3468 -3.702e+03 -3.709e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3701.780969) evaluations: (3468) 
#> Path [4] :Initial log joint density = -482361.778686 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.678e-03   2.996e-01    4.063e-01  1.000e+00      3278 -3.703e+03 -3.715e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3703.456210) evaluations: (3278) 
#> Path [5] :Initial log joint density = -481531.249410 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      5.101e-03   2.158e-01    1.000e+00  1.000e+00      2704 -3.709e+03 -3.721e+03                   
#> Path [5] :Best Iter: [33] ELBO (-3709.381746) evaluations: (2704) 
#> Path [6] :Initial log joint density = -481617.549582 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.831e-03   1.862e-01    5.122e-01  1.000e+00      3173 -3.711e+03 -3.720e+03                   
#> Path [6] :Best Iter: [42] ELBO (-3711.108676) evaluations: (3173) 
#> Path [7] :Initial log joint density = -481680.537667 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.166e-03   2.164e-01    6.981e-01  6.981e-01      2863 -3.705e+03 -3.724e+03                   
#> Path [7] :Best Iter: [51] ELBO (-3705.333508) evaluations: (2863) 
#> Path [8] :Initial log joint density = -481630.409563 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.447e-03   2.286e-01    1.000e+00  1.000e+00      2874 -3.709e+03 -3.716e+03                   
#> Path [8] :Best Iter: [48] ELBO (-3708.532817) evaluations: (2874) 
#> Path [9] :Initial log joint density = -483534.184390 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      2.068e-02   2.853e-01    1.000e+00  1.000e+00      3611 -3.697e+03 -3.701e+03                   
#> Path [9] :Best Iter: [60] ELBO (-3697.356292) evaluations: (3611) 
#> Path [10] :Initial log joint density = -481534.546050 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.248e-02   2.686e-01    1.000e+00  1.000e+00      3422 -3.698e+03 -3.703e+03                   
#> Path [10] :Best Iter: [57] ELBO (-3698.076720) evaluations: (3422) 
#> Path [11] :Initial log joint density = -481339.610799 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.546e-02   2.247e-01    1.000e+00  1.000e+00      3520 -3.697e+03 -3.702e+03                   
#> Path [11] :Best Iter: [59] ELBO (-3696.676688) evaluations: (3520) 
#> Path [12] :Initial log joint density = -481520.886449 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      8.162e-03   2.647e-01    1.000e+00  1.000e+00      2876 -3.709e+03 -3.706e+03                   
#> Path [12] :Best Iter: [52] ELBO (-3705.760576) evaluations: (2876) 
#> Path [13] :Initial log joint density = -483046.868786 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.828e-03   2.188e-01    1.000e+00  1.000e+00      3649 -3.699e+03 -3.700e+03                   
#> Path [13] :Best Iter: [57] ELBO (-3699.041194) evaluations: (3649) 
#> Path [14] :Initial log joint density = -481719.817104 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.586e-02   2.673e-01    1.000e+00  1.000e+00      3136 -3.698e+03 -3.703e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3697.662174) evaluations: (3136) 
#> Path [15] :Initial log joint density = -481443.459531 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.066e-02   3.599e-01    1.000e+00  1.000e+00      3026 -3.708e+03 -3.708e+03                   
#> Path [15] :Best Iter: [55] ELBO (-3707.703186) evaluations: (3026) 
#> Path [16] :Initial log joint density = -481707.040474 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.233e-02   3.941e-01    1.000e+00  1.000e+00      2944 -3.705e+03 -3.715e+03                   
#> Path [16] :Best Iter: [50] ELBO (-3704.913616) evaluations: (2944) 
#> Path [17] :Initial log joint density = -481912.411725 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      7.903e-03   2.280e-01    1.000e+00  1.000e+00      2975 -3.706e+03 -3.707e+03                   
#> Path [17] :Best Iter: [50] ELBO (-3705.880134) evaluations: (2975) 
#> Path [18] :Initial log joint density = -481648.510742 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.572e-03   1.971e-01    8.187e-01  8.187e-01      3485 -3.703e+03 -3.716e+03                   
#> Path [18] :Best Iter: [55] ELBO (-3703.195035) evaluations: (3485) 
#> Path [19] :Initial log joint density = -481575.899175 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.144e-02   2.817e-01    5.228e-01  1.000e+00      3136 -3.708e+03 -3.709e+03                   
#> Path [19] :Best Iter: [49] ELBO (-3707.995278) evaluations: (3136) 
#> Path [20] :Initial log joint density = -488440.609809 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.361e-02   2.405e-01    1.000e+00  1.000e+00      3238 -3.705e+03 -3.702e+03                   
#> Path [20] :Best Iter: [55] ELBO (-3701.824604) evaluations: (3238) 
#> Path [21] :Initial log joint density = -481760.379316 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.117e-02   2.482e-01    1.000e+00  1.000e+00      3220 -3.698e+03 -3.701e+03                   
#> Path [21] :Best Iter: [57] ELBO (-3698.464849) evaluations: (3220) 
#> Path [22] :Initial log joint density = -481779.188557 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.964e-03   2.349e-01    1.000e+00  1.000e+00      3136 -3.707e+03 -3.711e+03                   
#> Path [22] :Best Iter: [44] ELBO (-3706.525131) evaluations: (3136) 
#> Path [23] :Initial log joint density = -481562.467890 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.381e-02   1.951e-01    1.000e+00  1.000e+00      3324 -3.697e+03 -3.696e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3696.037430) evaluations: (3324) 
#> Path [24] :Initial log joint density = -481529.147109 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.666e-03   2.034e-01    1.000e+00  1.000e+00      3133 -3.704e+03 -3.704e+03                   
#> Path [24] :Best Iter: [56] ELBO (-3703.936553) evaluations: (3133) 
#> Path [25] :Initial log joint density = -481519.113375 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.784e-03   2.309e-01    7.739e-01  7.739e-01      3087 -3.708e+03 -3.721e+03                   
#> Path [25] :Best Iter: [47] ELBO (-3707.692217) evaluations: (3087) 
#> Path [26] :Initial log joint density = -481599.545941 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.470e-03   2.133e-01    1.000e+00  1.000e+00      3099 -3.705e+03 -3.712e+03                   
#> Path [26] :Best Iter: [42] ELBO (-3704.923755) evaluations: (3099) 
#> Path [27] :Initial log joint density = -481596.328389 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.916e-03   2.491e-01    7.395e-01  7.395e-01      3177 -3.704e+03 -3.718e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3704.495205) evaluations: (3177) 
#> Path [28] :Initial log joint density = -485278.241539 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      6.226e-03   1.879e-01    9.938e-01  9.938e-01      3689 -3.697e+03 -3.708e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3696.562339) evaluations: (3689) 
#> Path [29] :Initial log joint density = -481222.312707 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.074e-02   2.057e-01    1.000e+00  1.000e+00      2863 -3.708e+03 -3.710e+03                   
#> Path [29] :Best Iter: [50] ELBO (-3708.409093) evaluations: (2863) 
#> Path [30] :Initial log joint density = -484894.370585 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.801e-03   2.586e-01    1.000e+00  1.000e+00      3399 -3.703e+03 -3.700e+03                   
#> Path [30] :Best Iter: [59] ELBO (-3699.544693) evaluations: (3399) 
#> Path [31] :Initial log joint density = -481558.560173 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.732e-03   2.025e-01    1.000e+00  1.000e+00      2831 -3.707e+03 -3.713e+03                   
#> Path [31] :Best Iter: [48] ELBO (-3706.722480) evaluations: (2831) 
#> Path [32] :Initial log joint density = -481406.383872 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.086e-03   2.194e-01    1.000e+00  1.000e+00      3040 -3.708e+03 -3.713e+03                   
#> Path [32] :Best Iter: [45] ELBO (-3707.803001) evaluations: (3040) 
#> Path [33] :Initial log joint density = -481808.754665 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.496e-02   3.043e-01    1.000e+00  1.000e+00      3276 -3.703e+03 -3.703e+03                   
#> Path [33] :Best Iter: [56] ELBO (-3702.526089) evaluations: (3276) 
#> Path [34] :Initial log joint density = -481703.963497 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.327e-03   1.973e-01    7.020e-01  7.020e-01      3305 -3.702e+03 -3.714e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3701.742085) evaluations: (3305) 
#> Path [35] :Initial log joint density = -481634.139029 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.008e-02   2.906e-01    1.000e+00  1.000e+00      3570 -3.705e+03 -3.704e+03                   
#> Path [35] :Best Iter: [58] ELBO (-3704.425235) evaluations: (3570) 
#> Path [36] :Initial log joint density = -481564.583756 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      1.437e-02   3.077e-01    1.000e+00  1.000e+00      3072 -3.702e+03 -3.703e+03                   
#> Path [36] :Best Iter: [54] ELBO (-3701.742071) evaluations: (3072) 
#> Path [37] :Initial log joint density = -481852.902124 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.245e-03   2.977e-01    5.592e-01  5.592e-01      3078 -3.710e+03 -3.713e+03                   
#> Path [37] :Best Iter: [52] ELBO (-3709.697079) evaluations: (3078) 
#> Path [38] :Initial log joint density = -481492.447858 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      4.268e-03   2.171e-01    9.343e-01  9.343e-01      2790 -3.708e+03 -3.718e+03                   
#> Path [38] :Best Iter: [47] ELBO (-3707.583497) evaluations: (2790) 
#> Path [39] :Initial log joint density = -481917.567362 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.317e-03   2.286e-01    1.000e+00  1.000e+00      3328 -3.699e+03 -3.704e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3699.193197) evaluations: (3328) 
#> Path [40] :Initial log joint density = -482553.963470 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      8.003e-03   2.276e-01    4.433e-01  1.000e+00      3341 -3.697e+03 -3.710e+03                   
#> Path [40] :Best Iter: [56] ELBO (-3696.533466) evaluations: (3341) 
#> Path [41] :Initial log joint density = -482489.595121 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.041e-02   3.228e-01    1.000e+00  1.000e+00      3304 -3.700e+03 -3.708e+03                   
#> Path [41] :Best Iter: [56] ELBO (-3700.185166) evaluations: (3304) 
#> Path [42] :Initial log joint density = -481763.763155 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      9.933e-03   2.419e-01    1.000e+00  1.000e+00      2991 -3.707e+03 -3.708e+03                   
#> Path [42] :Best Iter: [47] ELBO (-3706.675045) evaluations: (2991) 
#> Path [43] :Initial log joint density = -482063.277097 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      5.944e-03   2.105e-01    8.308e-01  8.308e-01      3412 -3.697e+03 -3.710e+03                   
#> Path [43] :Best Iter: [57] ELBO (-3697.117810) evaluations: (3412) 
#> Path [44] :Initial log joint density = -481451.534993 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.365e-02   3.555e-01    1.000e+00  1.000e+00      2827 -3.705e+03 -3.717e+03                   
#> Path [44] :Best Iter: [47] ELBO (-3705.129501) evaluations: (2827) 
#> Path [45] :Initial log joint density = -481656.916258 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      2.921e-03   1.619e-01    7.108e-01  7.108e-01      3485 -3.698e+03 -3.714e+03                   
#> Path [45] :Best Iter: [58] ELBO (-3698.015919) evaluations: (3485) 
#> Path [46] :Initial log joint density = -482455.480268 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      9.557e-03   3.567e-01    5.109e-01  1.000e+00      3417 -3.702e+03 -3.711e+03                   
#> Path [46] :Best Iter: [58] ELBO (-3701.880309) evaluations: (3417) 
#> Path [47] :Initial log joint density = -481688.401390 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.087e-03   1.733e-01    1.000e+00  1.000e+00      2991 -3.707e+03 -3.714e+03                   
#> Path [47] :Best Iter: [41] ELBO (-3707.342228) evaluations: (2991) 
#> Path [48] :Initial log joint density = -481710.207369 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.594e-03   2.194e-01    9.536e-01  9.536e-01      2825 -3.710e+03 -3.719e+03                   
#> Path [48] :Best Iter: [39] ELBO (-3710.199231) evaluations: (2825) 
#> Path [49] :Initial log joint density = -481513.204411 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      9.122e-03   2.399e-01    8.344e-01  8.344e-01      3653 -3.698e+03 -3.707e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3698.180400) evaluations: (3653) 
#> Path [50] :Initial log joint density = -481940.265448 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.061e-02   1.965e-01    1.000e+00  1.000e+00      3219 -3.702e+03 -3.703e+03                   
#> Path [50] :Best Iter: [56] ELBO (-3701.575758) evaluations: (3219) 
#> Finished in  13.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
# }

```
