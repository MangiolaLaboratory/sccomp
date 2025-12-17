# Plot 2D Intervals for Mean-Variance Association

This function creates a 2D interval plot for mean-variance association,
highlighting significant differences based on a given significance
threshold.

## Usage

``` r
plot_2D_intervals(
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
  differences. Default is 0.025.

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

A ggplot object representing the 2D interval plot.

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
      ~type,
      "sample",
      "cell_group",
      "count",
      cores = 1
    ) |> 
    sccomp_test()
    
  # Example usage:
  my_plot = plot_2D_intervals(estimate)
    
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept), typecancer
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481597.479556 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.132e-02   5.749e+03    1.575e-02  1.575e-02      8435 -3.787e+03 -5.960e+05                   
#> Path [1] :Best Iter: [37] ELBO (-3786.811524) evaluations: (8435) 
#> Path [2] :Initial log joint density = -482946.252425 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.395e-02   6.545e+03    1.531e-02  1.531e-02      8207 -3.799e+03 -2.266e+06                   
#> Path [2] :Best Iter: [38] ELBO (-3799.252101) evaluations: (8207) 
#> Path [3] :Initial log joint density = -483899.219127 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.158e-01   5.257e+03    2.818e-02  2.818e-02      8322 -3.792e+03 -1.912e+11                   
#> Path [3] :Best Iter: [45] ELBO (-3791.835480) evaluations: (8322) 
#> Path [4] :Initial log joint density = -482593.869009 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.617e-02   2.054e+03    2.425e-02  2.425e-02      8156 -3.785e+03 -8.857e+05                   
#> Path [4] :Best Iter: [40] ELBO (-3785.117901) evaluations: (8156) 
#> Path [5] :Initial log joint density = -481655.705826 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.282e-02   1.112e+03    2.317e-02  2.317e-02      8081 -3.771e+03 -9.088e+03                   
#> Path [5] :Best Iter: [41] ELBO (-3770.830500) evaluations: (8081) 
#> Path [6] :Initial log joint density = -481553.763429 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.083e-02   1.987e+03    2.706e-02  2.706e-02      8276 -3.776e+03 -1.170e+04                   
#> Path [6] :Best Iter: [45] ELBO (-3775.832706) evaluations: (8276) 
#> Path [7] :Initial log joint density = -481417.084915 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.930e-02   2.243e+04    7.729e-03  7.729e-03      8335 -3.810e+03 -6.521e+04                   
#> Path [7] :Best Iter: [36] ELBO (-3810.032070) evaluations: (8335) 
#> Path [8] :Initial log joint density = -481684.603159 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.051e-02   9.150e+03    9.812e-03  9.812e-03      8247 -3.803e+03 -6.060e+11                   
#> Path [8] :Best Iter: [38] ELBO (-3802.759701) evaluations: (8247) 
#> Path [9] :Initial log joint density = -488330.727436 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.710e-02   1.240e+03    2.844e-02  2.844e-02      8214 -3.771e+03 -9.012e+03                   
#> Path [9] :Best Iter: [51] ELBO (-3771.317827) evaluations: (8214) 
#> Path [10] :Initial log joint density = -481868.537444 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.464e-02   3.402e+02    6.458e-02  6.458e-02      8050 -3.764e+03 -4.420e+03                   
#> Path [10] :Best Iter: [51] ELBO (-3763.545131) evaluations: (8050) 
#> Path [11] :Initial log joint density = -481883.239560 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.240e-02   6.046e+03    2.258e-02  2.258e-02      8168 -3.796e+03 -2.427e+04                   
#> Path [11] :Best Iter: [42] ELBO (-3795.967872) evaluations: (8168) 
#> Path [12] :Initial log joint density = -481865.584664 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.268e-02   5.477e+03    2.916e-02  2.916e-02      8550 -3.784e+03 -2.823e+05                   
#> Path [12] :Best Iter: [46] ELBO (-3784.078547) evaluations: (8550) 
#> Path [13] :Initial log joint density = -481709.849506 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.308e+00   1.008e+02    1.675e-01  1.675e-01      7865 -3.732e+03 -4.563e+03                   
#> Path [13] :Best Iter: [78] ELBO (-3732.353027) evaluations: (7865) 
#> Path [14] :Initial log joint density = -485549.132760 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.679e-01   1.970e+03    5.214e-02  5.214e-02      8050 -3.775e+03 -4.307e+05                   
#> Path [14] :Best Iter: [53] ELBO (-3774.830027) evaluations: (8050) 
#> Path [15] :Initial log joint density = -481973.542059 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.535e-02   6.299e+03    2.438e-02  2.438e-02      8554 -3.817e+03 -7.157e+03                   
#> Path [15] :Best Iter: [30] ELBO (-3817.047548) evaluations: (8554) 
#> Path [16] :Initial log joint density = -482213.169603 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.035e-02   1.662e-01    4.822e-01  1.000e+00      6235 -3.733e+03 -3.755e+03                   
#> Path [16] :Best Iter: [76] ELBO (-3733.234305) evaluations: (6235) 
#> Path [17] :Initial log joint density = -482057.782612 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      6.862e-03   2.391e-01    6.575e-01  6.575e-01      5834 -3.734e+03 -3.755e+03                   
#> Path [17] :Best Iter: [82] ELBO (-3733.677683) evaluations: (5834) 
#> Path [18] :Initial log joint density = -481553.704194 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      9.546e-02   6.227e+03    2.784e-02  2.784e-02      8441 -3.795e+03 -1.135e+05                   
#> Path [18] :Best Iter: [40] ELBO (-3794.829968) evaluations: (8441) 
#> Path [19] :Initial log joint density = -481906.935406 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      6.629e-02   4.645e+03    2.511e-02  6.291e-02      8387 -3.796e+03 -2.440e+08                   
#> Path [19] :Best Iter: [39] ELBO (-3796.474297) evaluations: (8387) 
#> Path [20] :Initial log joint density = -481488.509206 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.293e-02   4.786e+03    8.934e-03  2.608e-02      8123 -3.797e+03 -3.137e+07                   
#> Path [20] :Best Iter: [40] ELBO (-3797.209252) evaluations: (8123) 
#> Path [21] :Initial log joint density = -481624.983410 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.186e-01   1.967e+03    3.947e-02  3.947e-02      7905 -3.780e+03 -2.218e+04                   
#> Path [21] :Best Iter: [41] ELBO (-3779.999510) evaluations: (7905) 
#> Path [22] :Initial log joint density = -482032.303308 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.074e-02   1.006e+03    1.297e-02  2.351e-02      8204 -3.767e+03 -6.404e+04                   
#> Path [22] :Best Iter: [55] ELBO (-3767.222243) evaluations: (8204) 
#> Path [23] :Initial log joint density = -481558.098722 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      8.492e-02   4.713e+03    2.266e-02  2.266e-02      8251 -3.808e+03 -8.369e+04                   
#> Path [23] :Best Iter: [32] ELBO (-3807.984027) evaluations: (8251) 
#> Path [24] :Initial log joint density = -482061.158968 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.551e-01   1.207e+02    1.179e-01  2.657e-01      7992 -3.740e+03 -1.895e+04                   
#> Path [24] :Best Iter: [77] ELBO (-3740.388748) evaluations: (7992) 
#> Path [25] :Initial log joint density = -481721.246228 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.763e-02   2.225e+03    3.240e-02  3.240e-02      8182 -3.769e+03 -2.435e+05                   
#> Path [25] :Best Iter: [51] ELBO (-3768.995390) evaluations: (8182) 
#> Path [26] :Initial log joint density = -482078.266274 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.573e-02   1.915e-01    1.000e+00  1.000e+00      6170 -3.732e+03 -3.738e+03                   
#> Path [26] :Best Iter: [86] ELBO (-3731.542720) evaluations: (6170) 
#> Path [27] :Initial log joint density = -485063.451206 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.253e-01   6.083e+02    4.971e-02  4.971e-02      8053 -3.772e+03 -3.703e+06                   
#> Path [27] :Best Iter: [57] ELBO (-3771.915623) evaluations: (8053) 
#> Path [28] :Initial log joint density = -482378.323345 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.282e-01   1.420e+02    1.562e-01  1.562e-01      8073 -3.737e+03 -5.615e+03                   
#> Path [28] :Best Iter: [76] ELBO (-3737.110287) evaluations: (8073) 
#> Path [29] :Initial log joint density = -481355.569507 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.208e-02   1.997e+03    2.441e-02  2.441e-02      8385 -3.779e+03 -6.926e+05                   
#> Path [29] :Best Iter: [38] ELBO (-3779.244262) evaluations: (8385) 
#> Path [30] :Initial log joint density = -481910.612965 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.129e-01   2.641e+03    4.024e-02  4.024e-02      8200 -3.774e+03 -9.735e+04                   
#> Path [30] :Best Iter: [43] ELBO (-3774.285315) evaluations: (8200) 
#> Path [31] :Initial log joint density = -481735.284109 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.772e-02   1.151e+03    2.034e-02  2.034e-02      8137 -3.762e+03 -3.825e+07                   
#> Path [31] :Best Iter: [57] ELBO (-3762.484244) evaluations: (8137) 
#> Path [32] :Initial log joint density = -481827.798675 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      1.656e-02   2.102e-01    1.000e+00  1.000e+00      5428 -3.734e+03 -3.745e+03                   
#> Path [32] :Best Iter: [78] ELBO (-3733.885196) evaluations: (5428) 
#> Path [33] :Initial log joint density = -481784.738248 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.070e-02   3.122e+03    2.022e-02  3.930e-02      8351 -3.809e+03 -4.146e+04                   
#> Path [33] :Best Iter: [38] ELBO (-3808.859075) evaluations: (8351) 
#> Path [34] :Initial log joint density = -481697.446474 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      6.917e-02   3.750e+03    1.384e-02  2.650e-02      8493 -3.804e+03 -5.760e+09                   
#> Path [34] :Best Iter: [36] ELBO (-3803.743817) evaluations: (8493) 
#> Path [35] :Initial log joint density = -481776.510109 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      9.183e-02   2.710e+03    2.128e-02  6.484e-02      8237 -3.794e+03 -1.437e+05                   
#> Path [35] :Best Iter: [36] ELBO (-3793.873661) evaluations: (8237) 
#> Path [36] :Initial log joint density = -482011.186690 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.327e-02   3.419e+03    2.188e-02  2.188e-02      8149 -3.778e+03 -2.501e+05                   
#> Path [36] :Best Iter: [41] ELBO (-3778.071386) evaluations: (8149) 
#> Path [37] :Initial log joint density = -481873.669336 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.981e-02   3.356e+03    1.559e-02  1.559e-02      8652 -3.790e+03 -2.709e+05                   
#> Path [37] :Best Iter: [40] ELBO (-3790.032543) evaluations: (8652) 
#> Path [38] :Initial log joint density = -481759.165078 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      9.865e-02   1.462e+04    4.303e-02  4.303e-02      8529 -3.812e+03 -4.165e+04                   
#> Path [38] :Best Iter: [32] ELBO (-3811.587329) evaluations: (8529) 
#> Path [39] :Initial log joint density = -481694.096055 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.169e-01   1.862e+03    1.891e-02  3.407e-02      8360 -3.767e+03 -6.284e+05                   
#> Path [39] :Best Iter: [56] ELBO (-3767.185397) evaluations: (8360) 
#> Path [40] :Initial log joint density = -482693.319618 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.072e-01   5.129e+02    4.791e-02  4.791e-02      8084 -3.762e+03 -5.019e+04                   
#> Path [40] :Best Iter: [46] ELBO (-3761.734869) evaluations: (8084) 
#> Path [41] :Initial log joint density = -482230.240788 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      6.682e-03   1.687e-01    9.290e-01  9.290e-01      6094 -3.730e+03 -3.755e+03                   
#> Path [41] :Best Iter: [79] ELBO (-3729.957099) evaluations: (6094) 
#> Path [42] :Initial log joint density = -481619.358493 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.578e-01   1.001e+03    4.799e-02  4.799e-02      8294 -3.758e+03 -4.054e+04                   
#> Path [42] :Best Iter: [59] ELBO (-3758.317041) evaluations: (8294) 
#> Path [43] :Initial log joint density = -481864.066810 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      8.639e-02   5.176e+03    2.502e-02  2.502e-02      8288 -3.811e+03 -2.912e+04                   
#> Path [43] :Best Iter: [37] ELBO (-3811.168762) evaluations: (8288) 
#> Path [44] :Initial log joint density = -481423.973525 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.375e-02   2.603e+03    9.504e-03  9.504e-03      8263 -3.778e+03 -4.503e+10                   
#> Path [44] :Best Iter: [45] ELBO (-3777.726776) evaluations: (8263) 
#> Path [45] :Initial log joint density = -481987.071785 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.787e+05      7.410e-03   2.417e-01    1.000e+00  1.000e+00      6582 -3.733e+03 -3.747e+03                   
#> Path [45] :Best Iter: [85] ELBO (-3732.978050) evaluations: (6582) 
#> Path [46] :Initial log joint density = -482642.076068 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.986e-01   6.530e+02    8.160e-02  8.160e-02      7968 -3.766e+03 -7.866e+04                   
#> Path [46] :Best Iter: [59] ELBO (-3766.393718) evaluations: (7968) 
#> Path [47] :Initial log joint density = -481565.382487 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      8.047e-03   1.801e-01    7.973e-01  7.973e-01      5867 -3.732e+03 -3.759e+03                   
#> Path [47] :Best Iter: [80] ELBO (-3732.392389) evaluations: (5867) 
#> Path [48] :Initial log joint density = -483956.474593 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      9.448e-02   5.849e+03    3.183e-02  6.278e-02      8415 -3.799e+03 -6.484e+04                   
#> Path [48] :Best Iter: [40] ELBO (-3798.507927) evaluations: (8415) 
#> Path [49] :Initial log joint density = -481595.111520 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.566e-02   2.950e+03    2.219e-02  2.219e-02      8181 -3.779e+03 -5.741e+05                   
#> Path [49] :Best Iter: [52] ELBO (-3779.159335) evaluations: (8181) 
#> Path [50] :Initial log joint density = -481843.794082 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.295e-02   7.275e+03    1.967e-02  1.967e-02      8497 -3.823e+03 -6.512e+04                   
#> Path [50] :Best Iter: [32] ELBO (-3822.691190) evaluations: (8497) 
#> Finished in  25.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
