# Plot 2D Intervals for Mean-Variance Association

This function creates a 2D interval plot for mean-variance association,
handling both single and bimodal models. It highlights significant
differences based on a given significance threshold.

## Usage

``` r
sccomp_plot_intervals_2D(
  .data,
  factor = NULL,
  significance_threshold = 0.05,
  test_composition_above_CLR_fold_change = attr(.data,
    "test_composition_above_CLR_fold_change"),
  show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR"),
  add_marginal_density = TRUE,
  omit_ci = FALSE
)
```

## Arguments

- .data:

  Data frame containing the main data.

- factor:

  Optional character string selecting one model factor to plot. If
  provided, plots are restricted to that factor plus `(Intercept)`.

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences. Default is 0.05.

- test_composition_above_CLR_fold_change:

  A positive integer. It is the effect threshold used for the hypothesis
  test.

- show_fdr_message:

  Logical. Whether to show the Bayesian FDR interpretation message on
  the plot. Default is TRUE.

- significance_statistic:

  Character vector indicating which statistic to highlight. Default is
  "pH0".

- add_marginal_density:

  Logical. Whether to add marginal density plots on adjusted panels.
  Default is TRUE.

- omit_ci:

  Logical. Whether to omit credible interval error bars. Default is
  FALSE.

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
      cores = 1,
      bimodal_mean_variability_association = TRUE
    ) |>
    sccomp_test()

    # Example usage:
    my_plot = sccomp_plot_intervals_2D(estimate)

  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 2 parameters: (Intercept), typecancer
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481802.285461 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.065e-02   3.122e+02    8.401e-02  1.539e-01      8157 -3.795e+03 -4.435e+03                   
#> Path [1] :Best Iter: [36] ELBO (-3795.024614) evaluations: (8157) 
#> Path [2] :Initial log joint density = -481921.690470 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.591e-02   7.056e+02    5.114e-02  5.114e-02      8167 -3.784e+03 -5.460e+03                   
#> Path [2] :Best Iter: [38] ELBO (-3784.099396) evaluations: (8167) 
#> Path [3] :Initial log joint density = -482431.268002 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.531e-01   4.678e+02    1.818e-01  1.818e-01      8203 -3.777e+03 -4.270e+03                   
#> Path [3] :Best Iter: [45] ELBO (-3777.377556) evaluations: (8203) 
#> Path [4] :Initial log joint density = -481624.554643 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.496e-01   3.411e+02    3.812e-01  3.812e-01      8306 -3.828e+03 -4.160e+03                   
#> Path [4] :Best Iter: [37] ELBO (-3828.430713) evaluations: (8306) 
#> Path [5] :Initial log joint density = -482634.843527 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.851e-02   8.633e+02    3.438e-01  3.438e-01      8289 -3.814e+03 -4.804e+03                   
#> Path [5] :Best Iter: [36] ELBO (-3814.313816) evaluations: (8289) 
#> Path [6] :Initial log joint density = -483276.905500 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.097e-03   3.399e+02    5.460e-02  5.460e-02      8345 -3.780e+03 -4.399e+03                   
#> Path [6] :Best Iter: [41] ELBO (-3779.891027) evaluations: (8345) 
#> Path [7] :Initial log joint density = -482380.861191 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.544e-02   4.231e+02    5.671e-02  5.671e-02      8128 -3.781e+03 -4.727e+03                   
#> Path [7] :Best Iter: [40] ELBO (-3780.880108) evaluations: (8128) 
#> Path [8] :Initial log joint density = -481905.108558 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.960e-02   3.108e+02    5.114e-01  5.114e-01      8145 -3.800e+03 -4.170e+03                   
#> Path [8] :Best Iter: [40] ELBO (-3800.247246) evaluations: (8145) 
#> Path [9] :Initial log joint density = -483254.097045 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.926e-02   7.481e+02    4.029e-02  4.029e-02      8164 -3.823e+03 -4.483e+03                   
#> Path [9] :Best Iter: [38] ELBO (-3822.771652) evaluations: (8164) 
#> Path [10] :Initial log joint density = -481564.993439 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.756e-02   9.118e+02    2.343e-01  2.343e-01      8480 -3.805e+03 -4.312e+03                   
#> Path [10] :Best Iter: [35] ELBO (-3804.608272) evaluations: (8480) 
#> Path [11] :Initial log joint density = -481502.683304 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.916e-03   1.606e+02    1.216e-01  1.216e-01      8207 -3.807e+03 -4.244e+03                   
#> Path [11] :Best Iter: [39] ELBO (-3806.943346) evaluations: (8207) 
#> Path [12] :Initial log joint density = -481714.260363 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.010e-02   8.714e+02    2.870e-01  2.870e-01      8363 -3.799e+03 -4.172e+03                   
#> Path [12] :Best Iter: [34] ELBO (-3798.634398) evaluations: (8363) 
#> Path [13] :Initial log joint density = -481692.667945 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.758e-02   4.200e+02    1.240e-01  1.240e-01      8159 -3.842e+03 -4.221e+03                   
#> Path [13] :Best Iter: [35] ELBO (-3841.640459) evaluations: (8159) 
#> Path [14] :Initial log joint density = -482580.497482 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      6.543e-02   6.276e+02    1.197e-01  1.197e-01      8067 -3.794e+03 -4.346e+03                   
#> Path [14] :Best Iter: [45] ELBO (-3793.735927) evaluations: (8067) 
#> Path [15] :Initial log joint density = -481613.087897 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.687e-02   4.356e+02    2.848e-01  2.848e-01      8366 -3.802e+03 -4.218e+03                   
#> Path [15] :Best Iter: [39] ELBO (-3801.537748) evaluations: (8366) 
#> Path [16] :Initial log joint density = -481730.235019 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.970e-02   9.043e+02    2.648e-01  2.648e-01      8488 -3.845e+03 -4.396e+03                   
#> Path [16] :Best Iter: [34] ELBO (-3844.624022) evaluations: (8488) 
#> Path [17] :Initial log joint density = -482116.909150 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.988e-02   3.335e+02    1.907e-01  1.907e-01      8362 -3.781e+03 -4.214e+03                   
#> Path [17] :Best Iter: [41] ELBO (-3781.294393) evaluations: (8362) 
#> Path [18] :Initial log joint density = -482657.041721 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.595e-02   3.834e+02    2.374e-01  2.374e-01      8356 -3.781e+03 -4.271e+03                   
#> Path [18] :Best Iter: [39] ELBO (-3781.330345) evaluations: (8356) 
#> Path [19] :Initial log joint density = -481686.881837 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.154e-02   5.897e+02    5.512e-01  5.512e-01      8051 -3.802e+03 -4.127e+03                   
#> Path [19] :Best Iter: [35] ELBO (-3801.719285) evaluations: (8051) 
#> Path [20] :Initial log joint density = -482513.828711 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.615e-02   2.851e+02    3.490e-01  3.490e-01      8351 -3.819e+03 -4.464e+03                   
#> Path [20] :Best Iter: [37] ELBO (-3818.799557) evaluations: (8351) 
#> Path [21] :Initial log joint density = -482360.448817 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.426e-03   3.374e+02    4.684e-02  4.684e-02      7980 -3.808e+03 -4.423e+03                   
#> Path [21] :Best Iter: [42] ELBO (-3808.398990) evaluations: (7980) 
#> Path [22] :Initial log joint density = -481645.631479 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.864e-02   5.374e+02    3.793e-01  3.793e-01      8281 -3.793e+03 -4.133e+03                   
#> Path [22] :Best Iter: [39] ELBO (-3792.943402) evaluations: (8281) 
#> Path [23] :Initial log joint density = -482804.188928 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.311e-02   5.108e+02    1.191e-01  1.191e-01      8178 -3.821e+03 -4.607e+03                   
#> Path [23] :Best Iter: [37] ELBO (-3820.746522) evaluations: (8178) 
#> Path [24] :Initial log joint density = -483254.688356 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.414e-02   6.374e+02    8.531e-02  8.531e-02      8165 -3.830e+03 -4.489e+03                   
#> Path [24] :Best Iter: [39] ELBO (-3829.794740) evaluations: (8165) 
#> Path [25] :Initial log joint density = -481443.518480 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.644e-02   1.116e+03    1.459e-01  1.459e-01      8423 -3.821e+03 -4.338e+03                   
#> Path [25] :Best Iter: [34] ELBO (-3820.884003) evaluations: (8423) 
#> Path [26] :Initial log joint density = -481655.168628 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.099e-01   4.133e+02    5.765e-01  5.765e-01      8383 -3.808e+03 -4.140e+03                   
#> Path [26] :Best Iter: [33] ELBO (-3807.996130) evaluations: (8383) 
#> Path [27] :Initial log joint density = -481644.087945 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      9.425e-03   3.075e+02    5.794e-02  5.794e-02      8158 -3.778e+03 -4.235e+03                   
#> Path [27] :Best Iter: [40] ELBO (-3777.554811) evaluations: (8158) 
#> Path [28] :Initial log joint density = -481914.589058 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.194e-01   4.418e+02    2.673e-01  2.673e-01      8268 -3.850e+03 -4.276e+03                   
#> Path [28] :Best Iter: [35] ELBO (-3849.603199) evaluations: (8268) 
#> Path [29] :Initial log joint density = -484445.630584 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.614e-02   3.711e+02    2.390e-01  2.390e-01      8403 -3.884e+03 -4.212e+03                   
#> Path [29] :Best Iter: [35] ELBO (-3884.403528) evaluations: (8403) 
#> Path [30] :Initial log joint density = -481642.848983 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.045e-02   2.366e+02    2.303e-01  2.303e-01      8380 -3.821e+03 -4.237e+03                   
#> Path [30] :Best Iter: [36] ELBO (-3821.336888) evaluations: (8380) 
#> Path [31] :Initial log joint density = -481663.587381 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.411e-02   2.210e+02    4.907e-01  4.907e-01      8145 -3.787e+03 -4.234e+03                   
#> Path [31] :Best Iter: [38] ELBO (-3786.513542) evaluations: (8145) 
#> Path [32] :Initial log joint density = -481861.467824 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.555e-02   5.999e+02    2.579e-01  2.579e-01      8355 -3.794e+03 -4.213e+03                   
#> Path [32] :Best Iter: [38] ELBO (-3794.037679) evaluations: (8355) 
#> Path [33] :Initial log joint density = -481638.583494 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.630e-02   1.005e+03    1.092e-01  1.092e-01      8217 -3.819e+03 -4.670e+03                   
#> Path [33] :Best Iter: [39] ELBO (-3819.378907) evaluations: (8217) 
#> Path [34] :Initial log joint density = -482643.104091 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.465e-02   3.339e+02    6.375e-02  1.207e-01      8217 -3.808e+03 -4.401e+03                   
#> Path [34] :Best Iter: [46] ELBO (-3808.401381) evaluations: (8217) 
#> Path [35] :Initial log joint density = -485006.333366 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.977e-02   2.995e+02    5.720e-02  5.720e-02      8184 -3.805e+03 -4.611e+03                   
#> Path [35] :Best Iter: [42] ELBO (-3805.159555) evaluations: (8184) 
#> Path [36] :Initial log joint density = -482494.454084 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.752e-02   4.919e+02    1.035e-01  1.035e-01      7991 -3.784e+03 -4.531e+03                   
#> Path [36] :Best Iter: [46] ELBO (-3784.240172) evaluations: (7991) 
#> Path [37] :Initial log joint density = -482702.748577 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.673e-02   2.950e+02    4.033e-01  4.033e-01      8596 -3.814e+03 -4.244e+03                   
#> Path [37] :Best Iter: [40] ELBO (-3814.177120) evaluations: (8596) 
#> Path [38] :Initial log joint density = -482083.571987 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.944e-02   9.525e+02    8.318e-02  8.318e-02      8051 -3.804e+03 -4.477e+03                   
#> Path [38] :Best Iter: [35] ELBO (-3803.560014) evaluations: (8051) 
#> Path [39] :Initial log joint density = -482125.512457 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.423e-02   3.575e+02    1.268e-01  1.268e-01      8572 -3.789e+03 -4.457e+03                   
#> Path [39] :Best Iter: [38] ELBO (-3789.390741) evaluations: (8572) 
#> Path [40] :Initial log joint density = -481473.390677 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.481e-02   5.622e+02    2.370e-01  2.370e-01      8009 -3.798e+03 -4.255e+03                   
#> Path [40] :Best Iter: [39] ELBO (-3798.267762) evaluations: (8009) 
#> Path [41] :Initial log joint density = -485576.322909 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.763e-02   1.204e+03    1.264e-01  1.264e-01      8273 -3.794e+03 -4.283e+03                   
#> Path [41] :Best Iter: [44] ELBO (-3794.446283) evaluations: (8273) 
#> Path [42] :Initial log joint density = -481472.083452 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.965e-02   3.787e+02    2.185e-01  2.185e-01      8262 -3.805e+03 -4.166e+03                   
#> Path [42] :Best Iter: [31] ELBO (-3804.745538) evaluations: (8262) 
#> Path [43] :Initial log joint density = -481487.268826 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.116e-02   2.251e+02    3.982e-01  3.982e-01      8281 -3.793e+03 -4.144e+03                   
#> Path [43] :Best Iter: [35] ELBO (-3793.098793) evaluations: (8281) 
#> Path [44] :Initial log joint density = -482162.819395 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.983e-02   3.826e+02    1.857e-01  1.857e-01      8123 -3.819e+03 -4.247e+03                   
#> Path [44] :Best Iter: [38] ELBO (-3818.622940) evaluations: (8123) 
#> Path [45] :Initial log joint density = -481802.104086 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.063e-02   6.435e+02    1.008e-01  1.008e-01      8054 -3.792e+03 -4.419e+03                   
#> Path [45] :Best Iter: [38] ELBO (-3791.981126) evaluations: (8054) 
#> Path [46] :Initial log joint density = -481799.292419 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.850e-02   1.517e+02    2.669e-01  2.669e-01      8334 -3.832e+03 -4.118e+03                   
#> Path [46] :Best Iter: [39] ELBO (-3832.124370) evaluations: (8334) 
#> Path [47] :Initial log joint density = -481814.928485 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.861e-02   3.969e+02    3.178e-01  3.178e-01      8362 -3.796e+03 -4.263e+03                   
#> Path [47] :Best Iter: [44] ELBO (-3796.013110) evaluations: (8362) 
#> Path [48] :Initial log joint density = -482089.827674 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.630e-02   3.115e+02    1.585e-01  1.585e-01      8255 -3.783e+03 -4.161e+03                   
#> Path [48] :Best Iter: [42] ELBO (-3782.742540) evaluations: (8255) 
#> Path [49] :Initial log joint density = -481477.334223 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.320e-02   3.561e+02    2.097e-01  2.097e-01      8089 -3.798e+03 -4.228e+03                   
#> Path [49] :Best Iter: [37] ELBO (-3798.123372) evaluations: (8089) 
#> Path [50] :Initial log joint density = -481797.394449 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.875e-02   8.855e+02    5.196e-01  5.196e-01      8085 -3.814e+03 -4.137e+03                   
#> Path [50] :Best Iter: [35] ELBO (-3813.968967) evaluations: (8085) 
#> Finished in  24.1 seconds.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Joining with `by = join_by(cell_group, M, parameter)`
#> === Bimodal Model Parameters ===
#> 
#> (Intercept):
#>   Component 1: v = -(4.195 + -0.308 × c)
#>   Component 2: v = -(5.659 + -0.576 × c)
#> 
#> typecancer:
#>   Component 1: v = -(-0.197 + 0.584 × c)
#>   Component 2: v = -(0.301 + 0.078 × c)
#> 
# }
```
