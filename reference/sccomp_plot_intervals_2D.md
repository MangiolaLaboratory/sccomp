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
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
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

- test_composition_above_logit_fold_change:

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
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept), typecancer
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481818.013546 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.630e-02   5.988e+02    1.239e-01  1.239e-01      8279 -3.778e+03 -4.246e+03                   
#> Path [1] :Best Iter: [36] ELBO (-3778.075919) evaluations: (8279) 
#> Path [2] :Initial log joint density = -481939.973738 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.490e-02   5.327e+02    8.358e-02  8.358e-02      8098 -3.775e+03 -4.326e+03                   
#> Path [2] :Best Iter: [40] ELBO (-3775.190444) evaluations: (8098) 
#> Path [3] :Initial log joint density = -482444.581094 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.329e-02   4.346e+02    1.922e-01  1.922e-01      8202 -3.771e+03 -4.282e+03                   
#> Path [3] :Best Iter: [44] ELBO (-3771.441396) evaluations: (8202) 
#> Path [4] :Initial log joint density = -481639.803970 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.029e-02   2.782e+02    3.389e-01  3.389e-01      8489 -3.812e+03 -4.120e+03                   
#> Path [4] :Best Iter: [37] ELBO (-3811.972871) evaluations: (8489) 
#> Path [5] :Initial log joint density = -482647.948737 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.774e-02   8.556e+02    1.427e-01  1.427e-01      8476 -3.774e+03 -4.808e+03                   
#> Path [5] :Best Iter: [34] ELBO (-3773.520804) evaluations: (8476) 
#> Path [6] :Initial log joint density = -483294.104786 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.029e-02   6.509e+02    5.892e-02  5.892e-02      8351 -3.766e+03 -4.947e+03                   
#> Path [6] :Best Iter: [41] ELBO (-3765.559766) evaluations: (8351) 
#> Path [7] :Initial log joint density = -482398.014342 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.742e-02   6.539e+02    8.226e-02  8.226e-02      8124 -3.764e+03 -4.674e+03                   
#> Path [7] :Best Iter: [47] ELBO (-3764.482753) evaluations: (8124) 
#> Path [8] :Initial log joint density = -481917.076661 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.483e-02   2.143e+02    1.299e-01  1.299e-01      8219 -3.786e+03 -4.168e+03                   
#> Path [8] :Best Iter: [40] ELBO (-3786.272584) evaluations: (8219) 
#> Path [9] :Initial log joint density = -483270.503550 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.145e-02   3.115e+02    1.277e-01  1.277e-01      8144 -3.809e+03 -4.231e+03                   
#> Path [9] :Best Iter: [34] ELBO (-3809.420248) evaluations: (8144) 
#> Path [10] :Initial log joint density = -481578.476606 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.905e-02   6.666e+02    3.398e-01  3.398e-01      8119 -3.796e+03 -4.226e+03                   
#> Path [10] :Best Iter: [33] ELBO (-3796.175651) evaluations: (8119) 
#> Path [11] :Initial log joint density = -481518.971309 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.757e-02   2.380e+02    2.801e-01  2.801e-01      8357 -3.801e+03 -4.237e+03                   
#> Path [11] :Best Iter: [39] ELBO (-3801.188896) evaluations: (8357) 
#> Path [12] :Initial log joint density = -481728.388710 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.544e-02   9.176e+02    6.787e-02  6.787e-02      8285 -3.781e+03 -5.568e+03                   
#> Path [12] :Best Iter: [36] ELBO (-3781.257153) evaluations: (8285) 
#> Path [13] :Initial log joint density = -481706.540751 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.144e-01   4.982e+02    2.352e-01  2.352e-01      8196 -3.830e+03 -4.256e+03                   
#> Path [13] :Best Iter: [37] ELBO (-3829.580276) evaluations: (8196) 
#> Path [14] :Initial log joint density = -482597.775423 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.642e-02   6.004e+02    8.257e-02  8.257e-02      7700 -3.766e+03 -5.906e+03                   
#> Path [14] :Best Iter: [45] ELBO (-3766.173879) evaluations: (7700) 
#> Path [15] :Initial log joint density = -481625.732356 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.782e-02   1.885e+02    2.566e-01  2.566e-01      7902 -3.793e+03 -4.170e+03                   
#> Path [15] :Best Iter: [42] ELBO (-3793.280544) evaluations: (7902) 
#> Path [16] :Initial log joint density = -481745.325128 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.412e-02   9.491e+02    2.105e-01  2.105e-01      8286 -3.819e+03 -4.399e+03                   
#> Path [16] :Best Iter: [36] ELBO (-3819.457599) evaluations: (8286) 
#> Path [17] :Initial log joint density = -482133.452623 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.155e-01   6.306e+02    1.767e-01  1.767e-01      8302 -3.768e+03 -6.569e+03                   
#> Path [17] :Best Iter: [42] ELBO (-3768.003106) evaluations: (8302) 
#> Path [18] :Initial log joint density = -482668.322440 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.113e-02   3.056e+02    1.846e-01  1.846e-01      8140 -3.762e+03 -4.200e+03                   
#> Path [18] :Best Iter: [41] ELBO (-3762.433091) evaluations: (8140) 
#> Path [19] :Initial log joint density = -481698.903654 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.952e-02   5.296e+02    2.270e-01  2.270e-01      8115 -3.795e+03 -4.174e+03                   
#> Path [19] :Best Iter: [40] ELBO (-3794.548850) evaluations: (8115) 
#> Path [20] :Initial log joint density = -482527.785006 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.211e-02   3.472e+02    3.192e-01  3.192e-01      8320 -3.798e+03 -4.612e+03                   
#> Path [20] :Best Iter: [36] ELBO (-3798.000533) evaluations: (8320) 
#> Path [21] :Initial log joint density = -482378.269938 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.504e-02   3.851e+02    2.702e-02  1.028e-01      8225 -3.782e+03 -5.168e+03                   
#> Path [21] :Best Iter: [41] ELBO (-3782.185217) evaluations: (8225) 
#> Path [22] :Initial log joint density = -481660.819543 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.898e-02   8.721e+02    1.858e-01  3.479e-01      8250 -3.761e+03 -4.205e+03                   
#> Path [22] :Best Iter: [36] ELBO (-3761.456567) evaluations: (8250) 
#> Path [23] :Initial log joint density = -482818.133978 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.762e-02   4.760e+02    1.725e-01  3.217e-01      8178 -3.788e+03 -4.332e+03                   
#> Path [23] :Best Iter: [42] ELBO (-3788.063188) evaluations: (8178) 
#> Path [24] :Initial log joint density = -483272.582104 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.179e-01   4.423e+02    1.838e-01  1.838e-01      8253 -3.812e+03 -4.180e+03                   
#> Path [24] :Best Iter: [41] ELBO (-3811.635687) evaluations: (8253) 
#> Path [25] :Initial log joint density = -481458.865091 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.715e-02   2.889e+02    1.666e-01  1.666e-01      8149 -3.791e+03 -4.169e+03                   
#> Path [25] :Best Iter: [42] ELBO (-3791.140360) evaluations: (8149) 
#> Path [26] :Initial log joint density = -481671.327218 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.139e-02   3.314e+02    3.789e-01  3.789e-01      8449 -3.791e+03 -4.125e+03                   
#> Path [26] :Best Iter: [42] ELBO (-3791.088613) evaluations: (8449) 
#> Path [27] :Initial log joint density = -481662.036099 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.689e-02   3.033e+02    6.264e-02  2.306e-01      8258 -3.759e+03 -4.242e+03                   
#> Path [27] :Best Iter: [38] ELBO (-3759.111315) evaluations: (8258) 
#> Path [28] :Initial log joint density = -481930.851229 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.284e-02   2.961e+02    1.997e-01  1.997e-01      7928 -3.834e+03 -4.154e+03                   
#> Path [28] :Best Iter: [39] ELBO (-3834.017396) evaluations: (7928) 
#> Path [29] :Initial log joint density = -484455.889574 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.173e-02   3.905e+02    2.047e-01  2.047e-01      8421 -3.871e+03 -4.230e+03                   
#> Path [29] :Best Iter: [45] ELBO (-3870.869611) evaluations: (8421) 
#> Path [30] :Initial log joint density = -481653.245587 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.498e-02   4.319e+02    1.081e-01  2.962e-01      8319 -3.796e+03 -4.187e+03                   
#> Path [30] :Best Iter: [37] ELBO (-3795.643364) evaluations: (8319) 
#> Path [31] :Initial log joint density = -481681.424981 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.393e-02   3.146e+02    3.687e-01  3.687e-01      8210 -3.766e+03 -4.182e+03                   
#> Path [31] :Best Iter: [37] ELBO (-3766.466741) evaluations: (8210) 
#> Path [32] :Initial log joint density = -481876.417866 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.594e-02   3.615e+02    3.048e-01  9.288e-01      8482 -3.765e+03 -4.159e+03                   
#> Path [32] :Best Iter: [41] ELBO (-3765.448005) evaluations: (8482) 
#> Path [33] :Initial log joint density = -481652.803061 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.743e-02   6.689e+02    6.446e-02  1.724e-01      8278 -3.816e+03 -4.328e+03                   
#> Path [33] :Best Iter: [37] ELBO (-3816.370070) evaluations: (8278) 
#> Path [34] :Initial log joint density = -482661.929480 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.134e-02   6.070e+02    5.481e-02  5.481e-02      8210 -3.778e+03 -4.535e+03                   
#> Path [34] :Best Iter: [42] ELBO (-3778.448716) evaluations: (8210) 
#> Path [35] :Initial log joint density = -485025.870963 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.214e-01   1.012e+03    1.145e-01  1.145e-01      7986 -3.794e+03 -6.363e+03                   
#> Path [35] :Best Iter: [39] ELBO (-3794.366688) evaluations: (7986) 
#> Path [36] :Initial log joint density = -482511.877102 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.797e-03   4.249e+02    4.867e-02  4.867e-02      8023 -3.765e+03 -4.549e+03                   
#> Path [36] :Best Iter: [48] ELBO (-3765.295629) evaluations: (8023) 
#> Path [37] :Initial log joint density = -482720.460399 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.241e-02   2.713e+02    4.420e-01  4.420e-01      8473 -3.788e+03 -4.231e+03                   
#> Path [37] :Best Iter: [36] ELBO (-3788.250439) evaluations: (8473) 
#> Path [38] :Initial log joint density = -482101.026914 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      7.639e-02   1.120e+03    1.007e-01  1.007e-01      8225 -3.806e+03 -4.762e+03                   
#> Path [38] :Best Iter: [32] ELBO (-3805.869119) evaluations: (8225) 
#> Path [39] :Initial log joint density = -482141.658327 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.286e-02   3.540e+02    1.827e-01  1.827e-01      8355 -3.768e+03 -4.382e+03                   
#> Path [39] :Best Iter: [40] ELBO (-3768.430688) evaluations: (8355) 
#> Path [40] :Initial log joint density = -481487.158882 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.930e-02   2.816e+02    2.817e-01  2.817e-01      8214 -3.772e+03 -4.148e+03                   
#> Path [40] :Best Iter: [39] ELBO (-3771.717560) evaluations: (8214) 
#> Path [41] :Initial log joint density = -485591.025507 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.017e-02   6.357e+02    5.010e-02  5.010e-02      8153 -3.779e+03 -4.947e+03                   
#> Path [41] :Best Iter: [40] ELBO (-3778.881979) evaluations: (8153) 
#> Path [42] :Initial log joint density = -481486.840387 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.218e-02   4.754e+02    1.644e-01  1.644e-01      8074 -3.786e+03 -4.195e+03                   
#> Path [42] :Best Iter: [35] ELBO (-3785.893372) evaluations: (8074) 
#> Path [43] :Initial log joint density = -481503.842108 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.446e-01   6.859e+02    1.631e-01  1.631e-01      8114 -3.780e+03 -4.726e+03                   
#> Path [43] :Best Iter: [34] ELBO (-3779.824596) evaluations: (8114) 
#> Path [44] :Initial log joint density = -482176.782904 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.295e-02   4.738e+02    2.276e-01  2.276e-01      8104 -3.819e+03 -4.162e+03                   
#> Path [44] :Best Iter: [44] ELBO (-3819.046592) evaluations: (8104) 
#> Path [45] :Initial log joint density = -481818.282637 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.484e-02   6.167e+02    2.758e-01  2.758e-01      8050 -3.783e+03 -4.091e+03                   
#> Path [45] :Best Iter: [38] ELBO (-3782.716032) evaluations: (8050) 
#> Path [46] :Initial log joint density = -481815.935745 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.365e-02   5.913e+02    2.913e-01  2.913e-01      8317 -3.805e+03 -4.090e+03                   
#> Path [46] :Best Iter: [37] ELBO (-3804.929411) evaluations: (8317) 
#> Path [47] :Initial log joint density = -481827.547365 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.550e-02   5.807e+02    2.421e-01  2.421e-01      8001 -3.775e+03 -4.225e+03                   
#> Path [47] :Best Iter: [39] ELBO (-3775.298159) evaluations: (8001) 
#> Path [48] :Initial log joint density = -482103.055330 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.959e-02   5.123e+02    9.233e-02  9.233e-02      8132 -3.770e+03 -4.300e+03                   
#> Path [48] :Best Iter: [41] ELBO (-3769.741588) evaluations: (8132) 
#> Path [49] :Initial log joint density = -481489.271237 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.302e-02   7.024e+02    6.928e-02  6.928e-02      8251 -3.781e+03 -5.673e+03                   
#> Path [49] :Best Iter: [37] ELBO (-3780.772558) evaluations: (8251) 
#> Path [50] :Initial log joint density = -481809.858878 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.032e-01   7.250e+02    5.717e-01  5.717e-01      8218 -3.798e+03 -4.113e+03                   
#> Path [50] :Best Iter: [35] ELBO (-3797.706566) evaluations: (8218) 
#> Finished in  27.9 seconds.
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Warning: Dropping 'draws_df' class as required metadata was removed.
#> Joining with `by = join_by(cell_group, M, parameter)`
#> === Bimodal Model Parameters ===
#> 
#> (Intercept):
#>   Component 1: v = -(4.134 + -0.232 × c)
#>   Component 2: v = -(5.626 + -0.551 × c)
#> 
#> typecancer:
#>   Component 1: v = -(-0.186 + 0.628 × c)
#>   Component 2: v = -(0.325 + 0.032 × c)
#> 
# }
```
