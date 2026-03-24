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
  significance_statistic = c("pH0", "FDR")
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
  "pH0".

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
#> Path [1] :Initial log joint density = -487720.254774 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.394e-02   2.146e+03    2.330e-02  2.330e-02      8214 -3.754e+03 -1.079e+08                   
#> Path [1] :Best Iter: [55] ELBO (-3754.055173) evaluations: (8214) 
#> Path [2] :Initial log joint density = -483255.243152 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.874e-02   1.548e+03    2.009e-02  4.349e-02      8453 -3.770e+03 -2.419e+04                   
#> Path [2] :Best Iter: [36] ELBO (-3770.243563) evaluations: (8453) 
#> Path [3] :Initial log joint density = -481486.501872 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.196e-01   2.328e+03    4.306e-02  4.306e-02      8285 -3.779e+03 -1.922e+04                   
#> Path [3] :Best Iter: [39] ELBO (-3778.956833) evaluations: (8285) 
#> Path [4] :Initial log joint density = -481904.143288 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      9.836e-02   5.114e+03    3.539e-02  3.539e-02      8453 -3.795e+03 -2.582e+06                   
#> Path [4] :Best Iter: [33] ELBO (-3794.634464) evaluations: (8453) 
#> Path [5] :Initial log joint density = -483510.034020 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.525e-02   4.147e+03    9.969e-03  2.349e-02      8425 -3.769e+03 -2.110e+06                   
#> Path [5] :Best Iter: [38] ELBO (-3768.690572) evaluations: (8425) 
#> Path [6] :Initial log joint density = -482410.061556 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.096e-02   2.457e-01    6.636e-01  6.636e-01      5950 -3.716e+03 -3.730e+03                   
#> Path [6] :Best Iter: [82] ELBO (-3715.674810) evaluations: (5950) 
#> Path [7] :Initial log joint density = -481831.236976 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.085e-02   2.174e-01    4.200e-01  1.000e+00      6259 -3.714e+03 -3.740e+03                   
#> Path [7] :Best Iter: [80] ELBO (-3714.447849) evaluations: (6259) 
#> Path [8] :Initial log joint density = -481606.590984 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.158e-02   7.678e+03    8.648e-03  8.648e-03      8338 -3.775e+03 -2.676e+05                   
#> Path [8] :Best Iter: [38] ELBO (-3774.727861) evaluations: (8338) 
#> Path [9] :Initial log joint density = -481652.261322 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.625e-01   1.106e+03    5.530e-02  5.530e-02      8315 -3.749e+03 -3.293e+05                   
#> Path [9] :Best Iter: [50] ELBO (-3748.794336) evaluations: (8315) 
#> Path [10] :Initial log joint density = -481747.968665 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.841e-02   4.052e+03    1.488e-02  1.488e-02      8067 -3.761e+03 -3.843e+05                   
#> Path [10] :Best Iter: [52] ELBO (-3760.660360) evaluations: (8067) 
#> Path [11] :Initial log joint density = -482291.529104 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      3.305e-02   1.977e-01    1.000e+00  1.000e+00      6146 -3.713e+03 -3.715e+03                   
#> Path [11] :Best Iter: [82] ELBO (-3712.761240) evaluations: (6146) 
#> Path [12] :Initial log joint density = -481536.232770 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.059e-01   6.488e+03    7.534e-02  7.534e-02      8354 -3.832e+03 -5.460e+03                   
#> Path [12] :Best Iter: [28] ELBO (-3831.862971) evaluations: (8354) 
#> Path [13] :Initial log joint density = -481622.098271 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.378e-02   6.170e+03    6.518e-03  2.107e-02      8164 -3.778e+03 -3.816e+13                   
#> Path [13] :Best Iter: [37] ELBO (-3777.983718) evaluations: (8164) 
#> Path [14] :Initial log joint density = -482962.378266 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.722e-02   2.158e+03    3.347e-02  3.347e-02      8542 -3.767e+03 -2.273e+04                   
#> Path [14] :Best Iter: [39] ELBO (-3767.072550) evaluations: (8542) 
#> Path [15] :Initial log joint density = -481797.094700 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      9.565e-02   9.882e+02    2.901e-02  2.901e-02      8033 -3.749e+03 -2.401e+04                   
#> Path [15] :Best Iter: [56] ELBO (-3748.937375) evaluations: (8033) 
#> Path [16] :Initial log joint density = -481853.785592 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.147e-02   2.212e-01    1.000e+00  1.000e+00      5613 -3.715e+03 -3.716e+03                   
#> Path [16] :Best Iter: [78] ELBO (-3715.424951) evaluations: (5613) 
#> Path [17] :Initial log joint density = -481873.383054 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.202e-02   1.207e-01    1.000e+00  1.000e+00      6091 -3.714e+03 -3.732e+03                   
#> Path [17] :Best Iter: [83] ELBO (-3713.522827) evaluations: (6091) 
#> Path [18] :Initial log joint density = -482095.529142 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      4.734e-03   1.623e-01    4.470e-01  1.000e+00      6211 -3.714e+03 -3.738e+03                   
#> Path [18] :Best Iter: [83] ELBO (-3714.101234) evaluations: (6211) 
#> Path [19] :Initial log joint density = -482591.012068 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.295e-02   4.443e+03    2.984e-02  2.984e-02      8102 -3.759e+03 -6.586e+07                   
#> Path [19] :Best Iter: [42] ELBO (-3758.630982) evaluations: (8102) 
#> Path [20] :Initial log joint density = -481550.540691 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.062e-02   1.545e+03    1.969e-02  1.969e-02      8045 -3.750e+03 -5.029e+04                   
#> Path [20] :Best Iter: [51] ELBO (-3749.877388) evaluations: (8045) 
#> Path [21] :Initial log joint density = -483291.685919 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      1.197e-02   2.446e-01    1.000e+00  1.000e+00      6000 -3.711e+03 -3.730e+03                   
#> Path [21] :Best Iter: [79] ELBO (-3711.276906) evaluations: (6000) 
#> Path [22] :Initial log joint density = -481671.592245 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      6.272e-02   3.710e+03    1.688e-02  1.688e-02      8201 -3.764e+03 -2.520e+08                   
#> Path [22] :Best Iter: [45] ELBO (-3764.243487) evaluations: (8201) 
#> Path [23] :Initial log joint density = -481840.768634 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.559e-02   9.726e+02    1.836e-02  1.836e-02      8187 -3.749e+03 -2.176e+04                   
#> Path [23] :Best Iter: [55] ELBO (-3749.096171) evaluations: (8187) 
#> Path [24] :Initial log joint density = -481693.732190 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.257e-01   1.438e+02    9.565e-02  9.565e-02      7796 -3.719e+03 -4.159e+03                   
#> Path [24] :Best Iter: [81] ELBO (-3718.947403) evaluations: (7796) 
#> Path [25] :Initial log joint density = -481864.566949 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.385e-02   9.433e+03    4.862e-02  4.862e-02      8210 -3.783e+03 -3.473e+04                   
#> Path [25] :Best Iter: [37] ELBO (-3783.359446) evaluations: (8210) 
#> Path [26] :Initial log joint density = -481840.465854 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.667e-02   5.042e+02    5.188e-02  5.188e-02      8513 -3.746e+03 -8.367e+03                   
#> Path [26] :Best Iter: [56] ELBO (-3745.817672) evaluations: (8513) 
#> Path [27] :Initial log joint density = -483120.957867 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              96      -4.787e+05      2.785e-02   2.576e-01    8.895e-01  8.895e-01      7336 -3.707e+03 -3.721e+03                   
#> Path [27] :Best Iter: [89] ELBO (-3706.664769) evaluations: (7336) 
#> Path [28] :Initial log joint density = -481930.772568 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      1.264e-02   1.863e-01    1.000e+00  1.000e+00      5997 -3.718e+03 -3.719e+03                   
#> Path [28] :Best Iter: [82] ELBO (-3718.247201) evaluations: (5997) 
#> Path [29] :Initial log joint density = -482086.620780 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      2.690e-02   2.454e-01    1.000e+00  1.000e+00      5806 -3.719e+03 -3.717e+03                   
#> Path [29] :Best Iter: [83] ELBO (-3716.608120) evaluations: (5806) 
#> Path [30] :Initial log joint density = -482024.298770 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.037e-02   1.701e-01    8.057e-01  8.057e-01      6114 -3.716e+03 -3.743e+03                   
#> Path [30] :Best Iter: [83] ELBO (-3716.206010) evaluations: (6114) 
#> Path [31] :Initial log joint density = -481444.550214 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.494e-01   2.957e+02    1.256e-01  1.256e-01      8024 -3.733e+03 -8.765e+03                   
#> Path [31] :Best Iter: [76] ELBO (-3732.835460) evaluations: (8024) 
#> Path [32] :Initial log joint density = -481880.852177 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      8.092e-03   2.752e-01    6.953e-01  6.953e-01      5556 -3.715e+03 -3.744e+03                   
#> Path [32] :Best Iter: [75] ELBO (-3715.325551) evaluations: (5556) 
#> Path [33] :Initial log joint density = -481721.714034 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.412e-02   3.581e+03    1.547e-02  1.547e-02      8380 -3.782e+03 -1.068e+06                   
#> Path [33] :Best Iter: [37] ELBO (-3782.407400) evaluations: (8380) 
#> Path [34] :Initial log joint density = -481618.917783 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.794e-02   4.547e+03    9.812e-03  9.812e-03      8406 -3.788e+03 -6.429e+05                   
#> Path [34] :Best Iter: [34] ELBO (-3787.805170) evaluations: (8406) 
#> Path [35] :Initial log joint density = -481420.959890 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.346e-01   1.049e+04    3.132e-02  3.132e-02      8421 -3.780e+03 -3.315e+09                   
#> Path [35] :Best Iter: [40] ELBO (-3780.075570) evaluations: (8421) 
#> Path [36] :Initial log joint density = -481575.913035 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.639e-01   3.620e+03    4.139e-02  4.139e-02      8245 -3.757e+03 -5.823e+11                   
#> Path [36] :Best Iter: [49] ELBO (-3757.457237) evaluations: (8245) 
#> Path [37] :Initial log joint density = -481788.568944 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      2.017e-02   3.213e-01    1.000e+00  1.000e+00      5700 -3.713e+03 -3.727e+03                   
#> Path [37] :Best Iter: [82] ELBO (-3713.059845) evaluations: (5700) 
#> Path [38] :Initial log joint density = -481810.975540 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.980e-02   5.300e+03    1.278e-02  1.278e-02      8187 -3.780e+03 -1.194e+04                   
#> Path [38] :Best Iter: [40] ELBO (-3779.719772) evaluations: (8187) 
#> Path [39] :Initial log joint density = -481569.801399 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.737e-02   3.728e+03    1.440e-02  3.088e-02      8222 -3.780e+03 -7.837e+07                   
#> Path [39] :Best Iter: [38] ELBO (-3779.673906) evaluations: (8222) 
#> Path [40] :Initial log joint density = -485694.132216 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      1.099e-02   1.716e-01    1.000e+00  1.000e+00      5922 -3.715e+03 -3.730e+03                   
#> Path [40] :Best Iter: [80] ELBO (-3715.294351) evaluations: (5922) 
#> Path [41] :Initial log joint density = -481764.009190 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.200e-02   6.953e+03    1.273e-02  1.273e-02      8607 -3.798e+03 -2.630e+07                   
#> Path [41] :Best Iter: [36] ELBO (-3798.281192) evaluations: (8607) 
#> Path [42] :Initial log joint density = -481671.122171 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.346e-02   3.886e+03    2.046e-02  2.046e-02      8140 -3.777e+03 -3.573e+05                   
#> Path [42] :Best Iter: [44] ELBO (-3776.870000) evaluations: (8140) 
#> Path [43] :Initial log joint density = -481902.960625 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      5.199e-03   1.664e-01    1.000e+00  1.000e+00      6133 -3.716e+03 -3.730e+03                   
#> Path [43] :Best Iter: [83] ELBO (-3715.889303) evaluations: (6133) 
#> Path [44] :Initial log joint density = -481387.494888 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.869e-02   2.091e+03    1.226e-02  1.226e-02      8339 -3.767e+03 -8.588e+04                   
#> Path [44] :Best Iter: [43] ELBO (-3766.843087) evaluations: (8339) 
#> Path [45] :Initial log joint density = -481912.689228 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.104e-02   1.255e+03    2.483e-02  2.483e-02      8082 -3.759e+03 -1.658e+05                   
#> Path [45] :Best Iter: [51] ELBO (-3758.657849) evaluations: (8082) 
#> Path [46] :Initial log joint density = -481561.527414 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.031e-01   1.101e+03    4.642e-02  4.642e-02      8081 -3.749e+03 -1.519e+04                   
#> Path [46] :Best Iter: [51] ELBO (-3749.364917) evaluations: (8081) 
#> Path [47] :Initial log joint density = -482256.158294 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      7.232e-03   1.852e-01    1.000e+00  1.000e+00      5591 -3.715e+03 -3.724e+03                   
#> Path [47] :Best Iter: [77] ELBO (-3714.541327) evaluations: (5591) 
#> Path [48] :Initial log joint density = -481495.909585 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.218e-02   2.322e+03    1.314e-02  1.314e-02      8324 -3.773e+03 -2.426e+05                   
#> Path [48] :Best Iter: [44] ELBO (-3773.041354) evaluations: (8324) 
#> Path [49] :Initial log joint density = -481612.142308 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.438e-02   2.628e+03    1.452e-02  2.729e-02      8525 -3.771e+03 -9.449e+04                   
#> Path [49] :Best Iter: [40] ELBO (-3771.267320) evaluations: (8525) 
#> Path [50] :Initial log joint density = -482750.261273 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.153e-02   1.209e-01    1.000e+00  1.000e+00      6036 -3.720e+03 -3.731e+03                   
#> Path [50] :Best Iter: [82] ELBO (-3719.702972) evaluations: (6036) 
#> Finished in  24.9 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
