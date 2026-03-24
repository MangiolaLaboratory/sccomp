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
#> Path [1] :Initial log joint density = -487723.819127 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.941e-02   7.555e+02    3.190e-02  3.190e-02      8318 -3.769e+03 -8.170e+03                   
#> Path [1] :Best Iter: [54] ELBO (-3768.532568) evaluations: (8318) 
#> Path [2] :Initial log joint density = -483256.412805 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.805e-02   1.571e+03    1.365e-02  3.253e-02      8228 -3.781e+03 -5.620e+04                   
#> Path [2] :Best Iter: [40] ELBO (-3781.038020) evaluations: (8228) 
#> Path [3] :Initial log joint density = -481487.714339 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      8.522e-02   5.178e+03    2.248e-02  2.248e-02      8421 -3.809e+03 -2.378e+06                   
#> Path [3] :Best Iter: [34] ELBO (-3808.920526) evaluations: (8421) 
#> Path [4] :Initial log joint density = -481906.567470 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.104e-02   6.488e+03    2.633e-02  2.633e-02      8398 -3.809e+03 -5.724e+04                   
#> Path [4] :Best Iter: [34] ELBO (-3808.776890) evaluations: (8398) 
#> Path [5] :Initial log joint density = -483512.115888 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.962e-02   7.193e+03    7.076e-03  1.410e-02      8282 -3.792e+03 -2.787e+06                   
#> Path [5] :Best Iter: [37] ELBO (-3792.474065) evaluations: (8282) 
#> Path [6] :Initial log joint density = -482410.573240 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.726e-02   3.774e-01    1.000e+00  1.000e+00      6118 -3.735e+03 -3.748e+03                   
#> Path [6] :Best Iter: [79] ELBO (-3734.571094) evaluations: (6118) 
#> Path [7] :Initial log joint density = -481833.264680 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      4.439e-03   1.374e-01    7.613e-01  7.613e-01      6310 -3.732e+03 -3.761e+03                   
#> Path [7] :Best Iter: [79] ELBO (-3732.155288) evaluations: (6310) 
#> Path [8] :Initial log joint density = -481607.300829 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.441e-02   1.046e+04    1.406e-02  1.406e-02      8651 -3.802e+03 -8.181e+11                   
#> Path [8] :Best Iter: [34] ELBO (-3802.138707) evaluations: (8651) 
#> Path [9] :Initial log joint density = -481654.491631 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.489e-02   2.361e+03    1.689e-02  1.689e-02      8166 -3.771e+03 -2.495e+07                   
#> Path [9] :Best Iter: [45] ELBO (-3771.018214) evaluations: (8166) 
#> Path [10] :Initial log joint density = -481748.407462 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.054e-02   4.998e+03    9.362e-03  9.362e-03      8236 -3.776e+03 -5.488e+07                   
#> Path [10] :Best Iter: [43] ELBO (-3775.990959) evaluations: (8236) 
#> Path [11] :Initial log joint density = -482292.043369 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.787e+05      1.666e-02   2.618e-01    1.000e+00  1.000e+00      6524 -3.734e+03 -3.744e+03                   
#> Path [11] :Best Iter: [88] ELBO (-3733.603163) evaluations: (6524) 
#> Path [12] :Initial log joint density = -481537.610596 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.402e-02   8.338e+03    3.009e-02  3.009e-02      8403 -3.847e+03 -6.279e+03                   
#> Path [12] :Best Iter: [33] ELBO (-3846.903471) evaluations: (8403) 
#> Path [13] :Initial log joint density = -481622.753160 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      5.140e-02   7.266e+03    2.241e-02  2.241e-02      8384 -3.796e+03 -3.088e+05                   
#> Path [13] :Best Iter: [40] ELBO (-3796.456836) evaluations: (8384) 
#> Path [14] :Initial log joint density = -482964.325446 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.630e-02   8.374e+03    1.478e-02  1.478e-02      8421 -3.794e+03 -3.437e+06                   
#> Path [14] :Best Iter: [37] ELBO (-3793.950480) evaluations: (8421) 
#> Path [15] :Initial log joint density = -481800.167777 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.789e-02   1.945e+03    2.751e-02  2.751e-02      7970 -3.771e+03 -7.345e+04                   
#> Path [15] :Best Iter: [47] ELBO (-3770.977673) evaluations: (7970) 
#> Path [16] :Initial log joint density = -481855.512450 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.560e-01   1.154e+03    4.256e-02  4.256e-02      8071 -3.767e+03 -2.280e+05                   
#> Path [16] :Best Iter: [59] ELBO (-3767.320942) evaluations: (8071) 
#> Path [17] :Initial log joint density = -481873.536493 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      1.378e-02   2.494e-01    1.000e+00  1.000e+00      6074 -3.732e+03 -3.743e+03                   
#> Path [17] :Best Iter: [83] ELBO (-3732.081275) evaluations: (6074) 
#> Path [18] :Initial log joint density = -482096.707618 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      1.528e-02   2.166e-01    1.000e+00  1.000e+00      6115 -3.735e+03 -3.744e+03                   
#> Path [18] :Best Iter: [85] ELBO (-3734.850586) evaluations: (6115) 
#> Path [19] :Initial log joint density = -482592.389343 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.968e-02   1.801e+03    2.896e-02  2.896e-02      8139 -3.782e+03 -5.218e+06                   
#> Path [19] :Best Iter: [38] ELBO (-3782.185814) evaluations: (8139) 
#> Path [20] :Initial log joint density = -481550.632248 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.661e-02   3.601e+03    1.989e-02  1.989e-02      8145 -3.779e+03 -3.484e+06                   
#> Path [20] :Best Iter: [45] ELBO (-3779.192940) evaluations: (8145) 
#> Path [21] :Initial log joint density = -483292.262311 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      9.215e-03   1.720e-01    1.000e+00  1.000e+00      6114 -3.737e+03 -3.744e+03                   
#> Path [21] :Best Iter: [82] ELBO (-3737.146907) evaluations: (6114) 
#> Path [22] :Initial log joint density = -481672.214971 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.674e-02   2.647e+03    1.110e-02  2.501e-02      8255 -3.789e+03 -1.257e+05                   
#> Path [22] :Best Iter: [44] ELBO (-3788.652139) evaluations: (8255) 
#> Path [23] :Initial log joint density = -481843.220788 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.093e-02   2.336e+03    2.012e-02  2.012e-02      8006 -3.777e+03 -3.504e+04                   
#> Path [23] :Best Iter: [44] ELBO (-3777.445565) evaluations: (8006) 
#> Path [24] :Initial log joint density = -481694.133194 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.358e-02   1.118e+03    3.950e-02  3.950e-02      8045 -3.764e+03 -6.140e+03                   
#> Path [24] :Best Iter: [54] ELBO (-3763.568296) evaluations: (8045) 
#> Path [25] :Initial log joint density = -481865.366899 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      8.058e-02   1.064e+04    9.264e-03  4.495e-02      8411 -3.803e+03 -2.595e+05                   
#> Path [25] :Best Iter: [38] ELBO (-3803.414801) evaluations: (8411) 
#> Path [26] :Initial log joint density = -481840.699072 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.827e-02   1.916e+03    2.522e-02  2.522e-02      8343 -3.776e+03 -8.509e+05                   
#> Path [26] :Best Iter: [52] ELBO (-3775.712205) evaluations: (8343) 
#> Path [27] :Initial log joint density = -483122.507548 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.401e-01   4.007e+02    5.451e-02  1.696e-01      8097 -3.751e+03 -5.340e+03                   
#> Path [27] :Best Iter: [73] ELBO (-3751.086995) evaluations: (8097) 
#> Path [28] :Initial log joint density = -481932.338488 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      2.063e-02   2.102e-01    1.000e+00  1.000e+00      5912 -3.735e+03 -3.738e+03                   
#> Path [28] :Best Iter: [75] ELBO (-3735.479883) evaluations: (5912) 
#> Path [29] :Initial log joint density = -482086.966826 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      9.923e-03   1.854e-01    1.000e+00  1.000e+00      5654 -3.736e+03 -3.763e+03                   
#> Path [29] :Best Iter: [77] ELBO (-3736.287666) evaluations: (5654) 
#> Path [30] :Initial log joint density = -482025.281989 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.042e-02   1.802e-01    1.000e+00  1.000e+00      5558 -3.734e+03 -3.737e+03                   
#> Path [30] :Best Iter: [77] ELBO (-3734.467637) evaluations: (5558) 
#> Path [31] :Initial log joint density = -481445.651808 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.128e-02   1.891e+03    2.900e-02  2.900e-02      8190 -3.768e+03 -1.248e+05                   
#> Path [31] :Best Iter: [48] ELBO (-3768.091348) evaluations: (8190) 
#> Path [32] :Initial log joint density = -481881.883849 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      1.376e-02   1.918e-01    1.000e+00  1.000e+00      5713 -3.735e+03 -3.741e+03                   
#> Path [32] :Best Iter: [75] ELBO (-3734.702248) evaluations: (5713) 
#> Path [33] :Initial log joint density = -481723.112757 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.588e-02   9.575e+03    1.240e-02  1.240e-02      8474 -3.790e+03 -1.504e+08                   
#> Path [33] :Best Iter: [35] ELBO (-3789.978159) evaluations: (8474) 
#> Path [34] :Initial log joint density = -481618.999389 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      2.509e-02   3.896e+03    1.700e-02  1.700e-02      8400 -3.805e+03 -5.219e+05                   
#> Path [34] :Best Iter: [34] ELBO (-3804.951168) evaluations: (8400) 
#> Path [35] :Initial log joint density = -481421.139312 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.544e-02   8.775e+03    6.720e-03  6.720e-03      8294 -3.791e+03 -1.525e+10                   
#> Path [35] :Best Iter: [43] ELBO (-3790.959380) evaluations: (8294) 
#> Path [36] :Initial log joint density = -481576.982082 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.015e-02   3.439e+03    1.308e-02  1.308e-02      8401 -3.788e+03 -9.346e+05                   
#> Path [36] :Best Iter: [38] ELBO (-3787.657845) evaluations: (8401) 
#> Path [37] :Initial log joint density = -481788.701034 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.560e-01   1.801e+03    2.662e-02  2.662e-02      7941 -3.764e+03 -1.015e+06                   
#> Path [37] :Best Iter: [52] ELBO (-3764.116400) evaluations: (7941) 
#> Path [38] :Initial log joint density = -481810.978914 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.704e-02   6.143e+03    2.002e-02  2.002e-02      8334 -3.797e+03 -4.816e+06                   
#> Path [38] :Best Iter: [37] ELBO (-3796.746700) evaluations: (8334) 
#> Path [39] :Initial log joint density = -481571.490412 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      3.070e-02   3.721e+03    8.753e-03  2.303e-02      8356 -3.800e+03 -1.463e+05                   
#> Path [39] :Best Iter: [43] ELBO (-3800.210768) evaluations: (8356) 
#> Path [40] :Initial log joint density = -485696.982324 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.421e-01   4.650e+02    1.371e-01  1.371e-01      8236 -3.751e+03 -5.629e+03                   
#> Path [40] :Best Iter: [73] ELBO (-3750.724212) evaluations: (8236) 
#> Path [41] :Initial log joint density = -481765.074436 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      4.403e-02   1.534e+04    4.371e-02  4.371e-02      8751 -3.818e+03 -1.248e+04                   
#> Path [41] :Best Iter: [35] ELBO (-3817.807949) evaluations: (8751) 
#> Path [42] :Initial log joint density = -481671.163321 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.370e-02   4.120e+03    3.054e-02  3.054e-02      8386 -3.795e+03 -1.228e+04                   
#> Path [42] :Best Iter: [37] ELBO (-3794.697103) evaluations: (8386) 
#> Path [43] :Initial log joint density = -481904.875953 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      3.820e-03   1.879e-01    6.908e-01  6.908e-01      5832 -3.731e+03 -3.760e+03                   
#> Path [43] :Best Iter: [81] ELBO (-3730.508495) evaluations: (5832) 
#> Path [44] :Initial log joint density = -481390.503018 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      6.310e-02   4.344e+03    1.182e-02  3.971e-02      8276 -3.798e+03 -9.575e+07                   
#> Path [44] :Best Iter: [38] ELBO (-3797.874813) evaluations: (8276) 
#> Path [45] :Initial log joint density = -481913.949462 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      9.783e-02   2.002e+03    1.453e-02  1.453e-02      8079 -3.792e+03 -4.749e+06                   
#> Path [45] :Best Iter: [44] ELBO (-3792.044862) evaluations: (8079) 
#> Path [46] :Initial log joint density = -481561.545777 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.080e-01   8.544e+02    3.805e-02  3.805e-02      8074 -3.771e+03 -3.058e+05                   
#> Path [46] :Best Iter: [47] ELBO (-3770.686975) evaluations: (8074) 
#> Path [47] :Initial log joint density = -482258.754345 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.196e-02   2.435e-01    1.000e+00  1.000e+00      5707 -3.735e+03 -3.749e+03                   
#> Path [47] :Best Iter: [77] ELBO (-3735.110360) evaluations: (5707) 
#> Path [48] :Initial log joint density = -481496.432045 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.403e-02   1.332e+03    2.013e-02  2.013e-02      8107 -3.791e+03 -4.414e+05                   
#> Path [48] :Best Iter: [45] ELBO (-3791.350193) evaluations: (8107) 
#> Path [49] :Initial log joint density = -481614.172036 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.786e+05      1.465e-01   5.013e+03    2.152e-02  2.152e-02      8315 -3.798e+03 -1.621e+07                   
#> Path [49] :Best Iter: [41] ELBO (-3798.102597) evaluations: (8315) 
#> Path [50] :Initial log joint density = -482751.745350 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.272e-01   8.980e+02    9.153e-02  9.153e-02      8174 -3.767e+03 -9.836e+03                   
#> Path [50] :Best Iter: [54] ELBO (-3767.051853) evaluations: (8174) 
#> Finished in  25.2 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }
```
