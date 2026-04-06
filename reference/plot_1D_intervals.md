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
#> Path [1] :Initial log joint density = -481568.724340 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.205e-02   1.406e-01    1.000e+00  1.000e+00      3628 -3.679e+03 -3.684e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3679.320117) evaluations: (3628) 
#> Path [2] :Initial log joint density = -481457.793922 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.520e-02   1.775e-01    1.000e+00  1.000e+00      3788 -3.680e+03 -3.681e+03                   
#> Path [2] :Best Iter: [59] ELBO (-3679.772788) evaluations: (3788) 
#> Path [3] :Initial log joint density = -481783.387554 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.302e-02   2.576e-01    1.000e+00  1.000e+00      3443 -3.686e+03 -3.694e+03                   
#> Path [3] :Best Iter: [57] ELBO (-3686.008018) evaluations: (3443) 
#> Path [4] :Initial log joint density = -481849.653455 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      6.411e-03   2.332e-01    1.000e+00  1.000e+00      2916 -3.690e+03 -3.698e+03                   
#> Path [4] :Best Iter: [44] ELBO (-3689.790477) evaluations: (2916) 
#> Path [5] :Initial log joint density = -481877.087118 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.150e-03   1.963e-01    1.000e+00  1.000e+00      3117 -3.693e+03 -3.687e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3686.691939) evaluations: (3117) 
#> Path [6] :Initial log joint density = -482218.087875 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      6.676e-03   2.509e-01    1.000e+00  1.000e+00      3417 -3.687e+03 -3.693e+03                   
#> Path [6] :Best Iter: [55] ELBO (-3686.872761) evaluations: (3417) 
#> Path [7] :Initial log joint density = -481858.159099 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      8.475e-03   2.021e-01    1.000e+00  1.000e+00      3453 -3.684e+03 -3.689e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3683.781912) evaluations: (3453) 
#> Path [8] :Initial log joint density = -485045.478860 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.236e-03   2.162e-01    1.000e+00  1.000e+00      3731 -3.682e+03 -3.692e+03                   
#> Path [8] :Best Iter: [60] ELBO (-3681.549983) evaluations: (3731) 
#> Path [9] :Initial log joint density = -481442.020578 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.211e-02   2.140e-01    1.000e+00  1.000e+00      3375 -3.683e+03 -3.690e+03                   
#> Path [9] :Best Iter: [56] ELBO (-3683.337369) evaluations: (3375) 
#> Path [10] :Initial log joint density = -481813.229065 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      8.941e-03   2.101e-01    9.506e-01  9.506e-01      3286 -3.687e+03 -3.694e+03                   
#> Path [10] :Best Iter: [55] ELBO (-3686.792112) evaluations: (3286) 
#> Path [11] :Initial log joint density = -481795.933084 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      8.096e-03   2.427e-01    1.000e+00  1.000e+00      3363 -3.685e+03 -3.698e+03                   
#> Path [11] :Best Iter: [55] ELBO (-3684.899774) evaluations: (3363) 
#> Path [12] :Initial log joint density = -481585.713750 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.715e-03   1.748e-01    8.414e-01  8.414e-01      3358 -3.684e+03 -3.692e+03                   
#> Path [12] :Best Iter: [57] ELBO (-3683.858176) evaluations: (3358) 
#> Path [13] :Initial log joint density = -482050.121739 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.117e-02   1.898e-01    1.000e+00  1.000e+00      3412 -3.682e+03 -3.686e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3682.403496) evaluations: (3412) 
#> Path [14] :Initial log joint density = -481617.039105 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.536e-03   3.300e-01    1.000e+00  1.000e+00      2951 -3.691e+03 -3.700e+03                   
#> Path [14] :Best Iter: [51] ELBO (-3690.774633) evaluations: (2951) 
#> Path [15] :Initial log joint density = -481706.462598 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      7.540e-03   2.456e-01    8.807e-01  8.807e-01      3583 -3.686e+03 -3.697e+03                   
#> Path [15] :Best Iter: [59] ELBO (-3686.267591) evaluations: (3583) 
#> Path [16] :Initial log joint density = -481709.693981 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      5.036e-03   1.588e-01    1.000e+00  1.000e+00      3432 -3.685e+03 -3.695e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3684.689868) evaluations: (3432) 
#> Path [17] :Initial log joint density = -481722.594060 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.454e-02   1.836e-01    1.000e+00  1.000e+00      3384 -3.684e+03 -3.686e+03                   
#> Path [17] :Best Iter: [57] ELBO (-3684.473767) evaluations: (3384) 
#> Path [18] :Initial log joint density = -482012.226279 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.233e-02   3.309e-01    1.000e+00  1.000e+00      3447 -3.683e+03 -3.688e+03                   
#> Path [18] :Best Iter: [56] ELBO (-3683.441156) evaluations: (3447) 
#> Path [19] :Initial log joint density = -481551.056153 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      6.388e-03   2.008e-01    1.000e+00  1.000e+00      3286 -3.689e+03 -3.686e+03                   
#> Path [19] :Best Iter: [56] ELBO (-3685.950380) evaluations: (3286) 
#> Path [20] :Initial log joint density = -482332.915202 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.193e-02   4.002e-01    1.000e+00  1.000e+00      3870 -3.684e+03 -3.694e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3684.133187) evaluations: (3870) 
#> Path [21] :Initial log joint density = -481487.045471 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.698e-02   2.019e-01    1.000e+00  1.000e+00      3736 -3.680e+03 -3.685e+03                   
#> Path [21] :Best Iter: [60] ELBO (-3679.743607) evaluations: (3736) 
#> Path [22] :Initial log joint density = -481303.254084 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      7.327e-03   2.350e-01    1.000e+00  1.000e+00      2750 -3.691e+03 -3.704e+03                   
#> Path [22] :Best Iter: [43] ELBO (-3690.696270) evaluations: (2750) 
#> Path [23] :Initial log joint density = -481737.576936 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      4.490e-03   2.193e-01    8.589e-01  8.589e-01      3507 -3.682e+03 -3.693e+03                   
#> Path [23] :Best Iter: [58] ELBO (-3681.871148) evaluations: (3507) 
#> Path [24] :Initial log joint density = -481598.019324 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      4.983e-03   2.781e-01    7.514e-01  7.514e-01      3119 -3.691e+03 -3.696e+03                   
#> Path [24] :Best Iter: [55] ELBO (-3691.411477) evaluations: (3119) 
#> Path [25] :Initial log joint density = -481923.853383 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.034e-03   2.296e-01    4.350e-01  1.000e+00      2985 -3.692e+03 -3.704e+03                   
#> Path [25] :Best Iter: [43] ELBO (-3691.582791) evaluations: (2985) 
#> Path [26] :Initial log joint density = -481541.374325 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      5.294e-03   1.464e-01    1.000e+00  1.000e+00      2989 -3.689e+03 -3.703e+03                   
#> Path [26] :Best Iter: [50] ELBO (-3689.031777) evaluations: (2989) 
#> Path [27] :Initial log joint density = -482143.966369 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.500e-03   1.884e-01    1.000e+00  1.000e+00      3277 -3.684e+03 -3.684e+03                   
#> Path [27] :Best Iter: [56] ELBO (-3684.029465) evaluations: (3277) 
#> Path [28] :Initial log joint density = -482684.340466 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      1.642e-02   3.706e-01    1.000e+00  1.000e+00      3587 -3.680e+03 -3.686e+03                   
#> Path [28] :Best Iter: [59] ELBO (-3680.298761) evaluations: (3587) 
#> Path [29] :Initial log joint density = -481664.420491 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      6.008e-03   2.110e-01    1.000e+00  1.000e+00      2844 -3.689e+03 -3.690e+03                   
#> Path [29] :Best Iter: [38] ELBO (-3688.606300) evaluations: (2844) 
#> Path [30] :Initial log joint density = -481509.254177 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      9.827e-03   3.470e-01    1.000e+00  1.000e+00      2811 -3.691e+03 -3.700e+03                   
#> Path [30] :Best Iter: [50] ELBO (-3690.594995) evaluations: (2811) 
#> Path [31] :Initial log joint density = -481411.385448 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      6.144e-03   2.363e-01    1.000e+00  1.000e+00      2783 -3.692e+03 -3.709e+03                   
#> Path [31] :Best Iter: [46] ELBO (-3692.263557) evaluations: (2783) 
#> Path [32] :Initial log joint density = -482675.758667 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.843e-03   1.594e-01    4.666e-01  1.000e+00      3236 -3.684e+03 -3.698e+03                   
#> Path [32] :Best Iter: [56] ELBO (-3684.417506) evaluations: (3236) 
#> Path [33] :Initial log joint density = -481606.054741 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      6.340e-03   2.106e-01    1.000e+00  1.000e+00      3176 -3.690e+03 -3.698e+03                   
#> Path [33] :Best Iter: [48] ELBO (-3689.978444) evaluations: (3176) 
#> Path [34] :Initial log joint density = -484998.856338 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.811e-03   2.989e-01    5.857e-01  5.857e-01      3567 -3.682e+03 -3.696e+03                   
#> Path [34] :Best Iter: [57] ELBO (-3681.872805) evaluations: (3567) 
#> Path [35] :Initial log joint density = -481618.126138 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      8.992e-03   3.321e-01    1.000e+00  1.000e+00      3311 -3.684e+03 -3.693e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3683.801251) evaluations: (3311) 
#> Path [36] :Initial log joint density = -481747.520420 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.577e-03   2.207e-01    1.000e+00  1.000e+00      2830 -3.689e+03 -3.694e+03                   
#> Path [36] :Best Iter: [51] ELBO (-3689.376752) evaluations: (2830) 
#> Path [37] :Initial log joint density = -481394.044492 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      3.759e-03   2.014e-01    7.977e-01  7.977e-01      3463 -3.682e+03 -3.692e+03                   
#> Path [37] :Best Iter: [58] ELBO (-3681.600555) evaluations: (3463) 
#> Path [38] :Initial log joint density = -481640.992500 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.254e-02   3.371e-01    1.000e+00  1.000e+00      2939 -3.693e+03 -3.698e+03                   
#> Path [38] :Best Iter: [49] ELBO (-3692.645653) evaluations: (2939) 
#> Path [39] :Initial log joint density = -481660.182446 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.398e-02   3.777e-01    1.000e+00  1.000e+00      3185 -3.682e+03 -3.688e+03                   
#> Path [39] :Best Iter: [56] ELBO (-3681.961513) evaluations: (3185) 
#> Path [40] :Initial log joint density = -481796.680862 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      6.264e-03   1.749e-01    1.000e+00  1.000e+00      3495 -3.683e+03 -3.698e+03                   
#> Path [40] :Best Iter: [57] ELBO (-3682.502392) evaluations: (3495) 
#> Path [41] :Initial log joint density = -483017.334642 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      5.405e-03   1.758e-01    1.000e+00  1.000e+00      3072 -3.687e+03 -3.693e+03                   
#> Path [41] :Best Iter: [43] ELBO (-3687.375992) evaluations: (3072) 
#> Path [42] :Initial log joint density = -481807.850827 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.490e-03   2.283e-01    7.703e-01  7.703e-01      3007 -3.688e+03 -3.707e+03                   
#> Path [42] :Best Iter: [51] ELBO (-3688.052744) evaluations: (3007) 
#> Path [43] :Initial log joint density = -481441.347704 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.863e-03   1.721e-01    1.000e+00  1.000e+00      3145 -3.692e+03 -3.687e+03                   
#> Path [43] :Best Iter: [55] ELBO (-3687.495515) evaluations: (3145) 
#> Path [44] :Initial log joint density = -482572.741036 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      9.153e-03   2.130e-01    1.000e+00  1.000e+00      3112 -3.690e+03 -3.686e+03                   
#> Path [44] :Best Iter: [55] ELBO (-3685.789453) evaluations: (3112) 
#> Path [45] :Initial log joint density = -481733.135350 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      7.776e-03   2.034e-01    1.000e+00  1.000e+00      3246 -3.689e+03 -3.688e+03                   
#> Path [45] :Best Iter: [57] ELBO (-3687.908579) evaluations: (3246) 
#> Path [46] :Initial log joint density = -481796.042267 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      7.602e-03   2.312e-01    8.558e-01  8.558e-01      3330 -3.683e+03 -3.692e+03                   
#> Path [46] :Best Iter: [56] ELBO (-3682.944886) evaluations: (3330) 
#> Path [47] :Initial log joint density = -481818.209991 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      1.240e-02   3.210e-01    1.000e+00  1.000e+00      2830 -3.689e+03 -3.701e+03                   
#> Path [47] :Best Iter: [43] ELBO (-3689.232444) evaluations: (2830) 
#> Path [48] :Initial log joint density = -481552.921425 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.523e-03   2.108e-01    9.495e-01  9.495e-01      3081 -3.690e+03 -3.701e+03                   
#> Path [48] :Best Iter: [53] ELBO (-3690.106362) evaluations: (3081) 
#> Path [49] :Initial log joint density = -481959.054824 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      5.273e-03   1.835e-01    1.000e+00  1.000e+00      3682 -3.685e+03 -3.686e+03                   
#> Path [49] :Best Iter: [58] ELBO (-3684.922046) evaluations: (3682) 
#> Path [50] :Initial log joint density = -481713.130670 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      3.998e-03   1.414e-01    1.000e+00  1.000e+00      3137 -3.687e+03 -3.695e+03                   
#> Path [50] :Best Iter: [38] ELBO (-3686.687272) evaluations: (3137) 
#> Finished in  13.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
# }

```
