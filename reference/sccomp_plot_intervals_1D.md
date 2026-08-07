# Plot 1D Intervals for Cell-group Effects

This function creates a series of 1D interval plots for cell-group
effects, highlighting significant differences based on a given
significance threshold.

## Usage

``` r
sccomp_plot_intervals_1D(
  .data,
  factor = NULL,
  significance_threshold = 0.05,
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
  show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR"),
  sort_by = c("none", "effect", "significance", "alphabetical")
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

- sort_by:

  Character vector indicating how to sort taxa. Options are "none"
  (default), "effect" (by effect size), "significance" (by FDR/pH0), or
  "alphabetical".

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
  my_plot = sccomp_plot_intervals_1D(estimate, sort_by = "effect")

  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482162.087811 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.116e-02   1.893e-01    1.000e+00  1.000e+00      4340 -3.688e+03 -3.685e+03                   
#> Path [1] :Best Iter: [69] ELBO (-3684.705552) evaluations: (4340) 
#> Path [2] :Initial log joint density = -481882.810251 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      9.165e-03   2.182e-01    1.000e+00  1.000e+00      3832 -3.684e+03 -3.697e+03                   
#> Path [2] :Best Iter: [61] ELBO (-3684.430116) evaluations: (3832) 
#> Path [3] :Initial log joint density = -482144.538754 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.724e-03   1.689e-01    7.525e-01  7.525e-01      4467 -3.686e+03 -3.703e+03                   
#> Path [3] :Best Iter: [68] ELBO (-3685.968922) evaluations: (4467) 
#> Path [4] :Initial log joint density = -481690.850991 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      8.071e-03   2.111e-01    9.411e-01  9.411e-01      4163 -3.683e+03 -3.696e+03                   
#> Path [4] :Best Iter: [65] ELBO (-3682.995119) evaluations: (4163) 
#> Path [5] :Initial log joint density = -481908.807007 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.230e-02   2.055e-01    1.000e+00  1.000e+00      4420 -3.688e+03 -3.691e+03                   
#> Path [5] :Best Iter: [63] ELBO (-3687.721888) evaluations: (4420) 
#> Path [6] :Initial log joint density = -483238.373024 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.419e-02   1.568e-01    1.000e+00  1.000e+00      4447 -3.686e+03 -3.693e+03                   
#> Path [6] :Best Iter: [62] ELBO (-3686.163958) evaluations: (4447) 
#> Path [7] :Initial log joint density = -484811.849619 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      2.078e-03   2.244e-01    5.985e-01  5.985e-01      3436 -3.689e+03 -3.703e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3688.687093) evaluations: (3436) 
#> Path [8] :Initial log joint density = -481575.993334 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.403e-02   2.263e-01    1.000e+00  1.000e+00      3162 -3.692e+03 -3.691e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3690.812693) evaluations: (3162) 
#> Path [9] :Initial log joint density = -481603.430358 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.204e-02   2.451e-01    1.000e+00  1.000e+00      3960 -3.686e+03 -3.688e+03                   
#> Path [9] :Best Iter: [61] ELBO (-3686.375921) evaluations: (3960) 
#> Path [10] :Initial log joint density = -481865.773706 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      9.066e-03   2.133e-01    1.000e+00  1.000e+00      4590 -3.686e+03 -3.693e+03                   
#> Path [10] :Best Iter: [63] ELBO (-3686.228565) evaluations: (4590) 
#> Path [11] :Initial log joint density = -481580.051716 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.134e-02   3.370e-01    1.000e+00  1.000e+00      3026 -3.696e+03 -3.700e+03                   
#> Path [11] :Best Iter: [42] ELBO (-3695.742937) evaluations: (3026) 
#> Path [12] :Initial log joint density = -481572.245318 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      6.645e-03   1.962e-01    1.000e+00  1.000e+00      3368 -3.691e+03 -3.695e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3691.474748) evaluations: (3368) 
#> Path [13] :Initial log joint density = -482730.161681 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.485e-02   1.376e-01    1.000e+00  1.000e+00      4441 -3.681e+03 -3.690e+03                   
#> Path [13] :Best Iter: [66] ELBO (-3681.451059) evaluations: (4441) 
#> Path [14] :Initial log joint density = -481593.387025 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      1.498e-02   2.582e-01    1.000e+00  1.000e+00      2962 -3.696e+03 -3.699e+03                   
#> Path [14] :Best Iter: [49] ELBO (-3696.008672) evaluations: (2962) 
#> Path [15] :Initial log joint density = -481457.396248 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.258e-02   3.180e-01    1.000e+00  1.000e+00      4443 -3.683e+03 -3.687e+03                   
#> Path [15] :Best Iter: [68] ELBO (-3682.518570) evaluations: (4443) 
#> Path [16] :Initial log joint density = -481650.140359 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.256e-02   1.448e-01    1.000e+00  1.000e+00      4372 -3.685e+03 -3.687e+03                   
#> Path [16] :Best Iter: [61] ELBO (-3685.103353) evaluations: (4372) 
#> Path [17] :Initial log joint density = -481774.302230 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.620e-02   2.562e-01    6.690e-01  6.690e-01      4242 -3.686e+03 -3.693e+03                   
#> Path [17] :Best Iter: [58] ELBO (-3686.223776) evaluations: (4242) 
#> Path [18] :Initial log joint density = -481737.991455 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.047e-02   2.677e-01    1.000e+00  1.000e+00      4295 -3.685e+03 -3.692e+03                   
#> Path [18] :Best Iter: [61] ELBO (-3685.267596) evaluations: (4295) 
#> Path [19] :Initial log joint density = -481712.617961 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      1.087e-02   2.351e-01    1.000e+00  1.000e+00      3014 -3.695e+03 -3.699e+03                   
#> Path [19] :Best Iter: [49] ELBO (-3694.973358) evaluations: (3014) 
#> Path [20] :Initial log joint density = -481951.668512 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      6.096e-03   1.539e-01    1.000e+00  1.000e+00      5272 -3.686e+03 -3.689e+03                   
#> Path [20] :Best Iter: [75] ELBO (-3685.943193) evaluations: (5272) 
#> Path [21] :Initial log joint density = -481485.891020 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      2.400e-03   2.448e-01    6.907e-01  6.907e-01      3387 -3.687e+03 -3.701e+03                   
#> Path [21] :Best Iter: [56] ELBO (-3686.961850) evaluations: (3387) 
#> Path [22] :Initial log joint density = -481464.103140 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      5.939e-03   2.355e-01    1.000e+00  1.000e+00      3040 -3.692e+03 -3.701e+03                   
#> Path [22] :Best Iter: [49] ELBO (-3692.349164) evaluations: (3040) 
#> Path [23] :Initial log joint density = -481565.797818 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      4.204e-03   1.345e-01    7.986e-01  7.986e-01      3342 -3.693e+03 -3.705e+03                   
#> Path [23] :Best Iter: [56] ELBO (-3692.552339) evaluations: (3342) 
#> Path [24] :Initial log joint density = -481395.118982 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.528e-02   2.533e-01    1.000e+00  1.000e+00      4094 -3.683e+03 -3.684e+03                   
#> Path [24] :Best Iter: [64] ELBO (-3682.767857) evaluations: (4094) 
#> Path [25] :Initial log joint density = -481703.574993 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      7.769e-03   2.173e-01    1.000e+00  1.000e+00      4357 -3.687e+03 -3.690e+03                   
#> Path [25] :Best Iter: [66] ELBO (-3686.789352) evaluations: (4357) 
#> Path [26] :Initial log joint density = -481502.500419 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.964e-02   3.921e-01    1.000e+00  1.000e+00      4561 -3.684e+03 -3.692e+03                   
#> Path [26] :Best Iter: [70] ELBO (-3683.707968) evaluations: (4561) 
#> Path [27] :Initial log joint density = -481456.557828 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      9.060e-03   2.228e-01    8.361e-01  8.361e-01      4001 -3.686e+03 -3.700e+03                   
#> Path [27] :Best Iter: [65] ELBO (-3685.856526) evaluations: (4001) 
#> Path [28] :Initial log joint density = -483794.040309 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.144e-02   2.339e-01    1.000e+00  1.000e+00      4229 -3.685e+03 -3.693e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3685.442969) evaluations: (4229) 
#> Path [29] :Initial log joint density = -481565.251881 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      7.071e-03   2.529e-01    1.000e+00  1.000e+00      2944 -3.696e+03 -3.703e+03                   
#> Path [29] :Best Iter: [47] ELBO (-3695.640428) evaluations: (2944) 
#> Path [30] :Initial log joint density = -481633.948224 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.986e-02   3.026e-01    1.000e+00  1.000e+00      5240 -3.681e+03 -3.689e+03                   
#> Path [30] :Best Iter: [77] ELBO (-3681.294314) evaluations: (5240) 
#> Path [31] :Initial log joint density = -481283.698361 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      7.223e-03   1.868e-01    8.576e-01  8.576e-01      3399 -3.687e+03 -3.703e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3686.905471) evaluations: (3399) 
#> Path [32] :Initial log joint density = -482416.572196 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.491e-03   1.646e-01    1.000e+00  1.000e+00      4804 -3.686e+03 -3.695e+03                   
#> Path [32] :Best Iter: [69] ELBO (-3685.916491) evaluations: (4804) 
#> Path [33] :Initial log joint density = -482669.542325 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      4.464e-03   3.039e-01    7.051e-01  7.051e-01      3327 -3.687e+03 -3.701e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3687.363187) evaluations: (3327) 
#> Path [34] :Initial log joint density = -481199.737397 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.115e-02   1.719e-01    1.000e+00  1.000e+00      4084 -3.684e+03 -3.686e+03                   
#> Path [34] :Best Iter: [63] ELBO (-3684.103974) evaluations: (4084) 
#> Path [35] :Initial log joint density = -483900.624315 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      7.441e-03   2.179e-01    1.000e+00  1.000e+00      3298 -3.690e+03 -3.695e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3689.793469) evaluations: (3298) 
#> Path [36] :Initial log joint density = -482099.505025 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      6.677e-03   2.396e-01    1.000e+00  1.000e+00      3497 -3.690e+03 -3.700e+03                   
#> Path [36] :Best Iter: [56] ELBO (-3689.919640) evaluations: (3497) 
#> Path [37] :Initial log joint density = -481574.403834 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      3.524e-03   2.003e-01    7.482e-01  7.482e-01      3962 -3.682e+03 -3.698e+03                   
#> Path [37] :Best Iter: [63] ELBO (-3682.193553) evaluations: (3962) 
#> Path [38] :Initial log joint density = -481495.839736 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.812e-03   2.125e-01    1.000e+00  1.000e+00      3879 -3.690e+03 -3.700e+03                   
#> Path [38] :Best Iter: [57] ELBO (-3689.889135) evaluations: (3879) 
#> Path [39] :Initial log joint density = -481392.056952 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      8.315e-03   1.880e-01    1.000e+00  1.000e+00      4695 -3.684e+03 -3.689e+03                   
#> Path [39] :Best Iter: [69] ELBO (-3684.203351) evaluations: (4695) 
#> Path [40] :Initial log joint density = -481732.847862 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.747e-03   2.462e-01    5.633e-01  5.633e-01      4548 -3.684e+03 -3.701e+03                   
#> Path [40] :Best Iter: [64] ELBO (-3684.174031) evaluations: (4548) 
#> Path [41] :Initial log joint density = -481542.363579 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.697e-02   3.231e-01    1.000e+00  1.000e+00      4183 -3.684e+03 -3.689e+03                   
#> Path [41] :Best Iter: [66] ELBO (-3683.594284) evaluations: (4183) 
#> Path [42] :Initial log joint density = -481719.896579 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.069e-02   1.892e-01    1.000e+00  1.000e+00      5025 -3.686e+03 -3.686e+03                   
#> Path [42] :Best Iter: [75] ELBO (-3685.843221) evaluations: (5025) 
#> Path [43] :Initial log joint density = -481647.400759 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.107e-03   2.034e-01    8.030e-01  8.030e-01      3942 -3.690e+03 -3.700e+03                   
#> Path [43] :Best Iter: [61] ELBO (-3690.195563) evaluations: (3942) 
#> Path [44] :Initial log joint density = -481736.642925 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.361e-02   1.930e-01    4.036e-01  1.000e+00      4214 -3.686e+03 -3.698e+03                   
#> Path [44] :Best Iter: [60] ELBO (-3685.824135) evaluations: (4214) 
#> Path [45] :Initial log joint density = -481445.135159 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      7.687e-03   2.561e-01    1.000e+00  1.000e+00      4430 -3.686e+03 -3.694e+03                   
#> Path [45] :Best Iter: [59] ELBO (-3686.459842) evaluations: (4430) 
#> Path [46] :Initial log joint density = -481379.263203 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.765e-02   1.426e-01    1.000e+00  1.000e+00      4233 -3.688e+03 -3.688e+03                   
#> Path [46] :Best Iter: [67] ELBO (-3687.570386) evaluations: (4233) 
#> Path [47] :Initial log joint density = -482312.464233 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      1.162e-02   4.385e-01    1.000e+00  1.000e+00      3077 -3.697e+03 -3.704e+03                   
#> Path [47] :Best Iter: [53] ELBO (-3696.729754) evaluations: (3077) 
#> Path [48] :Initial log joint density = -481726.758339 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.593e-02   3.009e-01    9.812e-01  9.812e-01      4264 -3.690e+03 -3.697e+03                   
#> Path [48] :Best Iter: [65] ELBO (-3690.290101) evaluations: (4264) 
#> Path [49] :Initial log joint density = -481571.313945 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.411e-02   2.095e-01    1.000e+00  1.000e+00      4043 -3.683e+03 -3.685e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3682.671460) evaluations: (4043) 
#> Path [50] :Initial log joint density = -481520.353531 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      6.388e-03   1.668e-01    1.000e+00  1.000e+00      3026 -3.693e+03 -3.703e+03                   
#> Path [50] :Best Iter: [52] ELBO (-3693.399316) evaluations: (3026) 
#> Finished in  15.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Joining with `by = join_by(cell_group, M, parameter)`
# }

```
