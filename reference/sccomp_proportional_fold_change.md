# Calculate Proportional Fold Change from sccomp Estimated Effects

This function calculates the proportional fold change between two
conditions using the estimated effects from a sccomp model. The fold
changes are derived from the model's posterior predictions rather than
raw counts, providing a more robust estimate that accounts for the
model's uncertainty and covariate effects.

Note! This statistic is descriptive and should not be used to define
significance - use sccomp_test() for that. While fold changes in
proportions are easier to interpret than changes in logit space, they
are not linear (the same proportional change has different meaning for
rare vs abundant cell types). In contrast, the logit scale used
internally by sccomp provides linear effects that are more appropriate
for statistical inference.

## Usage

``` r
sccomp_proportional_fold_change(.data, formula_composition, from, to)
```

## Arguments

- .data:

  A sccomp estimate object (of class 'sccomp_tbl') obtained from running
  sccomp_estimate(). This object contains the fitted model and estimated
  effects.

- formula_composition:

  The formula specifying which model effects to use for calculating fold
  changes. This should match or be a subset of the formula used in the
  original sccomp_estimate() call.

- from:

  Character string specifying the reference/control condition (e.g.,
  "benign").

- to:

  Character string specifying the comparison condition (e.g., "cancer").

## Value

A tibble with the following columns:

- cell_group - The cell group identifier

- proportion_fold_change - The estimated fold change in proportions
  between conditions. Positive values indicate increases, negative
  values indicate decreases.

- average_uncertainty - The average uncertainty in the fold change
  estimate, derived from the credible intervals

- statement - A text description of the fold change, including the
  direction and the estimated proportions

## Examples

``` r
print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))


# \donttest{
if (instantiate::stan_cmdstan_exists()) {
  # Load example data
  data("counts_obj")

  # First estimate the composition effects
  estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      "sample",
      "cell_group",
      "count",
      cores = 1
  )
 
  # Calculate proportional fold changes from the estimated effects
  estimate |> 
  sccomp_proportional_fold_change(
    formula_composition = ~  type, 
    from = "benign", 
    to = "cancer"
  ) 
}
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481914.254441 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.788e+05      1.344e-02   2.298e-01    1.000e+00  1.000e+00      3482 -3.699e+03 -3.699e+03                   
#> Path [1] :Best Iter: [59] ELBO (-3698.630534) evaluations: (3482) 
#> Path [2] :Initial log joint density = -481818.586658 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.788e+05      1.260e-02   3.344e-01    1.000e+00  1.000e+00      3757 -3.700e+03 -3.704e+03                   
#> Path [2] :Best Iter: [62] ELBO (-3700.255395) evaluations: (3757) 
#> Path [3] :Initial log joint density = -481765.536808 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      7.149e-03   2.278e-01    1.000e+00  1.000e+00      2871 -3.707e+03 -3.716e+03                   
#> Path [3] :Best Iter: [50] ELBO (-3706.618768) evaluations: (2871) 
#> Path [4] :Initial log joint density = -481876.510291 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      7.221e-03   2.516e-01    7.698e-01  7.698e-01      3655 -3.701e+03 -3.713e+03                   
#> Path [4] :Best Iter: [56] ELBO (-3701.045034) evaluations: (3655) 
#> Path [5] :Initial log joint density = -481627.121563 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      5.092e-03   2.148e-01    1.000e+00  1.000e+00      3253 -3.707e+03 -3.707e+03                   
#> Path [5] :Best Iter: [56] ELBO (-3706.865309) evaluations: (3253) 
#> Path [6] :Initial log joint density = -481706.097935 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.296e-02   3.193e-01    1.000e+00  1.000e+00      3424 -3.700e+03 -3.710e+03                   
#> Path [6] :Best Iter: [58] ELBO (-3700.091095) evaluations: (3424) 
#> Path [7] :Initial log joint density = -481637.581738 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.606e-03   2.053e-01    9.947e-01  9.947e-01      3033 -3.706e+03 -3.706e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3705.897658) evaluations: (3033) 
#> Path [8] :Initial log joint density = -481710.239087 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.233e-02   3.274e-01    1.000e+00  1.000e+00      3300 -3.700e+03 -3.706e+03                   
#> Path [8] :Best Iter: [56] ELBO (-3699.600965) evaluations: (3300) 
#> Path [9] :Initial log joint density = -482129.725456 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.592e-03   2.459e-01    7.563e-01  7.563e-01      3498 -3.697e+03 -3.716e+03                   
#> Path [9] :Best Iter: [57] ELBO (-3697.322375) evaluations: (3498) 
#> Path [10] :Initial log joint density = -481486.271671 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.421e-02   2.878e-01    8.221e-01  8.221e-01      2917 -3.705e+03 -3.719e+03                   
#> Path [10] :Best Iter: [43] ELBO (-3705.113627) evaluations: (2917) 
#> Path [11] :Initial log joint density = -481642.924797 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      4.838e-03   2.438e-01    1.000e+00  1.000e+00      2829 -3.705e+03 -3.718e+03                   
#> Path [11] :Best Iter: [44] ELBO (-3704.779052) evaluations: (2829) 
#> Path [12] :Initial log joint density = -481384.333874 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      4.927e-03   2.350e-01    7.933e-01  7.933e-01      2971 -3.707e+03 -3.709e+03                   
#> Path [12] :Best Iter: [42] ELBO (-3707.102648) evaluations: (2971) 
#> Path [13] :Initial log joint density = -481612.584775 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      6.533e-03   2.284e-01    1.000e+00  1.000e+00      3126 -3.709e+03 -3.705e+03                   
#> Path [13] :Best Iter: [55] ELBO (-3704.726354) evaluations: (3126) 
#> Path [14] :Initial log joint density = -481478.396531 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      8.919e-03   1.968e-01    1.000e+00  1.000e+00      3086 -3.708e+03 -3.708e+03                   
#> Path [14] :Best Iter: [55] ELBO (-3707.674444) evaluations: (3086) 
#> Path [15] :Initial log joint density = -481466.598342 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      1.125e-02   1.953e-01    1.000e+00  1.000e+00      2985 -3.708e+03 -3.706e+03                   
#> Path [15] :Best Iter: [53] ELBO (-3706.092157) evaluations: (2985) 
#> Path [16] :Initial log joint density = -482113.566318 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.020e-02   2.547e-01    1.000e+00  1.000e+00      3414 -3.700e+03 -3.708e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3700.092143) evaluations: (3414) 
#> Path [17] :Initial log joint density = -482596.193615 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.788e+05      1.224e-02   3.849e-01    1.000e+00  1.000e+00      3808 -3.701e+03 -3.704e+03                   
#> Path [17] :Best Iter: [63] ELBO (-3701.232497) evaluations: (3808) 
#> Path [18] :Initial log joint density = -481590.640570 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.788e+05      1.901e-02   3.175e-01    1.000e+00  1.000e+00      3624 -3.703e+03 -3.709e+03                   
#> Path [18] :Best Iter: [57] ELBO (-3703.015069) evaluations: (3624) 
#> Path [19] :Initial log joint density = -481700.531952 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      1.454e-02   2.026e-01    1.000e+00  1.000e+00      3444 -3.698e+03 -3.698e+03                   
#> Path [19] :Best Iter: [58] ELBO (-3697.739017) evaluations: (3444) 
#> Path [20] :Initial log joint density = -481723.215826 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      5.410e-03   1.960e-01    8.769e-01  8.769e-01      3330 -3.702e+03 -3.709e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3701.621845) evaluations: (3330) 
#> Path [21] :Initial log joint density = -481343.770810 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      7.267e-03   1.748e-01    1.000e+00  1.000e+00      2989 -3.706e+03 -3.712e+03                   
#> Path [21] :Best Iter: [45] ELBO (-3706.230721) evaluations: (2989) 
#> Path [22] :Initial log joint density = -481442.615271 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.775e-03   2.301e-01    1.000e+00  1.000e+00      3148 -3.704e+03 -3.704e+03                   
#> Path [22] :Best Iter: [57] ELBO (-3703.518289) evaluations: (3148) 
#> Path [23] :Initial log joint density = -481451.035027 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.177e-02   4.659e-01    1.000e+00  1.000e+00      2758 -3.703e+03 -3.716e+03                   
#> Path [23] :Best Iter: [49] ELBO (-3703.455689) evaluations: (2758) 
#> Path [24] :Initial log joint density = -481500.394152 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      3.958e-03   2.010e-01    2.809e-01  1.000e+00      3043 -3.708e+03 -3.722e+03                   
#> Path [24] :Best Iter: [42] ELBO (-3708.245358) evaluations: (3043) 
#> Path [25] :Initial log joint density = -481281.261877 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      5.369e-03   2.311e-01    1.000e+00  1.000e+00      2863 -3.707e+03 -3.714e+03                   
#> Path [25] :Best Iter: [49] ELBO (-3707.312498) evaluations: (2863) 
#> Path [26] :Initial log joint density = -481543.411051 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.570e-02   3.821e-01    9.827e-01  9.827e-01      3026 -3.709e+03 -3.725e+03                   
#> Path [26] :Best Iter: [52] ELBO (-3709.385014) evaluations: (3026) 
#> Path [27] :Initial log joint density = -483718.641252 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.125e-02   3.207e-01    1.000e+00  1.000e+00      3122 -3.702e+03 -3.706e+03                   
#> Path [27] :Best Iter: [55] ELBO (-3701.568425) evaluations: (3122) 
#> Path [28] :Initial log joint density = -481814.356931 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      7.203e-03   2.412e-01    1.000e+00  1.000e+00      3100 -3.705e+03 -3.707e+03                   
#> Path [28] :Best Iter: [52] ELBO (-3704.926026) evaluations: (3100) 
#> Path [29] :Initial log joint density = -481623.938103 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.789e-03   2.186e-01    1.000e+00  1.000e+00      3092 -3.709e+03 -3.704e+03                   
#> Path [29] :Best Iter: [55] ELBO (-3703.645047) evaluations: (3092) 
#> Path [30] :Initial log joint density = -483701.346680 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.788e+05      7.184e-03   1.939e-01    1.000e+00  1.000e+00      3112 -3.707e+03 -3.702e+03                   
#> Path [30] :Best Iter: [55] ELBO (-3701.611034) evaluations: (3112) 
#> Path [31] :Initial log joint density = -481617.480201 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      5.622e-03   2.388e-01    1.000e+00  1.000e+00      3235 -3.706e+03 -3.702e+03                   
#> Path [31] :Best Iter: [57] ELBO (-3702.439474) evaluations: (3235) 
#> Path [32] :Initial log joint density = -481172.912101 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.946e-03   2.309e-01    1.000e+00  1.000e+00      2911 -3.707e+03 -3.709e+03                   
#> Path [32] :Best Iter: [49] ELBO (-3707.479718) evaluations: (2911) 
#> Path [33] :Initial log joint density = -481566.371381 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.745e-03   1.669e-01    8.217e-01  8.217e-01      2966 -3.708e+03 -3.721e+03                   
#> Path [33] :Best Iter: [41] ELBO (-3707.916646) evaluations: (2966) 
#> Path [34] :Initial log joint density = -481553.819334 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      8.936e-03   2.420e-01    9.085e-01  9.085e-01      2863 -3.707e+03 -3.718e+03                   
#> Path [34] :Best Iter: [46] ELBO (-3707.413107) evaluations: (2863) 
#> Path [35] :Initial log joint density = -481770.814883 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      9.410e-03   2.219e-01    8.703e-01  8.703e-01      2863 -3.709e+03 -3.721e+03                   
#> Path [35] :Best Iter: [37] ELBO (-3708.586423) evaluations: (2863) 
#> Path [36] :Initial log joint density = -482082.557371 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      1.281e-02   3.331e-01    1.000e+00  1.000e+00      3409 -3.701e+03 -3.705e+03                   
#> Path [36] :Best Iter: [55] ELBO (-3701.427288) evaluations: (3409) 
#> Path [37] :Initial log joint density = -481838.593237 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      1.148e-02   3.537e-01    1.000e+00  1.000e+00      3074 -3.706e+03 -3.716e+03                   
#> Path [37] :Best Iter: [49] ELBO (-3706.063398) evaluations: (3074) 
#> Path [38] :Initial log joint density = -481490.145935 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      3.793e-03   2.896e-01    6.487e-01  6.487e-01      2961 -3.707e+03 -3.712e+03                   
#> Path [38] :Best Iter: [51] ELBO (-3706.623309) evaluations: (2961) 
#> Path [39] :Initial log joint density = -481821.751824 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      3.172e-03   3.097e-01    6.432e-01  6.432e-01      3405 -3.700e+03 -3.711e+03                   
#> Path [39] :Best Iter: [57] ELBO (-3700.216004) evaluations: (3405) 
#> Path [40] :Initial log joint density = -481692.994182 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.788e+05      4.107e-03   1.525e-01    1.000e+00  1.000e+00      2957 -3.709e+03 -3.719e+03                   
#> Path [40] :Best Iter: [41] ELBO (-3708.777497) evaluations: (2957) 
#> Path [41] :Initial log joint density = -481929.216371 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.788e+05      8.835e-03   1.593e-01    1.000e+00  1.000e+00      3541 -3.700e+03 -3.703e+03                   
#> Path [41] :Best Iter: [58] ELBO (-3700.321096) evaluations: (3541) 
#> Path [42] :Initial log joint density = -481299.443730 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.788e+05      9.371e-03   2.402e-01    1.000e+00  1.000e+00      3261 -3.699e+03 -3.702e+03                   
#> Path [42] :Best Iter: [57] ELBO (-3699.077056) evaluations: (3261) 
#> Path [43] :Initial log joint density = -481764.461205 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.788e+05      3.977e-03   2.180e-01    8.589e-01  8.589e-01      2879 -3.705e+03 -3.718e+03                   
#> Path [43] :Best Iter: [50] ELBO (-3705.088919) evaluations: (2879) 
#> Path [44] :Initial log joint density = -482119.354569 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      2.299e-03   2.250e-01    6.226e-01  6.226e-01      3415 -3.699e+03 -3.710e+03                   
#> Path [44] :Best Iter: [57] ELBO (-3698.509230) evaluations: (3415) 
#> Path [45] :Initial log joint density = -481710.707689 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              51      -4.788e+05      1.536e-02   3.680e-01    1.000e+00  1.000e+00      2802 -3.709e+03 -3.715e+03                   
#> Path [45] :Best Iter: [46] ELBO (-3709.048476) evaluations: (2802) 
#> Path [46] :Initial log joint density = -481592.748758 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.788e+05      4.756e-03   2.438e-01    7.818e-01  7.818e-01      3338 -3.702e+03 -3.711e+03                   
#> Path [46] :Best Iter: [57] ELBO (-3701.896142) evaluations: (3338) 
#> Path [47] :Initial log joint density = -481589.875077 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.788e+05      6.214e-03   2.134e-01    1.000e+00  1.000e+00      2993 -3.708e+03 -3.713e+03                   
#> Path [47] :Best Iter: [43] ELBO (-3707.609880) evaluations: (2993) 
#> Path [48] :Initial log joint density = -481573.131542 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.788e+05      1.184e-02   2.650e-01    1.000e+00  1.000e+00      3213 -3.701e+03 -3.701e+03                   
#> Path [48] :Best Iter: [56] ELBO (-3701.204024) evaluations: (3213) 
#> Path [49] :Initial log joint density = -481910.635912 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      1.193e-02   3.140e-01    1.000e+00  1.000e+00      3192 -3.706e+03 -3.714e+03                   
#> Path [49] :Best Iter: [55] ELBO (-3705.598498) evaluations: (3192) 
#> Path [50] :Initial log joint density = -481824.294116 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.788e+05      4.722e-03   1.632e-01    1.000e+00  1.000e+00      3188 -3.708e+03 -3.705e+03                   
#> Path [50] :Best Iter: [57] ELBO (-3705.453666) evaluations: (3188) 
#> Finished in  13.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 4
#>    cell_group proportion_fold_change average_uncertainty statement              
#>    <chr>                       <dbl>               <dbl> <glue>                 
#>  1 B1                          -1.93             0.199   1.9-fold decrease (fro…
#>  2 B2                          -2.04             0.243   2-fold decrease (from …
#>  3 B3                          -1.38             0.126   1.4-fold decrease (fro…
#>  4 BM                          -1.39             0.178   1.4-fold decrease (fro…
#>  5 CD4 1                        1.16             0.00381 1.2-fold increase (fro…
#>  6 CD4 2                        1.44             0.0627  1.4-fold increase (fro…
#>  7 CD4 3                       -2.34             0.349   2.3-fold decrease (fro…
#>  8 CD4 4                       -1.06             0.0930  1.1-fold decrease (fro…
#>  9 CD4 5                        1.02             0.122   1-fold increase (from …
#> 10 CD8 1                        1.10             0.0330  1.1-fold increase (fro…
#> # ℹ 26 more rows
# }
```
