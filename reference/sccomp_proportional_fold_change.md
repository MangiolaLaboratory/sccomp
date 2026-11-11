# Calculate Proportional Fold Change from sccomp Estimated Effects

This function calculates the proportional fold change between two
conditions using the estimated effects from a sccomp model. The fold
changes are derived from the model's posterior predictions rather than
raw counts, providing a more robust estimate that accounts for the
model's uncertainty and covariate effects.

Note! This statistic is descriptive and should not be used to define
significance - use sccomp_test() for that. While fold changes in
proportions are easier to interpret than changes in centered log-ratio
space, they are not linear (the same proportional change has different
meaning for rare vs abundant cell types). In contrast, the centered
log-ratio scale used internally by sccomp provides linear effects that
are more appropriate for statistical inference.

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
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), typecancer
#>   variability Xa: 1 parameter: (Intercept)
#>   random effects: none
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481873.445600 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.929e-03   2.169e-01    1.000e+00  1.000e+00      4751 -3.703e+03 -3.711e+03                   
#> Path [1] :Best Iter: [70] ELBO (-3703.262509) evaluations: (4751) 
#> Path [2] :Initial log joint density = -482110.628704 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.457e-02   2.043e-01    1.000e+00  1.000e+00      5686 -3.701e+03 -3.700e+03                   
#> Path [2] :Best Iter: [82] ELBO (-3699.926496) evaluations: (5686) 
#> Path [3] :Initial log joint density = -482405.881619 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.311e-02   1.806e-01    1.000e+00  1.000e+00      4438 -3.701e+03 -3.701e+03                   
#> Path [3] :Best Iter: [66] ELBO (-3701.249644) evaluations: (4438) 
#> Path [4] :Initial log joint density = -481515.077042 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      8.704e-03   2.291e-01    1.000e+00  1.000e+00      3869 -3.700e+03 -3.710e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3699.903199) evaluations: (3869) 
#> Path [5] :Initial log joint density = -481559.541343 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      2.118e-02   2.449e-01    1.000e+00  1.000e+00      4241 -3.705e+03 -3.703e+03                   
#> Path [5] :Best Iter: [68] ELBO (-3702.554557) evaluations: (4241) 
#> Path [6] :Initial log joint density = -482034.695132 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.639e-02   3.144e-01    1.000e+00  1.000e+00      4040 -3.704e+03 -3.706e+03                   
#> Path [6] :Best Iter: [64] ELBO (-3703.731394) evaluations: (4040) 
#> Path [7] :Initial log joint density = -481506.472591 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      2.466e-02   1.911e-01    1.000e+00  1.000e+00      4503 -3.702e+03 -3.704e+03                   
#> Path [7] :Best Iter: [62] ELBO (-3702.317549) evaluations: (4503) 
#> Path [8] :Initial log joint density = -481529.619480 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      4.224e-02   2.224e-01    1.000e+00  1.000e+00      3917 -3.703e+03 -3.702e+03                   
#> Path [8] :Best Iter: [64] ELBO (-3701.535895) evaluations: (3917) 
#> Path [9] :Initial log joint density = -481488.397008 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.313e-02   3.080e-01    1.000e+00  1.000e+00      3610 -3.706e+03 -3.711e+03                   
#> Path [9] :Best Iter: [58] ELBO (-3706.249016) evaluations: (3610) 
#> Path [10] :Initial log joint density = -481474.601207 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      9.839e-03   1.516e-01    1.000e+00  1.000e+00      4844 -3.700e+03 -3.703e+03                   
#> Path [10] :Best Iter: [73] ELBO (-3700.175342) evaluations: (4844) 
#> Path [11] :Initial log joint density = -482298.439806 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.028e-02   2.130e-01    1.000e+00  1.000e+00      3691 -3.703e+03 -3.708e+03                   
#> Path [11] :Best Iter: [59] ELBO (-3703.009071) evaluations: (3691) 
#> Path [12] :Initial log joint density = -481417.083632 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      3.119e-03   3.533e-01    6.478e-01  6.478e-01      3441 -3.709e+03 -3.723e+03                   
#> Path [12] :Best Iter: [57] ELBO (-3709.241555) evaluations: (3441) 
#> Path [13] :Initial log joint density = -481521.178805 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      6.240e-03   2.294e-01    1.000e+00  1.000e+00      3136 -3.712e+03 -3.721e+03                   
#> Path [13] :Best Iter: [48] ELBO (-3711.535740) evaluations: (3136) 
#> Path [14] :Initial log joint density = -481605.592211 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      9.500e-03   2.026e-01    9.216e-01  9.216e-01      3960 -3.700e+03 -3.713e+03                   
#> Path [14] :Best Iter: [63] ELBO (-3699.705708) evaluations: (3960) 
#> Path [15] :Initial log joint density = -481434.633914 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      1.533e-02   3.093e-01    1.000e+00  1.000e+00      3775 -3.698e+03 -3.708e+03                   
#> Path [15] :Best Iter: [61] ELBO (-3698.195461) evaluations: (3775) 
#> Path [16] :Initial log joint density = -481549.812410 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.566e-02   2.077e-01    8.799e-01  8.799e-01      4375 -3.701e+03 -3.713e+03                   
#> Path [16] :Best Iter: [62] ELBO (-3701.496835) evaluations: (4375) 
#> Path [17] :Initial log joint density = -481873.904932 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.460e-02   3.308e-01    1.000e+00  1.000e+00      4040 -3.702e+03 -3.714e+03                   
#> Path [17] :Best Iter: [65] ELBO (-3701.500327) evaluations: (4040) 
#> Path [18] :Initial log joint density = -481851.930741 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      9.954e-03   2.025e-01    1.000e+00  1.000e+00      5202 -3.698e+03 -3.706e+03                   
#> Path [18] :Best Iter: [75] ELBO (-3697.831632) evaluations: (5202) 
#> Path [19] :Initial log joint density = -481606.433377 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      5.221e-03   1.603e-01    1.000e+00  1.000e+00      3166 -3.713e+03 -3.718e+03                   
#> Path [19] :Best Iter: [51] ELBO (-3712.901733) evaluations: (3166) 
#> Path [20] :Initial log joint density = -485264.280941 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.549e-03   2.731e-01    8.090e-01  8.090e-01      3715 -3.705e+03 -3.713e+03                   
#> Path [20] :Best Iter: [58] ELBO (-3705.070099) evaluations: (3715) 
#> Path [21] :Initial log joint density = -481596.989132 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      7.381e-03   1.774e-01    1.000e+00  1.000e+00      4342 -3.700e+03 -3.702e+03                   
#> Path [21] :Best Iter: [66] ELBO (-3700.388656) evaluations: (4342) 
#> Path [22] :Initial log joint density = -481732.313804 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      9.367e-03   2.922e-01    1.000e+00  1.000e+00      2908 -3.708e+03 -3.715e+03                   
#> Path [22] :Best Iter: [52] ELBO (-3708.472846) evaluations: (2908) 
#> Path [23] :Initial log joint density = -481600.932463 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.021e-02   2.033e-01    1.000e+00  1.000e+00      3595 -3.705e+03 -3.714e+03                   
#> Path [23] :Best Iter: [59] ELBO (-3705.119780) evaluations: (3595) 
#> Path [24] :Initial log joint density = -481741.562636 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.048e-02   2.257e-01    1.000e+00  1.000e+00      5252 -3.699e+03 -3.705e+03                   
#> Path [24] :Best Iter: [67] ELBO (-3698.565274) evaluations: (5252) 
#> Path [25] :Initial log joint density = -482576.647533 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.064e-02   2.639e-01    1.000e+00  1.000e+00      4340 -3.699e+03 -3.702e+03                   
#> Path [25] :Best Iter: [66] ELBO (-3699.124450) evaluations: (4340) 
#> Path [26] :Initial log joint density = -481793.636182 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      4.304e-03   2.274e-01    6.857e-01  6.857e-01      4118 -3.700e+03 -3.714e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3699.861519) evaluations: (4118) 
#> Path [27] :Initial log joint density = -481293.862765 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.505e-03   2.188e-01    6.967e-01  6.967e-01      3651 -3.704e+03 -3.707e+03                   
#> Path [27] :Best Iter: [59] ELBO (-3703.934839) evaluations: (3651) 
#> Path [28] :Initial log joint density = -482076.395253 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      2.858e-02   3.453e-01    1.000e+00  1.000e+00      5192 -3.702e+03 -3.709e+03                   
#> Path [28] :Best Iter: [76] ELBO (-3701.627697) evaluations: (5192) 
#> Path [29] :Initial log joint density = -481684.331519 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      9.093e-03   2.340e-01    1.000e+00  1.000e+00      4182 -3.700e+03 -3.704e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3699.602167) evaluations: (4182) 
#> Path [30] :Initial log joint density = -481165.259553 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      5.415e-03   1.235e-01    1.000e+00  1.000e+00      3687 -3.702e+03 -3.715e+03                   
#> Path [30] :Best Iter: [58] ELBO (-3701.801907) evaluations: (3687) 
#> Path [31] :Initial log joint density = -484616.632056 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.946e-02   3.792e-01    9.073e-01  9.073e-01      4384 -3.698e+03 -3.713e+03                   
#> Path [31] :Best Iter: [63] ELBO (-3698.173086) evaluations: (4384) 
#> Path [32] :Initial log joint density = -481650.517452 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.375e-02   1.969e-01    1.000e+00  1.000e+00      5326 -3.701e+03 -3.704e+03                   
#> Path [32] :Best Iter: [73] ELBO (-3700.633272) evaluations: (5326) 
#> Path [33] :Initial log joint density = -481265.372052 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      8.095e-03   1.900e-01    1.000e+00  1.000e+00      3884 -3.702e+03 -3.705e+03                   
#> Path [33] :Best Iter: [60] ELBO (-3701.799219) evaluations: (3884) 
#> Path [34] :Initial log joint density = -481530.305558 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.331e-02   2.088e-01    1.000e+00  1.000e+00      3883 -3.704e+03 -3.702e+03                   
#> Path [34] :Best Iter: [63] ELBO (-3701.937545) evaluations: (3883) 
#> Path [35] :Initial log joint density = -481531.216846 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.915e-02   1.367e-01    1.000e+00  1.000e+00      4228 -3.705e+03 -3.705e+03                   
#> Path [35] :Best Iter: [66] ELBO (-3705.127616) evaluations: (4228) 
#> Path [36] :Initial log joint density = -481433.794507 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      1.579e-02   2.338e-01    1.000e+00  1.000e+00      4035 -3.701e+03 -3.702e+03                   
#> Path [36] :Best Iter: [62] ELBO (-3701.053596) evaluations: (4035) 
#> Path [37] :Initial log joint density = -487749.994728 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.496e-02   2.230e-01    1.000e+00  1.000e+00      5091 -3.701e+03 -3.703e+03                   
#> Path [37] :Best Iter: [74] ELBO (-3700.810111) evaluations: (5091) 
#> Path [38] :Initial log joint density = -481531.890069 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      2.750e-02   1.719e-01    1.000e+00  1.000e+00      4041 -3.704e+03 -3.702e+03                   
#> Path [38] :Best Iter: [65] ELBO (-3702.184669) evaluations: (4041) 
#> Path [39] :Initial log joint density = -481614.596994 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.627e-02   2.037e-01    1.000e+00  1.000e+00      5044 -3.699e+03 -3.702e+03                   
#> Path [39] :Best Iter: [71] ELBO (-3698.691803) evaluations: (5044) 
#> Path [40] :Initial log joint density = -481637.025310 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      2.312e-02   2.999e-01    1.000e+00  1.000e+00      5047 -3.697e+03 -3.702e+03                   
#> Path [40] :Best Iter: [73] ELBO (-3697.269849) evaluations: (5047) 
#> Path [41] :Initial log joint density = -481698.564667 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      6.366e-03   1.948e-01    1.000e+00  1.000e+00      3333 -3.703e+03 -3.716e+03                   
#> Path [41] :Best Iter: [55] ELBO (-3702.829614) evaluations: (3333) 
#> Path [42] :Initial log joint density = -481376.078452 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      2.844e-03   3.222e-01    6.418e-01  6.418e-01      3333 -3.706e+03 -3.725e+03                   
#> Path [42] :Best Iter: [56] ELBO (-3706.010786) evaluations: (3333) 
#> Path [43] :Initial log joint density = -481795.239603 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.505e-02   2.196e-01    1.000e+00  1.000e+00      5404 -3.702e+03 -3.709e+03                   
#> Path [43] :Best Iter: [77] ELBO (-3701.902763) evaluations: (5404) 
#> Path [44] :Initial log joint density = -481735.358719 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      2.402e-02   1.977e-01    1.000e+00  1.000e+00      4703 -3.700e+03 -3.702e+03                   
#> Path [44] :Best Iter: [70] ELBO (-3699.640812) evaluations: (4703) 
#> Path [45] :Initial log joint density = -481833.474628 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      9.322e-03   1.640e-01    1.000e+00  1.000e+00      4342 -3.699e+03 -3.708e+03                   
#> Path [45] :Best Iter: [66] ELBO (-3698.971854) evaluations: (4342) 
#> Path [46] :Initial log joint density = -481821.774364 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.307e-02   1.838e-01    1.000e+00  1.000e+00      5110 -3.700e+03 -3.700e+03                   
#> Path [46] :Best Iter: [72] ELBO (-3700.394674) evaluations: (5110) 
#> Path [47] :Initial log joint density = -481412.705341 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.871e-02   2.595e-01    1.000e+00  1.000e+00      3809 -3.699e+03 -3.709e+03                   
#> Path [47] :Best Iter: [62] ELBO (-3698.807887) evaluations: (3809) 
#> Path [48] :Initial log joint density = -481464.929036 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.357e-03   1.754e-01    1.000e+00  1.000e+00      4170 -3.702e+03 -3.714e+03                   
#> Path [48] :Best Iter: [62] ELBO (-3701.585053) evaluations: (4170) 
#> Path [49] :Initial log joint density = -481924.407074 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.801e-02   4.097e-01    1.000e+00  1.000e+00      5325 -3.699e+03 -3.711e+03                   
#> Path [49] :Best Iter: [73] ELBO (-3698.987967) evaluations: (5325) 
#> Path [50] :Initial log joint density = -481767.549461 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.441e-02   1.742e-01    1.000e+00  1.000e+00      5552 -3.697e+03 -3.701e+03                   
#> Path [50] :Best Iter: [78] ELBO (-3697.192769) evaluations: (5552) 
#> Finished in  14.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.066 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 4
#>    cell_group proportion_fold_change average_uncertainty statement              
#>    <chr>                       <dbl>               <dbl> <glue>                 
#>  1 B1                          -1.92             0.248   1.9-fold decrease (fro…
#>  2 B2                          -2.01             0.325   2-fold decrease (from …
#>  3 B3                          -1.36             0.113   1.4-fold decrease (fro…
#>  4 BM                          -1.39             0.158   1.4-fold decrease (fro…
#>  5 CD4 1                        1.16             0.0246  1.2-fold increase (fro…
#>  6 CD4 2                        1.43             0.0345  1.4-fold increase (fro…
#>  7 CD4 3                       -2.30             0.329   2.3-fold decrease (fro…
#>  8 CD4 4                        1.01             0.0946  1-fold increase (from …
#>  9 CD4 5                        1.03             0.0636  1-fold increase (from …
#> 10 CD8 1                        1.09             0.00381 1.1-fold increase (fro…
#> # ℹ 26 more rows
# }
```
