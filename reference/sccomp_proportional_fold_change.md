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
#> Path [1] :Initial log joint density = -481888.859592 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      4.322e-03   2.398e-01    8.101e-01  8.101e-01      4426 -3.687e+03 -3.699e+03                   
#> Path [1] :Best Iter: [67] ELBO (-3686.995374) evaluations: (4426) 
#> Path [2] :Initial log joint density = -482125.142323 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      7.167e-03   1.959e-01    9.639e-01  9.639e-01      4845 -3.685e+03 -3.696e+03                   
#> Path [2] :Best Iter: [71] ELBO (-3684.600914) evaluations: (4845) 
#> Path [3] :Initial log joint density = -482422.086871 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.044e-02   2.327e-01    6.992e-01  6.992e-01      4537 -3.685e+03 -3.700e+03                   
#> Path [3] :Best Iter: [64] ELBO (-3685.499406) evaluations: (4537) 
#> Path [4] :Initial log joint density = -481528.222256 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      8.625e-03   2.008e-01    8.161e-01  8.161e-01      4056 -3.685e+03 -3.699e+03                   
#> Path [4] :Best Iter: [57] ELBO (-3685.179739) evaluations: (4056) 
#> Path [5] :Initial log joint density = -481571.658603 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      5.990e-03   1.910e-01    1.000e+00  1.000e+00      4238 -3.686e+03 -3.700e+03                   
#> Path [5] :Best Iter: [61] ELBO (-3686.366055) evaluations: (4238) 
#> Path [6] :Initial log joint density = -482053.042155 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      3.160e-03   3.022e-01    7.332e-01  7.332e-01      3190 -3.694e+03 -3.702e+03                   
#> Path [6] :Best Iter: [53] ELBO (-3693.589225) evaluations: (3190) 
#> Path [7] :Initial log joint density = -481520.308763 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.590e-02   2.061e-01    1.000e+00  1.000e+00      3417 -3.687e+03 -3.694e+03                   
#> Path [7] :Best Iter: [57] ELBO (-3687.012425) evaluations: (3417) 
#> Path [8] :Initial log joint density = -481545.253393 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      1.866e-02   2.332e-01    1.000e+00  1.000e+00      3809 -3.685e+03 -3.685e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3685.111987) evaluations: (3809) 
#> Path [9] :Initial log joint density = -481503.947447 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      5.037e-03   2.154e-01    1.000e+00  1.000e+00      3740 -3.687e+03 -3.694e+03                   
#> Path [9] :Best Iter: [58] ELBO (-3686.728905) evaluations: (3740) 
#> Path [10] :Initial log joint density = -481488.661433 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.324e-02   2.457e-01    1.000e+00  1.000e+00      4391 -3.686e+03 -3.688e+03                   
#> Path [10] :Best Iter: [68] ELBO (-3686.158241) evaluations: (4391) 
#> Path [11] :Initial log joint density = -482311.715382 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      2.204e-03   2.123e-01    6.429e-01  6.429e-01      3917 -3.687e+03 -3.702e+03                   
#> Path [11] :Best Iter: [56] ELBO (-3687.300040) evaluations: (3917) 
#> Path [12] :Initial log joint density = -481431.564359 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.162e-03   2.087e-01    1.000e+00  1.000e+00      3347 -3.692e+03 -3.694e+03                   
#> Path [12] :Best Iter: [55] ELBO (-3691.862129) evaluations: (3347) 
#> Path [13] :Initial log joint density = -481533.936465 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      1.109e-02   2.967e-01    1.000e+00  1.000e+00      3043 -3.694e+03 -3.700e+03                   
#> Path [13] :Best Iter: [51] ELBO (-3694.440855) evaluations: (3043) 
#> Path [14] :Initial log joint density = -481621.797150 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      2.741e-02   2.233e-01    1.000e+00  1.000e+00      4054 -3.685e+03 -3.684e+03                   
#> Path [14] :Best Iter: [66] ELBO (-3683.731388) evaluations: (4054) 
#> Path [15] :Initial log joint density = -481456.398841 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      4.845e-03   2.472e-01    7.376e-01  7.376e-01      4245 -3.683e+03 -3.692e+03                   
#> Path [15] :Best Iter: [66] ELBO (-3682.811729) evaluations: (4245) 
#> Path [16] :Initial log joint density = -481564.803842 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      8.511e-03   1.832e-01    1.000e+00  1.000e+00      3382 -3.688e+03 -3.691e+03                   
#> Path [16] :Best Iter: [56] ELBO (-3687.917595) evaluations: (3382) 
#> Path [17] :Initial log joint density = -481891.940906 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      5.933e-03   1.522e-01    8.150e-01  8.150e-01      3985 -3.687e+03 -3.699e+03                   
#> Path [17] :Best Iter: [63] ELBO (-3686.569382) evaluations: (3985) 
#> Path [18] :Initial log joint density = -481872.875124 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      7.605e-03   1.794e-01    9.611e-01  9.611e-01      4595 -3.686e+03 -3.692e+03                   
#> Path [18] :Best Iter: [66] ELBO (-3685.978281) evaluations: (4595) 
#> Path [19] :Initial log joint density = -481623.388937 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      1.135e-02   3.371e-01    1.000e+00  1.000e+00      3204 -3.688e+03 -3.700e+03                   
#> Path [19] :Best Iter: [55] ELBO (-3687.634886) evaluations: (3204) 
#> Path [20] :Initial log joint density = -485279.277508 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.076e-03   1.449e-01    1.000e+00  1.000e+00      3417 -3.694e+03 -3.691e+03                   
#> Path [20] :Best Iter: [56] ELBO (-3690.770074) evaluations: (3417) 
#> Path [21] :Initial log joint density = -481608.039767 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.059e-02   1.927e-01    8.643e-01  8.643e-01      4442 -3.683e+03 -3.692e+03                   
#> Path [21] :Best Iter: [68] ELBO (-3683.475234) evaluations: (4442) 
#> Path [22] :Initial log joint density = -481749.226609 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      3.222e-03   2.268e-01    7.221e-01  7.221e-01      3288 -3.689e+03 -3.701e+03                   
#> Path [22] :Best Iter: [55] ELBO (-3689.306251) evaluations: (3288) 
#> Path [23] :Initial log joint density = -481614.309883 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.012e-02   2.219e-01    1.000e+00  1.000e+00      3505 -3.684e+03 -3.687e+03                   
#> Path [23] :Best Iter: [55] ELBO (-3684.445921) evaluations: (3505) 
#> Path [24] :Initial log joint density = -481757.266670 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.246e-02   2.437e-01    1.000e+00  1.000e+00      5295 -3.683e+03 -3.682e+03                   
#> Path [24] :Best Iter: [78] ELBO (-3682.259757) evaluations: (5295) 
#> Path [25] :Initial log joint density = -482596.027849 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.154e-03   2.505e-01    1.000e+00  1.000e+00      3559 -3.690e+03 -3.700e+03                   
#> Path [25] :Best Iter: [56] ELBO (-3690.473489) evaluations: (3559) 
#> Path [26] :Initial log joint density = -481811.493952 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.301e-02   2.061e-01    8.110e-01  8.110e-01      4386 -3.685e+03 -3.697e+03                   
#> Path [26] :Best Iter: [59] ELBO (-3684.815902) evaluations: (4386) 
#> Path [27] :Initial log joint density = -481306.801429 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.773e-02   2.307e-01    1.000e+00  1.000e+00      3927 -3.686e+03 -3.689e+03                   
#> Path [27] :Best Iter: [62] ELBO (-3686.409096) evaluations: (3927) 
#> Path [28] :Initial log joint density = -482093.433457 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.066e-02   2.946e-01    1.000e+00  1.000e+00      4336 -3.687e+03 -3.696e+03                   
#> Path [28] :Best Iter: [64] ELBO (-3687.469650) evaluations: (4336) 
#> Path [29] :Initial log joint density = -481697.653901 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      5.226e-03   1.932e-01    9.965e-01  9.965e-01      4508 -3.683e+03 -3.686e+03                   
#> Path [29] :Best Iter: [60] ELBO (-3682.876823) evaluations: (4508) 
#> Path [30] :Initial log joint density = -481180.839333 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      9.517e-03   3.486e-01    4.294e-01  1.000e+00      2914 -3.699e+03 -3.706e+03                   
#> Path [30] :Best Iter: [51] ELBO (-3699.128759) evaluations: (2914) 
#> Path [31] :Initial log joint density = -484633.988703 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.853e-02   1.765e-01    1.000e+00  1.000e+00      4625 -3.682e+03 -3.686e+03                   
#> Path [31] :Best Iter: [68] ELBO (-3682.492185) evaluations: (4625) 
#> Path [32] :Initial log joint density = -481673.785505 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.616e-03   2.165e-01    1.000e+00  1.000e+00      4115 -3.687e+03 -3.692e+03                   
#> Path [32] :Best Iter: [65] ELBO (-3686.686131) evaluations: (4115) 
#> Path [33] :Initial log joint density = -481276.134041 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      3.721e-03   3.106e-01    5.495e-01  5.495e-01      4269 -3.684e+03 -3.700e+03                   
#> Path [33] :Best Iter: [61] ELBO (-3684.123095) evaluations: (4269) 
#> Path [34] :Initial log joint density = -481544.646374 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              53      -4.787e+05      8.353e-03   2.344e-01    1.000e+00  1.000e+00      2959 -3.694e+03 -3.709e+03                   
#> Path [34] :Best Iter: [41] ELBO (-3694.378102) evaluations: (2959) 
#> Path [35] :Initial log joint density = -481547.346591 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      6.884e-03   2.152e-01    1.000e+00  1.000e+00      3927 -3.687e+03 -3.698e+03                   
#> Path [35] :Best Iter: [57] ELBO (-3687.270624) evaluations: (3927) 
#> Path [36] :Initial log joint density = -481452.859991 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.872e-02   3.548e-01    1.000e+00  1.000e+00      3856 -3.687e+03 -3.693e+03                   
#> Path [36] :Best Iter: [57] ELBO (-3686.676151) evaluations: (3856) 
#> Path [37] :Initial log joint density = -487772.080417 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      4.809e-03   2.420e-01    7.570e-01  7.570e-01      4145 -3.687e+03 -3.702e+03                   
#> Path [37] :Best Iter: [58] ELBO (-3686.541435) evaluations: (4145) 
#> Path [38] :Initial log joint density = -481550.323440 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              54      -4.787e+05      7.063e-03   2.652e-01    1.000e+00  1.000e+00      3041 -3.694e+03 -3.699e+03                   
#> Path [38] :Best Iter: [51] ELBO (-3693.721500) evaluations: (3041) 
#> Path [39] :Initial log joint density = -481630.374188 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.242e-02   2.461e-01    1.000e+00  1.000e+00      4562 -3.688e+03 -3.688e+03                   
#> Path [39] :Best Iter: [70] ELBO (-3688.298817) evaluations: (4562) 
#> Path [40] :Initial log joint density = -481653.292158 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.944e-02   1.951e-01    1.000e+00  1.000e+00      4279 -3.686e+03 -3.688e+03                   
#> Path [40] :Best Iter: [63] ELBO (-3686.274842) evaluations: (4279) 
#> Path [41] :Initial log joint density = -481714.879484 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      9.164e-03   3.118e-01    1.000e+00  1.000e+00      2875 -3.694e+03 -3.701e+03                   
#> Path [41] :Best Iter: [50] ELBO (-3694.412243) evaluations: (2875) 
#> Path [42] :Initial log joint density = -481389.212649 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.821e-02   2.072e-01    1.000e+00  1.000e+00      3720 -3.688e+03 -3.693e+03                   
#> Path [42] :Best Iter: [61] ELBO (-3688.112298) evaluations: (3720) 
#> Path [43] :Initial log joint density = -481812.513771 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      5.983e-03   1.739e-01    1.000e+00  1.000e+00      5008 -3.687e+03 -3.692e+03                   
#> Path [43] :Best Iter: [73] ELBO (-3686.523883) evaluations: (5008) 
#> Path [44] :Initial log joint density = -481750.394341 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      7.113e-03   2.594e-01    1.000e+00  1.000e+00      3597 -3.690e+03 -3.701e+03                   
#> Path [44] :Best Iter: [58] ELBO (-3689.545766) evaluations: (3597) 
#> Path [45] :Initial log joint density = -481852.750849 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.066e-02   1.759e-01    1.000e+00  1.000e+00      4884 -3.681e+03 -3.688e+03                   
#> Path [45] :Best Iter: [67] ELBO (-3680.553686) evaluations: (4884) 
#> Path [46] :Initial log joint density = -481837.912059 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.359e-02   2.712e-01    1.000e+00  1.000e+00      4337 -3.685e+03 -3.689e+03                   
#> Path [46] :Best Iter: [69] ELBO (-3685.461642) evaluations: (4337) 
#> Path [47] :Initial log joint density = -481427.630914 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.969e-03   3.242e-01    8.521e-01  8.521e-01      3199 -3.692e+03 -3.703e+03                   
#> Path [47] :Best Iter: [55] ELBO (-3691.625051) evaluations: (3199) 
#> Path [48] :Initial log joint density = -481479.166642 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.707e-02   2.008e-01    1.000e+00  1.000e+00      4512 -3.683e+03 -3.687e+03                   
#> Path [48] :Best Iter: [62] ELBO (-3682.958907) evaluations: (4512) 
#> Path [49] :Initial log joint density = -481941.686372 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.223e-02   1.913e-01    1.000e+00  1.000e+00      4800 -3.685e+03 -3.682e+03                   
#> Path [49] :Best Iter: [73] ELBO (-3682.300354) evaluations: (4800) 
#> Path [50] :Initial log joint density = -481784.593365 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      9.221e-03   2.334e-01    1.000e+00  1.000e+00      4468 -3.688e+03 -3.698e+03                   
#> Path [50] :Best Iter: [67] ELBO (-3688.114000) evaluations: (4468) 
#> Finished in  15.5 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.072 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 4
#>    cell_group proportion_fold_change average_uncertainty statement              
#>    <chr>                       <dbl>               <dbl> <glue>                 
#>  1 B1                          -1.92             0.202   1.9-fold decrease (fro…
#>  2 B2                          -2.01             0.284   2-fold decrease (from …
#>  3 B3                          -1.38             0.163   1.4-fold decrease (fro…
#>  4 BM                          -1.41             0.114   1.4-fold decrease (fro…
#>  5 CD4 1                        1.17             0.00784 1.2-fold increase (fro…
#>  6 CD4 2                        1.42             0.0454  1.4-fold increase (fro…
#>  7 CD4 3                       -2.29             0.329   2.3-fold decrease (fro…
#>  8 CD4 4                       -1.06             0.0711  1.1-fold decrease (fro…
#>  9 CD4 5                        1.04             0.0972  1-fold increase (from …
#> 10 CD8 1                        1.09             0.00157 1.1-fold increase (fro…
#> # ℹ 26 more rows
# }
```
