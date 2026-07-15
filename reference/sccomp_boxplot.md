# sccomp_boxplot

Creates a boxplot visualization of the model results from `sccomp`. This
function plots the estimated cell proportions across samples,
highlighting significant changes in cell composition according to a
specified factor.

## Usage

``` r
sccomp_boxplot(
  .data,
  factor,
  significance_threshold = 0.05,
  significance_statistic = c("pH0", "FDR"),
  test_composition_above_logit_fold_change = attr(.data,
    "test_composition_above_logit_fold_change"),
  remove_unwanted_effects = FALSE,
  cache_stan_model = sccomp_stan_models_cache_dir
)
```

## Arguments

- .data:

  A tibble containing the results from `sccomp_estimate` and
  `sccomp_test`, including the columns: cell_group name, sample name,
  read counts, factor(s), p-values, and significance indicators.

- factor:

  A character string specifying the factor of interest included in the
  model for stratifying the boxplot.

- significance_threshold:

  A numeric value indicating the threshold for labeling significant
  cell-groups. Defaults to 0.05.

- significance_statistic:

  Character vector indicating which statistic is used to colour
  significant groups. Defaults to `c("pH0", "FDR")`.

- test_composition_above_logit_fold_change:

  A positive numeric value representing the effect size threshold used
  in the hypothesis test. A value of 0.2 corresponds to a change in cell
  proportion of approximately 10% for a cell type with a baseline
  proportion of 50% (e.g., from 45% to 55%). This threshold is
  consistent on the logit-unconstrained scale, even when the baseline
  proportion is close to 0 or 1.

- remove_unwanted_effects:

  A logical value indicating whether to remove unwanted variation from
  the data before plotting. Defaults to `FALSE`.

- cache_stan_model:

  A character string specifying the cache directory for compiled Stan
  models. Default is `sccomp_stan_models_cache_dir` which points to
  `~/.sccomp_models`. Use a custom path in restricted environments where
  the default is not writable.

## Value

A `ggplot` object representing the boxplot of cell proportions across
samples, stratified by the specified factor.

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## Examples

``` r

print("cmdstanr is needed to run this example.")
#> [1] "cmdstanr is needed to run this example."
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# \donttest{
if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate <- sccomp_estimate(
      counts_obj,
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1
    ) |>
    sccomp_test()

    # Plot the boxplot of estimated cell proportions
    sccomp_boxplot(
        .data = estimate,
        factor = "type",
        significance_threshold = 0.05
    )
}
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: the composition design matrix has columns: (Intercept), typecancer
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -482884.236146 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      7.922e-03   3.287e-01    1.000e+00  1.000e+00      3254 -3.690e+03 -3.695e+03                   
#> Path [1] :Best Iter: [55] ELBO (-3689.601926) evaluations: (3254) 
#> Path [2] :Initial log joint density = -482282.675312 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.022e-02   2.226e-01    1.000e+00  1.000e+00      4558 -3.684e+03 -3.687e+03                   
#> Path [2] :Best Iter: [67] ELBO (-3683.849832) evaluations: (4558) 
#> Path [3] :Initial log joint density = -481471.368844 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.074e-02   2.179e-01    1.000e+00  1.000e+00      3443 -3.688e+03 -3.690e+03                   
#> Path [3] :Best Iter: [56] ELBO (-3687.921794) evaluations: (3443) 
#> Path [4] :Initial log joint density = -484857.970242 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      5.676e-02   2.321e-01    1.000e+00  1.000e+00      4573 -3.686e+03 -3.686e+03                   
#> Path [4] :Best Iter: [66] ELBO (-3685.643617) evaluations: (4573) 
#> Path [5] :Initial log joint density = -481559.670007 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.740e-03   1.993e-01    9.291e-01  9.291e-01      3186 -3.691e+03 -3.702e+03                   
#> Path [5] :Best Iter: [55] ELBO (-3690.776805) evaluations: (3186) 
#> Path [6] :Initial log joint density = -482875.868619 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      9.735e-03   2.239e-01    9.331e-01  9.331e-01      3575 -3.691e+03 -3.699e+03                   
#> Path [6] :Best Iter: [59] ELBO (-3690.919142) evaluations: (3575) 
#> Path [7] :Initial log joint density = -481603.238451 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      6.452e-03   2.009e-01    1.000e+00  1.000e+00      3777 -3.686e+03 -3.694e+03                   
#> Path [7] :Best Iter: [55] ELBO (-3686.365684) evaluations: (3777) 
#> Path [8] :Initial log joint density = -481287.533598 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              64      -4.787e+05      2.204e-02   2.027e-01    1.000e+00  1.000e+00      3918 -3.683e+03 -3.685e+03                   
#> Path [8] :Best Iter: [61] ELBO (-3683.460499) evaluations: (3918) 
#> Path [9] :Initial log joint density = -481364.139605 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      8.587e-03   1.769e-01    1.000e+00  1.000e+00      4401 -3.685e+03 -3.690e+03                   
#> Path [9] :Best Iter: [65] ELBO (-3684.729361) evaluations: (4401) 
#> Path [10] :Initial log joint density = -481339.631442 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              67      -4.787e+05      1.519e-02   2.650e-01    1.000e+00  1.000e+00      4088 -3.686e+03 -3.690e+03                   
#> Path [10] :Best Iter: [66] ELBO (-3686.333369) evaluations: (4088) 
#> Path [11] :Initial log joint density = -482265.411111 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      2.253e-02   3.300e-01    1.000e+00  1.000e+00      5222 -3.683e+03 -3.689e+03                   
#> Path [11] :Best Iter: [74] ELBO (-3682.621924) evaluations: (5222) 
#> Path [12] :Initial log joint density = -481579.684868 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      3.408e-03   1.712e-01    7.450e-01  7.450e-01      4500 -3.685e+03 -3.698e+03                   
#> Path [12] :Best Iter: [64] ELBO (-3684.934495) evaluations: (4500) 
#> Path [13] :Initial log joint density = -481811.078367 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.084e-02   1.200e-01    1.000e+00  1.000e+00      4742 -3.684e+03 -3.685e+03                   
#> Path [13] :Best Iter: [69] ELBO (-3684.239036) evaluations: (4742) 
#> Path [14] :Initial log joint density = -481409.245481 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      4.047e-03   2.196e-01    6.787e-01  6.787e-01      4279 -3.683e+03 -3.696e+03                   
#> Path [14] :Best Iter: [67] ELBO (-3682.642243) evaluations: (4279) 
#> Path [15] :Initial log joint density = -481737.159081 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      8.336e-03   2.861e-01    1.000e+00  1.000e+00      4688 -3.684e+03 -3.694e+03                   
#> Path [15] :Best Iter: [71] ELBO (-3683.802017) evaluations: (4688) 
#> Path [16] :Initial log joint density = -481392.065927 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              58      -4.787e+05      1.906e-03   2.613e-01    5.603e-01  5.603e-01      3431 -3.692e+03 -3.703e+03                   
#> Path [16] :Best Iter: [55] ELBO (-3692.175287) evaluations: (3431) 
#> Path [17] :Initial log joint density = -487499.539625 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      9.098e-03   2.358e-01    1.000e+00  1.000e+00      3702 -3.688e+03 -3.686e+03                   
#> Path [17] :Best Iter: [62] ELBO (-3686.286919) evaluations: (3702) 
#> Path [18] :Initial log joint density = -481505.418671 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      3.819e-03   1.610e-01    4.361e-01  1.000e+00      3225 -3.693e+03 -3.702e+03                   
#> Path [18] :Best Iter: [54] ELBO (-3693.237307) evaluations: (3225) 
#> Path [19] :Initial log joint density = -482650.530727 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      1.798e-02   2.041e-01    1.000e+00  1.000e+00      4636 -3.685e+03 -3.685e+03                   
#> Path [19] :Best Iter: [72] ELBO (-3684.975491) evaluations: (4636) 
#> Path [20] :Initial log joint density = -481679.054999 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      9.330e-03   1.912e-01    1.000e+00  1.000e+00      3298 -3.693e+03 -3.690e+03                   
#> Path [20] :Best Iter: [57] ELBO (-3690.108417) evaluations: (3298) 
#> Path [21] :Initial log joint density = -484312.787875 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      8.246e-03   2.295e-01    8.669e-01  8.669e-01      3617 -3.690e+03 -3.699e+03                   
#> Path [21] :Best Iter: [58] ELBO (-3689.542184) evaluations: (3617) 
#> Path [22] :Initial log joint density = -482045.534947 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      8.834e-03   2.829e-01    9.021e-01  9.021e-01      4268 -3.690e+03 -3.703e+03                   
#> Path [22] :Best Iter: [66] ELBO (-3689.818691) evaluations: (4268) 
#> Path [23] :Initial log joint density = -481611.229887 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      5.940e-03   1.448e-01    1.000e+00  1.000e+00      3263 -3.693e+03 -3.699e+03                   
#> Path [23] :Best Iter: [43] ELBO (-3693.367850) evaluations: (3263) 
#> Path [24] :Initial log joint density = -483904.083265 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      4.704e-03   3.010e-01    5.961e-01  5.961e-01      3960 -3.684e+03 -3.700e+03                   
#> Path [24] :Best Iter: [63] ELBO (-3683.582303) evaluations: (3960) 
#> Path [25] :Initial log joint density = -483720.299957 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      4.782e-03   2.145e-01    6.630e-01  6.630e-01      3489 -3.690e+03 -3.698e+03                   
#> Path [25] :Best Iter: [57] ELBO (-3690.147023) evaluations: (3489) 
#> Path [26] :Initial log joint density = -481569.938105 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      8.381e-03   2.590e-01    1.000e+00  1.000e+00      4427 -3.688e+03 -3.698e+03                   
#> Path [26] :Best Iter: [58] ELBO (-3688.045850) evaluations: (4427) 
#> Path [27] :Initial log joint density = -482866.075748 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              57      -4.787e+05      6.972e-03   2.661e-01    1.000e+00  1.000e+00      3260 -3.697e+03 -3.694e+03                   
#> Path [27] :Best Iter: [57] ELBO (-3693.601859) evaluations: (3260) 
#> Path [28] :Initial log joint density = -482425.128699 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      3.957e-03   2.547e-01    6.801e-01  6.801e-01      4529 -3.685e+03 -3.703e+03                   
#> Path [28] :Best Iter: [70] ELBO (-3685.392793) evaluations: (4529) 
#> Path [29] :Initial log joint density = -482465.155968 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              52      -4.787e+05      8.009e-03   3.243e-01    8.884e-01  8.884e-01      2903 -3.696e+03 -3.707e+03                   
#> Path [29] :Best Iter: [49] ELBO (-3696.208484) evaluations: (2903) 
#> Path [30] :Initial log joint density = -483498.254775 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      2.529e-03   2.233e-01    7.029e-01  7.029e-01      3569 -3.688e+03 -3.701e+03                   
#> Path [30] :Best Iter: [59] ELBO (-3687.901567) evaluations: (3569) 
#> Path [31] :Initial log joint density = -481608.945017 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      8.084e-03   1.845e-01    1.000e+00  1.000e+00      4470 -3.684e+03 -3.690e+03                   
#> Path [31] :Best Iter: [66] ELBO (-3683.603788) evaluations: (4470) 
#> Path [32] :Initial log joint density = -481467.770427 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              59      -4.787e+05      1.394e-02   2.734e-01    1.000e+00  1.000e+00      3319 -3.688e+03 -3.693e+03                   
#> Path [32] :Best Iter: [58] ELBO (-3688.316046) evaluations: (3319) 
#> Path [33] :Initial log joint density = -481277.584419 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              55      -4.787e+05      7.390e-03   1.950e-01    1.000e+00  1.000e+00      3026 -3.695e+03 -3.690e+03                   
#> Path [33] :Best Iter: [55] ELBO (-3689.650986) evaluations: (3026) 
#> Path [34] :Initial log joint density = -481667.001464 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      1.092e-02   1.728e-01    1.000e+00  1.000e+00      4242 -3.684e+03 -3.684e+03                   
#> Path [34] :Best Iter: [65] ELBO (-3684.159634) evaluations: (4242) 
#> Path [35] :Initial log joint density = -481796.634549 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      1.376e-02   3.035e-01    1.000e+00  1.000e+00      3666 -3.686e+03 -3.692e+03                   
#> Path [35] :Best Iter: [55] ELBO (-3686.459477) evaluations: (3666) 
#> Path [36] :Initial log joint density = -482015.028319 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.075e-02   1.899e-01    1.000e+00  1.000e+00      4473 -3.684e+03 -3.690e+03                   
#> Path [36] :Best Iter: [68] ELBO (-3684.120837) evaluations: (4473) 
#> Path [37] :Initial log joint density = -481110.385164 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      6.788e-03   1.839e-01    1.000e+00  1.000e+00      4308 -3.686e+03 -3.698e+03                   
#> Path [37] :Best Iter: [64] ELBO (-3685.922880) evaluations: (4308) 
#> Path [38] :Initial log joint density = -481408.365397 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      1.630e-02   2.011e-01    1.000e+00  1.000e+00      4565 -3.683e+03 -3.684e+03                   
#> Path [38] :Best Iter: [68] ELBO (-3682.597434) evaluations: (4565) 
#> Path [39] :Initial log joint density = -481635.684920 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.278e-03   2.814e-01    4.074e-01  1.000e+00      3155 -3.691e+03 -3.700e+03                   
#> Path [39] :Best Iter: [55] ELBO (-3690.814600) evaluations: (3155) 
#> Path [40] :Initial log joint density = -481425.963944 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      1.551e-02   2.072e-01    1.000e+00  1.000e+00      4175 -3.686e+03 -3.692e+03                   
#> Path [40] :Best Iter: [64] ELBO (-3685.546962) evaluations: (4175) 
#> Path [41] :Initial log joint density = -481664.384315 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      1.229e-02   1.502e-01    1.000e+00  1.000e+00      3451 -3.688e+03 -3.691e+03                   
#> Path [41] :Best Iter: [56] ELBO (-3688.191182) evaluations: (3451) 
#> Path [42] :Initial log joint density = -481639.856277 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              62      -4.787e+05      1.139e-02   2.751e-01    1.000e+00  1.000e+00      3923 -3.688e+03 -3.699e+03                   
#> Path [42] :Best Iter: [58] ELBO (-3688.378484) evaluations: (3923) 
#> Path [43] :Initial log joint density = -481470.165009 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.485e-02   2.558e-01    1.000e+00  1.000e+00      5130 -3.682e+03 -3.689e+03                   
#> Path [43] :Best Iter: [73] ELBO (-3681.672354) evaluations: (5130) 
#> Path [44] :Initial log joint density = -481389.860602 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              60      -4.787e+05      4.632e-03   1.516e-01    1.000e+00  1.000e+00      3483 -3.689e+03 -3.699e+03                   
#> Path [44] :Best Iter: [56] ELBO (-3689.197240) evaluations: (3483) 
#> Path [45] :Initial log joint density = -481411.096970 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              63      -4.787e+05      7.235e-03   2.194e-01    8.280e-01  8.280e-01      3740 -3.689e+03 -3.699e+03                   
#> Path [45] :Best Iter: [60] ELBO (-3688.645805) evaluations: (3740) 
#> Path [46] :Initial log joint density = -481760.091704 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              56      -4.787e+05      7.269e-03   2.154e-01    1.000e+00  1.000e+00      3109 -3.696e+03 -3.696e+03                   
#> Path [46] :Best Iter: [48] ELBO (-3695.864511) evaluations: (3109) 
#> Path [47] :Initial log joint density = -482797.368267 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.652e-02   2.094e-01    1.000e+00  1.000e+00      4521 -3.683e+03 -3.685e+03                   
#> Path [47] :Best Iter: [68] ELBO (-3683.368846) evaluations: (4521) 
#> Path [48] :Initial log joint density = -481650.636114 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.972e-03   2.582e-01    6.606e-01  6.606e-01      4024 -3.686e+03 -3.698e+03                   
#> Path [48] :Best Iter: [62] ELBO (-3686.159314) evaluations: (4024) 
#> Path [49] :Initial log joint density = -482801.782611 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              61      -4.787e+05      7.187e-03   2.655e-01    7.732e-01  7.732e-01      3592 -3.686e+03 -3.703e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3685.969914) evaluations: (3592) 
#> Path [50] :Initial log joint density = -481439.950692 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              66      -4.787e+05      3.656e-03   2.647e-01    6.911e-01  6.911e-01      4140 -3.683e+03 -3.703e+03                   
#> Path [50] :Best Iter: [63] ELBO (-3682.656652) evaluations: (4140) 
#> Finished in  15.3 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Joining with `by = join_by(cell_group, M, parameter)`
#> sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.
#> sccomp says: from version 2.1.25, the default `significance_statistic` for boxplots is `pH0` (previously `FDR`). Set `significance_statistic = "FDR"` to use the previous default.
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.111 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> Joining with `by = join_by(cell_group, sample)`
#> Joining with `by = join_by(cell_group, type)`
#> Warning: Ignoring unknown parameters: `median.linewidth`

# }
```
