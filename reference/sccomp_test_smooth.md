# Test differences along a smooth term

Tests whether the mean cell-group composition differs between two
positions (or two intervals) of a continuous covariate modelled with a
smooth `s()` term. The test is a contrast of the model's mean prediction
on the centered log-ratio scale, evaluated at the two positions and
summarised with the same probability-of-null (`pH0`) and
false-discovery-rate (`FDR`) machinery as
[`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md).

## Usage

``` r
sccomp_test_smooth(
  fit,
  smooth = NULL,
  from,
  to,
  comparison = c("endpoints", "average"),
  at = NULL,
  resolution = 20L,
  test_composition_above_CLR_fold_change = 0.1,
  percent_false_positive = 5,
  number_of_draws = 500,
  mcmc_seed = sample_seed(),
  robust = FALSE
)
```

## Arguments

- fit:

  The result of
  [`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
  with a smooth term in `formula_composition`.

- smooth:

  Character. The continuous covariate inside the smooth to test (e.g.
  `"pseudotime"`). Can be `NULL` when the model has exactly one smooth
  term over one continuous covariate, in which case it is auto-detected.

- from:

  Numeric. The reference position. A single value for
  `comparison = "endpoints"`, or a length-2 interval `c(lo, hi)` for
  `comparison = "average"`.

- to:

  Numeric. The comparison position, same shape as `from`.

- comparison:

  One of `"endpoints"` or `"average"`. See details.

- at:

  Optional named list fixing the value of other covariates (e.g. the
  grouping factor of a factor smooth). Covariates not supplied are held
  at their first observed value (which cancels for non-grouping
  covariates).

- resolution:

  Integer. Number of grid points used to approximate each interval
  average when `comparison = "average"`.

- test_composition_above_CLR_fold_change:

  Positive numeric. Effect threshold for the hypothesis test, on the
  centered log-ratio scale — identical meaning to the argument of
  [`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md).

- percent_false_positive:

  Numeric in (0, 100). Used for the credible interval width, as in
  [`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md).

- number_of_draws:

  Integer. Number of posterior draws used for the prediction pass.

- mcmc_seed:

  Integer. Seed for the generated-quantities pass.

- robust:

  Logical. Currently unused placeholder for API parity with
  [`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md);
  the effect is always the posterior mean of the contrast.

## Value

A tibble with one row per cell group:

- `cell_group` — the cell group tested.

- `smooth` — the continuous covariate tested.

- `from`, `to` — the compared positions (as labels).

- `c_lower`, `c_effect`, `c_upper` — 95% CI and posterior mean of the
  centered log-ratio-scale contrast \\\mu(\code{to}) -
  \mu(\code{from})\\.

- `c_pH0` — probability the effect is within the null region.

- `c_FDR` — false-discovery rate across cell groups.

## Details

Because the smooth curve can be non-linear, "difference across an
interval" is ambiguous. Two comparison modes are provided:

- `"endpoints"` (default): `from` and `to` are single values; the
  contrast is \\\mu(\code{to}) - \mu(\code{from})\\.

- `"average"`: `from` and `to` are length-2 intervals `c(lo, hi)`; the
  contrast is the average of the curve over the `to` interval minus the
  average over the `from` interval (each interval sampled at
  `resolution` points).

Any covariate not being varied cancels out of the contrast because the
centered log-ratio-scale predictor is additive; it therefore does not
matter what value those covariates take, as long as it is held constant.
For factor smooths (`s(x, g, bs = "fs")`) or by-factor smooths
(`s(x, by = g)`) the grouping factor selects *which* curve is tested;
set it via `at`, e.g. `at = list(tissue = "tumor")`.

The test is about the **mean** composition: it uses the expected linear
predictor (`mu_unconstrained`), not the overdispersed beta-binomial
realisation. It therefore answers "does the expected composition differ
between these covariate positions?", not "would an individual sample
differ?".

## Examples

``` r
# \donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")
    # add a continuous covariate
    counts_obj$pseudotime <- as.numeric(factor(counts_obj$sample))

    fit <- sccomp_estimate(
      counts_obj,
      ~ s(pseudotime, k = 4), ~ 1, "sample", "cell_group", "count",
      cores = 1
    )

    fit |> sccomp_test_smooth(from = 2, to = 8)
  }
#> sccomp says: count column is an integer. The sum-constrained beta binomial model will be used
#> sccomp says: estimation
#> sccomp says: model design matrices
#>   composition X: 2 parameters: (Intercept), s(pseudotime, k = 4)__lin1
#>   variability Xa: 1 parameter: (Intercept)
#>   random effect X1 [smooth s(pseudotime, k = 4)]: 2 parameters: s(pseudotime, k = 4)___basis01, s(pseudotime, k = 4)___basis02
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481770.094907 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      2.054e-02   2.258e-01    1.000e+00  1.000e+00      4573 -3.873e+03 -3.875e+03                   
#> Path [1] :Best Iter: [70] ELBO (-3872.993101) evaluations: (4573) 
#> Path [2] :Initial log joint density = -481540.360428 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              92      -4.787e+05      1.074e-02   2.192e-01    1.000e+00  1.000e+00      6931 -3.869e+03 -3.876e+03                   
#> Path [2] :Best Iter: [40] ELBO (-3869.425960) evaluations: (6931) 
#> Path [3] :Initial log joint density = -481084.680629 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.012e-02   1.573e-01    1.000e+00  1.000e+00      6356 -3.878e+03 -4.005e+03                   
#> Path [3] :Best Iter: [76] ELBO (-3877.888093) evaluations: (6356) 
#> Path [4] :Initial log joint density = -481793.311282 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.787e+05      7.772e-03   1.755e-01    8.844e-01  8.844e-01      7446 -3.882e+03 -3.895e+03                   
#> Path [4] :Best Iter: [61] ELBO (-3881.949534) evaluations: (7446) 
#> Path [5] :Initial log joint density = -482050.031749 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      1.116e-02   1.725e-01    8.984e-01  8.984e-01      6257 -3.868e+03 -4.009e+03                   
#> Path [5] :Best Iter: [41] ELBO (-3867.645198) evaluations: (6257) 
#> Path [6] :Initial log joint density = -481632.737910 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.479e-02   1.601e-01    1.000e+00  1.000e+00      5951 -3.862e+03 -3.884e+03                   
#> Path [6] :Best Iter: [30] ELBO (-3862.355336) evaluations: (5951) 
#> Path [7] :Initial log joint density = -481340.397406 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      2.628e-02   3.452e-01    1.000e+00  1.000e+00      5371 -3.871e+03 -4.064e+03                   
#> Path [7] :Best Iter: [36] ELBO (-3871.342280) evaluations: (5371) 
#> Path [8] :Initial log joint density = -481276.236901 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.787e+05      1.977e-02   2.092e-01    8.276e-01  8.276e-01      6959 -3.875e+03 -3.926e+03                   
#> Path [8] :Best Iter: [65] ELBO (-3874.662637) evaluations: (6959) 
#> Path [9] :Initial log joint density = -482020.072409 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      8.281e-03   2.251e-01    1.000e+00  1.000e+00      4918 -3.869e+03 -3.888e+03                   
#> Path [9] :Best Iter: [58] ELBO (-3869.121160) evaluations: (4918) 
#> Path [10] :Initial log joint density = -481252.080528 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      1.288e-02   1.409e-01    1.000e+00  1.000e+00      7668 -3.864e+03 -3.914e+03                   
#> Path [10] :Best Iter: [41] ELBO (-3864.486207) evaluations: (7668) 
#> Path [11] :Initial log joint density = -481583.895304 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              90      -4.787e+05      2.004e-02   2.067e-01    1.000e+00  1.000e+00      6664 -3.872e+03 -3.949e+03                   
#> Path [11] :Best Iter: [89] ELBO (-3871.842637) evaluations: (6664) 
#> Path [12] :Initial log joint density = -481456.838687 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.787e+05      1.801e-02   2.734e-01    1.000e+00  1.000e+00      7398 -3.870e+03 -3.903e+03                   
#> Path [12] :Best Iter: [70] ELBO (-3870.204932) evaluations: (7398) 
#> Path [13] :Initial log joint density = -481308.376325 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.787e+05      1.953e-02   1.962e-01    1.000e+00  1.000e+00      6580 -3.876e+03 -3.922e+03                   
#> Path [13] :Best Iter: [46] ELBO (-3875.800272) evaluations: (6580) 
#> Path [14] :Initial log joint density = -481576.494468 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              95      -4.787e+05      1.428e-02   1.438e-01    1.000e+00  1.000e+00      7408 -3.867e+03 -3.962e+03                   
#> Path [14] :Best Iter: [61] ELBO (-3866.806490) evaluations: (7408) 
#> Path [15] :Initial log joint density = -481573.696405 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      1.682e-02   1.521e-01    1.000e+00  1.000e+00      7714 -3.870e+03 -3.883e+03                   
#> Path [15] :Best Iter: [88] ELBO (-3870.197460) evaluations: (7714) 
#> Path [16] :Initial log joint density = -481526.136498 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              99      -4.787e+05      1.982e-02   2.216e-01    1.000e+00  1.000e+00      7617 -3.871e+03 -4.112e+03                   
#> Path [16] :Best Iter: [65] ELBO (-3871.481984) evaluations: (7617) 
#> Path [17] :Initial log joint density = -481540.184543 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              96      -4.787e+05      1.448e-02   1.877e-01    1.000e+00  1.000e+00      7280 -3.876e+03 -3.978e+03                   
#> Path [17] :Best Iter: [89] ELBO (-3875.518754) evaluations: (7280) 
#> Path [18] :Initial log joint density = -482122.344826 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              85      -4.787e+05      1.004e-02   1.901e-01    1.000e+00  1.000e+00      6090 -3.871e+03 -3.904e+03                   
#> Path [18] :Best Iter: [68] ELBO (-3870.834123) evaluations: (6090) 
#> Path [19] :Initial log joint density = -481847.424407 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.285e-02   2.938e-01    1.000e+00  1.000e+00      7823 -3.865e+03 -3.887e+03                   
#> Path [19] :Best Iter: [44] ELBO (-3864.845342) evaluations: (7823) 
#> Path [20] :Initial log joint density = -481973.312474 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      3.867e-02   4.090e-01    1.000e+00  1.000e+00      7780 -3.869e+03 -4.063e+03                   
#> Path [20] :Best Iter: [69] ELBO (-3868.894407) evaluations: (7780) 
#> Path [21] :Initial log joint density = -482642.371028 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              81      -4.787e+05      7.510e-03   2.518e-01    1.000e+00  1.000e+00      5738 -3.868e+03 -4.362e+03                   
#> Path [21] :Best Iter: [64] ELBO (-3868.264054) evaluations: (5738) 
#> Path [22] :Initial log joint density = -481664.672900 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              96      -4.787e+05      1.894e-02   1.662e-01    1.000e+00  1.000e+00      7527 -3.872e+03 -3.895e+03                   
#> Path [22] :Best Iter: [73] ELBO (-3872.280541) evaluations: (7527) 
#> Path [23] :Initial log joint density = -482404.092208 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              82      -4.787e+05      1.454e-02   2.204e-01    1.000e+00  1.000e+00      5839 -3.873e+03 -3.911e+03                   
#> Path [23] :Best Iter: [39] ELBO (-3872.758442) evaluations: (5839) 
#> Path [24] :Initial log joint density = -481686.423754 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      7.870e-03   1.810e-01    1.000e+00  1.000e+00      5540 -3.873e+03 -3.894e+03                   
#> Path [24] :Best Iter: [71] ELBO (-3872.758781) evaluations: (5540) 
#> Path [25] :Initial log joint density = -483489.585345 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              92      -4.787e+05      4.727e-03   2.464e-01    5.549e-01  5.549e-01      6962 -3.879e+03 -3.937e+03                   
#> Path [25] :Best Iter: [67] ELBO (-3878.563690) evaluations: (6962) 
#> Path [26] :Initial log joint density = -481590.113512 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              91      -4.787e+05      1.344e-02   2.173e-01    1.000e+00  1.000e+00      6824 -3.865e+03 -4.171e+03                   
#> Path [26] :Best Iter: [67] ELBO (-3864.835414) evaluations: (6824) 
#> Path [27] :Initial log joint density = -481587.207906 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      1.275e-02   2.245e-01    1.000e+00  1.000e+00      4740 -3.870e+03 -3.883e+03                   
#> Path [27] :Best Iter: [67] ELBO (-3869.578644) evaluations: (4740) 
#> Path [28] :Initial log joint density = -481606.599004 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      3.382e-02   3.516e-01    1.000e+00  1.000e+00      7567 -3.869e+03 -3.911e+03                   
#> Path [28] :Best Iter: [80] ELBO (-3869.086553) evaluations: (7567) 
#> Path [29] :Initial log joint density = -481404.004604 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.787e+05      1.003e-02   1.657e-01    9.294e-01  9.294e-01      6956 -3.866e+03 -3.908e+03                   
#> Path [29] :Best Iter: [42] ELBO (-3866.009566) evaluations: (6956) 
#> Path [30] :Initial log joint density = -481335.273935 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.349e-02   2.656e-01    1.000e+00  1.000e+00      4977 -3.878e+03 -3.883e+03                   
#> Path [30] :Best Iter: [50] ELBO (-3878.031119) evaluations: (4977) 
#> Path [31] :Initial log joint density = -482133.925809 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              94      -4.787e+05      1.782e-02   1.732e-01    1.000e+00  1.000e+00      7015 -3.865e+03 -3.947e+03                   
#> Path [31] :Best Iter: [60] ELBO (-3865.086228) evaluations: (7015) 
#> Path [32] :Initial log joint density = -482674.451020 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      1.476e-02   2.028e-01    1.000e+00  1.000e+00      5542 -3.871e+03 -3.995e+03                   
#> Path [32] :Best Iter: [61] ELBO (-3870.704026) evaluations: (5542) 
#> Path [33] :Initial log joint density = -481499.756801 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      4.486e-02   5.167e-01    1.000e+00  1.000e+00      7666 -3.869e+03 -3.903e+03                   
#> Path [33] :Best Iter: [88] ELBO (-3868.956056) evaluations: (7666) 
#> Path [34] :Initial log joint density = -481800.450134 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      2.152e-02   2.238e-01    1.000e+00  1.000e+00      6169 -3.871e+03 -3.917e+03                   
#> Path [34] :Best Iter: [45] ELBO (-3871.311053) evaluations: (6169) 
#> Path [35] :Initial log joint density = -481497.119082 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      1.712e-02   1.809e-01    1.000e+00  1.000e+00      7828 -3.861e+03 -3.923e+03                   
#> Path [35] :Best Iter: [97] ELBO (-3861.466379) evaluations: (7828) 
#> Path [36] :Initial log joint density = -481495.319139 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              83      -4.787e+05      3.469e-02   2.284e-01    1.000e+00  1.000e+00      5858 -3.857e+03 -4.024e+03                   
#> Path [36] :Best Iter: [34] ELBO (-3856.916943) evaluations: (5858) 
#> Path [37] :Initial log joint density = -481343.909635 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.572e-02   2.050e-01    1.000e+00  1.000e+00      7862 -3.875e+03 -3.921e+03                   
#> Path [37] :Best Iter: [36] ELBO (-3874.812306) evaluations: (7862) 
#> Path [38] :Initial log joint density = -481904.318785 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      1.841e-02   2.679e-01    1.000e+00  1.000e+00      4957 -3.874e+03 -3.902e+03                   
#> Path [38] :Best Iter: [55] ELBO (-3874.100635) evaluations: (4957) 
#> Path [39] :Initial log joint density = -482154.271703 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              93      -4.787e+05      1.523e-02   2.000e-01    1.000e+00  1.000e+00      7039 -3.876e+03 -3.940e+03                   
#> Path [39] :Best Iter: [65] ELBO (-3876.382658) evaluations: (7039) 
#> Path [40] :Initial log joint density = -481785.787992 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      2.478e-02   2.294e-01    9.739e-01  9.739e-01      7579 -3.867e+03 -3.906e+03                   
#> Path [40] :Best Iter: [87] ELBO (-3867.021162) evaluations: (7579) 
#> Path [41] :Initial log joint density = -481852.303310 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.848e-02   2.977e-01    1.000e+00  1.000e+00      7814 -3.869e+03 -3.946e+03                   
#> Path [41] :Best Iter: [52] ELBO (-3869.420200) evaluations: (7814) 
#> Path [42] :Initial log joint density = -484344.757204 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              96      -4.787e+05      1.802e-02   2.045e-01    1.000e+00  1.000e+00      7805 -3.861e+03 -3.930e+03                   
#> Path [42] :Best Iter: [33] ELBO (-3860.805699) evaluations: (7805) 
#> Path [43] :Initial log joint density = -481587.299968 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              84      -4.787e+05      3.558e-03   1.729e-01    7.638e-01  7.638e-01      5902 -3.866e+03 -4.173e+03                   
#> Path [43] :Best Iter: [72] ELBO (-3865.745730) evaluations: (5902) 
#> Path [44] :Initial log joint density = -481902.090846 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              73      -4.787e+05      6.577e-03   1.592e-01    8.558e-01  8.558e-01      4672 -3.873e+03 -3.920e+03                   
#> Path [44] :Best Iter: [67] ELBO (-3873.477916) evaluations: (4672) 
#> Path [45] :Initial log joint density = -481754.840923 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      3.610e-03   2.941e-01    5.932e-01  5.932e-01      5215 -3.874e+03 -3.917e+03                   
#> Path [45] :Best Iter: [76] ELBO (-3874.295210) evaluations: (5215) 
#> Path [46] :Initial log joint density = -481698.817433 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      1.243e-02   1.998e-01    1.000e+00  1.000e+00      5181 -3.879e+03 -3.884e+03                   
#> Path [46] :Best Iter: [74] ELBO (-3879.450126) evaluations: (5181) 
#> Path [47] :Initial log joint density = -481643.989964 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              99      -4.787e+05      9.083e-03   1.304e-01    8.799e-01  8.799e-01      7692 -3.872e+03 -3.931e+03                   
#> Path [47] :Best Iter: [75] ELBO (-3872.213609) evaluations: (7692) 
#> Path [48] :Initial log joint density = -481615.801751 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      5.839e-02   6.051e-01    1.000e+00  1.000e+00      7767 -3.879e+03 -3.884e+03                   
#> Path [48] :Best Iter: [45] ELBO (-3878.668990) evaluations: (7767) 
#> Path [49] :Initial log joint density = -481517.096005 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              88      -4.787e+05      1.416e-02   1.589e-01    1.000e+00  1.000e+00      6615 -3.866e+03 -3.904e+03                   
#> Path [49] :Best Iter: [44] ELBO (-3866.150952) evaluations: (6615) 
#> Path [50] :Initial log joint density = -481533.695539 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              80      -4.787e+05      4.394e-03   1.635e-01    6.904e-01  6.904e-01      5479 -3.875e+03 -3.979e+03                   
#> Path [50] :Best Iter: [74] ELBO (-3875.181771) evaluations: (5479) 
#> Finished in  28.6 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_CLR_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.083 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 12
#>    cell_group smooth  from  to    c_lower c_effect c_upper  c_pH0   c_FDR c_rhat
#>    <chr>      <chr>   <chr> <chr>   <dbl>    <dbl>   <dbl>  <dbl>   <dbl>  <dbl>
#>  1 B1         pseudo… 2     8     -0.802   -0.400  -0.0391 0.0520 0.0268   1.000
#>  2 B2         pseudo… 2     8     -0.423   -0.0503  0.358  0.608  0.225    0.998
#>  3 B3         pseudo… 2     8     -0.251    0.0963  0.487  0.526  0.150    1.000
#>  4 BM         pseudo… 2     8     -0.821   -0.382  -0.0157 0.0680 0.0330   1.02 
#>  5 CD4 1      pseudo… 2     8     -0.364   -0.0530  0.256  0.628  0.238    1.00 
#>  6 CD4 2      pseudo… 2     8     -0.0950   0.248   0.601  0.202  0.0741   1.00 
#>  7 CD4 3      pseudo… 2     8     -1.06    -0.615  -0.169  0.0180 0.00550  0.998
#>  8 CD4 4      pseudo… 2     8     -0.356    0.0169  0.405  0.726  0.278    0.999
#>  9 CD4 5      pseudo… 2     8     -0.694   -0.276   0.0507 0.182  0.067    1.000
#> 10 CD8 1      pseudo… 2     8     -0.122    0.186   0.494  0.276  0.0842   0.998
#> # ℹ 26 more rows
#> # ℹ 2 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>
# }
```
