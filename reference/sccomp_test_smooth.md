# Test differences along a smooth term

Tests whether the mean cell-group composition differs between two
positions (or two intervals) of a continuous covariate modelled with a
smooth `s()` term. The test is a contrast of the model's mean prediction
on the logit scale, evaluated at the two positions and summarised with
the same probability-of-null (`pH0`) and false-discovery-rate (`FDR`)
machinery as
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
  test_composition_above_logit_fold_change = 0.1,
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

- test_composition_above_logit_fold_change:

  Positive numeric. Effect threshold for the hypothesis test, on the
  logit scale — identical meaning to the argument of
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
  logit-scale contrast \\\mu(\code{to}) - \mu(\code{from})\\.

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
logit-scale predictor is additive; it therefore does not matter what
value those covariates take, as long as it is held constant. For factor
smooths (`s(x, g, bs = "fs")`) or by-factor smooths (`s(x, by = g)`) the
grouping factor selects *which* curve is tested; set it via `at`, e.g.
`at = list(tissue = "tumor")`.

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
#> sccomp says: the composition design matrix has columns: (Intercept), s(pseudotime, k = 4)__lin1
#> sccomp says: the variability design matrix has columns: (Intercept)
#> Loading model from cache...
#> Path [1] :Initial log joint density = -481788.714266 
#> Path [1] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.304e-02   2.016e-01    1.000e+00  1.000e+00      4580 -3.845e+03 -3.872e+03                   
#> Path [1] :Best Iter: [34] ELBO (-3844.745417) evaluations: (4580) 
#> Path [2] :Initial log joint density = -481554.518719 
#> Path [2] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.078e-02   1.115e-01    1.000e+00  1.000e+00      7847 -3.847e+03 -4.040e+03                   
#> Path [2] :Best Iter: [39] ELBO (-3846.543997) evaluations: (7847) 
#> Path [3] :Initial log joint density = -481097.615381 
#> Path [3] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              92      -4.787e+05      7.264e-03   1.214e-01    1.000e+00  1.000e+00      6763 -3.844e+03 -3.901e+03                   
#> Path [3] :Best Iter: [37] ELBO (-3844.180391) evaluations: (6763) 
#> Path [4] :Initial log joint density = -481808.306196 
#> Path [4] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              99      -4.787e+05      9.372e-03   1.538e-01    1.000e+00  1.000e+00      7685 -3.849e+03 -5.836e+03                   
#> Path [4] :Best Iter: [74] ELBO (-3848.796695) evaluations: (7685) 
#> Path [5] :Initial log joint density = -482067.336965 
#> Path [5] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              92      -4.787e+05      1.487e-02   2.216e-01    8.476e-01  8.476e-01      6771 -3.842e+03 -3.887e+03                   
#> Path [5] :Best Iter: [77] ELBO (-3841.691835) evaluations: (6771) 
#> Path [6] :Initial log joint density = -481651.721959 
#> Path [6] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              87      -4.787e+05      2.646e-02   1.711e-01    1.000e+00  1.000e+00      6270 -3.844e+03 -3.869e+03                   
#> Path [6] :Best Iter: [66] ELBO (-3843.925547) evaluations: (6270) 
#> Path [7] :Initial log joint density = -481354.266412 
#> Path [7] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.398e-02   2.551e-01    1.000e+00  1.000e+00      4461 -3.850e+03 -3.864e+03                   
#> Path [7] :Best Iter: [68] ELBO (-3849.578769) evaluations: (4461) 
#> Path [8] :Initial log joint density = -481288.458914 
#> Path [8] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              94      -4.787e+05      5.779e-03   2.339e-01    6.093e-01  6.093e-01      7158 -3.842e+03 -3.892e+03                   
#> Path [8] :Best Iter: [71] ELBO (-3842.109338) evaluations: (7158) 
#> Path [9] :Initial log joint density = -482041.378764 
#> Path [9] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      1.423e-02   2.157e-01    1.000e+00  1.000e+00      4767 -3.855e+03 -3.861e+03                   
#> Path [9] :Best Iter: [69] ELBO (-3854.754221) evaluations: (4767) 
#> Path [10] :Initial log joint density = -481265.689352 
#> Path [10] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      2.023e-02   2.379e-01    9.432e-01  9.432e-01      7960 -3.839e+03 -3.891e+03                   
#> Path [10] :Best Iter: [41] ELBO (-3839.427632) evaluations: (7960) 
#> Path [11] :Initial log joint density = -481602.656305 
#> Path [11] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              90      -4.787e+05      1.375e-02   1.817e-01    1.000e+00  1.000e+00      6665 -3.846e+03 -4.017e+03                   
#> Path [11] :Best Iter: [38] ELBO (-3846.006849) evaluations: (6665) 
#> Path [12] :Initial log joint density = -481477.129273 
#> Path [12] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      2.100e-02   2.469e-01    1.000e+00  1.000e+00      4666 -3.846e+03 -3.858e+03                   
#> Path [12] :Best Iter: [70] ELBO (-3846.454549) evaluations: (4666) 
#> Path [13] :Initial log joint density = -481327.758865 
#> Path [13] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              91      -4.787e+05      1.125e-02   1.766e-01    1.000e+00  1.000e+00      6819 -3.847e+03 -3.986e+03                   
#> Path [13] :Best Iter: [73] ELBO (-3846.869105) evaluations: (6819) 
#> Path [14] :Initial log joint density = -481597.989717 
#> Path [14] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              94      -4.787e+05      1.625e-02   2.083e-01    1.000e+00  1.000e+00      7296 -3.850e+03 -3.903e+03                   
#> Path [14] :Best Iter: [63] ELBO (-3850.092073) evaluations: (7296) 
#> Path [15] :Initial log joint density = -481593.730115 
#> Path [15] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.562e-02   1.726e-01    1.000e+00  1.000e+00      8042 -3.845e+03 -3.890e+03                   
#> Path [15] :Best Iter: [71] ELBO (-3845.289266) evaluations: (8042) 
#> Path [16] :Initial log joint density = -481542.248913 
#> Path [16] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              68      -4.787e+05      9.125e-03   1.664e-01    1.000e+00  1.000e+00      4245 -3.848e+03 -3.872e+03                   
#> Path [16] :Best Iter: [65] ELBO (-3847.651526) evaluations: (4245) 
#> Path [17] :Initial log joint density = -481558.836728 
#> Path [17] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.422e-02   2.274e-01    1.000e+00  1.000e+00      7847 -3.846e+03 -3.868e+03                   
#> Path [17] :Best Iter: [73] ELBO (-3846.198769) evaluations: (7847) 
#> Path [18] :Initial log joint density = -482141.496446 
#> Path [18] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              88      -4.787e+05      3.682e-02   2.999e-01    1.000e+00  1.000e+00      6535 -3.846e+03 -3.983e+03                   
#> Path [18] :Best Iter: [68] ELBO (-3846.093475) evaluations: (6535) 
#> Path [19] :Initial log joint density = -481863.172536 
#> Path [19] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              75      -4.787e+05      6.106e-03   1.832e-01    1.000e+00  1.000e+00      4904 -3.853e+03 -3.881e+03                   
#> Path [19] :Best Iter: [50] ELBO (-3853.215230) evaluations: (4904) 
#> Path [20] :Initial log joint density = -481989.154423 
#> Path [20] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.557e-02   2.276e-01    1.000e+00  1.000e+00      4335 -3.844e+03 -3.857e+03                   
#> Path [20] :Best Iter: [68] ELBO (-3844.305513) evaluations: (4335) 
#> Path [21] :Initial log joint density = -482657.074974 
#> Path [21] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.047e-02   1.983e-01    1.000e+00  1.000e+00      4653 -3.847e+03 -3.879e+03                   
#> Path [21] :Best Iter: [42] ELBO (-3847.440668) evaluations: (4653) 
#> Path [22] :Initial log joint density = -481679.297399 
#> Path [22] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.787e+05      2.889e-02   2.869e-01    1.000e+00  1.000e+00      7629 -3.844e+03 -3.864e+03                   
#> Path [22] :Best Iter: [79] ELBO (-3844.427561) evaluations: (7629) 
#> Path [23] :Initial log joint density = -482420.397494 
#> Path [23] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      3.046e-02   2.477e-01    1.000e+00  1.000e+00      6298 -3.854e+03 -3.878e+03                   
#> Path [23] :Best Iter: [72] ELBO (-3853.828758) evaluations: (6298) 
#> Path [24] :Initial log joint density = -481704.512927 
#> Path [24] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              78      -4.787e+05      7.509e-03   2.232e-01    8.936e-01  8.936e-01      5260 -3.856e+03 -3.893e+03                   
#> Path [24] :Best Iter: [71] ELBO (-3856.092729) evaluations: (5260) 
#> Path [25] :Initial log joint density = -483508.327502 
#> Path [25] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              92      -4.787e+05      1.891e-02   2.100e-01    1.000e+00  1.000e+00      7214 -3.853e+03 -3.917e+03                   
#> Path [25] :Best Iter: [73] ELBO (-3853.426920) evaluations: (7214) 
#> Path [26] :Initial log joint density = -481607.811712 
#> Path [26] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              96      -4.787e+05      2.147e-02   1.729e-01    1.000e+00  1.000e+00      7396 -3.842e+03 -3.881e+03                   
#> Path [26] :Best Iter: [35] ELBO (-3842.268448) evaluations: (7396) 
#> Path [27] :Initial log joint density = -481606.538655 
#> Path [27] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      2.151e-02   2.250e-01    1.000e+00  1.000e+00      4401 -3.847e+03 -3.869e+03                   
#> Path [27] :Best Iter: [43] ELBO (-3846.610715) evaluations: (4401) 
#> Path [28] :Initial log joint density = -481620.020424 
#> Path [28] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              70      -4.787e+05      5.559e-03   1.996e-01    1.000e+00  1.000e+00      4468 -3.856e+03 -3.878e+03                   
#> Path [28] :Best Iter: [60] ELBO (-3856.070445) evaluations: (4468) 
#> Path [29] :Initial log joint density = -481421.290639 
#> Path [29] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      7.017e-02   5.695e-01    1.000e+00  1.000e+00      7680 -3.848e+03 -4.174e+03                   
#> Path [29] :Best Iter: [70] ELBO (-3847.744453) evaluations: (7680) 
#> Path [30] :Initial log joint density = -481345.234941 
#> Path [30] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      9.968e-03   1.947e-01    1.000e+00  1.000e+00      4933 -3.856e+03 -3.874e+03                   
#> Path [30] :Best Iter: [71] ELBO (-3855.771504) evaluations: (4933) 
#> Path [31] :Initial log joint density = -482146.169570 
#> Path [31] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              97      -4.787e+05      5.408e-03   2.554e-01    5.585e-01  5.585e-01      7372 -3.851e+03       -inf                   
#> Path [31] :Best Iter: [44] ELBO (-3851.205337) evaluations: (7372) 
#> Path [32] :Initial log joint density = -482693.264600 
#> Path [32] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              71      -4.787e+05      1.386e-02   1.922e-01    1.000e+00  1.000e+00      4643 -3.851e+03 -3.862e+03                   
#> Path [32] :Best Iter: [68] ELBO (-3851.489414) evaluations: (4643) 
#> Path [33] :Initial log joint density = -481515.285347 
#> Path [33] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.075e-01   1.125e+00    1.000e+00  1.000e+00      7747 -3.847e+03 -3.878e+03                   
#> Path [33] :Best Iter: [76] ELBO (-3846.851257) evaluations: (7747) 
#> Path [34] :Initial log joint density = -481819.264793 
#> Path [34] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              89      -4.787e+05      1.534e-02   2.082e-01    1.000e+00  1.000e+00      6544 -3.851e+03       -inf                   
#> Path [34] :Best Iter: [45] ELBO (-3851.482323) evaluations: (6544) 
#> Path [35] :Initial log joint density = -481516.380285 
#> Path [35] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              65      -4.787e+05      6.460e-03   2.337e-01    6.934e-01  6.934e-01      4138 -3.851e+03 -3.896e+03                   
#> Path [35] :Best Iter: [63] ELBO (-3850.766673) evaluations: (4138) 
#> Path [36] :Initial log joint density = -481511.797983 
#> Path [36] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      9.591e-03   1.902e-01    1.000e+00  1.000e+00      5331 -3.852e+03 -3.858e+03                   
#> Path [36] :Best Iter: [66] ELBO (-3852.141983) evaluations: (5331) 
#> Path [37] :Initial log joint density = -481358.017569 
#> Path [37] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.574e-02   1.729e-01    1.000e+00  1.000e+00      5093 -3.846e+03 -3.880e+03                   
#> Path [37] :Best Iter: [70] ELBO (-3846.398663) evaluations: (5093) 
#> Path [38] :Initial log joint density = -481921.183069 
#> Path [38] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      8.915e-03   2.058e-01    1.000e+00  1.000e+00      5226 -3.855e+03 -3.891e+03                   
#> Path [38] :Best Iter: [71] ELBO (-3854.522443) evaluations: (5226) 
#> Path [39] :Initial log joint density = -482167.836887 
#> Path [39] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              98      -4.787e+05      1.947e-02   1.631e-01    1.000e+00  1.000e+00      7590 -3.847e+03 -3.883e+03                   
#> Path [39] :Best Iter: [75] ELBO (-3846.856139) evaluations: (7590) 
#> Path [40] :Initial log joint density = -481802.159080 
#> Path [40] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      8.837e-03   2.286e-01    7.904e-01  7.904e-01      7897 -3.850e+03 -4.021e+03                   
#> Path [40] :Best Iter: [75] ELBO (-3849.659092) evaluations: (7897) 
#> Path [41] :Initial log joint density = -481870.846771 
#> Path [41] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              79      -4.787e+05      1.884e-02   1.778e-01    1.000e+00  1.000e+00      5330 -3.848e+03 -3.856e+03                   
#> Path [41] :Best Iter: [78] ELBO (-3847.678778) evaluations: (5330) 
#> Path [42] :Initial log joint density = -484358.914685 
#> Path [42] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              99      -4.787e+05      1.199e-02   1.791e-01    1.000e+00  1.000e+00      8285 -3.848e+03 -3.891e+03                   
#> Path [42] :Best Iter: [64] ELBO (-3848.346544) evaluations: (8285) 
#> Path [43] :Initial log joint density = -481604.371051 
#> Path [43] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      9.960e-03   1.933e-01    1.000e+00  1.000e+00      4979 -3.844e+03 -3.858e+03                   
#> Path [43] :Best Iter: [42] ELBO (-3844.403293) evaluations: (4979) 
#> Path [44] :Initial log joint density = -481917.230464 
#> Path [44] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              69      -4.787e+05      1.176e-02   2.422e-01    1.000e+00  1.000e+00      4278 -3.853e+03 -3.877e+03                   
#> Path [44] :Best Iter: [66] ELBO (-3852.788978) evaluations: (4278) 
#> Path [45] :Initial log joint density = -481771.912084 
#> Path [45] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              77      -4.787e+05      1.066e-02   2.500e-01    5.179e-01  1.000e+00      5117 -3.854e+03 -3.913e+03                   
#> Path [45] :Best Iter: [66] ELBO (-3854.049018) evaluations: (5117) 
#> Path [46] :Initial log joint density = -481715.028195 
#> Path [46] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              76      -4.787e+05      1.043e-02   1.730e-01    1.000e+00  1.000e+00      4903 -3.864e+03 -3.874e+03                   
#> Path [46] :Best Iter: [74] ELBO (-3864.391384) evaluations: (4903) 
#> Path [47] :Initial log joint density = -481659.323096 
#> Path [47] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>             100      -4.787e+05      1.790e-02   1.308e-01    1.000e+00  1.000e+00      7888 -3.854e+03 -3.870e+03                   
#> Path [47] :Best Iter: [77] ELBO (-3853.835389) evaluations: (7888) 
#> Path [48] :Initial log joint density = -481629.588401 
#> Path [48] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              74      -4.787e+05      5.219e-03   1.708e-01    8.272e-01  8.272e-01      4700 -3.857e+03 -3.892e+03                   
#> Path [48] :Best Iter: [47] ELBO (-3856.698748) evaluations: (4700) 
#> Path [49] :Initial log joint density = -481532.917306 
#> Path [49] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              86      -4.787e+05      8.517e-03   1.736e-01    1.000e+00  1.000e+00      6226 -3.847e+03 -3.923e+03                   
#> Path [49] :Best Iter: [59] ELBO (-3846.806768) evaluations: (6226) 
#> Path [50] :Initial log joint density = -481545.583637 
#> Path [50] : Iter      log prob        ||dx||      ||grad||     alpha      alpha0      # evals       ELBO    Best ELBO        Notes  
#>              72      -4.787e+05      3.573e-03   2.246e-01    6.625e-01  6.625e-01      4603 -3.856e+03 -3.895e+03                   
#> Path [50] :Best Iter: [46] ELBO (-3856.403656) evaluations: (4603) 
#> Finished in  33.4 seconds.
#> sccomp says: to do hypothesis testing run `sccomp_test()`,
#>   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
#>   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
#>   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).
#> sccomp says: auto-cleanup removed 1 draw files from 'sccomp_draws_files'
#> Loading model from cache...
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 0.093 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
#> # A tibble: 36 × 12
#>    cell_group smooth from  to    c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
#>    <chr>      <chr>  <chr> <chr>   <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
#>  1 B1         pseud… 2     8      -0.791 -0.398   -0.0289 0.0520  0.0185   1.00 
#>  2 B2         pseud… 2     8      -0.419 -0.0371   0.299  0.63    0.224    1.00 
#>  3 B3         pseud… 2     8      -0.260  0.0852   0.447  0.572   0.162    0.999
#>  4 BM         pseud… 2     8      -0.753 -0.382   -0.0242 0.0580  0.0258   1.01 
#>  5 CD4 1      pseud… 2     8      -0.419 -0.0547   0.268  0.624   0.209    1.01 
#>  6 CD4 2      pseud… 2     8      -0.122  0.226    0.542  0.222   0.0779   0.999
#>  7 CD4 3      pseud… 2     8      -0.984 -0.602   -0.234  0.00400 0.00150  1.01 
#>  8 CD4 4      pseud… 2     8      -0.337  0.00209  0.365  0.722   0.312    1.00 
#>  9 CD4 5      pseud… 2     8      -0.676 -0.267    0.118  0.17    0.0559   0.999
#> 10 CD8 1      pseud… 2     8      -0.109  0.178    0.498  0.32    0.1      1.00 
#> # ℹ 26 more rows
#> # ℹ 2 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>
# }
```
