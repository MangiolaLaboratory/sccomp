sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)

<!-- badges: end -->

<a href="https://www.youtube.com/watch?v=R_lt58We9nA&ab_channel=RConsortium" target="_blank">
<img src="https://img.youtube.com/vi/R_lt58We9nA/mqdefault.jpg" alt="Watch the video" width="280" height="180" border="10" />
</a>

# <img src="inst/logo-01.png" height="139px" width="120px"/>

Sccomp is a generalised method for differential composition and
variability analyses.

## Characteristics

- Complex linear models with continuous and categorical covariates
- Multilevel modelling, with population (i.e. fixed) and group (random)
  effects
- Modelling counts
- Modelling proportionality
- Modelling cell-type specific variability
- Cell-type information share for variability adaptive shrinkage
- Testing differential variability
- Probabilistic outlier identification
- Cross-dataset learning (hyperpriors).

# Installation

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
BiocManager::install("sccomp")
```

**Github**

``` r
devtools::install_github("stemangiola/sccomp")
```

| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                  | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`           | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                      | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`                 | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                   | Predicts proportions, based on the mode, or part of the model                                                               |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plors summary plots to asses significance                                                                                   |

# Analysis

`sccomp` can model changes in composition and variability. By default,
the formula for variability is either `~1`, which assumes that the
cell-group variability is independent of any covariate or
`~ factor_of_interest`, which assumes that the model is dependent on the
factor of interest only. The variability model must be a subset of the
model for composition.

## Binary factor

### From Seurat, SingleCellExperiment, metadata objects

``` r
single_cell_object |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
```

### From counts

``` r
counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000394 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.94 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 2.839 seconds (Warm-up)
    ## Chain 1:                21.137 seconds (Sampling)
    ## Chain 1:                23.976 seconds (Total)
    ## Chain 1:

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.884     1.12   1.32   0       0         4441.
    ##  2 B1         typecancer type   -1.17     -0.760 -0.381  1.75e-3 6.51e-4   2922.
    ##  3 B2         (Intercep… <NA>    0.423     0.703  0.974  2.50e-4 1.32e-5   6265.
    ##  4 B2         typecancer type   -1.23     -0.730 -0.263  1.25e-2 3.53e-3   3965.
    ##  5 B3         (Intercep… <NA>   -0.657    -0.389 -0.114  8.51e-2 6.86e-3   4537.
    ##  6 B3         typecancer type   -0.719    -0.317  0.0704 2.76e-1 6.71e-2   3631.
    ##  7 BM         (Intercep… <NA>   -1.32     -1.03  -0.767  0       0         3997.
    ##  8 BM         typecancer type   -0.736    -0.312  0.0916 2.94e-1 8.78e-2   3760.
    ##  9 CD4 1      (Intercep… <NA>    0.0853    0.298  0.499  1.79e-1 3.16e-2   4733.
    ## 10 CD4 1      typecancer type   -0.0963    0.186  0.467  5.41e-1 1.52e-1   3783.
    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

## Contrasts

``` r
seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(
      contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000371 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.71 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 2.127 seconds (Warm-up)
    ## Chain 1:                17.85 seconds (Sampling)
    ## Chain 1:                19.977 seconds (Total)
    ## Chain 1:

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -1.91    -1.42   -0.937 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.937    1.42    1.91  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.37    -1.74   -1.11  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.11     1.74    2.37  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.25    -0.858  -0.514 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.514    0.858   1.25  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.829    1.83    2.96  2.50e-4 4.17e-5      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.96    -1.83   -0.829 2.50e-4 4.17e-5      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.346    1.01    1.74  6.51e-3 1.25e-3      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.74    -1.01   -0.346 6.51e-3 1.25e-3      NA
    ## # ℹ 50 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

## Categorical factor (e.g. Bayesian ANOVA)

This is achieved through model comparison with `loo`. In the following
example, the model with association with factors better fits the data
compared to the baseline model with no factor association. For
comparisons `check_outliers` must be set to FALSE as the leave-one-out
must work with the same amount of data, while outlier elimination does
not guarantee it.

If `elpd_diff` is away from zero of \> 5 `se_diff` difference of 5, we
are confident that a model is better than the other
[reference](https://discourse.mc-stan.org/t/interpreting-elpd-diff-loo-package/1628/2?u=stemangiola).
In this case, -79.9 / 11.5 = -6.9, therefore we can conclude that model
one, the one with factor association, is better than model two.

``` r
library(loo)

# Fit first model
model_with_factor_association = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000394 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 3.94 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 2.234 seconds (Warm-up)
    ## Chain 1:                20.698 seconds (Sampling)
    ## Chain 1:                22.932 seconds (Total)
    ## Chain 1:

``` r
# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    enable_loo = TRUE
  )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000286 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.86 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 2.229 seconds (Warm-up)
    ## Chain 1:                23.82 seconds (Sampling)
    ## Chain 1:                26.049 seconds (Total)
    ## Chain 1:

``` r
# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -81.1      11.3

## Differential variability, binary factor

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.001618 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 16.18 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 4.185 seconds (Warm-up)
    ## Chain 1:                33.902 seconds (Sampling)
    ## Chain 1:                38.087 seconds (Total)
    ## Chain 1:

``` r
res
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  (Interce… <NA>     0.548    0.934  1.30   2.50e-4 2.28e-5   5948.
    ##  2 B immature  typeheal… type     0.809    1.33   1.89   0       0         4471.
    ##  3 B mem       (Interce… <NA>    -1.30    -0.759 -0.159  3.08e-2 2.21e-3   5687.
    ##  4 B mem       typeheal… type     1.06     1.81   2.52   2.50e-4 8.34e-5   4832.
    ##  5 CD4 cm S10… (Interce… <NA>     1.73     1.98   2.23   0       0         5733.
    ##  6 CD4 cm S10… typeheal… type     0.301    0.689  1.07   8.26e-3 3.42e-3   4382.
    ##  7 CD4 cm hig… (Interce… <NA>    -0.881   -0.413  0.113  1.92e-1 2.72e-2   4595.
    ##  8 CD4 cm hig… typeheal… type    -3.22    -1.49   1.03   1.52e-1 5.59e-2   2946.
    ##  9 CD4 cm rib… (Interce… <NA>     0.144    0.471  0.809  5.13e-2 4.67e-3   3843.
    ## 10 CD4 cm rib… typeheal… type    -1.90    -1.09  -0.0711 4.25e-2 9.48e-3   3818.
    ## # ℹ 50 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Suggested settings

## For single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

# Visualisation

## Summary plots

``` r
plots = plot(res) 
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

    ## Warning: Expected 2 pieces. Additional pieces discarded in 4 rows [6, 7, 13,
    ## 14].

A plot of group proportion, faceted by groups. The blue boxplots
represent the posterior predictive check. If the model is likely to be
descriptively adequate to the data, the blue box plot should roughly
overlay with the black box plot, which represents the observed data. The
outliers are coloured in red. A box plot will be returned for every
(discrete) covariate present in `formula_composition`. The colour coding
represents the significant associations for composition and/or
variability.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = plot(res)
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

    ## Warning: Expected 2 pieces. Additional pieces discarded in 4 rows [6, 7, 13,
    ## 14].

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-17-1.png)<!-- -->

Plot 2D significance plot. Data points are cell groups. Error bars are
the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis that
no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (besides the Intercept) is
provided for each category of interest. The for each category of
interest, the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-18-1.png)<!-- -->

# Multilevel modelling

`sccomp` is cabable of estimating population (i.e. fixed) and group
(i.e. random) effects. The formulation is analogous to the `lme4`
package and `brms`.

!! For now, only one grouping is allowed (e.g. group2\_\_).

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type + continuous_covariate + (type | group2__), 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000458 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 4.58 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 5.428 seconds (Warm-up)
    ## Chain 1:                39.833 seconds (Sampling)
    ## Chain 1:                45.261 seconds (Total)
    ## Chain 1:

``` r
res
```

    ## # A tibble: 210 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature (Intercep… <NA>     0.543   1.05     1.61  0       0         2521.
    ##  2 B immature typehealt… type     0.611   1.26     1.88  0.00175 8.76e-4   2604.
    ##  3 B immature continuou… conti…  -0.207   0.0390   0.334 0.873   6.26e-1   6341.
    ##  4 B immature (Intercep… <NA>    -0.607  -0.0870   0.392 0.681   5.86e-1     NA 
    ##  5 B immature typehealt… <NA>    -0.392   0.0870   0.607 0.681   5.86e-1     NA 
    ##  6 B immature (Intercep… <NA>    -0.576  -0.0465   0.435 0.733   6.50e-1     NA 
    ##  7 B immature typehealt… <NA>    -0.435   0.0465   0.576 0.733   6.50e-1     NA 
    ##  8 B mem      (Intercep… <NA>    -1.04   -0.379    0.459 0.305   5.49e-2   2005.
    ##  9 B mem      typehealt… type     0.506   1.41     2.22  0.00475 2.70e-3   2206.
    ## 10 B mem      continuou… conti…  -0.234   0.0635   0.399 0.805   5.97e-1   7345.
    ## # ℹ 200 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Removal of unwanted variation

After you model your dataset, you can remove the unwanted variation from
your input data, **for visualisation purposes**

We decide to just keep the type population (i.e. fixed) effect for
abundance, and do not keep it for variability.

``` r
res |> sccomp_remove_unwanted_variation(~type)
```

    ## sccomp says: calculating residuals

    ## sccomp says: regressing out unwanted factors

    ## # A tibble: 600 × 5
    ##    sample       cell_group adjusted_proportion adjusted_counts logit_residuals
    ##    <chr>        <chr>                    <dbl>           <dbl>           <dbl>
    ##  1 10x_6K       B immature              0.0554           260.          -0.702 
    ##  2 10x_8K       B immature              0.144           1080.           0.366 
    ##  3 GSE115189    B immature              0.114            268.           0.0767
    ##  4 SCP345_580   B immature              0.0903           520.          -0.161 
    ##  5 SCP345_860   B immature              0.151            969.           0.421 
    ##  6 SCP424_pbmc1 B immature              0.110            295.           0.0129
    ##  7 SCP424_pbmc2 B immature              0.200            596.           0.755 
    ##  8 SCP591       B immature              0.0251            14.3         -1.52  
    ##  9 SI-GA-E5     B immature              0.0272           114.          -0.692 
    ## 10 SI-GA-E7     B immature              0.105            772.           0.747 
    ## # ℹ 590 more rows
