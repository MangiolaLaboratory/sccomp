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
| `sccomp_boxplot`                   | Plots the data across factors, overlaying model-generated data, for asses the goodness of fit                               |
| `sccomp_plot_summary`              | Plors summary plots to asses significance                                                                                   |

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
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.868     1.12   1.35   0       0         4803.
    ##  2 B1         typecancer type   -1.05     -0.643 -0.241  1.42e-2 4.90e-3   3590.
    ##  3 B2         (Intercep… <NA>    0.408     0.701  0.972  7.50e-4 6.58e-5   4625.
    ##  4 B2         typecancer type   -1.19     -0.724 -0.264  1.40e-2 3.86e-3   2444.
    ##  5 B3         (Intercep… <NA>   -0.678    -0.389 -0.111  8.98e-2 8.9 e-3   4102.
    ##  6 B3         typecancer type   -0.758    -0.319  0.0904 2.88e-1 7.64e-2   3031.
    ##  7 BM         (Intercep… <NA>   -1.33     -1.04  -0.741  0       0         4696.
    ##  8 BM         typecancer type   -0.770    -0.318  0.122  2.95e-1 9.60e-2   4254.
    ##  9 CD4 1      (Intercep… <NA>    0.0791    0.297  0.517  1.70e-1 2.35e-2   3819.
    ## 10 CD4 1      typecancer type   -0.119     0.182  0.480  5.50e-1 1.60e-1   3761.
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
    sccomp_test(
      contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -1.90    -1.39   -0.880 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.880    1.39    1.90  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.32    -1.71   -1.08  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.08     1.71    2.32  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.48    -1.04   -0.607 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.607    1.04    1.48  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.843    1.76    2.84  2.50e-4 5.00e-5      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.84    -1.76   -0.843 2.50e-4 5.00e-5      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.307    0.983   1.68  1.12e-2 3.35e-3      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.68    -0.983  -0.307 1.12e-2 3.35e-3      NA
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

# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -80.8      11.3

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
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )

res
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  (Interce… <NA>    0.343     0.762  1.16   5.00e-3 5.19e-4   6144.
    ##  2 B immature  typeheal… type    0.878     1.44   1.99   2.50e-4 8.33e-5   5459.
    ##  3 B mem       (Interce… <NA>   -1.46     -0.871 -0.178  2.68e-2 4.76e-3   4473.
    ##  4 B mem       typeheal… type    1.07      1.85   2.63   0       0         4004.
    ##  5 CD4 cm S10… (Interce… <NA>    1.31      1.66   1.99   0       0         6321.
    ##  6 CD4 cm S10… typeheal… type    0.472     0.941  1.41   1.00e-3 3.12e-4   4707.
    ##  7 CD4 cm hig… (Interce… <NA>   -1.08     -0.541  0.0370 1.05e-1 1.92e-2   4181.
    ##  8 CD4 cm hig… typeheal… type   -3.06     -1.26   1.01   2.01e-1 5.07e-2   3636.
    ##  9 CD4 cm rib… (Interce… <NA>   -0.0630    0.311  0.695  2.76e-1 3.34e-2   4496.
    ## 10 CD4 cm rib… typeheal… type   -1.80     -0.960  0.124  7.55e-2 2.43e-2   3299.
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
plots = plot_summary(res) 
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
plots = plot_summary(res)
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
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )

res
```

    ## # A tibble: 210 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature (Intercep… <NA>     0.338   0.871    1.49  1.10e-2 3.07e-3   2668.
    ##  2 B immature typehealt… type     0.625   1.32     1.98  5.00e-4 2.50e-4   2930.
    ##  3 B immature continuou… conti…  -0.241   0.0316   0.325 8.70e-1 6.43e-1   8525.
    ##  4 B immature (Intercep… <NA>    -0.568  -0.0497   0.419 7.35e-1 6.35e-1     NA 
    ##  5 B immature typehealt… <NA>    -0.419   0.0497   0.568 7.35e-1 6.35e-1     NA 
    ##  6 B immature (Intercep… <NA>    -0.630  -0.0329   0.475 7.39e-1 6.72e-1     NA 
    ##  7 B immature typehealt… <NA>    -0.475   0.0329   0.630 7.39e-1 6.72e-1     NA 
    ##  8 B mem      (Intercep… <NA>    -1.16   -0.472    0.388 2.44e-1 4.93e-2   1904.
    ##  9 B mem      typehealt… type     0.490   1.41     2.21  7.00e-3 2.56e-3   1984.
    ## 10 B mem      continuou… conti…  -0.257   0.0542   0.401 8.13e-1 6.18e-1   8974.
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
    ##  1 10x_6K       B immature              0.0553           260.          -0.675 
    ##  2 10x_8K       B immature              0.143           1079.           0.398 
    ##  3 GSE115189    B immature              0.114            267.           0.102 
    ##  4 SCP345_580   B immature              0.0902           520.          -0.128 
    ##  5 SCP345_860   B immature              0.151            967.           0.455 
    ##  6 SCP424_pbmc1 B immature              0.111            297.           0.0579
    ##  7 SCP424_pbmc2 B immature              0.200            598.           0.791 
    ##  8 SCP591       B immature              0.0248            14.1         -1.49  
    ##  9 SI-GA-E5     B immature              0.0261           109.          -0.724 
    ## 10 SI-GA-E7     B immature              0.104            763.           0.763 
    ## # ℹ 590 more rows
