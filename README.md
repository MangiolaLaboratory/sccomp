sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)

<!-- badges: end -->

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
    ##  1 B1         (Intercep… <NA>    0.873     1.12   1.35   0       0         6210.
    ##  2 B1         typecancer type   -1.06     -0.642 -0.262  1.18e-2 2.72e-3   4662.
    ##  3 B2         (Intercep… <NA>    0.398     0.695  0.996  7.50e-4 6.58e-5   5458.
    ##  4 B2         typecancer type   -1.21     -0.721 -0.233  1.80e-2 6.34e-3   3621.
    ##  5 B3         (Intercep… <NA>   -0.680    -0.391 -0.121  7.90e-2 8.27e-3   5356.
    ##  6 B3         typecancer type   -0.728    -0.310  0.0880 2.89e-1 8.44e-2   3963.
    ##  7 BM         (Intercep… <NA>   -1.32     -1.04  -0.758  0       0         6759.
    ##  8 BM         typecancer type   -0.737    -0.312  0.0969 2.94e-1 9.39e-2   4720.
    ##  9 CD4 1      (Intercep… <NA>    0.0807    0.298  0.503  1.83e-1 2.86e-2   4736.
    ## 10 CD4 1      typecancer type   -0.0988    0.188  0.483  5.30e-1 1.57e-1   4729.
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
    ##  1 B immature  typecanc… <NA>    -1.89    -1.40   -0.901 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.901    1.40    1.89  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.33    -1.72   -1.07  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.07     1.72    2.33  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.47    -1.04   -0.617 2.50e-4 6.25e-5      NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.617    1.04    1.47  2.50e-4 6.25e-5      NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.837    1.77    2.91  5.00e-4 1.50e-4      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.91    -1.77   -0.837 5.00e-4 1.50e-4      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.290    0.997   1.71  1.30e-2 3.57e-3      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.71    -0.997  -0.290 1.30e-2 3.57e-3      NA
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
    ## model2 -80.2      11.4

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
    ##  1 B immature  (Interce… <NA>    0.371     0.769  1.20   0.00550 5.19e-4   5537.
    ##  2 B immature  typeheal… type    0.851     1.43   2.00   0       0         5064.
    ##  3 B mem       (Interce… <NA>   -1.50     -0.878 -0.188  0.0270  4.69e-3   3986.
    ##  4 B mem       typeheal… type    1.07      1.87   2.65   0       0         3711.
    ##  5 CD4 cm S10… (Interce… <NA>    1.32      1.66   2.00   0       0         7824.
    ##  6 CD4 cm S10… typeheal… type    0.485     0.938  1.39   0.00225 5.62e-4   5599.
    ##  7 CD4 cm hig… (Interce… <NA>   -1.04     -0.544 -0.0184 0.101   1.93e-2   4656.
    ##  8 CD4 cm hig… typeheal… type   -3.08     -1.31   1.17   0.183   5.03e-2   3137.
    ##  9 CD4 cm rib… (Interce… <NA>   -0.0416    0.312  0.691  0.265   3.32e-2   4634.
    ## 10 CD4 cm rib… typeheal… type   -1.80     -0.970  0.0583 0.0695  2.08e-2   3588.
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

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

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

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

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

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->

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
    ##  1 B immature (Intercep… <NA>     0.327   0.868    1.45  0.00700 2.13e-3   2746.
    ##  2 B immature typehealt… type     0.658   1.33     1.98  0.00150 7.50e-4   2898.
    ##  3 B immature continuou… conti…  -0.222   0.0317   0.329 0.878   6.51e-1   7737.
    ##  4 B immature (Intercep… <NA>    -0.574  -0.0432   0.434 0.736   6.40e-1     NA 
    ##  5 B immature typehealt… <NA>    -0.434   0.0432   0.574 0.736   6.40e-1     NA 
    ##  6 B immature (Intercep… <NA>    -0.595  -0.0413   0.469 0.745   6.83e-1     NA 
    ##  7 B immature typehealt… <NA>    -0.469   0.0413   0.595 0.745   6.83e-1     NA 
    ##  8 B mem      (Intercep… <NA>    -1.20   -0.487    0.352 0.235   4.75e-2   2960.
    ##  9 B mem      typehealt… type     0.529   1.43     2.23  0.00700 3.00e-3   3028.
    ## 10 B mem      continuou… conti…  -0.255   0.0520   0.402 0.803   6.16e-1   7688.
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
    ##  1 10x_6K       B immature              0.0558           262.          -0.673 
    ##  2 10x_8K       B immature              0.144           1087.           0.403 
    ##  3 GSE115189    B immature              0.115            269.           0.105 
    ##  4 SCP345_580   B immature              0.0909           524.          -0.123 
    ##  5 SCP345_860   B immature              0.152            976.           0.461 
    ##  6 SCP424_pbmc1 B immature              0.113            300.           0.0707
    ##  7 SCP424_pbmc2 B immature              0.202            603.           0.799 
    ##  8 SCP591       B immature              0.0249            14.1         -1.49  
    ##  9 SI-GA-E5     B immature              0.0257           107.          -0.734 
    ## 10 SI-GA-E7     B immature              0.103            759.           0.759 
    ## # ℹ 590 more rows
