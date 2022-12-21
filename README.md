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

-   Modelling counts
-   Modelling proportionality
-   Modelling cell-type specific variability
-   Cell-type information share for variability shrinkage
-   Testing differential variability
-   Probabilistic outlier identification
-   Cross-dataset learning (hyperpriors).

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
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
```

### From counts

``` r
counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter   factor c_lower c_eff…¹ c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>       <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>     0.869   1.06   1.24   0       0         4529.
    ##  2 B1         typecancer  type    -1.21   -0.878 -0.538  0       0         4296.
    ##  3 B2         (Intercept) <NA>     0.423   0.697  0.959  5.00e-4 3.75e-5   3358.
    ##  4 B2         typecancer  type    -1.21   -0.784 -0.363  3.75e-3 6.11e-4   4084.
    ##  5 B3         (Intercept) <NA>    -0.643  -0.383 -0.136  7.17e-2 4.01e-3   5204.
    ##  6 B3         typecancer  type    -0.616  -0.224  0.137  4.54e-1 1.37e-1   4037.
    ##  7 BM         (Intercept) <NA>    -1.27   -1.02  -0.754  0       0         5698.
    ##  8 BM         typecancer  type    -0.764  -0.347  0.0320 2.24e-1 3.62e-2   4352.
    ##  9 CD4 1      (Intercept) <NA>     0.121   0.324  0.522  1.03e-1 1.15e-2   3652.
    ## 10 CD4 1      typecancer  type    -0.110   0.157  0.425  6.36e-1 2.19e-1   2916.
    ## # … with 62 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable name ¹​c_effect

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

## Contrasts

``` r
seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 0 + type, 
    contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )
```

    ## # A tibble: 60 × 18
    ##    cell_group     param…¹ factor c_lower c_eff…² c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>          <chr>   <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature     typeca… type    -2.01   -1.57   -1.12  0       0            NA
    ##  2 B immature     typehe… type     1.12    1.57    2.01  0       0            NA
    ##  3 B mem          typeca… type    -2.31   -1.71   -1.11  0       0            NA
    ##  4 B mem          typehe… type     1.11    1.71    2.31  0       0            NA
    ##  5 CD4 cm high c… typeca… type     0.940   1.86    3.00  0       0            NA
    ##  6 CD4 cm high c… typehe… type    -3.00   -1.86   -0.940 0       0            NA
    ##  7 CD4 cm riboso… typeca… type     0.571   1.24    1.97  1.50e-3 3.89e-4      NA
    ##  8 CD4 cm riboso… typehe… type    -1.97   -1.24   -0.571 1.50e-3 3.89e-4      NA
    ##  9 CD4 cm S100A4  typeca… type    -1.24   -0.882  -0.540 2.50e-4 5.00e-5      NA
    ## 10 CD4 cm S100A4  typehe… type     0.540   0.882   1.24  2.50e-4 5.00e-5      NA
    ## # … with 50 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​parameter, ²​c_effect

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
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    check_outliers = FALSE, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )

# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    check_outliers = FALSE, 
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

## Differential variability, binary factor

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
res = 
  seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  )

res
```

    ## # A tibble: 60 × 18
    ##    cell_group     param…¹ factor c_lower c_eff…² c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>          <chr>   <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature     (Inter… <NA>     0.582   0.965   1.33  5.00e-4 5.36e-5   5635.
    ##  2 B immature     typehe… type     0.967   1.49    2.00  0       0         3507.
    ##  3 B mem          (Inter… <NA>    -1.65   -1.08   -0.467 3.00e-3 4.03e-4   3873.
    ##  4 B mem          typehe… type     1.24    2.00    2.75  2.50e-4 8.33e-5   3301.
    ##  5 CD4 cm high c… (Inter… <NA>    -0.831  -0.390   0.103 2.08e-1 3.12e-2   4146.
    ##  6 CD4 cm high c… typehe… type    -3.19   -1.39    1.51  1.86e-1 4.81e-2   2131.
    ##  7 CD4 cm riboso… (Inter… <NA>     0.166   0.491   0.850 4.10e-2 2.82e-3   3992.
    ##  8 CD4 cm riboso… typehe… type    -2.45   -1.74   -0.841 2.25e-3 6.25e-4   4160.
    ##  9 CD4 cm S100A4  (Inter… <NA>     1.76    2.01    2.26  0       0         6267.
    ## 10 CD4 cm S100A4  typehe… type     0.337   0.743   1.16  2.50e-3 1.00e-3   3289.
    ## # … with 50 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​parameter, ²​c_effect

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

    ## Joining, by = c("cell_group", "sample")
    ## Joining, by = c("cell_group", "type")

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

    ## Joining, by = c("cell_group", "sample")
    ## Joining, by = c("cell_group", "type")

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
