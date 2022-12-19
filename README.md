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

Modelling counts Modelling proportionality Modelling cell-type specific
variability Cell-type information share for variability shrinkage
Performing differential variability Probabilistic outlier identification
Cross-dataset learning (hyperpriors).

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
   .cell_group = cell_group 
  )
```

### From counts

``` r
res =
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3

    ## sccomp says: outlier identification second pass - step 2/3

    ## sccomp says: outlier-free model fitting - step 3/3

    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
res
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter   factor c_lower c_eff…¹ c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>       <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>    0.967    1.15   1.32   0       0         6461.
    ##  2 B1         typecancer  type   -1.16    -0.855 -0.537  7.51e-4 2.00e-4   3732.
    ##  3 B2         (Intercept) <NA>    0.505    0.791  1.07   0       0         4535.
    ##  4 B2         typecancer  type   -1.16    -0.726 -0.274  1.23e-2 1.30e-3   3683.
    ##  5 B3         (Intercept) <NA>   -0.567   -0.312 -0.0582 1.84e-1 1.58e-2   5212.
    ##  6 B3         typecancer  type   -0.556   -0.184  0.191  5.38e-1 1.57e-1   4926.
    ##  7 BM         (Intercept) <NA>   -1.24    -0.943 -0.682  0       0         4496.
    ##  8 BM         typecancer  type   -0.710   -0.311  0.0996 2.89e-1 8.08e-2   3733.
    ##  9 CD4 1      (Intercept) <NA>    0.211    0.402  0.591  2.03e-2 1.77e-3   4034.
    ## 10 CD4 1      typecancer  type   -0.0671   0.189  0.455  5.30e-1 1.43e-1   3403.
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
    contrasts =  c("typecancer - typebenign", "typebenign - typecancer"),
    .sample = sample,
    .cell_group = cell_group
  )
```

    ## sccomp says: outlier identification first pass - step 1/3

    ## sccomp says: outlier identification second pass - step 2/3

    ## sccomp says: outlier-free model fitting - step 3/3

    ## sccomp says: the composition design matrix has columns: typecancer, typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## # A tibble: 30 × 11
    ##    cell_group param…¹ factor v_lower v_eff…² v_upper v_pH0 v_FDR v_n_eff v_R_k…³
    ##    <chr>      <chr>   <chr>    <dbl>   <dbl>   <dbl> <dbl> <dbl>   <dbl>   <dbl>
    ##  1 B immature (Inter… <NA>     -5.10   -4.55   -3.94     0     0   5740.   1.00 
    ##  2 B mem      (Inter… <NA>     -4.82   -4.31   -3.72     0     0   4922.   1.00 
    ##  3 CD4 cm hi… (Inter… <NA>     -5.16   -4.60   -3.94     0     0   4556.   1.00 
    ##  4 CD4 cm ri… (Inter… <NA>     -5.87   -5.24   -4.56     0     0   5370.   1.00 
    ##  5 CD4 cm S1… (Inter… <NA>     -5.77   -5.12   -4.46     0     0   5318.   0.999
    ##  6 CD4 em hi… (Inter… <NA>     -5.55   -4.95   -4.30     0     0   4777.   1.00 
    ##  7 CD4 naive  (Inter… <NA>     -4.94   -4.38   -3.72     0     0   5196.   1.00 
    ##  8 CD4 ribos… (Inter… <NA>     -5.22   -4.70   -4.09     0     0   5149.   0.999
    ##  9 CD8 em 1   (Inter… <NA>     -5.41   -4.87   -4.23     0     0   5373.   1.00 
    ## 10 CD8 em 2   (Inter… <NA>     -5.35   -4.64   -3.85     0     0   4182.   0.999
    ## # … with 20 more rows, 1 more variable: count_data <list>, and abbreviated
    ## #   variable names ¹​parameter, ²​v_effect, ³​v_R_k_hat

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
```

    ## This is loo version 2.5.1

    ## - Online documentation and vignettes at mc-stan.org/loo

    ## - As of v2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session.

    ## 
    ## Attaching package: 'loo'

    ## The following object is masked from 'package:rstan':
    ## 
    ##     loo

``` r
# Fit first model
model_with_factor_association = 
  seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ type, 
   .sample =  sample, 
   .cell_group = cell_group, 
   check_outliers = FALSE
  )
```

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ 1, 
   .sample =  sample, 
   .cell_group = cell_group, 
   check_outliers = FALSE
  )
```

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept)

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
loo_with_factor_association <- loo(model_with_factor_association |> attr("fit"))
```

    ## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

``` r
loo_without_association <- loo(model_without_association |> attr("fit"))
```

    ## Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

``` r
# Compare models
loo_compare(loo_with_factor_association, loo_without_association)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -79.3      11.3

## Differential variability, binary factor

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
  )
```

    ## sccomp says: outlier identification first pass - step 1/3

    ## sccomp says: outlier identification second pass - step 2/3

    ## sccomp says: outlier-free model fitting - step 3/3

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept), typehealthy

    ## # A tibble: 60 × 18
    ##    cell_group     param…¹ factor c_lower c_eff…² c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>          <chr>   <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature     (Inter… <NA>     0.596   0.967   1.35  0       0         5904.
    ##  2 B immature     typehe… type     0.950   1.47    2.00  0       0         3846.
    ##  3 B mem          (Inter… <NA>    -1.64   -1.09   -0.508 0.00225 1.77e-4   4790.
    ##  4 B mem          typehe… type     1.29    2.00    2.71  0       0         4149.
    ##  5 CD4 cm high c… (Inter… <NA>    -0.823  -0.377   0.118 0.226   3.17e-2   4845.
    ##  6 CD4 cm high c… typehe… type    -3.28   -1.40    1.52  0.194   5.05e-2   2814.
    ##  7 CD4 cm riboso… (Inter… <NA>     0.167   0.501   0.847 0.0380  2.59e-3   4961.
    ##  8 CD4 cm riboso… typehe… type    -2.46   -1.75   -0.823 0.00100 2.50e-4   2949.
    ##  9 CD4 cm S100A4  (Inter… <NA>     1.77    2.02    2.25  0       0         8247.
    ## 10 CD4 cm S100A4  typehe… type     0.341   0.732   1.14  0.00450 1.25e-3   4657.
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
