sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

# <img src="inst/logo-01.png" height="139px" width="120px" />

Single-cell transcriptomics allows the unbiased characterisation of the
cellular composition of tissues. The cellular composition can be
compared between biological or clinical conditions to identify potential
cellular drivers. This strategy has been critical to unveil drivers of
immune response in cancer and pathogen infection from single-cell data.
Developing a robust statistical method for differential composition
analyses from single-cell data is crucial for driving discoveries. The
compositional data from single-cell experiments has four main
properties. The data is in count form; counts underlie inversely
correlated proportions that sum to one; larger cell groups are more
variable across samples than small groups; real-world data is rich in
outlier observation. A model that covers more than two of these
properties is currently lacking. **Here, we present a robust and
outlier-aware method for testing differential tissue composition from
single-cell data. This model can also transfer knowledge from a large
set of integrated datasets to increase accuracy further. We present how
this model can be applied to identify novel compositional and
heterogeneity changes in existing studies.**

# Installation

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) {
   install.packages("BiocManager")
 }
 BiocManager::install("sccomp")
```

**Github**

``` r
devtools::install_github("stemangiola/sccomp")
```

# Analysis

## From Seurat Object

``` r
res =
  seurat_obj |>
   formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

``` r
res =
  sce_obj |>
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From data.frame

``` r
res =
  seurat_obj[[]] |>
  sccomp_glm(
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From counts

``` r
res =
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ 1, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

``` r
res
```

    ## # A tibble: 72 × 8
    ##    cell_group parameter   c_lower c_effect c_upper    c_pH0     c_FDR count_data
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>    <dbl>     <dbl> <list>    
    ##  1 B1         (Intercept)  0.489     0.675  0.866  0        0         <tibble […
    ##  2 B1         typecancer  -1.29     -0.914 -0.536  0.000250 0.0000500 <tibble […
    ##  3 B2         (Intercept)  0.0399    0.281  0.511  0.245    0.0157    <tibble […
    ##  4 B2         typecancer  -1.04     -0.598 -0.136  0.0407   0.00507   <tibble […
    ##  5 B3         (Intercept) -0.731    -0.535 -0.315  0.00125  0.000152  <tibble […
    ##  6 B3         typecancer  -0.718    -0.312  0.0863 0.294    0.0803    <tibble […
    ##  7 BM         (Intercept) -1.39     -1.18  -0.958  0        0         <tibble […
    ##  8 BM         typecancer  -0.745    -0.336  0.0748 0.257    0.0597    <tibble […
    ##  9 CD4 1      (Intercept)  0.266     0.410  0.561  0.00275  0.000260  <tibble […
    ## 10 CD4 1      typecancer  -0.111     0.179  0.495  0.558    0.144     <tibble […
    ## # … with 62 more rows

## Visualise data + inference

``` r
plots = plot_summary(res) 
```

    ## Joining, by = c("sample", "cell_group")

    ## Joining, by = c("cell_group", "type")

Plot of group proportion, faceted by groups. The blue boxplots represent
the posterior predictive check. If the model is likely be descriptively
adequate to the data, the blue boxplot should roughly overlay with the
black boxplot, which represent the observed data. The outliers are
coloured in red.

``` r
plots$boxplot
```

![](inst/figures/unnamed-chunk-10-1.png)<!-- -->

Plot of estimates of differential composition (c\_) on the x axis, and
differential variability (v\_) on the y axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-11-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

## Differential variability

We can model the cell-group variability also dependent on type, and so
test differences in variability

``` r
res = 
  counts_obj |>
  sccomp_glm( 
    formula_composition = ~ type, 
    formula_variability = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

``` r
res
```

    ## # A tibble: 72 × 13
    ##    cell_group parameter   c_lower c_effect  c_upper    c_pH0     c_FDR v_lower
    ##    <chr>      <chr>         <dbl>    <dbl>    <dbl>    <dbl>     <dbl>   <dbl>
    ##  1 B1         (Intercept)  0.462     0.643  0.822   0        0           -4.80
    ##  2 B1         typecancer  -1.40     -1.04  -0.675   0.000250 0.0000833    1.12
    ##  3 B2         (Intercept)  0.0409    0.272  0.509   0.264    0.0250      -4.54
    ##  4 B2         typecancer  -1.14     -0.659 -0.194   0.0270   0.00428      1.49
    ##  5 B3         (Intercept) -0.733    -0.537 -0.319   0.00150  0.0000909   -5.57
    ##  6 B3         typecancer  -0.764    -0.333  0.104   0.263    0.0813       1.88
    ##  7 BM         (Intercept) -1.40     -1.19  -0.966   0        0           -6.26
    ##  8 BM         typecancer  -0.842    -0.427 -0.00467 0.149    0.0449       1.21
    ##  9 CD4 1      (Intercept)  0.246     0.399  0.553   0.00575  0.000437    -5.33
    ## 10 CD4 1      typecancer  -0.198     0.107  0.414   0.725    0.245        1.55
    ## # … with 62 more rows, and 5 more variables: v_effect <dbl>, v_upper <dbl>,
    ## #   v_pH0 <dbl>, v_FDR <dbl>, count_data <list>

Plot 1D significance plot

``` r
plots = plot_summary(res)
```

    ## Joining, by = c("sample", "cell_group")

    ## Joining, by = c("cell_group", "type")

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

Plot 2D significance plot. This is possible if only differential
variability has been tested

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->
