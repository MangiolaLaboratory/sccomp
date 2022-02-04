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
    ##    cell_group parameter   c_lower c_effect c_upper   c_pH0    c_FDR count_data  
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl> <list>      
    ##  1 B1         (Intercept)   0.465    0.639  0.815  0       0        <tibble [20…
    ##  2 B1         typecancer   -1.22    -0.885 -0.560  0       0        <tibble [20…
    ##  3 B2         (Intercept)   0.120    0.336  0.547  0.104   0.0176   <tibble [20…
    ##  4 B2         typecancer   -1.16    -0.746 -0.344  0.00450 0.000775 <tibble [20…
    ##  5 B3         (Intercept)  -0.668   -0.484 -0.281  0.00400 0.000375 <tibble [20…
    ##  6 B3         typecancer   -0.584   -0.217  0.131  0.462   0.122    <tibble [20…
    ##  7 BM         (Intercept)  -1.38    -1.17  -0.974  0       0        <tibble [20…
    ##  8 BM         typecancer   -0.750   -0.353  0.0410 0.214   0.0521   <tibble [20…
    ##  9 CD4 1      (Intercept)   0.276    0.423  0.576  0.00150 0.000136 <tibble [20…
    ## 10 CD4 1      typecancer   -0.122    0.168  0.460  0.590   0.191    <tibble [20…
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
    ##    cell_group parameter   c_lower c_effect c_upper   c_pH0    c_FDR v_lower
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
    ##  1 B1         (Intercept)   0.452   0.612   0.782  0       0          -5.17
    ##  2 B1         typecancer   -1.27   -0.941  -0.586  0       0           1.25
    ##  3 B2         (Intercept)   0.171   0.366   0.572  0.0513  0.00246    -4.96
    ##  4 B2         typecancer   -1.06   -0.660  -0.252  0.0128  0.00156     1.65
    ##  5 B3         (Intercept)  -0.682  -0.495  -0.301  0.00250 0.000130   -5.86
    ##  6 B3         typecancer   -0.620  -0.232   0.178  0.435   0.148       1.44
    ##  7 BM         (Intercept)  -1.40   -1.18   -0.974  0       0          -6.40
    ##  8 BM         typecancer   -0.822  -0.411   0.0101 0.154   0.0384      1.05
    ##  9 CD4 1      (Intercept)   0.230   0.379   0.531  0.00725 0.000427   -5.48
    ## 10 CD4 1      typecancer   -0.220   0.0762  0.393  0.778   0.273       1.13
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
