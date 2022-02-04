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
    ##  1 B1         (Intercept)  0.565     0.759  0.957  0       0        <tibble [20…
    ##  2 B1         typecancer  -1.14     -0.758 -0.382  0.00275 0.00112  <tibble [20…
    ##  3 B2         (Intercept)  0.130     0.372  0.631  0.0735  0.00539  <tibble [20…
    ##  4 B2         typecancer  -1.13     -0.669 -0.198  0.0260  0.00339  <tibble [20…
    ##  5 B3         (Intercept) -0.744    -0.528 -0.316  0.00175 0.000228 <tibble [20…
    ##  6 B3         typecancer  -0.708    -0.308  0.0761 0.292   0.0959   <tibble [20…
    ##  7 BM         (Intercept) -1.39     -1.17  -0.965  0       0        <tibble [20…
    ##  8 BM         typecancer  -0.742    -0.326  0.0882 0.265   0.0774   <tibble [20…
    ##  9 CD4 1      (Intercept)  0.271     0.417  0.568  0.00200 0.000302 <tibble [20…
    ## 10 CD4 1      typecancer  -0.0930    0.182  0.476  0.551   0.134    <tibble [20…
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
plots$plot_associations
```

    ## NULL

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

    ## # A tibble: 72 × 13
    ##    cell_group parameter   c_lower c_effect c_upper   c_pH0    c_FDR v_lower
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
    ##  1 B1         (Intercept)  0.460     0.646  0.825  0       0           4.83
    ##  2 B1         typecancer  -1.39     -1.03  -0.672  0       0          -1.07
    ##  3 B2         (Intercept)  0.0476    0.274  0.502  0.253   0.0232      4.57
    ##  4 B2         typecancer  -1.12     -0.668 -0.202  0.0248  0.00411    -1.52
    ##  5 B3         (Intercept) -0.740    -0.537 -0.317  0.00150 0.000141    5.60
    ##  6 B3         typecancer  -0.773    -0.339  0.0868 0.251   0.0871     -1.84
    ##  7 BM         (Intercept) -1.40     -1.18  -0.969  0       0           6.27
    ##  8 BM         typecancer  -0.859    -0.433 -0.0111 0.148   0.0399     -1.21
    ##  9 CD4 1      (Intercept)  0.255     0.403  0.559  0.00600 0.000385    5.31
    ## 10 CD4 1      typecancer  -0.200     0.107  0.424  0.728   0.263      -1.51
    ## # … with 62 more rows, and 5 more variables: v_effect <dbl>, v_upper <dbl>,
    ## #   v_pH0 <dbl>, v_FDR <dbl>, count_data <list>
