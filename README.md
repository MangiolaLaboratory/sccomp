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

## (simple) Suggested for single-cell and CyTOF analyses

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

## (more complex and efficient, until further optimisation of the default installation) Suggested for microbiomics

**Github**

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
# Then, check the correct cmdstanr installation here
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# Then install sccomp with the cmdstanr branch
devtools::install_github("stemangiola/sccomp@cmdstanr")
```

# Analysis

## From Seurat Object

``` r
res =
  seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ type, 
    formula_variability = ~ 1, 
    sample, 
    cell_group 
  )
```

## From SingleCellExperiment Object

``` r
res =
  sce_obj |>
  sccomp_glm( 
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

    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
res
```

    ## # A tibble: 72 × 9
    ##    cell_group parameter   covariate c_lower c_effect c_upper   c_pH0   c_FDR
    ##    <chr>      <chr>       <chr>       <dbl>    <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>       0.557     0.710  0.866  0       0      
    ##  2 B1         typecancer  type      -1.19     -0.888 -0.587  0       0      
    ##  3 B2         (Intercept) <NA>       0.177     0.414  0.657  0.0368  0.00286
    ##  4 B2         typecancer  type      -1.18     -0.748 -0.312  0.00976 0.00177
    ##  5 B3         (Intercept) <NA>      -0.610    -0.409 -0.205  0.0225  0.00165
    ##  6 B3         typecancer  type      -0.593    -0.219  0.151  0.461   0.113  
    ##  7 BM         (Intercept) <NA>      -1.31     -1.10  -0.884  0       0      
    ##  8 BM         typecancer  type      -0.736    -0.344  0.0415 0.234   0.0643 
    ##  9 CD4 1      (Intercept) <NA>       0.362     0.490  0.627  0       0      
    ## 10 CD4 1      typecancer  type      -0.0882    0.161  0.420  0.622   0.163  
    ## # … with 62 more rows, and 1 more variable: count_data <list>

## Visualise data + inference

``` r
plots = plot_summary(res) 
```

    ## Joining, by = c("sample", "cell_group")
    ## Joining, by = c("cell_group", "type")

    ## Warning: Ignoring unknown aesthetics: label

Plot of group proportion, faceted by groups. The blue boxplots represent
the posterior predictive check. If the model is likely be descriptively
adequate to the data, the blue boxplot should roughly overlay with the
black boxplot, which represent the observed data. The outliers are
coloured in red.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-11-1.png)<!-- -->

Plot of estimates of differential composition (c\_) on the x axis, and
differential variability (v\_) on the y axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

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

    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept), typecancer

``` r
res
```

    ## # A tibble: 72 × 14
    ##    cell_group parameter   covariate c_lower c_effect c_upper    c_pH0     c_FDR
    ##    <chr>      <chr>       <chr>       <dbl>    <dbl>   <dbl>    <dbl>     <dbl>
    ##  1 B1         (Intercept) <NA>        0.627    0.796  0.972  0        0        
    ##  2 B1         typecancer  type       -1.26    -0.930 -0.602  0.000250 0.0000501
    ##  3 B2         (Intercept) <NA>        0.190    0.417  0.667  0.0315   0.00150  
    ##  4 B2         typecancer  type       -1.04    -0.561 -0.0636 0.0673   0.0108   
    ##  5 B3         (Intercept) <NA>       -0.557   -0.360 -0.150  0.0618   0.00526  
    ##  6 B3         typecancer  type       -0.518   -0.124  0.296  0.648    0.203    
    ##  7 BM         (Intercept) <NA>       -1.26    -1.05  -0.822  0        0        
    ##  8 BM         typecancer  type       -0.741   -0.298  0.163  0.323    0.0956   
    ##  9 CD4 1      (Intercept) <NA>        0.370    0.510  0.655  0        0        
    ## 10 CD4 1      typecancer  type       -0.117    0.176  0.470  0.569    0.139    
    ## # … with 62 more rows, and 6 more variables: v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, count_data <list>

Plot 1D significance plot

``` r
plots = plot_summary(res)
```

    ## Joining, by = c("sample", "cell_group")
    ## Joining, by = c("cell_group", "type")

    ## Warning: Ignoring unknown aesthetics: label

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

Plot 2D significance plot. This is possible if only differential
variability has been tested

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->
