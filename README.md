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
    ##  1 B1         (Intercept)   0.573    0.753  0.940  0        0         <tibble>  
    ##  2 B1         typecancer   -1.29    -0.943 -0.599  0        0         <tibble>  
    ##  3 B2         (Intercept)   0.119    0.367  0.627  0.0948   0.00509   <tibble>  
    ##  4 B2         typecancer   -1.08    -0.633 -0.160  0.0353   0.00432   <tibble>  
    ##  5 B3         (Intercept)  -0.609   -0.420 -0.211  0.0213   0.00189   <tibble>  
    ##  6 B3         typecancer   -0.600   -0.209  0.162  0.480    0.153     <tibble>  
    ##  7 BM         (Intercept)  -1.32    -1.11  -0.898  0        0         <tibble>  
    ##  8 BM         typecancer   -0.734   -0.342  0.0521 0.249    0.0649    <tibble>  
    ##  9 CD4 1      (Intercept)   0.328    0.456  0.583  0.000250 0.0000119 <tibble>  
    ## 10 CD4 1      typecancer   -0.117    0.113  0.351  0.759    0.243     <tibble>  
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

``` r
res
```

    ## # A tibble: 72 × 13
    ##    cell_group parameter   c_lower c_effect c_upper   c_pH0    c_FDR v_lower
    ##    <chr>      <chr>         <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
    ##  1 B1         (Intercept)  0.613     0.788   0.962 0       0         -4.89 
    ##  2 B1         typecancer  -1.27     -0.937  -0.591 0       0          0.973
    ##  3 B2         (Intercept)  0.239     0.460   0.688 0.00901 0.000635  -4.68 
    ##  4 B2         typecancer  -1.11     -0.660  -0.198 0.0253  0.00405    2.04 
    ##  5 B3         (Intercept) -0.618    -0.408  -0.191 0.0290  0.00239   -5.55 
    ##  6 B3         typecancer  -0.660    -0.222   0.232 0.464   0.142      2.00 
    ##  7 BM         (Intercept) -1.28     -1.05   -0.824 0       0         -6.10 
    ##  8 BM         typecancer  -0.776    -0.311   0.155 0.305   0.0713     1.22 
    ##  9 CD4 1      (Intercept)  0.393     0.529   0.679 0       0         -5.47 
    ## 10 CD4 1      typecancer  -0.0704    0.223   0.509 0.436   0.130      1.74 
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

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

Plot 2D significance plot. This is possible if only differential
variability has been tested

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->
