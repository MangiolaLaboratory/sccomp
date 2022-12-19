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
Cross-dataset learning (hyperpriors)

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

## From Seurat, SingleCellExperiment, metadata objects

``` r
single_cell_object |>
  sccomp_glm( 
   formula_composition = ~ type, 
   .sample =  sample, 
   .cell_group = cell_group 
  )
```

## From counts

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

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

    ## sccomp says: the composition design matrix has columns: (Intercept), typecancer

    ## sccomp says: the variability design matrix has columns: (Intercept)

``` r
res
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter   covar…¹ c_lower c_eff…² c_upper  c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>       <chr>     <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>     0.969    1.15   1.31   0      0         5023.
    ##  2 B1         typecancer  type    -1.18    -0.859 -0.531  0      0         4508.
    ##  3 B2         (Intercept) <NA>     0.500    0.784  1.06   0      0         5341.
    ##  4 B2         typecancer  type    -1.17    -0.735 -0.280  0.0110 0.00141   3916.
    ##  5 B3         (Intercept) <NA>    -0.570   -0.315 -0.0648 0.183  0.0172    5835.
    ##  6 B3         typecancer  type    -0.577   -0.184  0.190  0.533  0.148     5340.
    ##  7 BM         (Intercept) <NA>    -1.22    -0.941 -0.679  0      0         5737.
    ##  8 BM         typecancer  type    -0.740   -0.310  0.0910 0.288  0.0856    5186.
    ##  9 CD4 1      (Intercept) <NA>     0.203    0.406  0.582  0.0223 0.00167   4057.
    ## 10 CD4 1      typecancer  type    -0.0607   0.180  0.456  0.562  0.163     3437.
    ## # … with 62 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​covariate, ²​c_effect

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

## Differential variability

``` r
seurat_obj |>
  sccomp_glm( 
   formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept), typehealthy

    ## # A tibble: 60 × 18
    ##    cell_group    param…¹ covar…² c_lower c_eff…³ c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>         <chr>   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature    (Inter… <NA>      0.579   0.955   1.34  0       0         5406.
    ##  2 B immature    typehe… type      0.975   1.49    2.02  0       0         4021.
    ##  3 B mem         (Inter… <NA>     -1.28   -0.744  -0.132 0.0380  2.59e-3   4782.
    ##  4 B mem         typehe… type      0.937   1.66    2.38  0       0         4190.
    ##  5 CD4 cm high … (Inter… <NA>     -0.856  -0.401   0.105 0.207   2.49e-2   4442.
    ##  6 CD4 cm high … typehe… type     -3.18   -1.39    1.34  0.188   4.27e-2   3269.
    ##  7 CD4 cm ribos… (Inter… <NA>      0.153   0.484   0.823 0.0473  6.84e-3   4347.
    ##  8 CD4 cm ribos… typehe… type     -2.47   -1.72   -0.808 0.00225 8.01e-4   3980.
    ##  9 CD4 cm S100A4 (Inter… <NA>      1.74    2.00    2.24  0       0         6252.
    ## 10 CD4 cm S100A4 typehe… type      0.361   0.759   1.16  0.00250 1.08e-3   4833.
    ## # … with 50 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​parameter, ²​covariate, ³​c_effect

## With contrasts

``` r
seurat_obj |>
  sccomp_glm( 
    formula_composition = ~ 0 + type, 
    contrasts =  c("typecancer - typebenign", "typebenign - typecancer"),
    .sample = sample,
    .cell_group = cell_group
  )
```

    ## sccomp says: outlier identification first pass - step 1/3 [ETA: ~20s]

    ## sccomp says: outlier identification second pass - step 2/3 [ETA: ~60s]

    ## sccomp says: outlier-free model fitting - step 3/3 [ETA: ~20s]

    ## sccomp says: the composition design matrix has columns: typecancer, typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## # A tibble: 30 × 11
    ##    cell_gr…¹ param…² covar…³ v_lower v_eff…⁴ v_upper v_pH0 v_FDR v_n_eff v_R_k…⁵
    ##    <chr>     <chr>   <chr>     <dbl>   <dbl>   <dbl> <dbl> <dbl>   <dbl>   <dbl>
    ##  1 B immatu… (Inter… <NA>      -5.08   -4.56   -3.95     0     0   5559.   0.999
    ##  2 B mem     (Inter… <NA>      -4.82   -4.30   -3.73     0     0   5014.   1.00 
    ##  3 CD4 cm h… (Inter… <NA>      -5.17   -4.59   -3.94     0     0   4397.   0.999
    ##  4 CD4 cm r… (Inter… <NA>      -5.86   -5.24   -4.55     0     0   5086.   0.999
    ##  5 CD4 cm S… (Inter… <NA>      -5.77   -5.13   -4.46     0     0   5423.   0.999
    ##  6 CD4 em h… (Inter… <NA>      -5.54   -4.94   -4.30     0     0   4753.   0.999
    ##  7 CD4 naive (Inter… <NA>      -4.93   -4.38   -3.73     0     0   5521.   0.999
    ##  8 CD4 ribo… (Inter… <NA>      -5.20   -4.69   -4.10     0     0   4992.   0.999
    ##  9 CD8 em 1  (Inter… <NA>      -5.41   -4.87   -4.24     0     0   5639.   1.00 
    ## 10 CD8 em 2  (Inter… <NA>      -5.26   -4.61   -3.85     0     0   3693.   1.00 
    ## # … with 20 more rows, 1 more variable: count_data <list>, and abbreviated
    ## #   variable names ¹​cell_group, ²​parameter, ³​covariate, ⁴​v_effect, ⁵​v_R_k_hat

## Model comparison (e.g. Bayesian ANOVA) with `loo`

In the following example, the model with association with factors better
fits the data compared to the baseline model with no factor association.
For comparisons `check_outliers` must be set to FALSE as the
leave-one-out must work with the same amount of data, while outlier
elimination does not guarantee it.

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

    ## sccomp says: estimation [ETA: ~20s]

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

    ## sccomp says: estimation [ETA: ~20s]

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
    ## model2 -80.7      11.4

## Suggested settings for single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## Suggested settings for CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

## Visualise Data + inference

``` r
plots = plot_summary(res) 
```

    ## Joining, by = c("sample", "cell_group")
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

![](inst/figures/unnamed-chunk-11-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
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
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

## Differential variability

We can model the cell-group variability also dependent on the type, and
so test differences in variability

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

    ## # A tibble: 72 × 18
    ##    cell_group parameter   covar…¹ c_lower c_eff…² c_upper  c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>       <chr>     <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercept) <NA>      0.975  1.16    1.35   0      0         5071.
    ##  2 B1         typecancer  type     -1.36  -1.02   -0.682  0      0         3424.
    ##  3 B2         (Intercept) <NA>      0.523  0.774   1.02   0      0         5474.
    ##  4 B2         typecancer  type     -1.28  -0.805  -0.282  0.0128 0.00267   3351.
    ##  5 B3         (Intercept) <NA>     -0.542 -0.308  -0.0571 0.187  0.0228    4382.
    ##  6 B3         typecancer  type     -0.712 -0.307   0.125  0.299  0.116     3273.
    ##  7 BM         (Intercept) <NA>     -1.21  -0.912  -0.604  0      0         4838.
    ##  8 BM         typecancer  type     -0.893 -0.475  -0.0520 0.0988 0.0177    4144.
    ##  9 CD4 1      (Intercept) <NA>      0.230  0.413   0.604  0.0133 0.00114   5290.
    ## 10 CD4 1      typecancer  type     -0.216  0.0433  0.330  0.872  0.346     2632.
    ## # … with 62 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​covariate, ²​c_effect

Plot 1D significance plot

``` r
plots = plot_summary(res)
```

    ## Joining, by = c("sample", "cell_group")
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
