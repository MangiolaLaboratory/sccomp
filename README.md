sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

# <img src="inst/logo-01.png" height="139px" width="120px" />

Cell omics such as single-cell genomics, proteomics and microbiomics
allow the characterisation of tissue and microbial community
composition, which can be compared between conditions to identify
biological drivers. This strategy has been critical to unveiling markers
of disease progression such as cancer and pathogen infection. For cell
omic data, no method for differential variability analysis exists, and
methods for differential composition analysis only take a few
fundamental data properties into account. Here we introduce sccomp, a
generalised method for differential composition and variability analyses
able to jointly model data count distribution, compositionality,
group-specific variability and proportion mean-variability association,
with awareness against outliers. Sccomp is an extensive analysis
framework that allows realistic data simulation and cross-study
knowledge transfer. Here, we demonstrate that mean-variability
association is ubiquitous across technologies showing the inadequacy of
the very popular Dirichlet-multinomial modelling and provide mandatory
principles for differential variability analysis. We show that sccomp
accurately fits experimental data, with a 50% incremental improvement
over state-of-the-art algorithms. Using sccomp, we identified novel
differential constraints and composition in the microenvironment of
primary breast cancer.

# Installation

## (simple) Suggested for single-cell and CyTOF analyses

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")
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

`sccomp` can model changes in composition and variability. Normally the
furmula for variability is either `~1`, which assumes that the
cell-group variability is independent on any covariate, or
`~ factor_of_interest`, which assumes that the model is dependent on the
factor of interest only. However, more complex models for variability
are possible, is the sample size is large. In any case the model for
variability must be a subset of the model for composition.

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
    ##    cell_group parameter  covar…¹ c_lower c_eff…² c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>     0.965    1.15   1.31   0       0         4987.
    ##  2 B1         typecancer type    -1.18    -0.865 -0.551  2.50e-4 4.17e-5   3625.
    ##  3 B2         (Intercep… <NA>     0.494    0.785  1.05   5.01e-4 4.77e-5   5662.
    ##  4 B2         typecancer type    -1.18    -0.738 -0.303  9.01e-3 1.41e-3   4150.
    ##  5 B3         (Intercep… <NA>    -0.567   -0.312 -0.0684 1.80e-1 1.74e-2   4114.
    ##  6 B3         typecancer type    -0.560   -0.195  0.183  5.08e-1 1.57e-1   4593.
    ##  7 BM         (Intercep… <NA>    -1.22    -0.945 -0.671  0       0         4182.
    ##  8 BM         typecancer type    -0.721   -0.311  0.0905 2.92e-1 7.83e-2   4421.
    ##  9 CD4 1      (Intercep… <NA>     0.205    0.401  0.587  2.23e-2 2.38e-3   4064.
    ## 10 CD4 1      typecancer type    -0.0764   0.179  0.442  5.64e-1 1.72e-1   4342.
    ## # … with 62 more rows, 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>,
    ## #   v_R_k_hat <dbl>, count_data <list>, and abbreviated variable names
    ## #   ¹​covariate, ²​c_effect

Of the output table, the estimate columns star twith the prefix `c_`
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
    ##  1 B immature    (Inter… <NA>      0.584   0.950   1.31  0       0         6357.
    ##  2 B immature    typehe… type      0.966   1.49    2.02  0       0         3777.
    ##  3 B mem         (Inter… <NA>     -1.27   -0.749  -0.147 0.0368  2.59e-3   5002.
    ##  4 B mem         typehe… type      0.934   1.67    2.38  0       0         4563.
    ##  5 CD4 cm high … (Inter… <NA>     -0.847  -0.397   0.103 0.203   2.42e-2   3671.
    ##  6 CD4 cm high … typehe… type     -3.14   -1.38    1.38  0.191   4.69e-2   2966.
    ##  7 CD4 cm ribos… (Inter… <NA>      0.155   0.485   0.832 0.0450  6.54e-3   4872.
    ##  8 CD4 cm ribos… typehe… type     -2.43   -1.72   -0.783 0.00400 1.21e-3   1938.
    ##  9 CD4 cm S100A4 (Inter… <NA>      1.74    2.00    2.24  0       0         5896.
    ## 10 CD4 cm S100A4 typehe… type      0.369   0.752   1.15  0.00225 6.51e-4   4150.
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
    ##  1 B immatu… (Inter… <NA>      -5.05   -4.55   -3.96     0     0   5518.   0.999
    ##  2 B mem     (Inter… <NA>      -4.85   -4.33   -3.77     0     0   4581.   1.00 
    ##  3 CD4 cm h… (Inter… <NA>      -5.16   -4.60   -3.96     0     0   4133.   0.999
    ##  4 CD4 cm r… (Inter… <NA>      -5.84   -5.23   -4.53     0     0   5845.   1.00 
    ##  5 CD4 cm S… (Inter… <NA>      -5.74   -5.08   -4.42     0     0   5603.   0.999
    ##  6 CD4 em h… (Inter… <NA>      -5.53   -4.94   -4.31     0     0   4522.   1.00 
    ##  7 CD4 naive (Inter… <NA>      -5.19   -4.60   -3.96     0     0   5404.   0.999
    ##  8 CD4 ribo… (Inter… <NA>      -5.21   -4.69   -4.07     0     0   4907.   1.00 
    ##  9 CD8 em 1  (Inter… <NA>      -5.41   -4.85   -4.23     0     0   5386.   1.00 
    ## 10 CD8 em 2  (Inter… <NA>      -5.34   -4.64   -3.88     0     0   3577.   1.00 
    ## # … with 20 more rows, 1 more variable: count_data <list>, and abbreviated
    ## #   variable names ¹​cell_group, ²​parameter, ³​covariate, ⁴​v_effect, ⁵​v_R_k_hat

## Model comparison (e.g. Bayesian ANOVA) with `loo`

In the following example the model with association with factors better
fits the data compares to the baseline model with no factor association.
For comparisons `check_outliers` must be set to FALSE as the
leave-one-out must work with the same amount of data, while outlier
elimination does not guarantee it.

If `elpd_diff` is away from zero of \> 5 `se_diff` difference of 5 we
are confident that a model is better than the other
[reference](https://discourse.mc-stan.org/t/interpreting-elpd-diff-loo-package/1628/2?u=stemangiola).
In this case -79.9 / 11.5 = -6.9, therefore we can conclude that model
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
    ## model2 -79.4      11.6

## Suggested settings for single-cell RNA sequencing

We recommend to set `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## Suggested settings for CyTOF and microbiome data

We recommend to set `bimodal_mean_variability_association  = FALSE`
(Default).

## Visualise data + inference

``` r
plots = plot_summary(res) 
```

    ## Joining, by = c("sample", "cell_group")
    ## Joining, by = c("cell_group", "type")

    ## Warning: Expected 2 pieces. Additional pieces discarded in 4 rows [6, 7, 13,
    ## 14].

Plot of group proportion, faceted by groups. The blue boxplots represent
the posterior predictive check. If the model is likely be descriptively
adequate to the data, the blue box plot should roughly overlay with the
black box plot, which represent the observed data. The outliers are
coloured in red. A box plot will be returned for every (discrete)
covariates present in `formula_composition`. The color coding represent
the significant associations for composition and/or variability.

``` r
plots$boxplot
```

    ## [[1]]

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

Plot of estimates of differential composition (c\_) on the x axis, and
differential variability (v\_) on the y axis. The error bars represent
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
example we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

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

    ## # A tibble: 72 × 18
    ##    cell_group parameter  covar…¹ c_lower c_eff…² c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>      0.969  1.15    1.34   0       0         4854.
    ##  2 B1         typecancer type     -1.35  -1.02   -0.685  0       0         3819.
    ##  3 B2         (Intercep… <NA>      0.538  0.768   1.01   2.50e-4 1.19e-5   4958.
    ##  4 B2         typecancer type     -1.31  -0.800  -0.241  2.08e-2 5.31e-3   3335.
    ##  5 B3         (Intercep… <NA>     -0.558 -0.318  -0.0616 1.67e-1 1.97e-2   4432.
    ##  6 B3         typecancer type     -0.720 -0.304   0.119  3.10e-1 1.18e-1   3458.
    ##  7 BM         (Intercep… <NA>     -1.21  -0.925  -0.628  0       0         5407.
    ##  8 BM         typecancer type     -0.892 -0.472  -0.0640 9.41e-2 1.98e-2   3765.
    ##  9 CD4 1      (Intercep… <NA>      0.213  0.398   0.583  1.93e-2 1.66e-3   4963.
    ## 10 CD4 1      typecancer type     -0.216  0.0535  0.340  8.58e-1 3.40e-1   3075.
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

![](inst/figures/unnamed-chunk-16-1.png)<!-- -->

Plot 2D significance plot. Data points are cell groups. Error bars are
the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis,
that no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (beside the Intercept) is
provided for each of the categories of interest. The for each category
of interest, the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-17-1.png)<!-- -->
