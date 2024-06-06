sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)

<!-- badges: end -->

<a href="https://www.youtube.com/watch?v=R_lt58We9nA&ab_channel=RConsortium" target="_blank">
<img src="https://img.youtube.com/vi/R_lt58We9nA/mqdefault.jpg" alt="Watch the video" width="280" height="180" border="10" />
</a>

# <img src="inst/logo-01.png" height="139px" width="120px"/>

`sccomp` tests differences in cell type proportions from single-cell
data. It is robust against outliers, it models continuous and discrete
factors, and capable of random-effect/intercept modelling.

Please cite [PNAS - sccomp: Robust differential composition and
variability analysis for single-cell
data](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)

## Characteristics

- Complex linear models with continuous and categorical covariates
- Multilevel modelling, with population fixed and random
  effects/intercept
- Modelling data from counts
- Testing differences in cell-type proportionality
- Testing differences in cell-type specific variability
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

| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                  | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`           | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                      | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`                 | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                   | Predicts proportions, based on the mode, or part of the model                                                               |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plors summary plots to asses significance                                                                                   |

# Analysis

`sccomp` can model changes in composition and variability. By default,
the formula for variability is either `~1`, which assumes that the
cell-group variability is independent of any covariate or
`~ factor_of_interest`, which assumes that the model is dependent on the
factor of interest only. The variability model must be a subset of the
model for composition.

## Binary factor

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

### From Seurat, SingleCellExperiment, metadata objects

``` r
sccomp_result = 
  single_cell_object |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
```

### From counts

``` r
sccomp_result = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_remove_outliers(cores = 1, verbose = FALSE) |> # Optional
  sccomp_test()
```

Here you see the results of the fit, the effects of the factor on
composition and variability. You also can see the uncertainty around
those effects.

``` r
sccomp_result
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.842     1.08    1.32  0       0           NaN
    ##  2 B1         typecancer type   -0.966    -0.611  -0.207 0.00800 8.18e-4     NaN
    ##  3 B2         (Intercep… <NA>    0.499     0.816   1.12  0       0           NaN
    ##  4 B2         typecancer type   -1.21     -0.823  -0.413 0.00100 1.00e-4     NaN
    ##  5 B3         (Intercep… <NA>   -0.872    -0.525  -0.235 0       0           NaN
    ##  6 B3         typecancer type   -0.684    -0.249   0.220 0.268   5.12e-2     NaN
    ##  7 BM         (Intercep… <NA>   -1.33     -1.03   -0.716 0       0           NaN
    ##  8 BM         typecancer type   -0.660    -0.221   0.206 0.295   6.06e-2     NaN
    ##  9 CD4 1      (Intercep… <NA>    0.107     0.298   0.485 0.0200  1.52e-3     NaN
    ## 10 CD4 1      typecancer type   -0.0204    0.238   0.494 0.155   3.63e-2     NaN
    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

## An aid to result interpretation and communication

The estimated effects are expressed in the unconstrained space of the
parameters. Similarly, to differential expression analysis that express
change in terms of log fold change. However, for differences, in
proportion, logit foold change must be used. This measure is harder to
interpret and understand.

Therefore, we provide a more intuitive proportion, full change, that can
be easier understood. However, these cannot be used to infer
significance (use sccomp_test() instead), and a lot of care must be
taken given the nonlinearity of these measure (1 fold increase from
0.0001 to 0.0002 carried a different weight that 1 fold increase from
0.4 to 0.8).

From your estimates, you can state which effects you are interested
about (this can be a part of the full model, in case you want to not
consider unwanted effects), and the two points you would like to
compare.

In case of a chategorical variable, the starting and ending points are
categories.

``` r
sccomp_result |> 
   sccomp_proportional_fold_change(
     formula_composition = ~  type,
     from =  "healthy", 
     to = "cancer"
    ) |> 
  select(cell_group, statement)
```

    ## # A tibble: 36 × 2
    ##    cell_group statement                                
    ##    <chr>      <glue>                                   
    ##  1 B1         1.9-fold decrease (from 0.054 to 0.0288) 
    ##  2 B2         2.4-fold decrease (from 0.0421 to 0.0179)
    ##  3 B3         1.3-fold decrease (from 0.0109 to 0.0082)
    ##  4 BM         1.3-fold decrease (from 0.0066 to 0.0051)
    ##  5 CD4 1      1.2-fold increase (from 0.0247 to 0.0304)
    ##  6 CD4 2      1.4-fold increase (from 0.0533 to 0.0751)
    ##  7 CD4 3      2.6-fold decrease (from 0.0848 to 0.0321)
    ##  8 CD4 4      1-fold increase (from 0.0018 to 0.0018)  
    ##  9 CD4 5      1.1-fold increase (from 0.0292 to 0.0323)
    ## 10 CD8 1      1.3-fold increase (from 0.0999 to 0.1339)
    ## # ℹ 26 more rows

## Summary plots

``` r
plots = sccomp_result |> plot() 
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

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

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated, you also appreciate that
this relationship is by model for single cell RNA sequencing data.

`sccomp` models, these relationship to obtain a shrinkage effect on the
estimates of both the abundance and the variability. This shrinkage is
adaptive as it is modelled jointly, thanks for Bayesian inference.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

## Contrasts

``` r
seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test( contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"))
```

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -1.67    -1.18   -0.672 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.672    1.18    1.67  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.25    -1.63   -1.05  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.05     1.63    2.25  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.55    -1.11   -0.707 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.707    1.11    1.55  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.694    1.68    2.70  0.00100 1.43e-4      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.70    -1.68   -0.694 0.00100 1.43e-4      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.212    0.913   1.66  0.0120  2.33e-3      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.66    -0.913  -0.212 0.0120  2.33e-3      NA
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
    cores = 1, verbose = FALSE
  )

res
```

    ## # A tibble: 60 × 14
    ##    cell_group        parameter factor c_lower c_effect c_upper c_n_eff c_R_k_hat
    ##    <chr>             <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>     <dbl>
    ##  1 B immature        (Interce… <NA>     0.559    0.898  1.25       NaN      3.08
    ##  2 B immature        typeheal… type     0.885    1.29   1.74       NaN      2.93
    ##  3 B mem             (Interce… <NA>    -1.49    -1.02  -0.502      NaN      3.12
    ##  4 B mem             typeheal… type     1.53     2.13   2.77       NaN      3.11
    ##  5 CD4 cm S100A4     (Interce… <NA>     1.49     1.75   2.03       NaN      2.99
    ##  6 CD4 cm S100A4     typeheal… type     0.558    0.942  1.34       NaN      3.03
    ##  7 CD4 cm high cyto… (Interce… <NA>    -0.900   -0.412  0.0815     NaN      2.99
    ##  8 CD4 cm high cyto… typeheal… type    -2.20    -1.48  -0.796      NaN      2.96
    ##  9 CD4 cm ribosome   (Interce… <NA>    -0.111    0.356  0.800      NaN      3.00
    ## 10 CD4 cm ribosome   typeheal… type    -1.65    -1.09  -0.523      NaN      3.13
    ## # ℹ 50 more rows
    ## # ℹ 6 more variables: v_lower <dbl>, v_effect <dbl>, v_upper <dbl>,
    ## #   v_n_eff <dbl>, v_R_k_hat <dbl>, count_data <list>

# Suggested settings

## For single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

## Visualisation of the MCMC chains from the posterior distribution

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
res %>% attr("fit") %>% rstan::traceplot("beta[2,1]")
```

![](inst/figures/unnamed-chunk-19-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = res |> sccomp_test() |> plot()
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-20-1.png)<!-- -->

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

![](inst/figures/unnamed-chunk-21-1.png)<!-- -->

## The old framework

The new tidy framework was introduced in 2024, two, understand the
differences and improvements. Compared to the old framework, please read
this [blog
post](https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/).
