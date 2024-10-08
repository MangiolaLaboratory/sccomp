sccomp - Outlier-aware and count-based compositional analysis of
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)

<!-- badges: end -->

Cell omics such as single-cell genomics, proteomics, and microbiomics
allow the characterisation of tissue and microbial community
composition, which can be compared between conditions to identify
biological drivers. This strategy has been critical to unveiling markers
of disease progression, such as cancer and pathogen infection.

For cell omic data, no method for differential variability analysis
exists, and methods for differential composition analysis only take a
few fundamental data properties into account. Here we introduce sccomp,
a generalised method for differential composition and variability
analyses capable of jointly modelling data count distribution,
compositionality, group-specific variability, and proportion
mean-variability association, with awareness against outliers.

Sccomp is an extensive analysis framework that allows realistic data
simulation and cross-study knowledge transfer. We demonstrate that
mean-variability association is ubiquitous across technologies,
highlighting the inadequacy of the very popular Dirichlet-multinomial
modelling and providing essential principles for differential
variability analysis.

We show that sccomp accurately fits experimental data, with a 50%
incremental improvement over state-of-the-art algorithms. Using sccomp,
we identified novel differential constraints and composition in the
microenvironment of primary breast cancer.

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

`sccomp` is based on `cmdstanr` which provides the latest version of
`cmdstan` the Bayesian modelling tool. `cmdstanr` is not on CRAN, so we
need to have 3 simple step process (that will be prompted to the user is
forgot).

1.  R installation of `sccomp`
2.  R installation of `cmdstanr`
3.  `cmdstanr` call to `cmdstan` installation

**Bioconductor**

``` r
if (!requireNamespace("BiocManager")) install.packages("BiocManager")

# Step 1
BiocManager::install("sccomp")

# Step 2
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# Step 3
cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
cmdstanr::install_cmdstan()
```

**Github**

``` r
# Step 1
devtools::install_github("stemangiola/sccomp")

# Step 2
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# Step 3
cmdstanr::check_cmdstan_toolchain(fix = TRUE) # Just checking system setting
cmdstanr::install_cmdstan()
```

| Function                           | Description                                                                                                                 |
|------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                  | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`           | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                      | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`                 | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                   | Predicts proportions, based on the model, or part of the model                                                               |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plots summary plots to asses significance                                                                                   |

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

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.
    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

Here you see the results of the fit, the effects of the factor on
composition and variability. You also can see the uncertainty around
those effects.

``` r
sccomp_result
```

    ## # A tibble: 72 × 14
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR v_lower
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>     0.886    1.05   1.23   0       0         -6.13
    ##  2 B1         typecancer type    -1.16    -0.884 -0.623  0       0         NA   
    ##  3 B2         (Intercep… <NA>     0.422    0.702  0.969  0       0         -5.78
    ##  4 B2         typecancer type    -1.18    -0.810 -0.429  2.50e-4 3.12e-5   NA   
    ##  5 B3         (Intercep… <NA>    -0.638   -0.377 -0.130  1.55e-2 1.06e-3   -6.78
    ##  6 B3         typecancer type    -0.606   -0.248  0.134  2.17e-1 4.81e-2   NA   
    ##  7 BM         (Intercep… <NA>    -1.27    -1.01  -0.744  0       0         -7.36
    ##  8 BM         typecancer type    -0.715   -0.350  0.0280 9.88e-2 1.87e-2   NA   
    ##  9 CD4 1      (Intercep… <NA>     0.146    0.318  0.503  8.25e-3 4.35e-4   -6.59
    ## 10 CD4 1      typecancer type    -0.113    0.132  0.376  3.99e-1 1.26e-1   NA   
    ## # ℹ 62 more rows
    ## # ℹ 5 more variables: v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>,
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

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## # A tibble: 36 × 2
    ##    cell_group statement                                
    ##    <chr>      <glue>                                   
    ##  1 B1         2.4-fold decrease (from 0.0528 to 0.0225)
    ##  2 B2         2.2-fold decrease (from 0.0373 to 0.0171)
    ##  3 B3         1.2-fold decrease (from 0.0126 to 0.0103)
    ##  4 BM         1.4-fold decrease (from 0.0068 to 0.005) 
    ##  5 CD4 1      1.2-fold increase (from 0.0253 to 0.0298)
    ##  6 CD4 2      1.7-fold increase (from 0.0476 to 0.0827)
    ##  7 CD4 3      3.3-fold decrease (from 0.1019 to 0.0307)
    ##  8 CD4 4      1.2-fold increase (from 0.0016 to 0.0019)
    ##  9 CD4 5      1-fold increase (from 0.0302 to 0.0312)  
    ## 10 CD8 1      1.2-fold increase (from 0.1027 to 0.1269)
    ## # ℹ 26 more rows

## Summary plots

A plot of group proportion, faceted by groups. The blue boxplots
represent the posterior predictive check. If the model is likely to be
descriptively adequate to the data, the blue box plot should roughly
overlay with the black box plot, which represents the observed data. The
outliers are coloured in red. A box plot will be returned for every
(discrete) covariate present in `formula_composition`. The colour coding
represents the significant associations for composition and/or
variability.

``` r
sccomp_result |> 
  sccomp_boxplot(factor = "type")
```

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Joining with `by = join_by(cell_group, sample)`

    ## Joining with `by = join_by(cell_group, type)`

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if bigger than the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
sccomp_result |> 
  plot_1D_intervals()
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated, you also appreciate that
this relationship is by model for single cell RNA sequencing data.

`sccomp` models, these relationship to obtain a shrinkage effect on the
estimates of both the abundance and the variability. This shrinkage is
adaptive as it is modelled jointly, thanks for Bayesian inference.

``` r
sccomp_result |> 
  plot_2D_intervals()
```

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

You can produce the series of plots calling the `plot` method.

``` r
sccomp_result |> plot() 
```

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

    ## # A tibble: 60 × 14
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR v_lower
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -1.89    -1.34   -0.780 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.780    1.34    1.89  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.25    -1.63   -0.988 0       0            NA
    ##  4 B mem       typeheal… <NA>     0.988    1.63    2.25  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.46    -0.986  -0.527 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.527    0.986   1.46  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.808    1.55    2.25  0       0            NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.25    -1.55   -0.808 0       0            NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.375    0.980   1.60  0.00250 5.31e-4      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.60    -0.980  -0.375 0.00250 5.31e-4      NA
    ## # ℹ 50 more rows
    ## # ℹ 5 more variables: v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>,
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
    inference_method = "hmc",
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
    inference_method = "hmc",
    enable_loo = TRUE
  )

# Compare models
loo_compare(
   attr(model_with_factor_association, "fit")$loo(),
   attr(model_without_association, "fit")$loo()
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

    ## # A tibble: 60 × 10
    ##    cell_group parameter factor c_lower c_effect c_upper v_lower v_effect v_upper
    ##    <chr>      <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
    ##  1 B immature (Interce… <NA>     0.350    0.759  1.17   -4.37     -3.94  -3.48  
    ##  2 B immature typeheal… type     0.768    1.34   1.87   -1.03     -0.289  0.330 
    ##  3 B mem      (Interce… <NA>    -1.35    -0.832 -0.351  -5.13     -4.58  -4.07  
    ##  4 B mem      typeheal… type     0.992    1.67   2.37   -1.69     -0.789 -0.0238
    ##  5 CD4 cm S1… (Interce… <NA>     1.35     1.68   2.01   -3.80     -3.36  -2.91  
    ##  6 CD4 cm S1… typeheal… type     0.364    0.811  1.25   -1.34     -0.752 -0.243 
    ##  7 CD4 cm hi… (Interce… <NA>    -1.05    -0.511  0.0150 -5.20     -4.67  -4.17  
    ##  8 CD4 cm hi… typeheal… type    -1.94    -0.936  0.0524  0.531     1.55   2.57  
    ##  9 CD4 cm ri… (Interce… <NA>    -0.158    0.301  0.773  -5.10     -4.57  -4.04  
    ## 10 CD4 cm ri… typeheal… type    -1.73    -1.03  -0.362  -0.0762    0.530  1.30  
    ## # ℹ 50 more rows
    ## # ℹ 1 more variable: count_data <list>

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
library(cmdstanr)
```

    ## This is cmdstanr version 0.8.1.9000

    ## - CmdStanR documentation and vignettes: mc-stan.org/cmdstanr

    ## - CmdStan path: /Users/a1234450/.cmdstan/cmdstan-2.35.0

    ## - CmdStan version: 2.35.0

``` r
library(posterior)
```

    ## This is posterior version 1.6.0

    ## 
    ## Attaching package: 'posterior'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     mad, sd, var

    ## The following objects are masked from 'package:base':
    ## 
    ##     %in%, match

``` r
library(bayesplot)
```

    ## This is bayesplot version 1.11.1

    ## - Online documentation and vignettes at mc-stan.org/bayesplot

    ## - bayesplot theme set to bayesplot::theme_default()

    ##    * Does _not_ affect other ggplot2 plots

    ##    * See ?bayesplot_theme_set for details on theme setting

    ## 
    ## Attaching package: 'bayesplot'

    ## The following object is masked from 'package:posterior':
    ## 
    ##     rhat

``` r
# Assuming res contains the fit object from cmdstanr
fit <- res %>% attr("fit")

# Extract draws for 'beta[2,1]'
draws <- as_draws_array(fit$draws("beta[2,1]"))

# Create a traceplot for 'beta[2,1]'
mcmc_trace(draws, pars = "beta[2,1]")
```

![](inst/figures/unnamed-chunk-19-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = res |> sccomp_test() |> plot()
```

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

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
provided for each category of interest. For each category of
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
