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

Sccomp is a generalised method for differential composition and
variability analyses.

## Characteristics

- Complex linear models with continuous and categorical covariates
- Multilevel modelling, with population (i.e. fixed) and group (random)
  effects
- Modelling counts
- Modelling proportionality
- Modelling cell-type specific variability
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

### From proportions

If counts are available we strongly discourage the use of proportions,
as an important source of uncertainty (i.e. for rare groups/cell-types)
is not modeled.

Proportions should be bigger than 0. Assuming that 0s derive from a
precision threshold (e.g. deconvolution), 0s are converted to the
smaller non 0 value.

``` r
res = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = proportion, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> 
    sccomp_test(test_composition_above_logit_fold_change = 0.2)
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000283 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.83 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 1.584 seconds (Warm-up)
    ## Chain 1:                13.742 seconds (Sampling)
    ## Chain 1:                15.326 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000155 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.55 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:     1 / 20299 [  0%]  (Warmup)
    ## Chain 1: Iteration:   301 / 20299 [  1%]  (Sampling)
    ## Chain 1: Iteration:  1300 / 20299 [  6%]  (Sampling)
    ## Chain 1: Iteration:  2300 / 20299 [ 11%]  (Sampling)
    ## Chain 1: Iteration:  3300 / 20299 [ 16%]  (Sampling)
    ## Chain 1: Iteration:  4300 / 20299 [ 21%]  (Sampling)
    ## Chain 1: Iteration:  5300 / 20299 [ 26%]  (Sampling)
    ## Chain 1: Iteration:  6300 / 20299 [ 31%]  (Sampling)
    ## Chain 1: Iteration:  7300 / 20299 [ 35%]  (Sampling)
    ## Chain 1: Iteration:  8300 / 20299 [ 40%]  (Sampling)
    ## Chain 1: Iteration:  9300 / 20299 [ 45%]  (Sampling)
    ## Chain 1: Iteration: 10300 / 20299 [ 50%]  (Sampling)
    ## Chain 1: Iteration: 11300 / 20299 [ 55%]  (Sampling)
    ## Chain 1: Iteration: 12300 / 20299 [ 60%]  (Sampling)
    ## Chain 1: Iteration: 13300 / 20299 [ 65%]  (Sampling)
    ## Chain 1: Iteration: 14300 / 20299 [ 70%]  (Sampling)
    ## Chain 1: Iteration: 15300 / 20299 [ 75%]  (Sampling)
    ## Chain 1: Iteration: 16300 / 20299 [ 80%]  (Sampling)
    ## Chain 1: Iteration: 17300 / 20299 [ 85%]  (Sampling)
    ## Chain 1: Iteration: 18300 / 20299 [ 90%]  (Sampling)
    ## Chain 1: Iteration: 19300 / 20299 [ 95%]  (Sampling)
    ## Chain 1: Iteration: 20299 / 20299 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 3.306 seconds (Warm-up)
    ## Chain 1:                164.122 seconds (Sampling)
    ## Chain 1:                167.428 seconds (Total)
    ## Chain 1: 
    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000223 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.23 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 1.446 seconds (Warm-up)
    ## Chain 1:                13.614 seconds (Sampling)
    ## Chain 1:                15.06 seconds (Total)
    ## Chain 1:

``` r
res
```

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>     0.967    1.21   1.45   0       0         7272.
    ##  2 B1         typecancer type    -1.14    -0.725 -0.333  5.00e-3 1.34e-3   7992.
    ##  3 B2         (Intercep… <NA>     0.474    0.780  1.08   5.00e-4 4.76e-5   5502.
    ##  4 B2         typecancer type    -1.29    -0.784 -0.320  8.  e-3 2.08e-3   9060.
    ##  5 B3         (Intercep… <NA>    -0.600   -0.313 -0.0284 2.13e-1 1.93e-2   6199.
    ##  6 B3         typecancer type    -0.831   -0.398  0.0233 1.74e-1 4.02e-2   7531.
    ##  7 BM         (Intercep… <NA>    -1.28    -0.988 -0.703  0       0         7054.
    ##  8 BM         typecancer type    -0.781   -0.354  0.0370 2.18e-1 6.52e-2   6956.
    ##  9 CD4 1      (Intercep… <NA>     0.158    0.362  0.577  5.50e-2 3.42e-3   6656.
    ## 10 CD4 1      typecancer type    -0.150    0.139  0.412  6.60e-1 1.79e-1   7488.
    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Visualisation

## Summary plots

``` r
plots = plot(res) 
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

### Visualise the descriptive adequacy of the model

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

![](inst/figures/unnamed-chunk-10-1.png)<!-- -->

### Mean/variability association

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

![](inst/figures/unnamed-chunk-11-1.png)<!-- -->

### Effects significance

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

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

## Contrasts

``` r
counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> 
    sccomp_test(
      contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"),
      test_composition_above_logit_fold_change = 0.2
    )
```

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
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1, 
    enable_loo = TRUE
  )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.00022 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.2 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 1.651 seconds (Warm-up)
    ## Chain 1:                13.689 seconds (Sampling)
    ## Chain 1:                15.34 seconds (Total)
    ## Chain 1:

``` r
# Fit second model
model_without_association = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 , 
    enable_loo = TRUE
  )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000214 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.14 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 1.399 seconds (Warm-up)
    ## Chain 1:                13.272 seconds (Sampling)
    ## Chain 1:                14.671 seconds (Total)
    ## Chain 1:

``` r
# Compare models
loo_compare(
  model_with_factor_association |> attr("fit") |> loo(),
  model_without_association |> attr("fit") |> loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -54.7      11.9

## Differential variability, binary factor

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    .count = proportion,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )
```

    ## 
    ## SAMPLING FOR MODEL 'glm_multi_beta_binomial' NOW (CHAIN 1).
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.000245 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 2.45 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Iteration:    1 / 4300 [  0%]  (Warmup)
    ## Chain 1: Iteration:  301 / 4300 [  7%]  (Sampling)
    ## Chain 1: Iteration: 1300 / 4300 [ 30%]  (Sampling)
    ## Chain 1: Iteration: 2300 / 4300 [ 53%]  (Sampling)
    ## Chain 1: Iteration: 3300 / 4300 [ 76%]  (Sampling)
    ## Chain 1: Iteration: 4300 / 4300 [100%]  (Sampling)
    ## Chain 1: 
    ## Chain 1:  Elapsed Time: 3.673 seconds (Warm-up)
    ## Chain 1:                19.874 seconds (Sampling)
    ## Chain 1:                23.547 seconds (Total)
    ## Chain 1:

    ## # A tibble: 72 × 18
    ##    cell_group parameter  factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>      <chr>      <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B1         (Intercep… <NA>    0.948     1.18   1.42   0       0         4457.
    ##  2 B1         typecancer type   -1.03     -0.588 -0.123  4.35e-2 1.15e-2   3975.
    ##  3 B2         (Intercep… <NA>    0.454     0.784  1.12   7.51e-4 7.51e-5   5010.
    ##  4 B2         typecancer type   -1.25     -0.714 -0.154  3.68e-2 5.38e-3   4371.
    ##  5 B3         (Intercep… <NA>   -0.585    -0.343 -0.0762 1.28e-1 1.73e-2   4599.
    ##  6 B3         typecancer type   -0.707    -0.231  0.251  4.41e-1 1.21e-1   3998.
    ##  7 BM         (Intercep… <NA>   -1.25     -0.967 -0.650  0       0         4372.
    ##  8 BM         typecancer type   -0.741    -0.343  0.0617 2.33e-1 6.66e-2   4847.
    ##  9 CD4 1      (Intercep… <NA>    0.164     0.357  0.562  5.86e-2 4.09e-3   3887.
    ## 10 CD4 1      typecancer type   -0.0558    0.223  0.502  4.29e-1 1.07e-1   4770.
    ## # ℹ 62 more rows
    ## # ℹ 9 more variables: c_R_k_hat <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_n_eff <dbl>, v_R_k_hat <dbl>,
    ## #   count_data <list>

# Suggested settings

## For single-cell RNA sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and microbiome data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

# Multilevel modelling

`sccomp` is cabable of estimating population (i.e. fixed) and group
(i.e. random) effects. The formulation is analogous to the `lme4`
package and `brms`.

!! For now, only one grouping is allowed (e.g. group2\_\_).

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type + continuous_covariate + (type | group2__), 
    formula_variability = ~ type,
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_remove_outliers() |> 
    sccomp_test(
      test_composition_above_logit_fold_change = 0.2
    )

res
```

# Removal of unwanted variation

After you model your dataset, you can remove the unwanted variation from
your input data, **for visualisation purposes**

We decide to just keep the type population (i.e. fixed) effect for
abundance, and do not keep it for variability.

``` r
res |> sccomp_remove_unwanted_variation(~type)
```
