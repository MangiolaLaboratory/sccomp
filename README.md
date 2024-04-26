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
    cores = 1 
  ) |> 
  sccomp_remove_outliers(cores = 1) |> # Optional
  sccomp_test()
```

    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.006252 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 62.52 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 1] earlier than expected.
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100        -4003.375             1.000            1.000
    ## Chain 1:    200        -3725.256             0.537            1.000
    ## Chain 1:    300        -3709.548             0.360            0.075
    ## Chain 1:    400        -3702.257             0.270            0.075
    ## Chain 1:    500        -3703.945             0.216            0.004   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.005503 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 55.03 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Iteration: 250 / 250 [100%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 0.1].
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100        -4409.857             1.000            1.000
    ## Chain 1:    200        -3998.205             0.551            1.000
    ## Chain 1:    300        -3799.150             0.385            0.103
    ## Chain 1:    400        -3688.101             0.296            0.103
    ## Chain 1:    500        -3624.214             0.241            0.052
    ## Chain 1:    600        -3570.193             0.203            0.052
    ## Chain 1:    700        -3520.858             0.176            0.030
    ## Chain 1:    800        -3480.524             0.155            0.030
    ## Chain 1:    900        -3456.364             0.139            0.018
    ## Chain 1:   1000        -3429.880             0.126            0.018
    ## Chain 1:   1100        -3411.811             0.026            0.015
    ## Chain 1:   1200        -3393.713             0.017            0.014
    ## Chain 1:   1300        -3386.192             0.012            0.012
    ## Chain 1:   1400        -3371.879             0.009            0.008   MEAN ELBO CONVERGED   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.005745 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 57.45 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Iteration: 250 / 250 [100%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 0.1].
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100        -4733.528             1.000            1.000
    ## Chain 1:    200        -4296.589             0.551            1.000
    ## Chain 1:    300        -4091.552             0.384            0.102
    ## Chain 1:    400        -3983.170             0.295            0.102
    ## Chain 1:    500        -3922.063             0.239            0.050
    ## Chain 1:    600        -3871.097             0.201            0.050
    ## Chain 1:    700        -3826.858             0.174            0.027
    ## Chain 1:    800        -3791.614             0.154            0.027
    ## Chain 1:    900        -3771.938             0.137            0.016
    ## Chain 1:   1000        -3750.617             0.124            0.016
    ## Chain 1:   1100        -3737.067             0.024            0.013
    ## Chain 1:   1200        -3722.730             0.015            0.012
    ## Chain 1:   1300        -3718.192             0.010            0.009   MEAN ELBO CONVERGED   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.

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

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated, you also appreciate that
this relationship is by model for single cell RNA sequencing data.

`sccomp` models, these relationship to obtain a shrinkage effect on the
estimates of both the abundance and the variability. This shrinkage is
adaptive as it is modelled jointly, thanks for Bayesian inference.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

## Contrasts

``` r
seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    bimodal_mean_variability_association = TRUE,
    cores = 1 
  ) |> 
  sccomp_test( contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"))
```

    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.004847 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 48.47 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 1] earlier than expected.
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100        -3213.050             1.000            1.000
    ## Chain 1:    200        -3161.439             0.508            1.000
    ## Chain 1:    300        -3163.485             0.339            0.016
    ## Chain 1:    400        -3156.074             0.255            0.016
    ## Chain 1:    500        -3158.478             0.204            0.002   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.

    ## # A tibble: 60 × 18
    ##    cell_group  parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_n_eff
    ##    <chr>       <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ##  1 B immature  typecanc… <NA>    -1.87    -1.41   -0.953 0       0            NA
    ##  2 B immature  typeheal… <NA>     0.953    1.41    1.87  0       0            NA
    ##  3 B mem       typecanc… <NA>    -2.31    -1.77   -1.22  0       0            NA
    ##  4 B mem       typeheal… <NA>     1.22     1.77    2.31  0       0            NA
    ##  5 CD4 cm S10… typecanc… <NA>    -1.40    -0.950  -0.493 0       0            NA
    ##  6 CD4 cm S10… typeheal… <NA>     0.493    0.950   1.40  0       0            NA
    ##  7 CD4 cm hig… typecanc… <NA>     0.803    1.80    2.75  0.00100 1.43e-4      NA
    ##  8 CD4 cm hig… typeheal… <NA>    -2.75    -1.80   -0.803 0.00100 1.43e-4      NA
    ##  9 CD4 cm rib… typecanc… <NA>     0.398    0.999   1.62  0.00200 4.44e-4      NA
    ## 10 CD4 cm rib… typeheal… <NA>    -1.62    -0.999  -0.398 0.00200 4.44e-4      NA
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
    cores = 1 
  )
```

    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.005985 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 59.85 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 1] earlier than expected.
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100        -3415.586             1.000            1.000
    ## Chain 1:    200        -3180.228             0.537            1.000
    ## Chain 1:    300        -3171.852             0.359            0.074
    ## Chain 1:    400        -3175.715             0.269            0.074
    ## Chain 1:    500        -3167.118             0.216            0.003   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.

``` r
res
```

    ## # A tibble: 60 × 14
    ##    cell_group        parameter factor c_lower c_effect c_upper c_n_eff c_R_k_hat
    ##    <chr>             <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>     <dbl>
    ##  1 B immature        (Interce… <NA>     0.500    0.824   1.15      NaN      3.17
    ##  2 B immature        typeheal… type     0.949    1.43    1.92      NaN      3.09
    ##  3 B mem             (Interce… <NA>    -1.14    -0.725  -0.321     NaN      3.15
    ##  4 B mem             typeheal… type     1.08     1.65    2.19      NaN      3.21
    ##  5 CD4 cm S100A4     (Interce… <NA>     1.42     1.69    1.95      NaN      3.22
    ##  6 CD4 cm S100A4     typeheal… type     0.692    1.07    1.48      NaN      3.17
    ##  7 CD4 cm high cyto… (Interce… <NA>    -1.54    -0.938  -0.297     NaN      3.10
    ##  8 CD4 cm high cyto… typeheal… type    -1.88    -1.03   -0.128     NaN      3.16
    ##  9 CD4 cm ribosome   (Interce… <NA>    -0.159    0.240   0.656     NaN      3.17
    ## 10 CD4 cm ribosome   typeheal… type    -1.48    -0.893  -0.308     NaN      3.09
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

![](inst/figures/unnamed-chunk-17-1.png)<!-- -->

Plot 1D significance plot

``` r
plots = res |> sccomp_test() |> plot()
```

    ## Joining with `by = join_by(cell_group, sample)`
    ## Joining with `by = join_by(cell_group, type)`

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-18-1.png)<!-- -->

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

![](inst/figures/unnamed-chunk-19-1.png)<!-- -->

## The old framework

The new tidy framework was introduced in 2024, two, understand the
differences and improvements. Compared to the old framework, please read
this [blog
post](https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/).
