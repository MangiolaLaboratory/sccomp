sccomp - Tests differences in cell type proportions and variability from
single-cell data
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/MangiolaLaboratory/sccomp/workflows/R-CMD-check/badge.svg)](https://github.com/MangiolaLaboratory/sccomp/actions/)

<!-- badges: end -->

Cellular omics such as single-cell genomics, proteomics, and
microbiomics allow the characterization of tissue and microbial
community composition, which can be compared between conditions to
identify biological drivers. This strategy has been critical to
unveiling markers of disease progression in conditions such as cancer
and pathogen infections.

For cellular omic data, no method for differential variability analysis
exists, and methods for differential composition analysis only take a
few fundamental data properties into account. Here we introduce
**sccomp**, a generalised method for differential composition and
variability analyses capable of jointly modelling data count
distribution, compositionality, group-specific variability, and
proportion mean-variability association, while being robust to outliers.

**sccomp** is an extensive analysis framework that allows realistic data
simulation and cross-study knowledge transfer. We demonstrate that
mean-variability association is ubiquitous across technologies,
highlighting the inadequacy of the very popular Dirichlet-multinomial
modeling and providing essential principles for differential variability
analysis.

<img src="inst/cartoon_methods.png" width="100%"/>

### Comparison with other methods

- **I**: Data are modelled as counts.
- **II**: Group proportions are modelled as compositional.
- **III**: The proportion variability is modelled as cell-type specific.
- **IV**: Information sharing across cell types, mean–variability
  association.
- **V**: Outlier detection or robustness.
- **VI**: Differential variability analysis.

| Method           | Year | Model                         | I   | II  | III | IV  | V   | VI  |
|------------------|------|-------------------------------|-----|-----|-----|-----|-----|-----|
| **sccomp**       | 2023 | Sum-constrained Beta-binomial | ●   | ●   | ●   | ●   | ●   | ●   |
| **scCODA**       | 2021 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |
| **quasi-binom.** | 2021 | Quasi-binomial                | ●   |     | ●   |     |     |     |
| **rlm**          | 2021 | Robust-log-linear             |     | ●   |     |     | ●   |     |
| **propeller**    | 2021 | Logit-linear + limma          |     | ●   | ●   | ●   |     |     |
| **ANCOM-BC**     | 2020 | Log-linear                    |     | ●   | ●   |     |     |     |
| **corncob**      | 2020 | Beta-binomial                 | ●   |     | ●   |     |     |     |
| **scDC**         | 2019 | Log-linear                    |     | ●   | ●   |     |     |     |
| **dmbvs**        | 2017 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |
| **MixMC**        | 2016 | Zero-inflated Log-linear      |     | ●   | ●   |     |     |     |
| **ALDEx2**       | 2014 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |

### Cite

Mangiola, Stefano, Alexandra J. Roth-Schulze, Marie Trussart, Enrique
Zozaya-Valdés, Mengyao Ma, Zijie Gao, Alan F. Rubin, Terence P. Speed,
Heejung Shim, and Anthony T. Papenfuss. 2023. “Sccomp: Robust
Differential Composition and Variability Analysis for Single-Cell Data.”
Proceedings of the National Academy of Sciences of the United States of
America 120 (33): e2203828120. <https://doi.org/10.1073/pnas.2203828120>
[PNAS - sccomp: Robust differential composition and variability analysis
for single-cell
data](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)

### Talk

<a href="https://www.youtube.com/watch?v=R_lt58We9nA&ab_channel=RConsortium" target="_blank">
<img src="https://img.youtube.com/vi/R_lt58We9nA/mqdefault.jpg" alt="Watch the video" width="280" height="180" border="10" />
</a>

# <img src="inst/logo-01.png" height="139px" width="120px"/>

`sccomp` tests differences in cell type proportions from single-cell
data. It is robust against outliers, it models continuous and discrete
factors, and capable of random-effect/intercept modelling.

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
devtools::install_github("MangiolaLaboratory/sccomp")

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
| `sccomp_predict`                   | Predicts proportions, based on the model, or part of the model                                                              |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors                                                                                |
| `plot`                             | Plots summary plots to asses significance                                                                                   |

# Analysis

``` r
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)
data("seurat_obj")
data("sce_obj")
data("counts_obj")
```

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
  sce_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample =  sample, 
    .cell_group = cell_group, 
    cores = 1 
  ) |> 
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
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
```

Here you see the results of the fit, the effects of the factor on
composition and variability. You also can see the uncertainty around
those effects.

The output is a tibble containing the **Following columns**

- `cell_group` - The cell groups being tested.
- `parameter` - The parameter being estimated from the design matrix
  described by the input `formula_composition` and
  `formula_variability`.
- `factor` - The covariate factor in the formula, if applicable (e.g.,
  not present for Intercept or contrasts).
- `c_lower` - Lower (2.5%) quantile of the posterior distribution for a
  composition (c) parameter.
- `c_effect` - Mean of the posterior distribution for a composition (c)
  parameter.
- `c_upper` - Upper (97.5%) quantile of the posterior distribution for a
  composition (c) parameter.
- `c_pH0` - Probability of the null hypothesis (no difference) for a
  composition (c). This is not a p-value.
- `c_FDR` - False-discovery rate of the null hypothesis for a
  composition (c).
- `v_lower` - Lower (2.5%) quantile of the posterior distribution for a
  variability (v) parameter.
- `v_effect` - Mean of the posterior distribution for a variability (v)
  parameter.
- `v_upper` - Upper (97.5%) quantile of the posterior distribution for a
  variability (v) parameter.
- `v_pH0` - Probability of the null hypothesis for a variability (v).
- `v_FDR` - False-discovery rate of the null hypothesis for a
  variability (v).
- `count_data` - Nested input count data.

``` r
sccomp_result
```

    ## # A tibble: 72 × 20
    ##    cell_group parameter   factor c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>      <chr>       <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B1         (Intercept) <NA>    0.890     1.17   1.46   0       0         1.00
    ##  2 B1         typecancer  type   -1.08     -0.668 -0.256  3.75e-3 6.75e-4   1.00
    ##  3 B2         (Intercept) <NA>    0.420     0.749  1.07   2.50e-4 1.19e-5   1.00
    ##  4 B2         typecancer  type   -1.19     -0.729 -0.254  3.75e-3 9.55e-4   1.00
    ##  5 B3         (Intercept) <NA>   -0.677    -0.352 -0.0364 5.53e-2 4.78e-3   1.00
    ##  6 B3         typecancer  type   -0.757    -0.314  0.162  1.82e-1 4.96e-2   1.00
    ##  7 BM         (Intercept) <NA>   -1.31     -0.989 -0.666  0       0         1.00
    ##  8 BM         typecancer  type   -0.771    -0.308  0.152  1.88e-1 5.56e-2   1.00
    ##  9 CD4 1      (Intercept) <NA>    0.154     0.338  0.532  6.75e-3 1.00e-3   1.00
    ## 10 CD4 1      typecancer  type   -0.0763    0.168  0.414  2.88e-1 7.20e-2   1.00
    ## # ℹ 62 more rows
    ## # ℹ 11 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_rhat <dbl>,
    ## #   v_ess_bulk <dbl>, v_ess_tail <dbl>, count_data <list>

## Outlier identification

`sccomp` can identify outliers probabilistically and exclude them from
the estimation.

``` r
sccomp_result = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = count, 
    cores = 1, verbose = FALSE
  ) |> 
  
  # max_sampling_iterations is used her to not exceed the maximum file size for GitHub action of 100Mb
  sccomp_remove_outliers(cores = 1, verbose = FALSE, max_sampling_iterations = 2000) |> # Optional
  
  sccomp_test()
```

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

## Summary plots

A plot of group proportions, faceted by groups. The blue boxplots
represent the posterior predictive check. If the model is descriptively
adequate for the data, the blue boxplots should roughly overlay the
black boxplots, which represent the observed data. The outliers are
coloured in red. A boxplot will be returned for every (discrete)
covariate present in formula_composition. The colour coding represents
the significant associations for composition and/or variability.

``` r
sccomp_result |> 
  sccomp_boxplot(factor = "type")
```

    ## sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Joining with `by = join_by(cell_group, sample)`

    ## Joining with `by = join_by(cell_group, type)`

![](inst/figures/unnamed-chunk-12-1.png)<!-- -->

You can plot proportions adjusted for unwanted effects. This is helpful
especially for complex models, where multiple factors can significantly
impact the proportions.

``` r
sccomp_result |> 
  sccomp_boxplot(factor = "type", remove_unwanted_effects = TRUE)
```

    ## sccomp says: calculating residuals

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## sccomp says: regressing out unwanted factors
    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Joining with `by = join_by(cell_group, sample)`

    ## Joining with `by = join_by(cell_group, type)`

![](inst/figures/unnamed-chunk-13-1.png)<!-- -->

A plot of estimates of differential composition (c\_) on the x-axis and
differential variability (v\_) on the y-axis. The error bars represent
95% credible intervals. The dashed lines represent the minimal effect
that the hypothesis test is based on. An effect is labelled as
significant if it exceeds the minimal effect according to the 95%
credible interval. Facets represent the covariates in the model.

``` r
sccomp_result |> 
  plot_1D_intervals()
```

    ## sccomp says: Some FDR-significant populations may cross the fold change threshold. 
    ##  This because, as sccomp is a Bayesian method, the FDR is calculated according to Stephens (doi: 10.1093/biostatistics/kxw041), 
    ##  by sorting the probability of the null hypothesis in ascending order and calculating the cumulative average.

![](inst/figures/unnamed-chunk-14-1.png)<!-- -->

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated. sccomp models this
relationship to obtain a shrinkage effect on the estimates of both the
abundance and the variability. This shrinkage is adaptive as it is
modelled jointly, thanks to Bayesian inference.

``` r
sccomp_result |> 
  plot_2D_intervals()
```

![](inst/figures/unnamed-chunk-15-1.png)<!-- -->

You can produce the series of plots calling the `plot` method.

``` r
sccomp_result |> plot() 
```

## Model proportions directly (e.g. from deconvolution)

**Note:** If counts are available, we strongly discourage the use of
proportions, as an important source of uncertainty (i.e., for rare
groups/cell types) is not modeled.

The use of proportions is better suited for modelling deconvolution
results (e.g., of bulk RNA data), in which case counts are not
available.

Proportions should be greater than 0. Assuming that zeros derive from a
precision threshold (e.g., deconvolution), zeros are converted to the
smallest non-zero value.

``` r
sccomp_result = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    .sample = sample,
    .cell_group = cell_group,
    .count = proportion, 
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
```

## Continuous factor

`sccomp` is able to fit erbitrary complex models. In this example we
have a continuous and binary covariate.

``` r
res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type + continuous_covariate, 
      .sample = sample, .cell_group = cell_group,
      cores = 1, verbose=FALSE
    )
```

    ## Loading required package: SeuratObject

    ## Loading required package: sp

    ## 'SeuratObject' was built with package 'Matrix' 1.7.0 but the current
    ## version is 1.7.3; it is recomended that you reinstall 'SeuratObject' as
    ## the ABI for 'Matrix' may have changed

    ## 
    ## Attaching package: 'SeuratObject'

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, t

    ## sccomp says: count column is an integer. The sum-constrained beta binomial model will be used

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy, continuous_covariate

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## Loading model from cache...

    ## sccomp says: to do hypothesis testing run `sccomp_test()`,
    ##   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
    ##   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
    ##   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).

``` r
res
```

    ## # A tibble: 90 × 16
    ##    cell_group       parameter factor c_lower c_effect  c_upper c_rhat c_ess_bulk
    ##    <chr>            <chr>     <chr>    <dbl>    <dbl>    <dbl>  <dbl>      <dbl>
    ##  1 B immature       (Interce… <NA>    0.390    0.753   1.12      1.00      4271.
    ##  2 B immature       typeheal… type    0.869    1.36    1.88      1.00      4118.
    ##  3 B immature       continuo… conti… -0.232    0.0614  0.371     1.00      3977.
    ##  4 B mem            (Interce… <NA>   -1.26    -0.815  -0.370     1.00      3889.
    ##  5 B mem            typeheal… type    1.05     1.68    2.29      1.00      3923.
    ##  6 B mem            continuo… conti… -0.216    0.0890  0.408     1.00      3938.
    ##  7 CD4 cm S100A4    (Interce… <NA>    1.18     1.50    1.81      1.00      3825.
    ##  8 CD4 cm S100A4    typeheal… type    0.690    1.12    1.55      1.00      3664.
    ##  9 CD4 cm S100A4    continuo… conti… -0.0722   0.185   0.443     1.00      4024.
    ## 10 CD4 cm high cyt… (Interce… <NA>   -0.909   -0.458  -0.00205   1.00      3541.
    ## # ℹ 80 more rows
    ## # ℹ 8 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>,
    ## #   count_data <list>

## Random Effect Modeling (mixed-effect modeling, multilevel-modeling, hierarchical modeling)

`sccomp` supports multilevel modeling by allowing the inclusion of
random effects in the compositional and variability formulas. This is
particularly useful when your data has hierarchical or grouped
structures, such as measurements nested within subjects, batches, or
experimental units. By incorporating random effects, sccomp can account
for variability at different levels of your data, improving model fit
and inference accuracy.

### Random Intercept Model

In this example, we demonstrate how to fit a random intercept model
using sccomp. We’ll model the cell-type proportions with both fixed
effects (e.g., treatment) and random effects (e.g., subject-specific
variability).

Here is the input data

``` r
seurat_obj[[]] |> as_tibble()
```

    ## # A tibble: 106,297 × 9
    ##    cell_group nCount_RNA nFeature_RNA group__ group__wrong sample type  group2__
    ##    <chr>           <dbl>        <int> <chr>   <chr>        <chr>  <chr> <chr>   
    ##  1 CD4 naive           0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  2 Mono clas…          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  3 CD4 cm S1…          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  4 B immature          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  5 CD8 naive           0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  6 CD4 naive           0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  7 Mono clas…          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  8 CD4 cm S1…          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ##  9 CD4 cm hi…          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ## 10 B immature          0            0 GROUP1  1            SI-GA… canc… GROUP21 
    ## # ℹ 106,287 more rows
    ## # ℹ 1 more variable: continuous_covariate <dbl>

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type + (1 | group__), 
    .sample = sample,
    .cell_group = cell_group,
    bimodal_mean_variability_association = TRUE,
    cores = 1, verbose = FALSE
  ) 
```

    ## sccomp says: count column is an integer. The sum-constrained beta binomial model will be used

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## Loading model from cache...

    ## sccomp says: to do hypothesis testing run `sccomp_test()`,
    ##   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
    ##   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
    ##   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).

``` r
res
```

    ## # A tibble: 180 × 16
    ##    cell_group parameter       factor c_lower c_effect  c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>           <chr>    <dbl>    <dbl>    <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)     <NA>    0.552     0.853  1.16      1.03      116. 
    ##  2 B immature typehealthy     type    0.551     1.03   1.45      1.00       91.4
    ##  3 B immature (Intercept)___… <NA>   -0.245     0.171  0.851    NA          NA  
    ##  4 B immature (Intercept)___… <NA>   -0.0440    0.315  0.851    NA          NA  
    ##  5 B immature (Intercept)___… <NA>   -0.101     0.243  0.705    NA          NA  
    ##  6 B immature (Intercept)___… <NA>   -0.797    -0.348 -0.00387  NA          NA  
    ##  7 B mem      (Intercept)     <NA>   -0.772    -0.292  0.0913    1.01      135. 
    ##  8 B mem      typehealthy     type    0.264     0.930  1.52      1.02       98.5
    ##  9 B mem      (Intercept)___… <NA>   -0.275     0.127  0.845    NA          NA  
    ## 10 B mem      (Intercept)___… <NA>   -0.0196    0.366  0.948    NA          NA  
    ## # ℹ 170 more rows
    ## # ℹ 8 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>,
    ## #   count_data <list>

### Random Effect Model (random slopes)

`sccomp` can model random slopes. We providean example below.

``` r
res = 
  seurat_obj |>
  sccomp_estimate(
      formula_composition = ~ type + (type | group__),
      .sample = sample,
      .cell_group = cell_group,
      bimodal_mean_variability_association = TRUE,
      cores = 1, verbose = FALSE
    )
```

    ## sccomp says: count column is an integer. The sum-constrained beta binomial model will be used

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## Loading model from cache...

    ## sccomp says: to do hypothesis testing run `sccomp_test()`,
    ##   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
    ##   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
    ##   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).

``` r
res
```

    ## # A tibble: 240 × 16
    ##    cell_group parameter        factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>            <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)      <NA>    0.417    0.831   1.20     1.02       77.7
    ##  2 B immature typehealthy      type    0.484    1.00    1.54     1.00       70.3
    ##  3 B immature (Intercept)___G… <NA>   -0.149    0.116   0.581   NA          NA  
    ##  4 B immature typehealthy___G… <NA>   -0.168    0.0909  0.621   NA          NA  
    ##  5 B immature (Intercept)___G… <NA>   -0.0589   0.192   0.639   NA          NA  
    ##  6 B immature typehealthy___G… <NA>   -0.0795   0.162   0.638   NA          NA  
    ##  7 B immature (Intercept)___G… <NA>   -0.0763   0.199   0.604   NA          NA  
    ##  8 B immature (Intercept)___G… <NA>   -0.727   -0.266   0.0548  NA          NA  
    ##  9 B mem      (Intercept)      <NA>   -0.868   -0.386   0.0415   1.00      117. 
    ## 10 B mem      typehealthy      type    0.250    1.00    1.63     1.01       80.0
    ## # ℹ 230 more rows
    ## # ℹ 8 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>,
    ## #   count_data <list>

### Nested Random Effects

If you have a more complex hierarchy, such as measurements nested within
subjects and subjects nested within batches, you can include multiple
grouping variables. Here `group2__` is nested within `group__`.

``` r
res = 
  seurat_obj |>
  sccomp_estimate(
      formula_composition = ~ type + (type | group__) + (1 | group2__),
      .sample = sample,
      .cell_group = cell_group,
      bimodal_mean_variability_association = TRUE,
      cores = 1, verbose = FALSE
    )
```

    ## sccomp says: count column is an integer. The sum-constrained beta binomial model will be used

    ## sccomp says: estimation

    ## sccomp says: the composition design matrix has columns: (Intercept), typehealthy

    ## sccomp says: the variability design matrix has columns: (Intercept)

    ## Loading model from cache...

    ## sccomp says: to do hypothesis testing run `sccomp_test()`,
    ##   the `test_composition_above_logit_fold_change` = 0.1 equates to a change of ~10%, and
    ##   0.7 equates to ~100% increase, if the baseline is ~0.1 proportion.
    ##   Use `sccomp_proportional_fold_change` to convert c_effect (linear) to proportion difference (non-linear).

``` r
res
```

    ## # A tibble: 300 × 16
    ##    cell_group parameter        factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>            <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)      <NA>    0.441    0.807   1.26     1.00      106. 
    ##  2 B immature typehealthy      type    0.714    1.17    1.62     1.06       63.2
    ##  3 B immature (Intercept)___G… <NA>   -0.182    0.0159  0.386   NA          NA  
    ##  4 B immature typehealthy___G… <NA>   -0.158    0.0200  0.343   NA          NA  
    ##  5 B immature (Intercept)___G… <NA>   -0.233    0.0634  0.392   NA          NA  
    ##  6 B immature typehealthy___G… <NA>   -0.116    0.0595  0.314   NA          NA  
    ##  7 B immature (Intercept)___G… <NA>   -0.0516   0.119   0.588   NA          NA  
    ##  8 B immature (Intercept)___G… <NA>   -0.651   -0.149   0.0318  NA          NA  
    ##  9 B immature (Intercept)___G… <NA>   -0.364   -0.0654  0.0941  NA          NA  
    ## 10 B immature (Intercept)___G… <NA>   -0.0522   0.132   0.508   NA          NA  
    ## # ℹ 290 more rows
    ## # ℹ 8 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>,
    ## #   count_data <list>

## An aid to result interpretation and communication

The estimated effects are expressed in the unconstrained space of the
parameters, similar to differential expression analysis that expresses
changes in terms of log fold change. However, for differences in
proportion, logit fold change must be used, which is harder to interpret
and understand.

Therefore, we provide a more intuitive proportional fold change that can
be more easily understood. However, these cannot be used to infer
significance (use sccomp_test() instead), and a lot of care must be
taken given the nonlinearity of these measures (a 1-fold increase from
0.0001 to 0.0002 carries a different weight than a 1-fold increase from
0.4 to 0.8).

From your estimates, you can specify which effects you are interested in
(this can be a subset of the full model if you wish to exclude unwanted
effects), and the two points you would like to compare.

In the case of a categorical variable, the starting and ending points
are categories.

``` r
sccomp_result |> 
   sccomp_proportional_fold_change(
     formula_composition = ~  type,
     from =  "benign", 
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
    ##  1 B1         2-fold decrease (from 0.0591 to 0.0293)  
    ##  2 B2         2.2-fold decrease (from 0.0388 to 0.0181)
    ##  3 B3         1.4-fold decrease (from 0.013 to 0.0091) 
    ##  4 BM         1.4-fold decrease (from 0.0066 to 0.0048)
    ##  5 CD4 1      1.2-fold increase (from 0.0248 to 0.0299)
    ##  6 CD4 2      1.4-fold increase (from 0.0516 to 0.0738)
    ##  7 CD4 3      2.6-fold decrease (from 0.0849 to 0.0327)
    ##  8 CD4 4      1-fold increase (from 0.0017 to 0.0018)  
    ##  9 CD4 5      1.1-fold increase (from 0.03 to 0.0316)  
    ## 10 CD8 1      1.1-fold increase (from 0.111 to 0.1207) 
    ## # ℹ 26 more rows

## Contrasts

``` r
seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    .sample = sample,
    .cell_group = cell_group, 
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test( contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"))
```

    ## # A tibble: 60 × 12
    ##    cell_group   parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>        <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B immature   typecanc… <NA>    -1.89    -1.34   -0.799 0       0           NA
    ##  2 B immature   typeheal… <NA>     0.799    1.34    1.89  0       0           NA
    ##  3 B mem        typecanc… <NA>    -2.26    -1.65   -1.02  0       0           NA
    ##  4 B mem        typeheal… <NA>     1.02     1.65    2.26  0       0           NA
    ##  5 CD4 cm S100… typecanc… <NA>    -1.46    -0.983  -0.513 5.00e-4 1.25e-4     NA
    ##  6 CD4 cm S100… typeheal… <NA>     0.513    0.983   1.46  5.00e-4 1.25e-4     NA
    ##  7 CD4 cm high… typecanc… <NA>     0.842    1.53    2.26  0       0           NA
    ##  8 CD4 cm high… typeheal… <NA>    -2.26    -1.53   -0.842 0       0           NA
    ##  9 CD4 cm ribo… typecanc… <NA>     0.281    0.943   1.57  6.50e-3 1.80e-3     NA
    ## 10 CD4 cm ribo… typeheal… <NA>    -1.57    -0.943  -0.281 6.50e-3 1.80e-3     NA
    ## # ℹ 50 more rows
    ## # ℹ 3 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>, count_data <list>

## Categorical factor (e.g. Bayesian ANOVA)

This is achieved through model comparison with `loo`. In the following
example, the model with association with factors better fits the data
compared to the baseline model with no factor association. For model
comparisons `sccomp_remove_outliers()` must not be executed as the
leave-one-out must work with the same amount of data, while outlier
elimination does not guarantee it.

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
    inference_method = "hmc",
    enable_loo = TRUE
  )
```

    ## Running MCMC with 6 parallel chains, with 3 thread(s) per chain...
    ## 
    ## Chain 1 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 1 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 2 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 2 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 3 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 3 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 4 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 5 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 6 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 1 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 4 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 1 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 3 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 6 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 2 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 1 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 3 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 6 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 4 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 4 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 5 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 5 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 1 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 2 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 6 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 6 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 5 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 2 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 3 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 6 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 1 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 2 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 5 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 3 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 4 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 5 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 1 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 2 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 3 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 6 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 5 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 2 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 3 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 6 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 1 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 4 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 2 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 3 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 6 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 finished in 2.5 seconds.
    ## Chain 2 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 3 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 4 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 2 finished in 2.5 seconds.
    ## Chain 3 finished in 2.5 seconds.
    ## Chain 4 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 5 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 4 finished in 2.5 seconds.
    ## Chain 5 finished in 2.5 seconds.
    ## Chain 6 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 finished in 2.5 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 2.5 seconds.
    ## Total execution time: 3.0 seconds.

``` r
# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    .sample =  sample, 
    .cell_group = cell_group, 
    inference_method = "hmc",
    enable_loo = TRUE
  )
```

    ## Running MCMC with 6 parallel chains, with 3 thread(s) per chain...
    ## 
    ## Chain 1 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 1 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 2 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 2 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 3 Iteration:   1 / 966 [  0%]  (Warmup) 
    ## Chain 3 Iteration: 100 / 966 [ 10%]  (Warmup)

    ## Chain 4 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 5 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 6 Iteration:   1 / 966 [  0%]  (Warmup)

    ## Chain 1 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 2 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 4 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 1 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 1 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 3 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 5 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 2 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 2 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 4 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 6 Iteration: 100 / 966 [ 10%]  (Warmup) 
    ## Chain 3 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 3 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 5 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 1 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 4 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 4 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 6 Iteration: 200 / 966 [ 20%]  (Warmup) 
    ## Chain 2 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 5 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 5 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 1 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 3 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 4 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 6 Iteration: 300 / 966 [ 31%]  (Warmup) 
    ## Chain 6 Iteration: 301 / 966 [ 31%]  (Sampling) 
    ## Chain 2 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 5 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 3 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 4 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 400 / 966 [ 41%]  (Sampling) 
    ## Chain 2 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 5 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 1 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 3 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 4 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 6 Iteration: 500 / 966 [ 51%]  (Sampling) 
    ## Chain 2 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 5 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 1 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 3 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 4 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 6 Iteration: 600 / 966 [ 62%]  (Sampling) 
    ## Chain 2 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 4 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 5 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 1 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 3 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 6 Iteration: 700 / 966 [ 72%]  (Sampling) 
    ## Chain 1 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 2 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 1 finished in 2.4 seconds.
    ## Chain 3 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 4 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 2 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 4 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 5 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 6 Iteration: 800 / 966 [ 82%]  (Sampling) 
    ## Chain 2 finished in 2.5 seconds.
    ## Chain 4 finished in 2.4 seconds.
    ## Chain 3 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 3 finished in 2.5 seconds.
    ## Chain 5 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 Iteration: 900 / 966 [ 93%]  (Sampling) 
    ## Chain 5 finished in 2.5 seconds.
    ## Chain 6 Iteration: 966 / 966 [100%]  (Sampling) 
    ## Chain 6 finished in 2.5 seconds.
    ## 
    ## All 6 chains finished successfully.
    ## Mean chain execution time: 2.5 seconds.
    ## Total execution time: 3.1 seconds.

``` r
# Compare models
loo_compare(
   attr(model_with_factor_association, "fit")$loo(),
   attr(model_without_association, "fit")$loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -81.2      10.7

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
    cores = 1, verbose = FALSE
  )

res
```

    ## # A tibble: 60 × 16
    ##    cell_group        parameter factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>             <chr>     <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature        (Interce… <NA>     0.341    0.752  1.16     1.00     1372. 
    ##  2 B immature        typeheal… type     0.829    1.36   1.90     1.00     1759. 
    ##  3 B mem             (Interce… <NA>    -1.33    -0.831 -0.293    1.01      215. 
    ##  4 B mem             typeheal… type     1.01     1.68   2.34     1.01      275. 
    ##  5 CD4 cm S100A4     (Interce… <NA>     1.32     1.66   2.01     1.00      674. 
    ##  6 CD4 cm S100A4     typeheal… type     0.406    0.854  1.28     1.00      799. 
    ##  7 CD4 cm high cyto… (Interce… <NA>    -1.02    -0.537 -0.0142   1.00     1572. 
    ##  8 CD4 cm high cyto… typeheal… type    -2.01    -1.15  -0.207    1.00       88.6
    ##  9 CD4 cm ribosome   (Interce… <NA>    -0.155    0.337  0.820    1.00      601. 
    ## 10 CD4 cm ribosome   typeheal… type    -1.72    -1.05  -0.379    1.00     1182. 
    ## # ℹ 50 more rows
    ## # ℹ 8 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>,
    ## #   count_data <list>

**Plot 1D significance plot**

``` r
plots = res |> sccomp_test() |> plot()
```

    ## sccomp says: When visualising proportions, especially for complex models, consider setting `remove_unwanted_effects=TRUE`. This will adjust the proportions, preserving only the observed effect.

    ## Loading model from cache...

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Joining with `by = join_by(cell_group, sample)`

    ## Joining with `by = join_by(cell_group, type)`
    ## sccomp says: Some FDR-significant populations may cross the fold change
    ## threshold.  This because, as sccomp is a Bayesian method, the FDR is calculated
    ## according to Stephens (doi: 10.1093/biostatistics/kxw041), by sorting the
    ## probability of the null hypothesis in ascending order and calculating the
    ## cumulative average.

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-27-1.png)<!-- -->

**Plot 2D significance plot** Data points are cell groups. Error bars
are the 95% credible interval. The dashed lines represent the default
threshold fold change for which the probabilities (c_pH0, v_pH0) are
calculated. pH0 of 0 represent the rejection of the null hypothesis that
no effect is observed.

This plot is provided only if differential variability has been tested.
The differential variability estimates are reliable only if the linear
association between mean and variability for `(intercept)` (left-hand
side facet) is satisfied. A scatterplot (besides the Intercept) is
provided for each category of interest. For each category of interest,
the composition and variability effects should be generally
uncorrelated.

``` r
plots$credible_intervals_2D
```

![](inst/figures/unnamed-chunk-28-1.png)<!-- -->

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

    ## 
    ## A newer version of CmdStan is available. See ?install_cmdstan() to install it.
    ## To disable this check set option or environment variable cmdstanr_no_ver_check=TRUE.

``` r
library(posterior)
```

    ## This is posterior version 1.6.1

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
fit <- res |> attr("fit")

# Extract draws for 'beta[2,1]'
draws <- as_draws_array(fit$draws("beta[2,1]"))

# Create a traceplot for 'beta[2,1]'
mcmc_trace(draws, pars = "beta[2,1]")
```

![](inst/figures/unnamed-chunk-29-1.png)<!-- -->

## The old framework

The new tidy framework was introduced in 2024. To understand the
differences and improvements compared to the old framework, please read
this [blog
post](https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/).
