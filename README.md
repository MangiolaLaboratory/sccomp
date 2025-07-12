sccomp: Differential Composition and Variability Analysis for
Single-Cell Data
================
Stefano Mangiola
2025-07-13

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/tidyseurat/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/tidyseurat/actions/)
<!-- badges: end -->

# sccomp: Advanced Differential Composition and Variability Analysis for Single-Cell Data

**sccomp** is a powerful R package designed for comprehensive
differential composition and variability analysis in single-cell
genomics, proteomics, and microbiomics data. This cutting-edge tool
provides robust Bayesian modeling capabilities for analyzing cell type
proportions and variability across different experimental conditions.

## What is sccomp?

**sccomp** (Single-Cell Composition Analysis) is a specialized R package
that addresses the critical need for advanced statistical methods in
single-cell data analysis. It enables researchers to:

- **Analyze cell type proportions** with high statistical rigor
- **Detect differential variability** between experimental conditions  
- **Identify outliers** using probabilistic methods
- **Model complex experimental designs** with random effects
- **Handle multiple data types** including scRNA-seq, CyTOF, and
  microbiome data

## Key Applications

- Single-Cell RNA Sequencing (scRNA-seq)
- CyTOF (Cytometry by Time-of-Flight)
- Microbiome Analysis

## Why sccomp?

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

<img src="../inst/cartoon_methods.png" width="100%"/>

### Comprehensive Method Comparison

- **I**: Data are modelled as counts.
- **II**: Group proportions are modelled as compositional.
- **III**: The proportion variability is modelled as cell-type specific.
- **IV**: Information sharing across cell types, mean–variability
  association.
- **V**: Outlier detection or robustness.
- **VI**: Differential variability analysis.
- **VII** Mixed effect modelling
- **VIII** Removal unwanted effects

| Method | Year | Model | I | II | III | IV | V | VI | VII | VIII |
|----|----|----|----|----|----|----|----|----|----|----|
| **sccomp** | 2023 | Sum-constrained Beta-binomial | ● | ● | ● | ● | ● | ● | ● | ● |
| **scCODA** | 2021 | Dirichlet-multinomial | ● | ● |  |  |  |  |  |  |
| **quasi-binom.** | 2021 | Quasi-binomial | ● |  | ● |  |  |  |  |  |
| **rlm** | 2021 | Robust-log-linear |  | ● |  |  | ● |  |  |  |
| **propeller** | 2021 | Logit-linear + limma |  | ● | ● | ● |  |  |  |  |
| **ANCOM-BC** | 2020 | Log-linear |  | ● | ● |  |  |  |  |  |
| **corncob** | 2020 | Beta-binomial | ● |  | ● |  |  |  |  |  |
| **scDC** | 2019 | Log-linear |  | ● | ● |  |  |  |  |  |
| **dmbvs** | 2017 | Dirichlet-multinomial | ● | ● |  |  |  |  |  |  |
| **MixMC** | 2016 | Zero-inflated Log-linear |  | ● | ● |  |  |  |  |  |
| **ALDEx2** | 2014 | Dirichlet-multinomial | ● | ● |  |  |  |  |  |  |

### Scientific Citation

Mangiola, Stefano, Alexandra J. Roth-Schulze, Marie Trussart, Enrique
Zozaya-Valdés, Mengyao Ma, Zijie Gao, Alan F. Rubin, Terence P. Speed,
Heejung Shim, and Anthony T. Papenfuss. 2023. “Sccomp: Robust
Differential Composition and Variability Analysis for Single-Cell Data.”
Proceedings of the National Academy of Sciences of the United States of
America 120 (33): e2203828120. <https://doi.org/10.1073/pnas.2203828120>
[PNAS - sccomp: Robust differential composition and variability analysis
for single-cell
data](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)

# Installation Guide

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

## Core Functions

| Function | Description |
|----|----|
| `sccomp_estimate` | Fit the model onto the data, and estimate the coefficients |
| `sccomp_remove_outliers` | Identify outliers probabilistically based on the model fit, and exclude them from the estimation |
| `sccomp_test` | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate` | Simulate data from the model, or part of the model |
| `sccomp_predict` | Predicts proportions, based on the model, or part of the model |
| `sccomp_remove_unwanted_effects` | Removes the variability for unwanted factors |
| `plot` | Plots summary plots to assess significance |

# Analysis Tutorial

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

## Binary Factor Analysis

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

### From Seurat, SingleCellExperiment, metadata objects

``` r
sccomp_result = 
  sce_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    sample = "sample", 
    cell_group = "cell_group", 
    cores = 1,
    verbose = FALSE
  ) |> 
  sccomp_test()
```

### From counts

``` r
sccomp_result = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count", 
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

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 36 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 72 × 19
    ##    cell_group parameter factor  c_lower c_effect  c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>      <chr>     <chr>     <dbl>    <dbl>    <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B1         (Interce… <NA>    0.956      1.20   1.45    0       0        1.000
    ##  2 B1         typecanc… type   -0.925     -0.617 -0.297   7.50e-4 9.09e-5  1.00 
    ##  3 B2         (Interce… <NA>    0.509      0.771  1.03    0       0        1.00 
    ##  4 B2         typecanc… type   -0.965     -0.668 -0.359   2.50e-4 2.50e-5  1.00 
    ##  5 B3         (Interce… <NA>   -0.576     -0.330 -0.0747  3.83e-2 3.21e-3  1.00 
    ##  6 B3         typecanc… type   -0.575     -0.282  0.00542 1.16e-1 2.36e-2  1.000
    ##  7 BM         (Interce… <NA>   -1.23      -0.969 -0.712   0       0        1.00 
    ##  8 BM         typecanc… type   -0.608     -0.288  0.0201  1.14e-1 1.53e-2  1.000
    ##  9 CD4 1      (Interce… <NA>    0.202      0.365  0.529   1.00e-3 7.69e-5  1.00 
    ## 10 CD4 1      typecanc… type   -0.00989    0.208  0.420   1.58e-1 2.90e-2  1.00 
    ## # ℹ 62 more rows
    ## # ℹ 10 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_rhat <dbl>,
    ## #   v_ess_bulk <dbl>, v_ess_tail <dbl>

## Outlier Identification

`sccomp` can identify outliers probabilistically and exclude them from
the estimation.

``` r
sccomp_result = 
  counts_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count", 
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

## Visualization and Summary Plots

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

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Warning in stat_summary(aes(!!as.symbol(factor_of_interest), (generated_proportions)), : Ignoring unknown parameters: `outlier.shape`, `outlier.colour`, and
    ## `outlier.size`

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

You can plot proportions adjusted for unwanted effects. This is helpful
especially for complex models, where multiple factors can significantly
impact the proportions.

``` r
sccomp_result |> 
  sccomp_boxplot(factor = "type", remove_unwanted_effects = TRUE)
```

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Warning in stat_summary(aes(!!as.symbol(factor_of_interest), (generated_proportions)), : Ignoring unknown parameters: `outlier.shape`, `outlier.colour`, and
    ## `outlier.size`

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

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

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated. sccomp models this
relationship to obtain a shrinkage effect on the estimates of both the
abundance and the variability. This shrinkage is adaptive as it is
modelled jointly, thanks to Bayesian inference.

``` r
sccomp_result |> 
  plot_2D_intervals()
```

    ## Warning: `geom_hline()`: Ignoring `mapping` because `yintercept` was provided.

    ## Warning: `geom_hline()`: Ignoring `data` because `yintercept` was provided.

    ## Warning in geom_hline(yintercept = mean_corrected, data = data.frame(parameter
    ## = "(Intercept, corrected)"), : Ignoring unknown parameters: `inherit.aes`

    ## Warning: 1 unknown level in `f`: (Intercept, corrected)
    ## 1 unknown level in `f`: (Intercept, corrected)

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

You can produce the series of plots calling the `plot` method.

``` r
sccomp_result |> plot() 
```

## Model Proportions Directly (e.g. from deconvolution)

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
    sample = "sample",
    cell_group = "cell_group",
    abundance = "proportion", 
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
```

## Continuous Factor Analysis

`sccomp` is able to fit arbitrary complex models. In this example we
have a continuous and binary covariate.

``` r
res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type + continuous_covariate, 
      sample = "sample", 
      cell_group = "cell_group",
      cores = 1, verbose=FALSE
    )

res
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type + continuous_covariate 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 90 × 15
    ##    cell_group        parameter factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>             <chr>     <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature        (Interce… <NA>     0.556   0.835    1.12   1.00       4067.
    ##  2 B immature        typeheal… type     1.03    1.35     1.66   1.000      3963.
    ##  3 B immature        continuo… conti…  -0.274   0.0491   0.378  1.00       3993.
    ##  4 B mem             (Interce… <NA>    -0.981  -0.671   -0.360  1.00       3959.
    ##  5 B mem             typeheal… type     1.20    1.58     1.94   1.00       3951.
    ##  6 B mem             continuo… conti…  -0.267   0.0689   0.403  1.00       4079.
    ##  7 CD4 cm S100A4     (Interce… <NA>     1.31    1.56     1.79   1.00       3935.
    ##  8 CD4 cm S100A4     typeheal… type     0.851   1.13     1.41   1.000      3904.
    ##  9 CD4 cm S100A4     continuo… conti…  -0.116   0.173    0.458  1.00       4006.
    ## 10 CD4 cm high cyto… (Interce… <NA>    -0.920  -0.575   -0.236  1.000      3531.
    ## # ℹ 80 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

## Random Effect Modeling (Mixed-Effect Modeling)

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
    sample = "sample",
    cell_group = "cell_group",
    bimodal_mean_variability_association = TRUE,
    cores = 1, verbose = FALSE
  ) 

res
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type + (1 | group__) 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 180 × 15
    ##    cell_group parameter        factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>            <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)      <NA>    0.505    0.843   1.24    1.000      150. 
    ##  2 B immature typehealthy      type    0.736    1.20    1.62    1.00       113. 
    ##  3 B immature (Intercept)___G… <NA>   -0.329    0.0385  0.490   1.01        95.2
    ##  4 B immature (Intercept)___G… <NA>   -0.0716   0.232   0.683   1.01        94.6
    ##  5 B immature (Intercept)___G… <NA>   -0.181    0.202   0.634   1.01       128. 
    ##  6 B immature (Intercept)___G… <NA>   -0.734   -0.273   0.0441  1.01       146. 
    ##  7 B mem      (Intercept)      <NA>   -0.859   -0.436  -0.0417  1.00       111. 
    ##  8 B mem      typehealthy      type    0.591    1.21    1.77    1.01        61.5
    ##  9 B mem      (Intercept)___G… <NA>   -0.325    0.0132  0.524   1.01       137. 
    ## 10 B mem      (Intercept)___G… <NA>   -0.0531   0.305   0.812   1.01        81.0
    ## # ℹ 170 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

### Random Effect Model (random slopes)

`sccomp` can model random slopes. We provide an example below.

``` r
res = 
  seurat_obj |>
  sccomp_estimate(
      formula_composition = ~ type + (type | group__),
      sample = "sample",
      cell_group = "cell_group",
      bimodal_mean_variability_association = TRUE,
      cores = 1, verbose = FALSE
    )

res
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type + (type | group__) 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 240 × 15
    ##    cell_group parameter        factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>            <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)      <NA>    0.438    0.834   1.17     1.00      152. 
    ##  2 B immature typehealthy      type    0.771    1.19    1.75     1.05      100. 
    ##  3 B immature (Intercept)___G… <NA>   -0.388    0.0250  0.368    1.01      150. 
    ##  4 B immature typehealthy___G… <NA>   -0.282    0.0116  0.341    1.00      193. 
    ##  5 B immature (Intercept)___G… <NA>   -0.333    0.119   0.506    1.02      109. 
    ##  6 B immature typehealthy___G… <NA>   -0.134    0.116   0.538    1.01      210. 
    ##  7 B immature (Intercept)___G… <NA>   -0.0860   0.214   0.629    1.01      135. 
    ##  8 B immature (Intercept)___G… <NA>   -0.709   -0.247   0.0294   1.02       92.6
    ##  9 B mem      (Intercept)      <NA>   -0.907   -0.538  -0.103    1.01      172. 
    ## 10 B mem      typehealthy      type    0.751    1.29    1.78     1.01       99.1
    ## # ℹ 230 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

### Nested Random Effects

If you have a more complex hierarchy, such as measurements nested within
subjects and subjects nested within batches, you can include multiple
grouping variables. Here `group2__` is nested within `group__`.

``` r
res = 
  seurat_obj |>
  sccomp_estimate(
      formula_composition = ~ type + (type | group__) + (1 | group2__),
      sample = "sample",
      cell_group = "cell_group",
      bimodal_mean_variability_association = TRUE,
      cores = 1, verbose = FALSE
    )

res
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type + (type | group__) + (1 | group2__) 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 300 × 15
    ##    cell_group parameter        factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>      <chr>            <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature (Intercept)      <NA>    0.450    0.790   1.50     1.05       56.5
    ##  2 B immature typehealthy      type    0.798    1.24    1.74     1.00       87.9
    ##  3 B immature (Intercept)___G… <NA>   -0.282    0.0832  0.496    1.01       84.8
    ##  4 B immature typehealthy___G… <NA>   -0.174    0.0726  0.436    1.04       63.4
    ##  5 B immature (Intercept)___G… <NA>   -0.309    0.0459  0.327    1.01      142. 
    ##  6 B immature typehealthy___G… <NA>   -0.213    0.0382  0.310    1.02      152. 
    ##  7 B immature (Intercept)___G… <NA>   -0.0480   0.214   0.631    1.01       63.4
    ##  8 B immature (Intercept)___G… <NA>   -0.723   -0.270   0.0482   1.01       84.9
    ##  9 B immature (Intercept)___G… <NA>   -0.562   -0.122   0.163    1.02      124. 
    ## 10 B immature (Intercept)___G… <NA>   -0.160    0.245   0.609    1.02       46.3
    ## # ℹ 290 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

## Result Interpretation and Communication

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
res |> 
   sccomp_proportional_fold_change(
     formula_composition = ~  type,
     from =  "healthy", 
     to = "cancer"
    ) |> 
  select(cell_group, statement)
```

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## # A tibble: 30 × 2
    ##    cell_group           statement                                
    ##    <chr>                <glue>                                   
    ##  1 B immature           2.1-fold decrease (from 0.1083 to 0.0518)
    ##  2 B mem                2.4-fold decrease (from 0.0335 to 0.0141)
    ##  3 CD4 cm high cytokine 7.9-fold increase (from 0.0015 to 0.0121)
    ##  4 CD4 cm ribosome      3.5-fold increase (from 0.0076 to 0.0261)
    ##  5 CD4 cm S100A4        1.6-fold decrease (from 0.1509 to 0.0969)
    ##  6 CD4 em high cytokine 4.7-fold increase (from 0.0022 to 0.0105)
    ##  7 CD4 naive            1.4-fold decrease (from 0.1099 to 0.0771)
    ##  8 CD4 ribosome         3.1-fold decrease (from 0.0877 to 0.0285)
    ##  9 CD8 em 1             1.2-fold increase (from 0.0495 to 0.0591)
    ## 10 CD8 em 2             3.8-fold increase (from 0.0056 to 0.0215)
    ## # ℹ 20 more rows

## Contrasts Analysis

``` r
seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 0 + type, 
    sample = "sample",
    cell_group = "cell_group", 
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test( contrasts =  c("typecancer - typehealthy", "typehealthy - typecancer"))
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~0 + type 
    ##   Variability formula: ~1 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 60 × 11
    ##    cell_group   parameter factor c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>        <chr>     <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B immature   typecanc… <NA>    -1.89    -1.35   -0.806 0       0           NA
    ##  2 B immature   typeheal… <NA>     0.806    1.35    1.89  0       0           NA
    ##  3 B mem        typecanc… <NA>    -2.22    -1.65   -1.12  0       0           NA
    ##  4 B mem        typeheal… <NA>     1.12     1.65    2.22  0       0           NA
    ##  5 CD4 cm S100… typecanc… <NA>    -1.47    -0.991  -0.514 0       0           NA
    ##  6 CD4 cm S100… typeheal… <NA>     0.514    0.991   1.47  0       0           NA
    ##  7 CD4 cm high… typecanc… <NA>     0.957    1.57    2.18  0       0           NA
    ##  8 CD4 cm high… typeheal… <NA>    -2.18    -1.57   -0.957 0       0           NA
    ##  9 CD4 cm ribo… typecanc… <NA>     0.359    0.945   1.53  0.00275 7.95e-4     NA
    ## 10 CD4 cm ribo… typeheal… <NA>    -1.53    -0.945  -0.359 0.00275 7.95e-4     NA
    ## # ℹ 50 more rows
    ## # ℹ 2 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>

## Categorical Factor Analysis (Bayesian ANOVA)

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
    sample = "sample", 
    cell_group = "cell_group", 
    inference_method = "hmc",
    enable_loo = TRUE,
    verbose = FALSE
  )

# Fit second model
model_without_association = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ 1, 
    sample = "sample", 
    cell_group = "cell_group", 
    inference_method = "hmc",
    enable_loo = TRUE,
    verbose = FALSE
  )

# Compare models
loo_compare(
   attr(model_with_factor_association, "fit")$loo(),
   attr(model_without_association, "fit")$loo()
)
```

    ##        elpd_diff se_diff
    ## model1   0.0       0.0  
    ## model2 -80.6      10.8

## Differential Variability Analysis

We can model the cell-group variability also dependent on the type, and
so test differences in variability

``` r
res = 
  seurat_obj |>
  sccomp_estimate( 
    formula_composition = ~ type, 
    formula_variability = ~ type,
    sample = "sample",
    cell_group = "cell_group",
    cores = 1, verbose = FALSE
  )

res
```

    ## sccomp model
    ## ============
    ## 
    ## Model specifications:
    ##   Family: multi_beta_binomial 
    ##   Composition formula: ~type 
    ##   Variability formula: ~type 
    ##   Inference method: pathfinder 
    ## 
    ## Data: Samples: 20   Cell groups: 30 
    ## 
    ## Column prefixes: c_ -> composition parameters  v_ -> variability parameters
    ## 
    ## Convergence diagnostics:
    ##   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
    ##   scale reduction factor on split chains (at convergence, R_k_hat = 1).
    ## 
    ## # A tibble: 60 × 15
    ##    cell_group        parameter factor c_lower c_effect c_upper c_rhat c_ess_bulk
    ##    <chr>             <chr>     <chr>    <dbl>    <dbl>   <dbl>  <dbl>      <dbl>
    ##  1 B immature        (Interce… <NA>    0.523     0.823   1.12   1.000     2942. 
    ##  2 B immature        typeheal… type    1.03      1.38    1.74   1.01       687. 
    ##  3 B mem             (Interce… <NA>   -1.06     -0.692  -0.313  1.00       290. 
    ##  4 B mem             typeheal… type    1.20      1.61    2.02   1.01       221. 
    ##  5 CD4 cm S100A4     (Interce… <NA>    1.45      1.71    1.96   1.00      1968. 
    ##  6 CD4 cm S100A4     typeheal… type    0.604     0.889   1.17   1.01       515. 
    ##  7 CD4 cm high cyto… (Interce… <NA>   -0.967    -0.584  -0.189  1.00      2508. 
    ##  8 CD4 cm high cyto… typeheal… type   -1.94     -1.27   -0.654  1.09        10.3
    ##  9 CD4 cm ribosome   (Interce… <NA>   -0.0768    0.302   0.654  1.00      3948. 
    ## 10 CD4 cm ribosome   typeheal… type   -1.31     -0.917  -0.496  1.01       743. 
    ## # ℹ 50 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

**Plot 1D significance plot**

``` r
plots = res |> sccomp_test() |> plot()
```

    ## Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
    ## 
    ## Chain 1 finished in 0.0 seconds.

    ## Warning in stat_summary(aes(!!as.symbol(factor_of_interest), (generated_proportions)), : Ignoring unknown parameters: `outlier.shape`, `outlier.colour`, and
    ## `outlier.size`

    ## Warning: `geom_hline()`: Ignoring `mapping` because `yintercept` was provided.

    ## Warning: `geom_hline()`: Ignoring `data` because `yintercept` was provided.

    ## Warning in geom_hline(yintercept = mean_corrected, data = data.frame(parameter
    ## = "(Intercept, corrected)"), : Ignoring unknown parameters: `inherit.aes`

``` r
plots$credible_intervals_1D
```

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

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

    ## Warning: 1 unknown level in `f`: (Intercept, corrected)
    ## 1 unknown level in `f`: (Intercept, corrected)

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

# Recommended Settings for Different Data Types

## For Single-Cell RNA Sequencing

We recommend setting `bimodal_mean_variability_association  = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

## For CyTOF and Microbiome Data

We recommend setting `bimodal_mean_variability_association  = FALSE`
(Default).

## MCMC Chain Visualization

It is possible to directly evaluate the posterior distribution. In this
example, we plot the Monte Carlo chain for the slope parameter of the
first cell type. We can see that it has converged and is negative with
probability 1.

``` r
library(cmdstanr)
library(posterior)
library(bayesplot)

# Assuming res contains the fit object from cmdstanr
fit <- res |> attr("fit")

# Extract draws for 'beta[2,1]'
draws <- as_draws_array(fit$draws("beta[2,1]"))

# Create a traceplot for 'beta[2,1]'
mcmc_trace(draws, pars = "beta[2,1]")
```

![](/Users/a1234450/Documents/GitHub/sccomp/README_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Sonoma 14.6.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Australia/Adelaide
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] bayesplot_1.12.0   posterior_1.6.1    cmdstanr_0.9.0     loo_2.8.0         
    ## [5] tidyr_1.3.1        forcats_1.0.0      ggplot2_3.5.2.9001 sccomp_2.1.14     
    ## [9] dplyr_1.1.4       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rlang_1.1.6                 magrittr_2.0.3             
    ##  [3] matrixStats_1.5.0           compiler_4.5.0             
    ##  [5] callr_3.7.6                 vctrs_0.6.5                
    ##  [7] reshape2_1.4.4              stringr_1.5.1              
    ##  [9] pkgconfig_2.0.3             crayon_1.5.3               
    ## [11] fastmap_1.2.0               backports_1.5.0            
    ## [13] XVector_0.48.0              labeling_0.4.3             
    ## [15] utf8_1.2.6                  rmarkdown_2.29             
    ## [17] tzdb_0.5.0                  UCSC.utils_1.4.0           
    ## [19] ps_1.9.1                    purrr_1.0.4                
    ## [21] xfun_0.52                   cachem_1.1.0               
    ## [23] GenomeInfoDb_1.44.0         jsonlite_2.0.0             
    ## [25] DelayedArray_0.34.1         parallel_4.5.0             
    ## [27] R6_2.6.1                    bslib_0.9.0                
    ## [29] stringi_1.8.7               RColorBrewer_1.1-3         
    ## [31] parallelly_1.45.0           GenomicRanges_1.60.0       
    ## [33] jquerylib_0.1.4             Rcpp_1.1.0                 
    ## [35] SummarizedExperiment_1.38.1 knitr_1.50                 
    ## [37] future.apply_1.20.0         readr_2.1.5                
    ## [39] IRanges_2.42.0              Matrix_1.7-3               
    ## [41] tidyselect_1.2.1            rstudioapi_0.17.1          
    ## [43] abind_1.4-8                 yaml_2.3.10                
    ## [45] codetools_0.2-20            processx_3.8.6             
    ## [47] listenv_0.9.1               lattice_0.22-7             
    ## [49] tibble_3.3.0                plyr_1.8.9                 
    ## [51] Biobase_2.68.0              withr_3.0.2                
    ## [53] S7_0.2.0                    evaluate_1.0.4             
    ## [55] future_1.58.0               pillar_1.11.0              
    ## [57] MatrixGenerics_1.20.0       tensorA_0.36.2.1           
    ## [59] checkmate_2.3.2             stats4_4.5.0               
    ## [61] distributional_0.5.0        generics_0.1.4             
    ## [63] rprojroot_2.0.4             sp_2.2-0                   
    ## [65] S4Vectors_0.46.0            hms_1.1.3                  
    ## [67] scales_1.4.0                globals_0.18.0             
    ## [69] glue_1.8.0                  tools_4.5.0                
    ## [71] data.table_1.17.8           fs_1.6.6                   
    ## [73] dotCall64_1.2               grid_4.5.0                 
    ## [75] instantiate_0.2.3           SingleCellExperiment_1.30.1
    ## [77] GenomeInfoDbData_1.2.14     patchwork_1.3.1            
    ## [79] cli_3.6.5                   spam_2.11-1                
    ## [81] fansi_1.0.6                 S4Arrays_1.8.1             
    ## [83] gtable_0.3.6                sass_0.4.10                
    ## [85] digest_0.6.37               progressr_0.15.1           
    ## [87] BiocGenerics_0.54.0         SparseArray_1.8.0          
    ## [89] ggrepel_0.9.6               SeuratObject_5.1.0         
    ## [91] farver_2.1.2                htmltools_0.5.8.1          
    ## [93] lifecycle_1.0.4             httr_1.4.7
