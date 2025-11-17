# sccomp: Differential Composition and Variability Analysis for Single-Cell Data

Abstract

Sccomp is a comprehensive R package for differential composition and
variability analysis in single-cell RNA sequencing, CyTOF, and
microbiome data. It provides robust Bayesian modeling with outlier
detection, random effects, and advanced statistical methods for cell
type proportion analysis. Perfect for cancer research, immunology, and
single-cell genomics.

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![R build
status](https://github.com/stemangiola/sccomp/workflows/R-CMD-check/badge.svg)](https://github.com/stemangiola/sccomp/actions/)

## ![](../inst/logo-01.png)

**sccomp** is a powerful R package designed for comprehensive
differential composition and variability analysis in single-cell
genomics, proteomics, and microbiomics data.

### Why sccomp?

For cellular omic data, no method for differential variability analysis
exists, and methods for differential composition analysis only take a
few fundamental data properties into account. Here we introduce
**sccomp**, a generalised method for differential composition and
variability analyses capable of jointly modelling data count
distribution, compositionality, group-specific variability, and
proportion mean-variability association, while being robust to outliers.

![](../inst/cartoon_methods.png)

#### Comprehensive Method Comparison

- **I**: Data are modelled as counts.
- **II**: Group proportions are modelled as compositional.
- **III**: The proportion variability is modelled as cell-type specific.
- **IV**: Information sharing across cell types, mean–variability
  association.
- **V**: Outlier detection or robustness.
- **VI**: Differential variability analysis.
- **VII** Mixed effect modelling
- **VIII** Removal unwanted effects

| Method           | Year | Model                         | I   | II  | III | IV  | V   | VI  | VII | VIII |
|------------------|------|-------------------------------|-----|-----|-----|-----|-----|-----|-----|------|
| **sccomp**       | 2023 | Sum-constrained Beta-binomial | ●   | ●   | ●   | ●   | ●   | ●   | ●   | ●    |
| **scCODA**       | 2021 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |     |      |
| **quasi-binom.** | 2021 | Quasi-binomial                | ●   |     | ●   |     |     |     |     |      |
| **rlm**          | 2021 | Robust-log-linear             |     | ●   |     |     | ●   |     |     |      |
| **propeller**    | 2021 | Logit-linear + limma          |     | ●   | ●   | ●   |     |     |     |      |
| **ANCOM-BC**     | 2020 | Log-linear                    |     | ●   | ●   |     |     |     |     |      |
| **corncob**      | 2020 | Beta-binomial                 | ●   |     | ●   |     |     |     |     |      |
| **scDC**         | 2019 | Log-linear                    |     | ●   | ●   |     |     |     |     |      |
| **dmbvs**        | 2017 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |     |      |
| **MixMC**        | 2016 | Zero-inflated Log-linear      |     | ●   | ●   |     |     |     |     |      |
| **ALDEx2**       | 2014 | Dirichlet-multinomial         | ●   | ●   |     |     |     |     |     |      |

#### Scientific Citation

Mangiola, Stefano, Alexandra J. Roth-Schulze, Marie Trussart, Enrique
Zozaya-Valdés, Mengyao Ma, Zijie Gao, Alan F. Rubin, Terence P. Speed,
Heejung Shim, and Anthony T. Papenfuss. 2023. “Sccomp: Robust
Differential Composition and Variability Analysis for Single-Cell Data.”
Proceedings of the National Academy of Sciences of the United States of
America 120 (33): e2203828120. <https://doi.org/10.1073/pnas.2203828120>
[PNAS - sccomp: Robust differential composition and variability analysis
for single-cell
data](https://www.pnas.org/doi/full/10.1073/pnas.2203828120)

#### Talk

[![Watch the
video](https://img.youtube.com/vi/R_lt58We9nA/mqdefault.jpg)](https://www.youtube.com/watch?v=R_lt58We9nA&ab_channel=RConsortium)

## Installation Guide

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

### Core Functions

| Function                         | Description                                                                                                                 |
|----------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| `sccomp_estimate`                | Fit the model onto the data, and estimate the coefficients                                                                  |
| `sccomp_remove_outliers`         | Identify outliers probabilistically based on the model fit, and exclude them from the estimation                            |
| `sccomp_test`                    | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate`               | Simulate data from the model, or part of the model                                                                          |
| `sccomp_predict`                 | Predicts proportions, based on the model, or part of the model                                                              |
| `sccomp_remove_unwanted_effects` | Removes the variability for unwanted factors                                                                                |
| `plot`                           | Plots summary plots to assess significance                                                                                  |

## Analysis Tutorial

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

### Binary Factor Analysis

Of the output table, the estimate columns start with the prefix `c_`
indicate `composition`, or with `v_` indicate `variability` (when
formula_variability is set).

#### From Seurat, SingleCellExperiment, metadata objects

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

#### From counts

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
    ##    cell_group parameter  factor  c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>      <chr>      <chr>     <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B1         (Intercep… NA      0.949      1.20   1.46   0       0        1.000
    ##  2 B1         typecancer type   -0.926     -0.613 -0.308  0       0        1.000
    ##  3 B2         (Intercep… NA      0.509      0.770  1.03   0       0        1.00 
    ##  4 B2         typecancer type   -0.968     -0.669 -0.357  0.00100 9.09e-5  1.00 
    ##  5 B3         (Intercep… NA     -0.591     -0.333 -0.0663 0.0413  3.55e-3  1.00 
    ##  6 B3         typecancer type   -0.573     -0.277  0.0289 0.118   2.39e-2  1.000
    ##  7 BM         (Intercep… NA     -1.23      -0.968 -0.717  0       0        1.00 
    ##  8 BM         typecancer type   -0.595     -0.287  0.0199 0.109   1.98e-2  1.00 
    ##  9 CD4 1      (Intercep… NA      0.205      0.368  0.533  0       0        1.00 
    ## 10 CD4 1      typecancer type   -0.00772    0.209  0.422  0.162   2.94e-2  1.000
    ## # ℹ 62 more rows
    ## # ℹ 10 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_rhat <dbl>,
    ## #   v_ess_bulk <dbl>, v_ess_tail <dbl>

### Outlier Identification

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

### Visualization and Summary Plots

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

![](inst/figures/unnamed-chunk-8-1.png)

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

![](inst/figures/unnamed-chunk-9-1.png)

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

![](inst/figures/unnamed-chunk-10-1.png)

We can plot the relationship between abundance and variability. As we
can see below, they are positively correlated. sccomp models this
relationship to obtain a shrinkage effect on the estimates of both the
abundance and the variability. This shrinkage is adaptive as it is
modelled jointly, thanks to Bayesian inference.

``` r
sccomp_result |> 
  plot_2D_intervals()
```

![](inst/figures/unnamed-chunk-11-1.png)

You can produce the series of plots calling the `plot` method.

``` r
sccomp_result |> plot() 
```

### Model Proportions Directly (e.g. from deconvolution)

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

### Continuous Factor Analysis

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
    ##  1 B immature        (Interce… NA       0.548   0.835    1.13   1.000      3595.
    ##  2 B immature        typeheal… type     1.02    1.35     1.68   1.00       3617.
    ##  3 B immature        continuo… conti…  -0.262   0.0541   0.374  1.000      3718.
    ##  4 B mem             (Interce… NA      -0.997  -0.674   -0.352  1.000      3932.
    ##  5 B mem             typeheal… type     1.22    1.57     1.94   1.00       4304.
    ##  6 B mem             continuo… conti…  -0.269   0.0675   0.413  1.00       4165.
    ##  7 CD4 cm S100A4     (Interce… NA       1.32    1.55     1.80   1.00       3714.
    ##  8 CD4 cm S100A4     typeheal… type     0.859   1.13     1.41   1.00       4075.
    ##  9 CD4 cm S100A4     continuo… conti…  -0.109   0.169    0.466  1.00       3746.
    ## 10 CD4 cm high cyto… (Interce… NA      -0.932  -0.582   -0.244  1.00       3829.
    ## # ℹ 80 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

### Random Effect Modeling (Mixed-Effect Modeling)

`sccomp` supports multilevel modeling by allowing the inclusion of
random effects in the compositional and variability formulas. This is
particularly useful when your data has hierarchical or grouped
structures, such as measurements nested within subjects, batches, or
experimental units. By incorporating random effects, sccomp can account
for variability at different levels of your data, improving model fit
and inference accuracy.

#### Random Intercept Model

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
    ##  1 B immature (Intercept)      NA      0.535   0.836    1.22     1.01      149. 
    ##  2 B immature typehealthy      type    0.778   1.21     1.61     1.01      127. 
    ##  3 B immature (Intercept)___G… NA     -0.309   0.0223   0.407    1.00      165. 
    ##  4 B immature (Intercept)___G… NA     -0.0872  0.232    0.653    1.00      122. 
    ##  5 B immature (Intercept)___G… NA     -0.134   0.231    0.617    1.01      136. 
    ##  6 B immature (Intercept)___G… NA     -0.802  -0.288    0.0357   1.02      111. 
    ##  7 B mem      (Intercept)      NA     -0.858  -0.489   -0.0126   1.01       68.9
    ##  8 B mem      typehealthy      type    0.584   1.32     1.81     1.02       61.5
    ##  9 B mem      (Intercept)___G… NA     -0.488  -0.00999  0.477    1.03      105. 
    ## 10 B mem      (Intercept)___G… NA     -0.151   0.254    0.780    1.00      128. 
    ## # ℹ 170 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

#### Random Effect Model (random slopes)

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
    ##  1 B immature (Intercept)      NA       0.491   0.861   1.28     1.00      131. 
    ##  2 B immature typehealthy      type     0.697   1.18    1.66     1.00      114. 
    ##  3 B immature (Intercept)___G… NA      -0.294   0.0331  0.402    1.00      121. 
    ##  4 B immature typehealthy___G… NA      -0.277   0.0131  0.348    1.03      190. 
    ##  5 B immature (Intercept)___G… NA      -0.177   0.132   0.559    1.00      136. 
    ##  6 B immature typehealthy___G… NA      -0.172   0.0975  0.578    1.00      145. 
    ##  7 B immature (Intercept)___G… NA      -0.227   0.172   0.574    1.01      122. 
    ##  8 B immature (Intercept)___G… NA      -0.785  -0.224   0.0407   1.00       97.7
    ##  9 B mem      (Intercept)      NA      -1.06   -0.546  -0.101    1.01      114. 
    ## 10 B mem      typehealthy      type     0.723   1.35    1.97     1.01       81.8
    ## # ℹ 230 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

#### Nested Random Effects

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
    ##  1 B immature (Intercept)      NA      0.459    0.803   1.29     1.01      101. 
    ##  2 B immature typehealthy      type    0.617    1.18    1.71     1.03       97.9
    ##  3 B immature (Intercept)___G… NA     -0.226    0.133   0.538    1.00       65.6
    ##  4 B immature typehealthy___G… NA     -0.239    0.0784  0.449    1.02       72.9
    ##  5 B immature (Intercept)___G… NA     -0.311    0.0613  0.374    1.02       99.1
    ##  6 B immature typehealthy___G… NA     -0.179    0.0680  0.361    1.00      214. 
    ##  7 B immature (Intercept)___G… NA     -0.0129   0.254   0.654    1.00      114. 
    ##  8 B immature (Intercept)___G… NA     -0.734   -0.313   0.0604   1.03       48.4
    ##  9 B immature (Intercept)___G… NA     -0.493   -0.151   0.0807   1.02       82.5
    ## 10 B immature (Intercept)___G… NA      0.0101   0.241   0.638    1.01       90.8
    ## # ℹ 290 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

### Result Interpretation and Communication

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
    ##  1 B immature           2-fold decrease (from 0.1068 to 0.0527)  
    ##  2 B mem                2.3-fold decrease (from 0.033 to 0.0143) 
    ##  3 CD4 cm high cytokine 7.8-fold increase (from 0.0016 to 0.0123)
    ##  4 CD4 cm ribosome      3.5-fold increase (from 0.0076 to 0.0264)
    ##  5 CD4 cm S100A4        1.5-fold decrease (from 0.1416 to 0.0951)
    ##  6 CD4 em high cytokine 4.5-fold increase (from 0.0023 to 0.0105)
    ##  7 CD4 naive            1.5-fold decrease (from 0.116 to 0.078)  
    ##  8 CD4 ribosome         2.7-fold decrease (from 0.0832 to 0.0312)
    ##  9 CD8 em 1             1.2-fold increase (from 0.0497 to 0.059) 
    ## 10 CD8 em 2             3.7-fold increase (from 0.006 to 0.0219) 
    ## # ℹ 20 more rows

### Contrasts Analysis

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
    ##  1 B immature   typecanc… NA      -1.90    -1.35   -0.806 0       0           NA
    ##  2 B immature   typeheal… NA       0.806    1.35    1.90  0       0           NA
    ##  3 B mem        typecanc… NA      -2.22    -1.65   -1.07  0       0           NA
    ##  4 B mem        typeheal… NA       1.07     1.65    2.22  0       0           NA
    ##  5 CD4 cm S100… typecanc… NA      -1.46    -1.00   -0.509 0       0           NA
    ##  6 CD4 cm S100… typeheal… NA       0.509    1.00    1.46  0       0           NA
    ##  7 CD4 cm high… typecanc… NA       0.939    1.57    2.20  0       0           NA
    ##  8 CD4 cm high… typeheal… NA      -2.20    -1.57   -0.939 0       0           NA
    ##  9 CD4 cm ribo… typecanc… NA       0.342    0.951   1.55  0.00525 0.00129     NA
    ## 10 CD4 cm ribo… typeheal… NA      -1.55    -0.951  -0.342 0.00525 0.00129     NA
    ## # ℹ 50 more rows
    ## # ℹ 2 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>

### Categorical Factor Analysis (Bayesian ANOVA)

This is achieved through model comparison with `loo`. In the following
example, the model with association with factors better fits the data
compared to the baseline model with no factor association. For model
comparisons
[`sccomp_remove_outliers()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_remove_outliers.md)
must not be executed as the leave-one-out must work with the same amount
of data, while outlier elimination does not guarantee it.

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
    ## model2 -83.9      10.5

### Differential Variability Analysis

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
    ##  1 B immature        (Interce… NA      0.513     0.824   1.14    1.00     1718. 
    ##  2 B immature        typeheal… type    1.02      1.37    1.73    1.00     1257. 
    ##  3 B mem             (Interce… NA     -1.06     -0.677  -0.284   1.01      249. 
    ##  4 B mem             typeheal… type    1.13      1.58    2.02    1.00      236. 
    ##  5 CD4 cm S100A4     (Interce… NA      1.44      1.71    1.96    1.00     2082. 
    ##  6 CD4 cm S100A4     typeheal… type    0.589     0.889   1.19    1.00      504. 
    ##  7 CD4 cm high cyto… (Interce… NA     -0.972    -0.594  -0.208   1.00     1430. 
    ##  8 CD4 cm high cyto… typeheal… type   -1.77     -1.18   -0.444   1.01       56.0
    ##  9 CD4 cm ribosome   (Interce… NA     -0.0513    0.313   0.696   1.00     3474. 
    ## 10 CD4 cm ribosome   typeheal… type   -1.35     -0.931  -0.481   1.00      519. 
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

``` r
plots$credible_intervals_1D
```

![](inst/figures/unnamed-chunk-23-1.png)

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

![](inst/figures/unnamed-chunk-24-1.png)

## Recommended Settings for Different Data Types

### For Single-Cell RNA Sequencing

We recommend setting `bimodal_mean_variability_association = TRUE`. The
bimodality of the mean-variability association can be confirmed from the
plots\$credible_intervals_2D (see below).

### For CyTOF and Microbiome Data

We recommend setting `bimodal_mean_variability_association = FALSE`
(Default).

### MCMC Chain Visualization

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
mcmc_trace(draws, pars = "beta[2,1]") + theme_bw()
```

![](inst/figures/unnamed-chunk-25-1.png)

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] bayesplot_1.14.0  posterior_1.6.1   cmdstanr_0.9.0    loo_2.8.0        
    ##  [5] tidyr_1.3.1       forcats_1.0.1     ggplot2_4.0.0     sccomp_2.1.20    
    ##  [9] instantiate_0.2.3 dplyr_1.1.4      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1            farver_2.1.2               
    ##  [3] S7_0.2.0                    fastmap_1.2.0              
    ##  [5] SingleCellExperiment_1.32.0 tensorA_0.36.2.1           
    ##  [7] dotCall64_1.2               digest_0.6.38              
    ##  [9] lifecycle_1.0.4             SeuratObject_5.2.0         
    ## [11] processx_3.8.6              magrittr_2.0.4             
    ## [13] compiler_4.5.2              rlang_1.1.6                
    ## [15] sass_0.4.10                 tools_4.5.2                
    ## [17] utf8_1.2.6                  yaml_2.3.10                
    ## [19] data.table_1.17.8           knitr_1.50                 
    ## [21] labeling_0.4.3              S4Arrays_1.10.0            
    ## [23] htmlwidgets_1.6.4           sp_2.2-0                   
    ## [25] DelayedArray_0.36.0         plyr_1.8.9                 
    ## [27] RColorBrewer_1.1-3          abind_1.4-8                
    ## [29] withr_3.0.2                 purrr_1.2.0                
    ## [31] BiocGenerics_0.56.0         desc_1.4.3                 
    ## [33] grid_4.5.2                  stats4_4.5.2               
    ## [35] future_1.67.0               progressr_0.18.0           
    ## [37] globals_0.18.0              scales_1.4.0               
    ## [39] SummarizedExperiment_1.40.0 cli_3.6.5                  
    ## [41] rmarkdown_2.30              crayon_1.5.3               
    ## [43] ragg_1.5.0                  generics_0.1.4             
    ## [45] future.apply_1.20.0         reshape2_1.4.5             
    ## [47] tzdb_0.5.0                  cachem_1.1.0               
    ## [49] stringr_1.6.0               parallel_4.5.2             
    ## [51] XVector_0.50.0              matrixStats_1.5.0          
    ## [53] vctrs_0.6.5                 Matrix_1.7-4               
    ## [55] jsonlite_2.0.0              callr_3.7.6                
    ## [57] IRanges_2.44.0              hms_1.1.4                  
    ## [59] patchwork_1.3.2             S4Vectors_0.48.0           
    ## [61] ggrepel_0.9.6               listenv_0.10.0             
    ## [63] systemfonts_1.3.1           jquerylib_0.1.4            
    ## [65] spam_2.11-1                 parallelly_1.45.1          
    ## [67] glue_1.8.0                  pkgdown_2.2.0              
    ## [69] codetools_0.2-20            ps_1.9.1                   
    ## [71] distributional_0.5.0        stringi_1.8.7              
    ## [73] gtable_0.3.6                GenomicRanges_1.62.0       
    ## [75] tibble_3.3.0                pillar_1.11.1              
    ## [77] htmltools_0.5.8.1           Seqinfo_1.0.0              
    ## [79] R6_2.6.1                    textshaping_1.0.4          
    ## [81] evaluate_1.0.5              lattice_0.22-7             
    ## [83] Biobase_2.70.0              readr_2.1.5                
    ## [85] backports_1.5.0             bslib_0.9.0                
    ## [87] Rcpp_1.1.0                  SparseArray_1.10.1         
    ## [89] checkmate_2.3.3             xfun_0.54                  
    ## [91] fs_1.6.6                    MatrixGenerics_1.22.0      
    ## [93] prettydoc_0.4.1             pkgconfig_2.0.3
