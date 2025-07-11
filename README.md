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

| Function | Description |
|----|----|
| `sccomp_estimate` | Fit the model onto the data, and estimate the coefficients |
| `sccomp_remove_outliers` | Identify outliers probabilistically based on the model fit, and exclude them from the estimation |
| `sccomp_test` | Calculate the probability that the coefficients are outside the H0 interval (i.e. test_composition_above_logit_fold_change) |
| `sccomp_replicate` | Simulate data from the model, or part of the model |
| `sccomp_predict` | Predicts proportions, based on the model, or part of the model |
| `sccomp_remove_unwanted_variation` | Removes the variability for unwanted factors |
| `plot` | Plots summary plots to assess significance |

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
    sample = "sample", 
    cell_group = "cell_group", 
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
    ##    cell_group parameter   factor c_lower c_effect c_upper   c_pH0   c_FDR c_rhat
    ##    <chr>      <chr>       <chr>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>
    ##  1 B1         (Intercept) <NA>    0.904     1.18   1.44   0       0        1.00 
    ##  2 B1         typecancer  type   -1.07     -0.666 -0.272  0.00175 3.50e-4  1.00 
    ##  3 B2         (Intercept) <NA>    0.432     0.753  1.07   0       0        1.00 
    ##  4 B2         typecancer  type   -1.18     -0.723 -0.263  0.00200 5.00e-4  1.000
    ##  5 B3         (Intercept) <NA>   -0.659    -0.348 -0.0346 0.0597  7.38e-3  1.00 
    ##  6 B3         typecancer  type   -0.774    -0.315  0.135  0.18    5.37e-2  1.000
    ##  7 BM         (Intercept) <NA>   -1.30     -0.988 -0.678  0       0        1.00 
    ##  8 BM         typecancer  type   -0.734    -0.304  0.144  0.171   4.80e-2  1.00 
    ##  9 CD4 1      (Intercept) <NA>    0.155     0.340  0.519  0.00525 4.13e-4  1.00 
    ## 10 CD4 1      typecancer  type   -0.0769    0.170  0.412  0.288   6.92e-2  1.00 
    ## # ℹ 62 more rows
    ## # ℹ 10 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>, v_lower <dbl>,
    ## #   v_effect <dbl>, v_upper <dbl>, v_pH0 <dbl>, v_FDR <dbl>, v_rhat <dbl>,
    ## #   v_ess_bulk <dbl>, v_ess_tail <dbl>

## Outlier identification

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

    ## Warning in stat_summary(aes(!!as.symbol(factor_of_interest), (generated_proportions)), : Ignoring unknown parameters: `outlier.shape`, `outlier.colour`, and
    ## `outlier.size`

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

    ## Warning in stat_summary(aes(!!as.symbol(factor_of_interest), (generated_proportions)), : Ignoring unknown parameters: `outlier.shape`, `outlier.colour`, and
    ## `outlier.size`

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

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

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

    ## Warning: annotation$theme is not a valid theme.
    ## Please use `theme()` to construct themes.

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
    sample = "sample",
    cell_group = "cell_group",
    abundance = "proportion", 
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
```

## Continuous factor

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
```

    ## Loading required namespace: SeuratObject

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
    ##  1 B immature        (Interce… <NA>    0.379    0.758   1.12    1.00       4213.
    ##  2 B immature        typeheal… type    0.847    1.36    1.86    1.000      4153.
    ##  3 B immature        continuo… conti… -0.236    0.0547  0.345   1.00       3737.
    ##  4 B mem             (Interce… <NA>   -1.24    -0.804  -0.357   1.00       3904.
    ##  5 B mem             typeheal… type    1.06     1.66    2.27    1.00       3900.
    ##  6 B mem             continuo… conti… -0.228    0.0899  0.421   1.00       4116.
    ##  7 CD4 cm S100A4     (Interce… <NA>    1.19     1.50    1.82    1.00       4184.
    ##  8 CD4 cm S100A4     typeheal… type    0.685    1.12    1.56    1.000      4076.
    ##  9 CD4 cm S100A4     continuo… conti… -0.0726   0.189   0.445   1.000      3759.
    ## 10 CD4 cm high cyto… (Interce… <NA>   -0.941   -0.460   0.0213  1.000      4035.
    ## # ℹ 80 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

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
    sample = "sample",
    cell_group = "cell_group",
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
    ##  1 B immature (Intercept)      <NA>    0.510    0.827  1.17      1.00      128. 
    ##  2 B immature typehealthy      type    0.585    1.07   1.49      1.01       96.9
    ##  3 B immature (Intercept)___G… <NA>   -0.356    0.125  0.711     1.01       95.0
    ##  4 B immature (Intercept)___G… <NA>   -0.0298   0.301  0.835     1.00      121. 
    ##  5 B immature (Intercept)___G… <NA>   -0.0678   0.243  0.656     1.00      187. 
    ##  6 B immature (Intercept)___G… <NA>   -0.797   -0.327  0.00890   1.02      165. 
    ##  7 B mem      (Intercept)      <NA>   -0.741   -0.327  0.100     1.00      115. 
    ##  8 B mem      typehealthy      type    0.210    0.982  1.54      1.01       85.7
    ##  9 B mem      (Intercept)___G… <NA>   -0.400    0.0851 0.752     1.00      113. 
    ## 10 B mem      (Intercept)___G… <NA>   -0.0228   0.341  0.953     1.00      143. 
    ## # ℹ 170 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

### Random Effect Model (random slopes)

`sccomp` can model random slopes. We providean example below.

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
    ##  1 B immature (Intercept)      <NA>    0.457    0.861   1.33     1.00       96.4
    ##  2 B immature typehealthy      type    0.524    1.00    1.52     1.01       96.8
    ##  3 B immature (Intercept)___G… <NA>   -0.302    0.0696  0.475    1.01       88.4
    ##  4 B immature typehealthy___G… <NA>   -0.234    0.0549  0.437    1.01      115. 
    ##  5 B immature (Intercept)___G… <NA>   -0.106    0.154   0.523    1.00      134. 
    ##  6 B immature typehealthy___G… <NA>   -0.103    0.140   0.535    1.03      101. 
    ##  7 B immature (Intercept)___G… <NA>   -0.0764   0.176   0.579    1.01       98.2
    ##  8 B immature (Intercept)___G… <NA>   -0.661   -0.233   0.0545   1.02       65.0
    ##  9 B mem      (Intercept)      <NA>   -0.827   -0.404  -0.0354   1.00      209. 
    ## 10 B mem      typehealthy      type    0.429    1.01    1.64     1.02       59.1
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
    ##  1 B immature (Intercept)      <NA>    0.463    0.814   1.39     1.00       55.5
    ##  2 B immature typehealthy      type    0.566    1.03    1.52     1.02       74.3
    ##  3 B immature (Intercept)___G… <NA>   -0.261    0.0980  0.515    1.04       87.0
    ##  4 B immature typehealthy___G… <NA>   -0.214    0.0994  0.554    1.04       59.3
    ##  5 B immature (Intercept)___G… <NA>   -0.189    0.0832  0.426    1.02      113. 
    ##  6 B immature typehealthy___G… <NA>   -0.147    0.0902  0.385    1.00       77.6
    ##  7 B immature (Intercept)___G… <NA>   -0.0822   0.178   0.584    1.01       53.1
    ##  8 B immature (Intercept)___G… <NA>   -0.786   -0.273   0.0332   1.02       40.9
    ##  9 B immature (Intercept)___G… <NA>   -0.366   -0.0550  0.260    1.02       83.5
    ## 10 B immature (Intercept)___G… <NA>   -0.0169   0.197   0.621    1.00       38.3
    ## # ℹ 290 more rows
    ## # ℹ 7 more variables: c_ess_tail <dbl>, v_lower <dbl>, v_effect <dbl>,
    ## #   v_upper <dbl>, v_rhat <dbl>, v_ess_bulk <dbl>, v_ess_tail <dbl>

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
    ##  1 B1         2-fold decrease (from 0.0585 to 0.0295)  
    ##  2 B2         2.1-fold decrease (from 0.0383 to 0.0182)
    ##  3 B3         1.4-fold decrease (from 0.0129 to 0.0091)
    ##  4 BM         1.4-fold decrease (from 0.0066 to 0.0048)
    ##  5 CD4 1      1.2-fold increase (from 0.0249 to 0.0297)
    ##  6 CD4 2      1.4-fold increase (from 0.0511 to 0.0732)
    ##  7 CD4 3      2.6-fold decrease (from 0.0856 to 0.0332)
    ##  8 CD4 4      1.1-fold increase (from 0.0017 to 0.0018)
    ##  9 CD4 5      1-fold increase (from 0.0305 to 0.0315)  
    ## 10 CD8 1      1.1-fold increase (from 0.11 to 0.1205)  
    ## # ℹ 26 more rows

## Contrasts

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
    ##  1 B immature   typecanc… <NA>    -1.90    -1.34   -0.793 2.50e-4 4.17e-5     NA
    ##  2 B immature   typeheal… <NA>     0.793    1.34    1.90  2.50e-4 4.17e-5     NA
    ##  3 B mem        typecanc… <NA>    -2.25    -1.65   -1.03  0       0           NA
    ##  4 B mem        typeheal… <NA>     1.03     1.65    2.25  0       0           NA
    ##  5 CD4 cm S100… typecanc… <NA>    -1.45    -0.990  -0.535 0       0           NA
    ##  6 CD4 cm S100… typeheal… <NA>     0.535    0.990   1.45  0       0           NA
    ##  7 CD4 cm high… typecanc… <NA>     0.822    1.54    2.23  0       0           NA
    ##  8 CD4 cm high… typeheal… <NA>    -2.23    -1.54   -0.822 0       0           NA
    ##  9 CD4 cm ribo… typecanc… <NA>     0.296    0.944   1.60  7.25e-3 2.10e-3     NA
    ## 10 CD4 cm ribo… typeheal… <NA>    -1.60    -0.944  -0.296 7.25e-3 2.10e-3     NA
    ## # ℹ 50 more rows
    ## # ℹ 2 more variables: c_ess_bulk <dbl>, c_ess_tail <dbl>

## Categorical factor (e.g. Bayesian ANOVA)

This is achieved through model comparison with `loo`. In the following
example, the model with association with factors better fits the data
compared to the baseline model with no factor association. For model
comparisons `sccomp_remove_outliers()` must not be executed as the
leave-one-out must work with the same amount of data, while outlier
elimination does not guarantee it.

If `elpd_diff`

## The old framework

The new tidy framework was introduced in 2024. To understand the
differences and improvements compared to the old framework, please read
this [blog
post](https://tidyomics.github.io/tidyomicsBlog/post/2023-12-07-tidy-sccomp/).
