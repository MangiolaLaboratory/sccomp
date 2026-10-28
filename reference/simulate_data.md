# DEPRECATED: simulate_data

This function is DEPRECATED. Please use
[`sccomp_simulate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_simulate.md)
instead.

## Usage

``` r
simulate_data(
  .data,
  .estimate_object,
  formula_composition,
  formula_variability = NULL,
  .sample = NULL,
  .cell_group = NULL,
  .coefficients = NULL,
  variability_multiplier = 1,
  number_of_draws = 1,
  mcmc_seed = sample_seed(),
  cores = detectCores(),
  sig_figs = 9,
  cache_stan_model = sccomp_stan_models_cache_dir
)
```

## Arguments

- .data:

  A tibble including a cell_group name column \| sample name column \|
  read counts column \| factor columns \| Pvalue column \| a
  significance column

- .estimate_object:

  The result of sccomp_estimate execution. This is used for sampling
  from real-data properties.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment

- formula_variability:

  A formula. The formula describing the model for differential
  variability, for example ~treatment

- .sample:

  A column name as symbol. The sample identifier

- .cell_group:

  A column name as symbol. The cell_group identifier

- .coefficients:

  The column names for coefficients, for example, c(b_0, b_1). Use the
  `coefficients` argument of
  [`sccomp_simulate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_simulate.md)
  instead.

- variability_multiplier:

  A real scalar. This can be used for artificially increasing the
  variability of the simulation for benchmarking purposes.

- number_of_draws:

  An integer. How may copies of the data you want to draw from the model
  joint posterior distribution.

- mcmc_seed:

  An integer. Used for Markov-chain Monte Carlo reproducibility. By
  default a random number is sampled from 1 to 999999. This itself can
  be controlled by set.seed()

- cores:

  Integer, the number of cores to be used for parallel calculations.

- sig_figs:

  Number of significant figures to use for Stan model output. Default is
  9.

- cache_stan_model:

  A character string specifying the cache directory for compiled Stan
  models. The sccomp version will be automatically appended to ensure
  version isolation. Default is `sccomp_stan_models_cache_dir` which
  points to `~/.sccomp_models`.

## Value

See
[`sccomp_simulate`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_simulate.md).
