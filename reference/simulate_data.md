# simulate_data

This function simulates data from a fitted model.

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
  variability_multiplier = 5,
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

  The column names for coefficients, for example, c(b_0, b_1)

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

A tibble (`tbl`) with the following columns:

- **sample** - A character column representing the sample name.

- **type** - A factor column representing the type of the sample.

- **phenotype** - A factor column representing the phenotype in the
  data.

- **count** - An integer column representing the original cell counts.

- **cell_group** - A character column representing the cell group
  identifier.

- **b_0** - A numeric column representing the first coefficient used for
  simulation.

- **b_1** - A numeric column representing the second coefficient used
  for simulation.

- **generated_proportions** - A numeric column representing the
  generated proportions from the simulation.

- **generated_counts** - An integer column representing the generated
  cell counts from the simulation.

- **replicate** - An integer column representing the replicate number
  for each draw from the posterior distribution.

## References

S. Mangiola, A.J. Roth-Schulze, M. Trussart, E. Zozaya-Valdés, M. Ma, Z.
Gao, A.F. Rubin, T.P. Speed, H. Shim, & A.T. Papenfuss, sccomp: Robust
differential composition and variability analysis for single-cell data,
Proc. Natl. Acad. Sci. U.S.A. 120 (33) e2203828120,
https://doi.org/10.1073/pnas.2203828120 (2023).

## Examples

``` r
# print("cmdstanr is needed to run this example.")
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

# \donttest{
#   if (instantiate::stan_cmdstan_exists()) {
#     data("counts_obj")
#     library(dplyr)

#     estimate = sccomp_estimate(
#       counts_obj,
#       ~ type, ~1, "sample", "cell_group", "count",
#       cores = 1
#     )

#     # Set coefficients for cell_groups. In this case all coefficients are 0 for simplicity.
#     counts_obj = counts_obj |> mutate(b_0 = 0, b_1 = 0)

#     # Simulate data
#     simulate_data(counts_obj, estimate, ~type, ~1, sample, cell_group, c(b_0, b_1))
#   }
# }
```
