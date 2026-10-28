# sccomp_simulate

Simulate counts from a fitted sccomp model, optionally replacing the
abundance coefficients and the mean-variability association with
user-chosen values. Variability is not taken from the posterior: for
every variability term it is derived from the abundance coefficient it
is associated with, `alpha = beta * slope + intercept`, plus the fitted
residual scatter around that trend (see `variability_scatter`), so
changing the abundance coefficients changes the simulated variability as
the model expects. Random effects and the remaining parameters come from
the posterior draws.

## Usage

``` r
sccomp_simulate(
  fit,
  formula_composition,
  formula_variability = NULL,
  new_data = NULL,
  coefficients = NULL,
  mean_dispersion_slope = NULL,
  mean_dispersion_intercept = NULL,
  sample = NULL,
  cell_group = NULL,
  variability_multiplier = 1,
  number_of_draws = 1,
  mcmc_seed = sample_seed(),
  cores = detectCores(),
  sig_figs = 9,
  cache_stan_model = sccomp_stan_models_cache_dir,
  variability_scatter = TRUE
)
```

## Arguments

- fit:

  The result of sccomp_estimate execution. This is used for sampling
  from real-data properties.

- formula_composition:

  A formula. The formula describing the model for differential
  abundance, for example ~treatment. It can be a sub-formula of the
  fitted one.

- formula_variability:

  A formula. The formula describing the model for differential
  variability, for example ~treatment. If NULL, the fitted variability
  formula is used.

- new_data:

  A tibble including sample-specific columns (sample identifier and
  factor columns from formula). If coefficients are provided separately,
  cell_group column is optional. Otherwise, should include cell_group
  column. If NULL, the data used for the fit is used.

- coefficients:

  A data frame/tibble with cell-type specific coefficients. Must contain
  a column matching the cell_group column name with every fitted cell
  group, and columns matching the design matrix column names (e.g.,
  "(Intercept)", "typecancer"). Design columns that are absent are set
  to zero. For each design column, the coefficients must sum to zero
  across cell groups. If NULL, the posterior coefficients are used.

- mean_dispersion_slope:

  Optional slope of the mean-variability association. Either a single
  number, applied to every variability term, or a vector named by
  variability design column (e.g. `c("(Intercept)" = -0.7)`), in which
  case unnamed terms keep their posterior slope. If NULL, the posterior
  slopes are used, or 0 for a model fitted with
  `exclude_mean_variability_association = TRUE`.

- mean_dispersion_intercept:

  Optional intercept of the mean-variability association, in the same
  forms as `mean_dispersion_slope`. If NULL, the posterior intercepts
  are used, except for terms whose slope is overridden, which take their
  posterior-mean intercept so that only the slope changes.

- sample:

  A column name as a character string for the sample identifier. If
  NULL, the one used for the fit.

- cell_group:

  A column name as a character string for the cell-group identifier. If
  NULL, the one used for the fit.

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

- variability_scatter:

  Logical. If `TRUE` (default), each cell group's variability is drawn
  around the mean-variability trend with the fitted residual standard
  deviation (`prec_sd`), as fitted cell groups scatter around it. If
  `FALSE`, every cell group sits exactly on the trend.

## Value

The rows of `new_data` (expanded to every cell group if `new_data` has
no cell_group column), with the following columns added:

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

#     # Coefficients for every cell group. All zero for simplicity.
#     coefficients = counts_obj |>
#       distinct(cell_group) |>
#       mutate(`(Intercept)` = 0, typecancer = 0)

#     # Simulate data
#     sccomp_simulate(
#       estimate, ~type, ~1,
#       new_data = counts_obj |> distinct(sample, type),
#       coefficients = coefficients
#     )
#   }
# }
```
