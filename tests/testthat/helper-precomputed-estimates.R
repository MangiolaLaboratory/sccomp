# Expensive pre-fitted sccomp objects shared by plot and association tests.
# Helpers run before any `test_*.R` file so objects exist regardless of
# alphabetical test order.

if (requireNamespace("instantiate", quietly = TRUE) && instantiate::stan_cmdstan_exists()) {
  data("seurat_obj", package = "sccomp", envir = environment())

  n_iterations <- 1000L
  set.seed(42)

  my_estimate <-
    seurat_obj |>
    sccomp::sccomp_estimate(
      formula_composition = ~ continuous_covariate * type,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )

  my_estimate_with_variance <-
    seurat_obj |>
    sccomp::sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ type,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )

  my_estimate_exclude_mean_variability_association <-
    seurat_obj |>
    sccomp::sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ type,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      exclude_mean_variability_association = TRUE,
      verbose = FALSE
    )

  my_estimate_with_variance_bimodal <-
    seurat_obj |>
    sccomp::sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ type,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      bimodal_mean_variability_association = TRUE,
      verbose = FALSE
    )

  my_estimate_intercept_only <-
    seurat_obj |>
    sccomp::sccomp_estimate(
      formula_composition = ~ 1,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )
}
