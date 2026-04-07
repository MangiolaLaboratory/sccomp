#' Unit tests for restricted environment cache handling
#'
#' Tests both approaches for using custom Stan model cache in restricted environments:
#' 1. assignInNamespace - set global cache before any sccomp call
#' 2. cache_stan_model - pass explicitly to sccomp_estimate, sccomp_remove_outliers, sccomp_boxplot
#'
#' @seealso https://github.com/MangiolaLaboratory/sccomp/issues/250
library(dplyr)
library(sccomp)
data("counts_obj")

test_that("assignInNamespace: custom cache dir works for sccomp_estimate, sccomp_remove_outliers, sccomp_boxplot", {
  skip_cmdstan()

  # Use temp dir as custom cache (simulates /opt/sccomp_models on a server)
  custom_cache <- file.path(tempdir(), "sccomp_test_cache_assignInNamespace")
  dir.create(custom_cache, showWarnings = FALSE, recursive = TRUE)

  # Save original and restore on exit to avoid affecting other tests
  orig_cache <- get("sccomp_stan_models_cache_dir", envir = asNamespace("sccomp"))
  on.exit(
    utils::assignInNamespace("sccomp_stan_models_cache_dir", orig_cache, ns = "sccomp"),
    add = TRUE
  )

  # Set custom cache before any sccomp call
  utils::assignInNamespace("sccomp_stan_models_cache_dir", custom_cache, ns = "sccomp")

  # Run full pipeline - estimate, remove_outliers, test (isolate output dir per test)
  output_dir <- file.path(tempdir(), "sccomp_draws_assignInNamespace")
  result <- expect_no_error({
    counts_obj |>
      sccomp_estimate(
        formula_composition = ~ type,
        formula_variability = ~ 1,
        sample = "sample",
        cell_group = "cell_group",
        abundance = "count",
        cores = 1,
        inference_method = "pathfinder",
        max_sampling_iterations = 500,
        verbose = FALSE,
        output_directory = output_dir
      ) |>
      sccomp_remove_outliers(
        percent_false_positive = 80,
        cores = 1,
        inference_method = "pathfinder",
        max_sampling_iterations = 300,
        verbose = FALSE,
        output_directory = output_dir
      ) |>
      sccomp_test()
  })

  # Verify compiled models were written to the custom cache
  sccomp_version <- as.character(packageVersion("sccomp"))
  versioned_cache <- file.path(custom_cache, sccomp_version)
  expect_true(dir.exists(versioned_cache), info = "versioned cache dir should exist")
  cached_rds <- list.files(versioned_cache, pattern = "\\.rds$")
  expect_true(length(cached_rds) >= 2L, info = "cache should contain at least 2 compiled models (estimate + generate_data)")
  expect_true("glm_multi_beta_binomial.rds" %in% cached_rds, info = "main model should be in cache")
  expect_true("glm_multi_beta_binomial_generate_data.rds" %in% cached_rds, info = "generate_data model should be in cache")

  # Cleanup
  unlink(custom_cache, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
})

test_that("cache_stan_model argument: works for sccomp_estimate, sccomp_remove_outliers, sccomp_boxplot", {
  skip_cmdstan()

  custom_cache <- file.path(tempdir(), "sccomp_test_cache_explicit")
  dir.create(custom_cache, showWarnings = FALSE, recursive = TRUE)

  # Pass cache_stan_model explicitly to each function (isolate output dir per test)
  output_dir <- file.path(tempdir(), "sccomp_draws_explicit")
  result <- expect_no_error({
    counts_obj |>
      sccomp_estimate(
        formula_composition = ~ type,
        formula_variability = ~ 1,
        sample = "sample",
        cell_group = "cell_group",
        abundance = "count",
        cores = 1,
        inference_method = "pathfinder",
        max_sampling_iterations = 500,
        verbose = FALSE,
        cache_stan_model = custom_cache,
        output_directory = output_dir
      ) |>
      sccomp_remove_outliers(
        percent_false_positive = 80,
        cores = 1,
        inference_method = "pathfinder",
        max_sampling_iterations = 300,
        verbose = FALSE,
        cache_stan_model = custom_cache,
        output_directory = output_dir
      ) |>
      sccomp_test()
  })

  # Verify compiled models were written to the custom cache
  sccomp_version <- as.character(packageVersion("sccomp"))
  versioned_cache <- file.path(custom_cache, sccomp_version)
  expect_true(dir.exists(versioned_cache), info = "versioned cache dir should exist")
  cached_rds <- list.files(versioned_cache, pattern = "\\.rds$")
  expect_true(length(cached_rds) >= 2L, info = "cache should contain at least 2 compiled models (estimate + generate_data)")
  expect_true("glm_multi_beta_binomial.rds" %in% cached_rds, info = "main model should be in cache")
  expect_true("glm_multi_beta_binomial_generate_data.rds" %in% cached_rds, info = "generate_data model should be in cache")

  unlink(custom_cache, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
})

test_that("sccomp_boxplot has cache_stan_model parameter", {
  expect_true("cache_stan_model" %in% names(formals(sccomp_boxplot)))
  expect_equal(formals(sccomp_boxplot)$cache_stan_model, quote(sccomp_stan_models_cache_dir))
})
