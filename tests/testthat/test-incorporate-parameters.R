library(dplyr)
library(tidyr)
library(sccomp)

# Shared fast estimate for draw-file / portability checks (one factor, pathfinder).
estimate_for_draw_tests <- function(output_directory, portable) {
  data("counts_obj", package = "sccomp", envir = environment())
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
      output_directory = output_directory,
      portable = portable
    )
}

test_that("non-portable estimate: deleting Stan output files before sccomp_test() errors without incorporation", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_draws_no_incorp_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_draw_tests(test_output_dir, portable = FALSE)
  fit <- attr(result, "fit")

  paths <- fit$output_files(include_failed = TRUE)
  paths <- paths[file.exists(paths)]
  expect_gt(length(paths), 0L, label = "Stan output files on disk")

  ok <- file.remove(paths)
  expect_true(all(ok), label = "removing Stan chain/output files")
  expect_false(any(file.exists(paths)), label = "recorded Stan paths should not exist after deletion")

  expect_error(
    sccomp_test(result),
    "Stan output files for this fit are not on disk"
  )
})

test_that("incorporate_parameters_into_fit_object keeps draws after CSV deletion", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_draws_fit_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_draw_tests(test_output_dir, portable = FALSE)
  fit <- attr(result, "fit")

  csv_files <- list.files(test_output_dir, pattern = "\\.csv$", full.names = TRUE)
  expect_gt(length(csv_files), 0L)

  sccomp:::incorporate_parameters_into_fit_object(fit)

  expect_no_error({
    beta_draws <- fit$draws(variables = "beta", format = "draws_df")
  })
  expect_no_error({
    fit$draws(variables = "alpha", format = "draws_df")
  })
  
  expect_no_error({
    prec_intercept_draws <- fit$draws(variables = "prec_intercept_1", format = "draws_df")
  })
  
  # Now delete the CSV files to simulate cleanup

  file.remove(csv_files)

  expect_no_error({
    beta_after <- fit$draws(variables = "beta", format = "draws_df")
  })
  expect_equal(beta_draws, beta_after)
})

test_that("incorporate_parameters_into_sccomp_object keeps draws after CSV deletion", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_draws_sccomp_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_draw_tests(test_output_dir, portable = FALSE)
  result <- sccomp:::incorporate_parameters_into_sccomp_object(result)

  csv_files <- list.files(test_output_dir, pattern = "\\.csv$", full.names = TRUE)
  expect_gt(length(csv_files), 0L)

  expect_no_error({
    beta_draws <- attr(result, "fit")$draws(variables = "beta", format = "draws_df")
  })

  file.remove(csv_files)

  expect_no_error({
    beta_after <- attr(result, "fit")$draws(variables = "beta", format = "draws_df")
  })
  expect_equal(beta_draws, beta_after)
})

test_that("portable = TRUE keeps draws available after package removes draw CSV files", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_draws_portable_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_draw_tests(test_output_dir, portable = TRUE)
  fit <- attr(result, "fit")

  remaining <- fit$output_files(include_failed = TRUE)
  remaining <- remaining[file.exists(remaining)]
  expect_equal(length(remaining), 0L, label = "Stan output files still on disk after portable cleanup")

  expect_no_error({
    fit$draws(variables = "beta", format = "draws_df")
  })
})

test_that("incorporate_parameters_into_fit_object handles models without random effects", {
  skip_cmdstan()

  data("counts_obj")

  test_output_dir <- tempfile("sccomp_test_draws_intercept_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- counts_obj |>
    sccomp_estimate(
      formula_composition = ~ 1,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE,
      output_directory = test_output_dir,
      portable = FALSE
    )

  fit <- attr(result, "fit")

  expect_no_error({
    sccomp:::incorporate_parameters_into_fit_object(fit)
  })

  expect_no_error({
    fit$draws(variables = "beta", format = "draws_df")
  })
  expect_no_error({
    fit$draws(variables = "alpha", format = "draws_df")
  })
})

test_that("incorporate_parameters_into_sccomp_object errors without fit attribute", {
  bad <- tibble::tibble(x = 1L)
  expect_error(
    sccomp:::incorporate_parameters_into_sccomp_object(bad),
    "expected a \"fit\" attribute on the sccomp object"
  )
})

test_that("incorporate_parameters_into_sccomp_object forwards fit and writes back attribute", {
  fit_in <- list(seed = 42L)
  obj <- tibble::tibble(x = 1L)
  attr(obj, "fit") <- fit_in
  class(obj) <- c("sccomp_tbl", class(obj))

  local_mocked_bindings(
    incorporate_parameters_into_fit_object = function(fit) {
      expect_identical(fit, fit_in)
      fit$incorporated <- TRUE
      fit
    },
    .package = "sccomp"
  )

  out <- sccomp:::incorporate_parameters_into_sccomp_object(obj)

  expect_s3_class(out, "sccomp_tbl")
  expect_true(attr(out, "fit")$incorporated)
  expect_identical(attr(out, "fit")$seed, 42L)
})
