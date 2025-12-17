library(dplyr)
library(tidyr)
library(sccomp)

test_that("incorporate_parameters_into_fit_object loads all parameters", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Create a test output directory
  test_output_dir <- tempfile("sccomp_test_draws_")
  dir.create(test_output_dir)
  
  # Run sccomp_estimate with cleanup_draw_files = FALSE to keep CSV files temporarily
  result <- counts_obj |>
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
      output_directory = test_output_dir,
      cleanup_draw_files = FALSE  # Keep files to test manual incorporation
    )
  
  # Get the fit object
  fit <- attr(result, "fit")
  
  # Check that CSV files exist
  csv_files <- list.files(test_output_dir, pattern = "\\.csv$", full.names = TRUE)
  expect_true(length(csv_files) > 0, info = "CSV files should exist before cleanup")
  
  # Call the function to incorporate parameters
  expect_no_error({
    sccomp:::incorporate_parameters_into_fit_object(fit)
  })
  
  # Verify that key parameters can be accessed after incorporation
  expect_no_error({
    beta_draws <- fit$draws(variables = "beta", format = "draws_df")
  })
  
  expect_no_error({
    alpha_draws <- fit$draws(variables = "alpha", format = "draws_df")
  })
  
  expect_no_error({
    prec_coeff_draws <- fit$draws(variables = "prec_coeff", format = "draws_df")
  })
  
  # Now delete the CSV files to simulate cleanup
  file.remove(csv_files)
  
  # Verify CSV files are gone
  csv_files_after <- list.files(test_output_dir, pattern = "\\.csv$", full.names = TRUE)
  expect_equal(length(csv_files_after), 0, info = "CSV files should be deleted")
  
  # Parameters should still be accessible because they were incorporated
  expect_no_error({
    beta_draws_after <- fit$draws(variables = "beta", format = "draws_df")
  })
  
  expect_no_error({
    alpha_draws_after <- fit$draws(variables = "alpha", format = "draws_df")
  })
  
  # Verify the draws are the same before and after CSV deletion
  expect_equal(beta_draws, beta_draws_after)
  
  # Clean up test directory
  unlink(test_output_dir, recursive = TRUE)
})

test_that("incorporate_parameters_into_fit_object handles models without random effects", {
  skip_cmdstan()
  
  data("counts_obj")
  
  test_output_dir <- tempfile("sccomp_test_draws_")
  dir.create(test_output_dir)
  
  # Run a simple model without random effects
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
      cleanup_draw_files = FALSE
    )
  
  fit <- attr(result, "fit")
  
  # Should handle models without random effects gracefully
  expect_no_error({
    sccomp:::incorporate_parameters_into_fit_object(fit)
  })
  
  # Basic parameters should still be accessible
  expect_no_error({
    fit$draws(variables = "beta", format = "draws_df")
  })
  
  expect_no_error({
    fit$draws(variables = "alpha", format = "draws_df")
  })
  
  # Clean up
  unlink(test_output_dir, recursive = TRUE)
})
