# Unit tests for sccomp_simulate function
# Testing simulation functionality with various model configurations

library(testthat)
library(sccomp)
library(dplyr)

# Helper function to skip tests if cmdstan is not available
skip_cmdstan <- function() {
  if (!instantiate::stan_cmdstan_exists()) {
    skip("CmdStan not available")
  }
}

test_that("sccomp_simulate works with simple model", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table (all zeros for simplicity)
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate data
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that result has expected columns
  expect_true("sample" %in% colnames(result))
  expect_true("cell_group" %in% colnames(result))
  
  # Check that result has data
  expect_true(nrow(result) > 0)
})

test_that("sccomp_simulate works with variability formula", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model with variability formula
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~type, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate data with variability formula
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~type,  # Use variability formula
    new_data = new_data_samples,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that result has expected columns
  expect_true("sample" %in% colnames(result))
  expect_true("cell_group" %in% colnames(result))
})

test_that("sccomp_simulate validates dimensions correctly", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table with more cell groups than original (should fail)
  coeffs_table_expanded = counts_obj |>
    bind_rows(
      counts_obj |>
        mutate(cell_group = paste0(cell_group, "_new"))
    ) |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # This should fail with dimension validation error
  expect_error(
    sccomp_simulate(
      estimate,
      ~type, 
      ~1, 
      new_data = new_data_samples,
      coefficients = coeffs_table_expanded,
      cores = 1
    ),
    "M_simulated.*cannot be larger than M"
  )
})

test_that("sccomp_simulate works with number_of_draws parameter", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate with multiple draws
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    number_of_draws = 3,
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that we have data (may have replicate column if multiple draws)
  expect_true(nrow(result) > 0)
})

test_that("sccomp_simulate works with variability_multiplier", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate with different variability multipliers
  result1 = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    variability_multiplier = 1,
    cores = 1
  )
  
  result2 = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    variability_multiplier = 10,
    cores = 1
  )
  
  # Both should work
  expect_s3_class(result1, "tbl")
  expect_s3_class(result2, "tbl")
})

test_that("sccomp_simulate handles missing formula_variability gracefully", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model with variability formula
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~type, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate without specifying formula_variability (should use from estimate)
  result = sccomp_simulate(
    estimate,
    ~type, 
    NULL,  # Let it use the formula from estimate
    new_data = new_data_samples,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Should work
  expect_s3_class(result, "tbl")
})

test_that("sccomp_simulate works with subset of original data", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Use subset of data (fewer samples)
  new_data_samples_subset = counts_obj |>
    filter(sample %in% unique(sample)[1:3]) |>
    distinct(sample, type)
  
  # Simulate data
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples_subset,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Should work with subset
  expect_s3_class(result, "tbl")
  expect_true(nrow(result) >= nrow(new_data_samples_subset))
})

test_that("sccomp_simulate works with separate coefficients table", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create separate coefficients table (cell-type specific)
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create sample-specific new_data (no cell_group needed)
  new_data_samples = counts_obj |>
    distinct(sample, type)
  
  # Simulate data with separate coefficients table
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Should work
  expect_s3_class(result, "tbl")
  expect_true("sample" %in% colnames(result))
  expect_true("cell_group" %in% colnames(result))
  expect_true(nrow(result) > 0)
  
  # Check that all cell_groups from coefficients table are present
  expect_true(all(coeffs_table$cell_group %in% result$cell_group))
})

test_that("sccomp_simulate works with separate coefficients table and subset data", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create separate coefficients table
  coeffs_table = counts_obj |>
    distinct(cell_group) |>
    mutate(`(Intercept)` = 0, `typecancer` = 0)
  
  # Create subset of samples
  new_data_samples = counts_obj |>
    filter(sample %in% unique(sample)[1:3]) |>
    distinct(sample, type)
  
  # Simulate data
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = new_data_samples,
    coefficients = coeffs_table,
    cores = 1
  )
  
  # Should work
  expect_s3_class(result, "tbl")
  expect_true(nrow(result) > 0)
})

test_that("sccomp_simulate uses posterior beta when coefficients are NULL", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Simulate without providing coefficients (should use posterior beta_raw)
  result = sccomp_simulate(
    estimate,
    ~type, 
    ~1, 
    new_data = NULL,  # Use original data
    coefficients = NULL,  # No coefficients provided
    cores = 1
  )
  
  # Should work
  expect_s3_class(result, "tbl")
  expect_true(nrow(result) > 0)
})

