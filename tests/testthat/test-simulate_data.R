# Unit tests for simulate_data function
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

test_that("simulate_data works with simple model", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Set coefficients for cell_groups (all zeros for simplicity)
  counts_obj_with_coefs = counts_obj |> 
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate data
  result = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    ~1, 
    sample, 
    cell_group, 
    c(b_0, b_1),
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that result has expected columns
  expect_true("sample" %in% colnames(result))
  expect_true("cell_group" %in% colnames(result))
  
  # Check that result has same number of rows as input (or more if multiple draws)
  expect_true(nrow(result) >= nrow(counts_obj_with_coefs))
})

test_that("simulate_data works with variability formula", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model with variability formula
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~type, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Set coefficients
  counts_obj_with_coefs = counts_obj |> 
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate data with variability formula
  result = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    ~type,  # Use variability formula
    sample, 
    cell_group, 
    c(b_0, b_1),
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that result has expected columns
  expect_true("sample" %in% colnames(result))
  expect_true("cell_group" %in% colnames(result))
})

test_that("simulate_data validates dimensions correctly", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Create data with more cell groups than original (should fail)
  counts_obj_expanded = counts_obj |>
    bind_rows(
      counts_obj |>
        mutate(cell_group = paste0(cell_group, "_new"))
    ) |>
    mutate(b_0 = 0, b_1 = 0)
  
  # This should fail with dimension validation error
  expect_error(
    simulate_data(
      counts_obj_expanded, 
      estimate, 
      ~type, 
      ~1, 
      sample, 
      cell_group, 
      c(b_0, b_1),
      cores = 1
    ),
    "M_simulated.*cannot be larger than M"
  )
})

test_that("simulate_data works with number_of_draws parameter", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Set coefficients
  counts_obj_with_coefs = counts_obj |> 
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate with multiple draws
  result = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    ~1, 
    sample, 
    cell_group, 
    c(b_0, b_1),
    number_of_draws = 3,
    cores = 1
  )
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl")
  
  # Check that we have data (may have replicate column if multiple draws)
  expect_true(nrow(result) > 0)
})

test_that("simulate_data works with variability_multiplier", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a simple model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Set coefficients
  counts_obj_with_coefs = counts_obj |> 
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate with different variability multipliers
  result1 = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    ~1, 
    sample, 
    cell_group, 
    c(b_0, b_1),
    variability_multiplier = 1,
    cores = 1
  )
  
  result2 = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    ~1, 
    sample, 
    cell_group, 
    c(b_0, b_1),
    variability_multiplier = 10,
    cores = 1
  )
  
  # Both should work
  expect_s3_class(result1, "tbl")
  expect_s3_class(result2, "tbl")
})

test_that("simulate_data handles missing formula_variability gracefully", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model with variability formula
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~type, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Set coefficients
  counts_obj_with_coefs = counts_obj |> 
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate without specifying formula_variability (should use from estimate)
  result = simulate_data(
    counts_obj_with_coefs, 
    estimate, 
    ~type, 
    NULL,  # Let it use the formula from estimate
    sample, 
    cell_group, 
    c(b_0, b_1),
    cores = 1
  )
  
  # Should work
  expect_s3_class(result, "tbl")
})

test_that("simulate_data works with subset of original data", {
  skip_cmdstan()
  
  # Load test data
  data("counts_obj")
  
  # Fit a model
  estimate = sccomp_estimate(
    counts_obj,
    ~ type, ~1, "sample", "cell_group", "count",
    cores = 1
  )
  
  # Use subset of data (fewer samples)
  counts_obj_subset = counts_obj |>
    filter(sample %in% unique(sample)[1:3]) |>
    mutate(b_0 = 0, b_1 = 0)
  
  # Simulate data
  result = simulate_data(
    counts_obj_subset, 
    estimate, 
    ~type, 
    ~1, 
    sample, 
    cell_group, 
    c(b_0, b_1),
    cores = 1
  )
  
  # Should work with subset
  expect_s3_class(result, "tbl")
  expect_true(nrow(result) >= nrow(counts_obj_subset))
})

