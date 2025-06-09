library(testthat)
library(dplyr)
library(sccomp)
library(rlang)

test_that("warnings when symbols are provided for key columns", {
  
  skip_cmdstan()
  
  # Create test data
  test_data <- tibble(
    donor = rep(c("sample1", "sample2", "sample3"), each = 3),
    cell_group = rep(c("A", "B", "C"), times = 3),
    count = c(10L, 20L, 30L, 15L, 25L, 35L, 12L, 22L, 32L),
    type = rep(c("healthy", "healthy", "cancer"), each = 3)
  )
  
  # Test with symbol for sample column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = donor,
        cell_group_column = "cell_group",
        abundance_column = "count"
      ),
    "object 'donor' not found"
  )
  
  # Test with symbol for cell_group column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = cell_group,
        abundance_column = "count"
      ),
    "object 'cell_group' not found"
  )
  
  # Test with symbol for abundance column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = "cell_group",
        abundance_column = count
      ),
    "sccomp says: abundance_column must be of character type"
  )
  
  # Test with correct character strings (should not warn)
  expect_no_warning(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = "cell_group",
        abundance_column = "count"
      )
  )
})

test_that("errors when non-character values are provided for column arguments", {
  
  skip_cmdstan()
  
  # Create test data
  test_data <- tibble(
    donor = c("sample1", "sample2", "sample3"),
    cell_group = c("A", "B", "C"),
    count = c(10L, 20L, 30L),
    type = c("healthy", "healthy", "cancer")
  )
  
  # Test with non-character sample_column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = 123,
        cell_group_column = "cell_group",
        abundance_column = "count"
      ),
    "sccomp says: sample_column must be of character type"
  )
  
  # Test with non-character cell_group_column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = 456,
        abundance_column = "count"
      ),
    "sccomp says: cell_group_column must be of character type"
  )
  
  # Test with non-character abundance_column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = "cell_group",
        abundance_column = 789
      ),
    "sccomp says: abundance_column must be of character type"
  )
})

test_that("warnings when using deprecated column names", {
  
  skip_cmdstan()
  
  # Create test data
  test_data <- tibble(
    donor = rep(c("sample1", "sample2", "sample3"), each = 3),
    cell_group = rep(c("A", "B", "C"), times = 3),
    count = c(10L, 20L, 30L, 15L, 25L, 35L, 12L, 22L, 32L),
    type = rep(c("healthy", "healthy", "cancer"), each = 3)
  )
  
  # Test with deprecated .sample
  expect_warning(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        .sample = "donor",
        cell_group_column = "cell_group",
        abundance_column = "count"
      ),
    "The `.sample` argument of*"
  )
  
  # Test with deprecated .cell_group
  expect_warning(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        .cell_group = "cell_group",
        abundance_column = "count"
      ),
    ".cell_group argument.*have been deprecated in favour of cell_group_column"
  )
  
  # Test with deprecated .abundance
  expect_warning(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "donor",
        cell_group_column = "cell_group",
        .abundance = "count"
      ),
    ".abundance argument.*have been deprecated in favour of abundance_column"
  )
})


