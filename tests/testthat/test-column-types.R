library(testthat)
library(dplyr)
library(sccomp)
library(rlang)

test_that("warnings when symbols are provided for key columns", {
  # Create test data
  test_data <- tibble(
    sample = rep(c("sample1", "sample2", "sample3"), each = 3),
    cell_group = rep(c("A", "B", "C"), times = 3),
    count = c(10L, 20L, 30L, 15L, 25L, 35L, 12L, 22L, 32L),
    type = rep(c("healthy", "healthy", "cancer"), each = 3)
  )
  
  # Test with symbol for sample column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = sample,
        cell_group_column = "cell_group",
        abundance_column = "count"
      ),
    "sccomp says: sample_column must be of character type"
  )
  
  # Test with symbol for cell_group column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "sample",
        cell_group_column = cell_group,
        abundance_column = "count"
      ),
    "sccomp says: cell_group_column must be of character type"
  )
  
  # Test with symbol for abundance column
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "sample",
        cell_group_column = "cell_group",
        abundance_column = count
      ),
    "sccomp says: abundance_column must be of character type"
  )
  
  # Test with all symbols
  expect_error(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = sample,
        cell_group_column = cell_group,
        abundance_column = count
      ),
    "sccomp says: sample_column must be of character type"
  )
  
  # Test with correct character strings (should not warn)
  expect_no_warning(
    test_data |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample_column = "sample",
        cell_group_column = "cell_group",
        abundance_column = "count"
      )
  )
})

test_that("errors when non-character values are provided for column arguments", {
  # Create test data
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3"),
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
        sample_column = "sample",
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
        sample_column = "sample",
        cell_group_column = "cell_group",
        abundance_column = 789
      ),
    "sccomp says: abundance_column must be of character type"
  )
})

test_that("warnings when deprecated arguments are used", {
  test_data <- tibble(
    sample = c("s1", "s2"),
    cell_group = c("cg1", "cg2"),
    count = c(10, 20),
    type = c("A", "B")
  )
  
  # Test deprecated .sample argument
  expect_warning(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   .sample = sample,
                   cell_group_column = "cell_group",
                   abundance_column = "count"),
    "The argument '.sample' is deprecated"
  )
  
  # Test deprecated .cell_group argument
  expect_warning(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   sample_column = "sample",
                   .cell_group = cell_group,
                   abundance_column = "count"),
    "The argument '.cell_group' is deprecated"
  )
  
  # Test deprecated .abundance argument
  expect_warning(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   sample_column = "sample",
                   cell_group_column = "cell_group",
                   .abundance = count),
    "The argument '.abundance' is deprecated"
  )
})

test_that("errors when *_column arguments are provided with symbols", {
  test_data <- tibble(
    sample = c("s1", "s2"),
    cell_group = c("cg1", "cg2"),
    count = c(10, 20),
    type = c("A", "B")
  )
  
  # Test sample_column with symbol
  expect_error(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   sample_column = sample,
                   cell_group_column = "cell_group",
                   abundance_column = "count"),
    "sccomp says: sample_column must be of character type"
  )
  
  # Test cell_group_column with symbol
  expect_error(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   sample_column = "sample",
                   cell_group_column = cell_group,
                   abundance_column = "count"),
    "sccomp says: cell_group_column must be of character type"
  )
  
  # Test abundance_column with symbol
  expect_error(
    sccomp_estimate(test_data, 
                   formula_composition = ~type,
                   sample_column = "sample",
                   cell_group_column = "cell_group",
                   abundance_column = count),
    "sccomp says: abundance_column must be of character type"
  )
})

test_that("sample_column is a quosure with an expression that includes a double quote", {
  test_data <- tibble(
    sample = c("s1", "s2"),
    cell_group = c("cg1", "cg2"),
    count = c(10, 20),
    type = c("A", "B")
  )
  
  # Capture the quosure
  sample_quo <- rlang::enquo(sample)
  
  # Check if the expression is a symbol
  expect_true(rlang::is_symbol(rlang::quo_get_expr(sample_quo)))
  
  # Check if the symbol name is "sample"
  expect_identical(rlang::as_string(rlang::quo_get_expr(sample_quo)), "sample")
}) 