library(testthat)
library(tibble)
library(dplyr)
library(rlang)
library(sccomp)

test_that("get_random_effect_design2 creates correct design matrix", {
  # Create test input data
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    type = c("healthy", "healthy", "cancer", "cancer"),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2")  # Now we have both types in both groups
  )
  
  # Create formula
  formula <- ~ type + (type|group__)
  
  # Run function using internal function
  result <- sccomp:::get_random_effect_design2(test_data, sample, formula)
  
  # Extract and format the design matrix
  design_matrix <- result |> 
    pull(design) |> 
    _[[1]] |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  # Print actual output for debugging
  print("Actual output:")
  print(design_matrix)
  
  # Expected output
  expected_matrix <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    `(Intercept)___GROUP1` = c(1, 0, 1, 0),
    `typehealthy___GROUP1` = c(1, 0, 0, 0),
    `(Intercept)___GROUP2` = c(0, 1, 0, 1),
    `typehealthy___GROUP2` = c(0, 1, 0, 0)
  )
  
  print("Expected output:")
  print(expected_matrix)
  
  # Compare results
  expect_equal(design_matrix, expected_matrix)
}) 