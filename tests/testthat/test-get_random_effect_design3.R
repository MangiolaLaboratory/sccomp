library(testthat)
library(tibble)
library(dplyr)
library(tidyr)
library(rlang)
library(sccomp)

test_that("get_random_effect_design3 with new NA groups", {
  # Create test input data
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    type = c("healthy", "healthy", "cancer", NA),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2")
  )
  
  formula <- ~ type
  grouping <- "group__"
  
  # Expect an error because NAs are present and accept_NA_as_average_effect is FALSE
  expect_error(
    sccomp:::get_random_effect_design3(test_data, formula, grouping, sample),
    "sccomp says: NA values are present in the design matrix factors. Set accept_NA_as_average_effect = TRUE to handle NAs as average effects across factor levels."
  )
})

test_that("get_random_effect_design3 handles NA values in design matrix correctly", {
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
    type = c("healthy", "healthy", "cancer", "cancer", NA),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2", "GROUP1")
  )
  
  formula <- ~ type
  grouping <- "group__"
  
  result <- sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = TRUE)
  
  design_matrix <- result |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  print("Actual output:")
  print(design_matrix)
  
  expected_matrix <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
    `(Intercept)___GROUP1` = c(1, 0, 1, 0, 1),
    `typehealthy___GROUP1` = c(1, 0, 0, 0, 0.5),
    `(Intercept)___GROUP2` = c(0, 1, 0, 1, 0),
    `typehealthy___GROUP2` = c(0, 1, 0, 0, 0)
  )
  
  print("Expected output:")
  print(expected_matrix)
  
  expect_equal(design_matrix, expected_matrix)
})

test_that("get_random_effect_design3 throws error when NAs present and accept_NA_as_average_effect is FALSE", {
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
    type = c("healthy", "healthy", "cancer", "cancer", NA),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2", "GROUP1")
  )
  
  formula <- ~ type
  grouping <- "group__"
  
  # Test that an error is thrown with default accept_NA_as_average_effect = FALSE
  expect_error(
    sccomp:::get_random_effect_design3(test_data, formula, grouping, sample),
    "sccomp says: NA values are present in the design matrix factors. Set accept_NA_as_average_effect = TRUE to handle NAs as average effects across factor levels."
  )
  
  # Test that an error is thrown with explicit accept_NA_as_average_effect = FALSE
  expect_error(
    sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = FALSE),
    "sccomp says: NA values are present in the design matrix factors. Set accept_NA_as_average_effect = TRUE to handle NAs as average effects across factor levels."
  )
})

test_that("get_random_effect_design3 handles complex formulas with multiple factors", {
  # Create test data with multiple factors (scaled continuous variables)
  age_raw <- c(25, 30, NA, 35, 40, 45)
  age_scaled <- as.numeric(scale(age_raw, center = TRUE, scale = TRUE))
  sex <- c("M", "F", "M", NA, "F", "M")
  type <- c("healthy", "healthy", "cancer", "cancer", NA, "healthy")
  group__ <- c("GROUP1", "GROUP2", "GROUP1", "GROUP2", "GROUP1", "GROUP2")
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6"),
    type = type,
    age = age_scaled,
    sex = sex,
    group__ = group__
  )
  
  # Test with main effects and interactions
  formula <- ~ type + age + sex + type:age
  grouping <- "group__"
  
  result <- sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = TRUE)
  
  design_matrix <- result |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  print("Actual output (complex formulas with multiple factors):")
  print(design_matrix)
  
  # Verify the structure of the output
  expect_true(all(c("(Intercept)___GROUP1", "typehealthy___GROUP1", "age___GROUP1") %in% colnames(design_matrix)))
  expect_true(all(c("(Intercept)___GROUP2", "typehealthy___GROUP2", "age___GROUP2") %in% colnames(design_matrix)))
  
  # Verify specific values
  print("age___GROUP1 for sample1:")
  print(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"])
  print("typehealthy___GROUP1 for sample1:")
  print(design_matrix$`typehealthy___GROUP1`[design_matrix$sample == "sample1"])
  print("age___GROUP1 for sample5:")
  print(design_matrix$`age___GROUP1`[design_matrix$sample == "sample5"])
  expect_equal(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"], age_scaled[1])
  expect_equal(design_matrix$`typehealthy___GROUP1`[design_matrix$sample == "sample1"], 1)
  expect_equal(design_matrix$`age___GROUP1`[design_matrix$sample == "sample5"], age_scaled[5])
})

test_that("get_random_effect_design3 handles continuous variables with NAs", {
  age_raw <- c(25, NA, 35, 40)
  age_scaled <- as.numeric(scale(age_raw, center = TRUE, scale = TRUE))
  weight_raw <- c(70, 75, NA, 80)
  weight_scaled <- as.numeric(scale(weight_raw, center = TRUE, scale = TRUE))
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4"),
    age = age_scaled,
    weight = weight_scaled,
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2")
  )
  
  formula <- ~ age + weight
  grouping <- "group__"
  
  result <- sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = TRUE)
  
  design_matrix <- result |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  print("Actual output (continuous variables with NAs):")
  print(design_matrix)
  
  # Verify that NAs are handled correctly for continuous variables
  print("age___GROUP1 for sample1:")
  print(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"])
  print("age___GROUP2 for sample2:")
  print(design_matrix$`age___GROUP2`[design_matrix$sample == "sample2"])
  print("weight___GROUP1 for sample3:")
  print(design_matrix$`weight___GROUP1`[design_matrix$sample == "sample3"])
  expect_equal(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"], age_scaled[1])
  expect_equal(design_matrix$`age___GROUP2`[design_matrix$sample == "sample2"], 0) # NA replaced with 0
  expect_equal(design_matrix$`weight___GROUP1`[design_matrix$sample == "sample3"], 0) # NA replaced with 0
})

test_that("get_random_effect_design3 handles complex interactions with NAs", {
  age_raw <- c(25, NA, 35, 40, 45)
  age_scaled <- as.numeric(scale(age_raw, center = TRUE, scale = TRUE))
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
    type = c("healthy", "healthy", "cancer", NA, "healthy"),
    age = age_scaled,
    sex = c("M", "F", NA, "M", "F"),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2", "GROUP1")
  )
  
  # Test with three-way interaction
  formula <- ~ type * age * sex
  grouping <- "group__"
  
  result <- sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = TRUE)
  
  design_matrix <- result |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  print("Actual output (complex interactions with NAs):")
  print(design_matrix)
  
  # Verify that interaction terms are present
  print("Column names:")
  print(colnames(design_matrix))
  expect_true(any(grepl("typehealthy:age:sexM", colnames(design_matrix))))
  
  # Verify specific interaction values
  print("typehealthy:age:sexM___GROUP1 for sample1:")
  print(design_matrix$`typehealthy:age:sexM___GROUP1`[design_matrix$sample == "sample1"])
  # These values are now scaled
  # Use the scaled value for age in the interaction
  expect_equal(design_matrix$`typehealthy:age:sexM___GROUP1`[design_matrix$sample == "sample1"], age_scaled[1])
})

test_that("get_random_effect_design3 handles multiple groups with complex factors", {
  age_raw <- c(25, 30, 35, 40, 45, 50)
  age_scaled <- as.numeric(scale(age_raw, center = TRUE, scale = TRUE))
  test_data <- tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6"),
    type = c("healthy", "healthy", "cancer", "cancer", NA, "healthy"),
    age = age_scaled,
    group__ = c("GROUP1", "GROUP2", "GROUP3", "GROUP1", "GROUP2", "GROUP3")
  )
  
  formula <- ~ type + age
  grouping <- "group__"
  
  result <- sccomp:::get_random_effect_design3(test_data, formula, grouping, sample, accept_NA_as_average_effect = TRUE)
  
  design_matrix <- result |> 
    select(sample, group___label, value) |>
    pivot_wider(names_from = group___label, values_from = value) |>
    mutate(across(everything(), ~ .x |> replace_na(0)))
  
  print("Actual output (multiple groups with complex factors):")
  print(design_matrix)
  
  # Verify that all groups are present
  print("Groups present:")
  print(unique(gsub(".*___", "", colnames(design_matrix))))
  expect_true(all(c("GROUP1", "GROUP2", "GROUP3") %in% unique(gsub(".*___", "", colnames(design_matrix)))))
  
  # Verify specific values for each group
  print("age___GROUP1 for sample1:")
  print(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"])
  print("age___GROUP2 for sample2:")
  print(design_matrix$`age___GROUP2`[design_matrix$sample == "sample2"])
  print("age___GROUP3 for sample3:")
  print(design_matrix$`age___GROUP3`[design_matrix$sample == "sample3"])
  expect_equal(design_matrix$`age___GROUP1`[design_matrix$sample == "sample1"], age_scaled[1])
  expect_equal(design_matrix$`age___GROUP2`[design_matrix$sample == "sample2"], age_scaled[2])
  expect_equal(design_matrix$`age___GROUP3`[design_matrix$sample == "sample3"], age_scaled[3])
})
