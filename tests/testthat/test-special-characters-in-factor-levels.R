library(testthat)
library(sccomp)
library(dplyr)

test_that("sccomp handles factor levels with hyphens correctly", {
  skip_cmdstan()
  
  # Load the included dataset
  data("counts_obj")
  
  # Create test data with factor levels containing hyphens
  # Based on counts_obj but with renamed factor levels that include hyphens
  test_data <- counts_obj |>
    # First convert to character to avoid factor issues
    mutate(type_char = as.character(type)) |>
    mutate(
      # Rename the 'type' factor levels to include hyphens
      type_with_hyphen = case_when(
        type_char == "benign" ~ "pre-treatment",
        type_char == "cancer" ~ "post-treatment",
        TRUE ~ type_char
      )
    ) |>
    # Keep original count column and new factor
    select(sample, cell_group, count, type_with_hyphen)
  
  # Verify we have hyphens in our factor levels
  expect_true(any(grepl("-", unique(test_data$type_with_hyphen))))
  
  # Run sccomp_estimate with the hyphenated factor
  sccomp_result <- test_data |>
    sccomp_estimate(
      formula_composition = ~ type_with_hyphen,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Basic checks
  expect_s3_class(sccomp_result, "sccomp_tbl")
  expect_true("cell_group" %in% colnames(sccomp_result))
  expect_true("c_effect" %in% colnames(sccomp_result))
  
  # Test sccomp_proportional_fold_change with hyphenated factor levels
  fold_changes <- sccomp_proportional_fold_change(
    sccomp_result,
    formula_composition = ~ type_with_hyphen,
    from = "pre-treatment",
    to = "post-treatment"
  )
  
  # Verify the results structure
  expect_s3_class(fold_changes, "tbl_df")
  expect_true("cell_group" %in% colnames(fold_changes))
  expect_true("proportion_fold_change" %in% colnames(fold_changes))
  expect_true("average_uncertainty" %in% colnames(fold_changes))
  expect_true("statement" %in% colnames(fold_changes))
  
  # Verify we got results for multiple cell groups
  expect_gt(nrow(fold_changes), 0)
  
  # Verify fold changes are numeric and not all identical
  expect_true(is.numeric(fold_changes$proportion_fold_change))
  expect_false(all(is.na(fold_changes$proportion_fold_change)))
  
  # Key test: fold changes should not all be trivial (1.0 or -1.0)
  # With real biological differences, we expect variation
  non_trivial <- fold_changes |>
    filter(abs(abs(proportion_fold_change) - 1) > 0.05)
  
  expect_gt(nrow(non_trivial), 0,
            label = "Number of non-trivial fold changes should be > 0")
  
  # Verify that the fold changes are finite
  expect_true(all(is.finite(fold_changes$proportion_fold_change)))
})

test_that("sccomp handles factor levels with multiple special characters", {
  skip_cmdstan()
  
  data("counts_obj")
  
  # Create test data with various special characters in factor levels
  test_data <- counts_obj |>
    mutate(type_char = as.character(type)) |>
    mutate(
      # Test multiple types of special characters that might appear in real data
      complex_type = case_when(
        type_char == "benign" ~ "Day-0_baseline",
        type_char == "cancer" ~ "Day-7_post.treatment",
        TRUE ~ type_char
      )
    ) |>
    select(sample, cell_group, count, complex_type)
  
  # Run sccomp_estimate
  sccomp_result <- test_data |>
    sccomp_estimate(
      formula_composition = ~ complex_type,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Test proportional fold change
  fold_changes <- sccomp_proportional_fold_change(
    sccomp_result,
    formula_composition = ~ complex_type,
    from = "Day-0_baseline",
    to = "Day-7_post.treatment"
  )
  
  # Verify results
  expect_s3_class(fold_changes, "tbl_df")
  expect_gt(nrow(fold_changes), 0)
  expect_true(all(is.finite(fold_changes$proportion_fold_change)))
})

test_that("sccomp_predict handles factor levels with hyphens", {
  skip_cmdstan()
  
  data("counts_obj")
  
  # Create simple test data
  test_data <- counts_obj |>
    mutate(type_char = as.character(type)) |>
    mutate(
      type_hyphen = if_else(type_char == "benign", "pre-treatment", "post-treatment")
    ) |>
    select(sample, cell_group, count, type_hyphen)
  
  # Estimate
  sccomp_result <- test_data |>
    sccomp_estimate(
      formula_composition = ~ type_hyphen,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE
    )
  
  # Create new_data for prediction with hyphenated factor levels
  new_data <- tibble(
    sample = c("test_sample_1", "test_sample_2"),
    type_hyphen = c("pre-treatment", "post-treatment")
  )
  
  # Test that prediction works with hyphenated factor levels
  predictions <- sccomp_predict(
    sccomp_result,
    formula_composition = ~ type_hyphen,
    new_data = new_data
  )
  
  # Verify predictions
  expect_s3_class(predictions, "tbl_df")
  expect_true("proportion_mean" %in% colnames(predictions))
  expect_true(all(is.finite(predictions$proportion_mean)))
  expect_true(all(predictions$proportion_mean >= 0))
  expect_true(all(predictions$proportion_mean <= 1))
})

test_that("factor levels with hyphens match correctly in design matrix", {
  skip_cmdstan()
  
  data("counts_obj")
  
  # Simplified test focusing on factor level matching
  test_data <- counts_obj |>
    mutate(type_char = as.character(type)) |>
    mutate(treatment = if_else(type_char == "benign", "pre-surgery", "post-surgery")) |>
    select(sample, cell_group, count, treatment)
  
  # Estimate
  result <- test_data |>
    sccomp_estimate(
      formula_composition = ~ treatment,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Check that we have treatment-related parameters
  # The parameter names will include the treatment factor
  has_treatment_param <- any(grepl("treatment", result$parameter, fixed = TRUE))
  expect_true(has_treatment_param)
  
  # Test fold change calculation
  fc <- sccomp_proportional_fold_change(
    result,
    formula_composition = ~ treatment,
    from = "pre-surgery",
    to = "post-surgery"
  )
  
  # Verify output
  expect_s3_class(fc, "tbl_df")
  expect_true(all(!is.na(fc$proportion_fold_change)))
  expect_true(all(is.finite(fc$proportion_fold_change)))
  
  # Check that the statement includes meaningful values
  expect_true(all(nchar(fc$statement) > 10))
})

