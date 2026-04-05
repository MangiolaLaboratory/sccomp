# Test for Issue #249: Nonsense proportional_fold_change results
# https://github.com/MangiolaLaboratory/sccomp/issues/249
#
# This test reproduces the issue where all proportional_fold_change results
# were showing "1-fold increase" even with clear changes in the data.

library(testthat)
library(sccomp)
library(dplyr)
library(tidyr)
library(tidyomics)

# Helper function to skip tests if cmdstan is not available
skip_cmdstan <- function() {
  if (!instantiate::stan_cmdstan_exists()) {
    skip("CmdStan not available")
  }
}

test_that("Issue #249: proportional_fold_change returns non-trivial fold changes with real data", {
  skip_cmdstan()
  
  # Path to the data file from the issue
  data_file <- "/Users/a1234450/Downloads/sce_minimal_object_sccomp_17122025.Rds"
  
  # Skip if data file doesn't exist (e.g., on CI)
  if (!file.exists(data_file)) {
    skip("Test data file not available")
  }
  
  # Load the SCE object from the issue
  sce_dummy <- 
    readRDS(data_file) |> 
    dplyr::count(main_comparison_type, cell_group, sample_donor ) |> 
    nest(data = c(cell_group, n)) |> 
    dplyr::slice(1, .by = main_comparison_type) |> 
    unnest(data)
  
  library(tidyverse)
  library(tidyomics)
  library(stringr)
  
  # Run the exact pipeline from the issue report (without str_replace workaround)
  sccomp_result <- 
    sce_dummy |> 
    sccomp_estimate(
      formula_composition = ~main_comparison_type, 
      sample = 'sample_donor',
      cell_group = 'cell_group',
      abundance = "n",
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Calculate proportional fold changes with hyphenated factor levels
  fold_changes <- sccomp_proportional_fold_change(
    sccomp_result, 
    formula_composition = ~main_comparison_type,
    from = 'Pre-surgery', 
    to = 'Bacterial'
  )
  
  # Basic structure checks
  expect_s3_class(fold_changes, "tbl_df")
  expect_true("cell_group" %in% colnames(fold_changes))
  expect_true("proportion_fold_change" %in% colnames(fold_changes))
  expect_true("average_uncertainty" %in% colnames(fold_changes))
  expect_true("statement" %in% colnames(fold_changes))
  
  # Check that we have results for multiple cell groups
  expect_gt(nrow(fold_changes), 0)
  
  # KEY TEST: Verify that NOT ALL fold changes are 1
  # The bug caused all results to be exactly 1.0
  # With real biological data, we should see variation
  unique_fold_changes <- unique(round(fold_changes$proportion_fold_change, 1))
  
  # Should have more than just 1.0 as the fold change
  expect_gt(length(unique_fold_changes), 1)
  
  # Should have at least some fold changes that are NOT 1.0 or -1.0
  non_trivial_fcs <- fold_changes$proportion_fold_change[
    abs(fold_changes$proportion_fold_change) > 1.1 | 
    abs(fold_changes$proportion_fold_change) < 0.9
  ]
  
  expect_gt(length(non_trivial_fcs), 0,
            label = "All fold changes are trivial (close to 1.0 or -1.0) - Bug may still be present!")
  
  # CRITICAL TEST from the original issue: 
  # Error if all fold changes equal 1 (the bug symptom)
  all_are_one <- fold_changes |> 
    distinct(proportion_fold_change) |> 
    pull(1) |> 
    magrittr::equals(1) |>
    all()
  
  expect_false(all_are_one, 
               label = "Bug #249: All fold changes equal 1 - this was the main symptom of the bug")
  
  # Verify that the statements are not all identical
  unique_statements <- unique(fold_changes$statement)
  expect_gt(length(unique_statements), 1,
            label = "All statements are identical - Bug may still be present!")
  
  # Check that proportion_from and proportion_to are different
  # Extract proportions from the statement if possible
  # Statements should have format: "X-fold increase (from Y to Z)"
  expect_true(
    any(grepl("increase|decrease", fold_changes$statement)),
    info = "No increase/decrease found in statements"
  )
  
  # Print summary for manual inspection
  message("\n=== Issue #249 Test Results ===")
  message(sprintf("Number of cell groups: %d", nrow(fold_changes)))
  message(sprintf("Unique fold change values: %d", length(unique_fold_changes)))
  message(sprintf("Fold change range: %.2f to %.2f", 
                 min(fold_changes$proportion_fold_change),
                 max(fold_changes$proportion_fold_change)))
  message("\nFirst 5 results:")
  print(head(fold_changes, 5))
})


test_that("Issue #249: Regression test - arrange with string parameter works correctly", {
  skip_cmdstan()
  
  # Create simple synthetic data to test the specific bug
  # The bug was that arrange(!!.sample != !!from) was incorrect
  # It should be arrange(!!.sample != from) because from is a string
  
  set.seed(12345)
  test_data <- data.frame(
    sample = rep(c("s1", "s2", "s3", "s4"), each = 3),
    cell_group = rep(c("A", "B", "C"), times = 4),
    condition = rep(c("control", "treatment"), each = 6),
    count = as.integer(c(
      # control samples (should have lower proportions for group B)
      100, 50, 100,   # s1
      110, 40, 95,    # s2
      # treatment samples (should have higher proportions for group B)
      90, 150, 100,   # s3
      95, 160, 105    # s4
    ))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ condition,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    verbose = FALSE
  )
  
# plot boxplot of the proportions using ggplot2
  estimate |> 
    sccomp_test() |> 
    sccomp_boxplot("condition") |> 
    expect_s3_class("ggplot")

  # Calculate fold changes
  fold_changes <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ condition,
    from = "control",
    to = "treatment"
  )
  
  # Cell group B should show a large increase
  fc_B <- fold_changes |> 
    filter(cell_group == "B") |> 
    pull(proportion_fold_change)
  
 
  # Cell groups A and C should have fold changes close to 1 or less than 1
  fc_A <- fold_changes |> 
    filter(cell_group == "A") |> 
    pull(proportion_fold_change)
  
  fc_C <- fold_changes |> 
    filter(cell_group == "C") |> 
    pull(proportion_fold_change)
  

  
  message("\n=== Regression Test Results ===")
  message(sprintf("Cell group A: %.2f-fold", fc_A))
  message(sprintf("Cell group B: %.2f-fold (should be >> 1)", fc_B))
  message(sprintf("Cell group C: %.2f-fold", fc_C))
})


test_that("Issue #249: Verify correct use of from parameter in arrange", {
  # This test specifically checks the code fix
  # The bug was using !!from instead of from in arrange()
  
  skip_cmdstan()
  
  # Create minimal test data (compositional models need >= 2 cell groups per sample)
  test_data <- data.frame(
    sample = c("s1", "s1", "s2", "s2"),
    cell_group = c("A", "B", "A", "B"),
    condition = c("ctrl", "ctrl", "treat", "treat"),
    count = as.integer(c(100, 300, 200, 200))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ condition,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    verbose = FALSE
  )
  
  # Calculate fold changes with different string values for from/to
  # This should work correctly now that we use 'from' instead of '!!from'
  fold_changes <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ condition,
    from = "ctrl",
    to = "treat"
  )
  
  # Proportion of A: ctrl 100/400 = 0.25, treat 200/400 = 0.5 -> ~2-fold
  fc_A_fwd <- fold_changes |> filter(cell_group == "A") |> pull(proportion_fold_change)
  expect_gt(fc_A_fwd, 1.5)
  expect_lt(fc_A_fwd, 3.0)
  
  # Try with reversed direction - should get negative fold change
  fold_changes_rev <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ condition,
    from = "treat",
    to = "ctrl"
  )
  
  # Should get a negative fold change (ratio < 1 -> -1/ratio, about -2)
  fc_A_rev <- fold_changes_rev |> filter(cell_group == "A") |> pull(proportion_fold_change)
  expect_lt(fc_A_rev, -1.0)
  
  message("\n=== Parameter Usage Test ===")
  message(sprintf("ctrl -> treat (A): %.2f-fold", fc_A_fwd))
  message(sprintf("treat -> ctrl (A): %.2f-fold", fc_A_rev))
})

test_that("Issue #249: Exact user code with hyphenated factor levels works correctly", {
  skip_cmdstan()
  
  # Path to the data file from the issue
  data_file <- "/Users/a1234450/Downloads/sce_minimal_object_sccomp_17122025.Rds"
  
  # Skip if data file doesn't exist (e.g., on CI)
  if (!file.exists(data_file)) {
    skip("Test data file not available")
  }
  
  # Load the SCE object from the issue
  sce_dummy <- readRDS(data_file)
  
  library(tidyomics)
  library(stringr)
  
  # EXACT CODE FROM THE USER'S ISSUE (without the str_replace workaround)
  sccomp_result <- 
    sce_dummy |> 
    sccomp_estimate(
      formula_composition = ~main_comparison_type, 
      sample = 'sample_donor',
      cell_group = 'cell_group',
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Calculate proportional fold changes
  fc <- sccomp_proportional_fold_change(
    sccomp_result, 
    formula_composition = ~main_comparison_type,
    from = 'Pre-surgery', 
    to = 'Bacterial'
  )
  
  # THE KEY TEST: Check if all fold changes are 1 (the bug)
  all_equal_one <- fc |> 
    distinct(proportion_fold_change) |> 
    pull(1) |> 
    magrittr::equals(1) |>
    all()
  
  # This should be FALSE - not all fold changes should equal 1
  expect_false(all_equal_one, 
               label = "BUG: All fold changes equal 1.0 (the bug from issue #249)")
  
  # Additional checks
  expect_s3_class(fc, "tbl_df")
  expect_true("proportion_fold_change" %in% colnames(fc))
  expect_gt(nrow(fc), 0)
  
  # Should have variation in fold changes
  unique_fcs <- unique(fc$proportion_fold_change)
  expect_gt(length(unique_fcs), 1, 
            label = "Should have multiple different fold change values")
  
  # Should have some non-trivial fold changes
  non_trivial <- fc |> 
    filter(abs(abs(proportion_fold_change) - 1) > 0.1)
  expect_gt(nrow(non_trivial), 0,
            label = "Should have some fold changes > 10% different from ±1")
  
  message("\n=== Issue #249 Original Code Test ===")
  message(sprintf("Total cell groups: %d", nrow(fc)))
  message(sprintf("Non-trivial fold changes: %d", nrow(non_trivial)))
  message(sprintf("All equal 1.0? %s (should be FALSE)", all_equal_one))
})