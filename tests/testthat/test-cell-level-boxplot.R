test_that("sccomp_boxplot works with cell-level input", {
  skip_cmdstan()
  
  # Test that sccomp_boxplot works correctly when using cell-level input
  # This addresses the issue where boxplot failed with cell-level data from SingleCellExperiment
  
  # Use the standard test data that mimics the reported issue
  res <- sce_obj |>
    sccomp_estimate(
      ~ type,
      sample = "sample",
      cell_group = "cell_group",
      cores = 1,
      verbose = FALSE
    ) |>
    sccomp_remove_outliers(cores = 1, verbose = FALSE) |>
    sccomp_test()
  
  # This should not throw an error (was previously failing with cell-level input)
  plot_result <- res |> sccomp_boxplot(factor = "type")
  
  # Verify the result is a ggplot
  expect_s3_class(plot_result, "ggplot")
  
  # Test with different significance threshold
  plot_result2 <- res |> sccomp_boxplot(factor = "type", significance_threshold = 0.01)
  expect_s3_class(plot_result2, "ggplot")
  
  # Test with unwanted effects removal
  plot_result3 <- res |> sccomp_boxplot(factor = "type", remove_unwanted_effects = TRUE)
  expect_s3_class(plot_result3, "ggplot")
})

test_that("sccomp_boxplot handles missing factor in count_data gracefully", {
  skip_cmdstan()
  
  # Create a mock sccomp result where count_data doesn't have the factor column
  # This simulates the edge case that caused the original bug
  
  res <- seurat_obj |>
    sccomp_estimate(
      ~ type,
      sample = "sample",
      cell_group = "cell_group", 
      cores = 1,
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Modify the count_data attribute to not include the "type" column
  # This simulates the problematic scenario
  count_data_original <- attr(res, "count_data")
  count_data_missing_factor <- count_data_original |> 
    select(-type)  # Remove the factor column
  
  attr(res, "count_data") <- count_data_missing_factor
  
  # This should not fail - the function should fallback to getting factor values elsewhere
  expect_no_error({
    plot_result <- res |> sccomp_boxplot(factor = "type")
  })
  
  # Verify the result is still a ggplot
  plot_result <- res |> sccomp_boxplot(factor = "type")
  expect_s3_class(plot_result, "ggplot")
})

test_that("sccomp_boxplot reproduces maintainer's suggested test case", {
  skip_cmdstan()
  
  # This test reproduces exactly the test case suggested by the maintainer
  # in the GitHub issue comments
  
  res <- sce_obj |>
    SingleCellExperiment::colData() |>
    sccomp_estimate(
      ~ type,
      .sample = sample,
      .cell_group = cell_group,
      cores = 1
    ) |>
    sccomp_remove_outliers(cores = 1) |>
    sccomp_test()
  
  # This was failing before the fix
  expect_no_error({
    plot_result <- res |> sccomp_boxplot(factor = "type")
  })
  
  # Verify the result
  plot_result <- res |> sccomp_boxplot(factor = "type")
  expect_s3_class(plot_result, "ggplot")
})