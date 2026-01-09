library(dplyr)
library(sccomp)
data("seurat_obj")

n_iterations = 50

# Create test estimates if cmdstan is available
if (instantiate::stan_cmdstan_exists()) {
  test_estimate = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1, 
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations, 
      verbose = FALSE
    )
  
  test_estimate_random = 
    seurat_obj |>
    sccomp_estimate(
      formula_composition = ~ type + (1 | group__),
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )
}

test_that("sccomp_predict returns unconstrained predictors in summary mode", {
  skip_cmdstan()
  
  # Skip if estimates weren't created
  if (!exists("test_estimate")) {
    skip("Test estimate not created")
  }
  
  # Create simple new data for prediction
  new_data_tibble = 
    seurat_obj[[]] |> 
    distinct(sample, type) |>
    head(2)
  
  # Test with summary mode (default)
  prediction_summary = 
    test_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble,
      summary_instead_of_draws = TRUE,
      robust = FALSE,
      number_of_draws = 10
    )
  
  # Check that unconstrained predictor columns exist
  expect_true("unconstrained_mean" %in% colnames(prediction_summary))
  expect_true("unconstrained_lower" %in% colnames(prediction_summary))
  expect_true("unconstrained_upper" %in% colnames(prediction_summary))
  
  # Check that standard columns still exist
  expect_true("proportion_mean" %in% colnames(prediction_summary))
  expect_true("proportion_lower" %in% colnames(prediction_summary))
  expect_true("proportion_upper" %in% colnames(prediction_summary))
  
  # Check that unconstrained values are numeric
  expect_true(is.numeric(prediction_summary$unconstrained_mean))
  expect_true(is.numeric(prediction_summary$unconstrained_lower))
  expect_true(is.numeric(prediction_summary$unconstrained_upper))
  
  # Check that unconstrained values are reasonable (not NA, finite)
  expect_false(any(is.na(prediction_summary$unconstrained_mean)))
  expect_true(all(is.finite(prediction_summary$unconstrained_mean)))
  
  # Check that unconstrained lower < mean < upper
  expect_true(all(prediction_summary$unconstrained_lower <= prediction_summary$unconstrained_mean, na.rm = TRUE))
  expect_true(all(prediction_summary$unconstrained_mean <= prediction_summary$unconstrained_upper, na.rm = TRUE))
})

test_that("sccomp_predict returns unconstrained predictors with robust=TRUE", {
  skip_cmdstan()
  
  if (!exists("test_estimate")) {
    skip("Test estimate not created")
  }
  
  new_data_tibble = 
    seurat_obj[[]] |> 
    distinct(sample, type) |>
    head(2)
  
  # Test with robust statistics
  prediction_robust = 
    test_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble,
      summary_instead_of_draws = TRUE,
      robust = TRUE,
      number_of_draws = 10
    )
  
  # Check that unconstrained predictor columns exist
  expect_true("unconstrained_mean" %in% colnames(prediction_robust))
  expect_true("unconstrained_lower" %in% colnames(prediction_robust))
  expect_true("unconstrained_upper" %in% colnames(prediction_robust))
  
  # Check values are numeric and finite
  expect_true(is.numeric(prediction_robust$unconstrained_mean))
  expect_false(any(is.na(prediction_robust$unconstrained_mean)))
  expect_true(all(is.finite(prediction_robust$unconstrained_mean)))
})

test_that("sccomp_predict returns unconstrained predictors in draws mode", {
  skip_cmdstan()
  
  if (!exists("test_estimate")) {
    skip("Test estimate not created")
  }
  
  new_data_tibble = 
    seurat_obj[[]] |> 
    distinct(sample, type) |>
    head(1)
  
  # Test with draws mode (summary_instead_of_draws = FALSE)
  prediction_draws = 
    test_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble,
      summary_instead_of_draws = FALSE,
      number_of_draws = 10
    )
  
  # Check that unconstrained column exists in draws mode
  expect_true("unconstrained" %in% colnames(prediction_draws))
  expect_true(".draw" %in% colnames(prediction_draws))
  expect_true("proportion" %in% colnames(prediction_draws))
  
  # Check that unconstrained values are numeric
  expect_true(is.numeric(prediction_draws$unconstrained))
  
  # Check that unconstrained values are reasonable
  expect_false(any(is.na(prediction_draws$unconstrained)))
  expect_true(all(is.finite(prediction_draws$unconstrained)))
  
  # Check that we have the expected number of draws
  expect_equal(length(unique(prediction_draws$.draw)), 10)
})

test_that("sccomp_calculate_residuals includes unconstrained predictors", {
  skip_cmdstan()
  
  if (!exists("test_estimate_random")) {
    skip("Test estimate not created")
  }
  
  # Calculate residuals
  residuals_result = 
    test_estimate_random |> 
    sccomp_calculate_residuals()
  
  # Check that standard columns exist
  expect_true("sample" %in% colnames(residuals_result))
  expect_true("cell_group" %in% colnames(residuals_result))
  expect_true("residuals" %in% colnames(residuals_result))
  expect_true("exposure" %in% colnames(residuals_result))
  
  # Check that unconstrained predictor column exists
  expect_true("residuals_unconstrained" %in% colnames(residuals_result))
  
  # Check that unconstrained values are numeric
  expect_true(is.numeric(residuals_result$residuals_unconstrained))
  
  # Check that unconstrained values are reasonable
  expect_false(any(is.na(residuals_result$residuals_unconstrained)))
  expect_true(all(is.finite(residuals_result$residuals_unconstrained)))
})

test_that("unconstrained predictors are different from proportions (before vs after softmax)", {
  skip_cmdstan()
  
  if (!exists("test_estimate")) {
    skip("Test estimate not created")
  }
  
  new_data_tibble = 
    seurat_obj[[]] |> 
    distinct(sample, type) |>
    head(2)
  
  # Get predictions
  prediction = 
    test_estimate |>
    sccomp_predict(
      formula_composition = ~ type,
      new_data = new_data_tibble,
      summary_instead_of_draws = TRUE,
      number_of_draws = 10
    )
  
  # Verify that unconstrained values exist and are different from proportions
  # (since unconstrained are on logit scale before softmax, they should be different)
  # Note: We can't directly compare because softmax transforms, but we can check
  # that unconstrained values are not bounded [0,1] like proportions
  expect_true(any(prediction$unconstrained_mean > 1 | prediction$unconstrained_mean < -1))
  
  # Check that proportions are in [0, 1]
  expect_true(all(prediction$proportion_mean >= 0 & prediction$proportion_mean <= 1, na.rm = TRUE))
  
  # Check that proportions sum to approximately 1 per sample
  proportions_sum = 
    prediction |>
    with_groups(sample, ~ .x |> summarise(sum_prop = sum(proportion_mean)))
  
  expect_true(all(proportions_sum$sum_prop > 0.95 & proportions_sum$sum_prop < 1.05))
})

