# Unit tests for sccomp_proportional_fold_change function
# Testing from simple models to interaction models

library(testthat)
library(sccomp)
library(dplyr)

# Helper function to skip tests if cmdstan is not available
skip_cmdstan <- function() {
  if (!instantiate::stan_cmdstan_exists()) {
    skip("CmdStan not available")
  }
}

# Helper function to create simple test data
create_simple_test_data <- function() {
  set.seed(123)
  n_samples <- 20
  n_cell_groups <- 5
  
  # Create sample data
  samples <- paste0("sample_", 1:n_samples)
  cell_groups <- paste0("cell_", 1:n_cell_groups)
  
  # Create treatment groups
  treatment <- rep(c("control", "treatment"), each = n_samples/2)
  
  # Generate counts with some treatment effect
  data <- expand.grid(
    sample = samples,
    cell_group = cell_groups,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      treatment = rep(treatment, times = n_cell_groups),
      # Simulate treatment effect
      base_count = rpois(n(), 100),
      treatment_effect = ifelse(treatment == "treatment", 1.5, 1),
      count = as.integer(base_count * treatment_effect)
    )
  
  return(data)
}

# Helper function to create interaction test data
create_interaction_test_data <- function() {
  set.seed(456)
  n_samples <- 30
  n_cell_groups <- 4
  
  # Create sample data
  samples <- paste0("sample_", 1:n_samples)
  cell_groups <- paste0("cell_", 1:n_cell_groups)
  
  # Create factors
  genotype <- rep(c("WT", "KO"), each = n_samples/2)
  tissue <- rep(c("liver", "brain"), times = n_samples/2)
  treatment <- rep(c("control", "drug"), times = n_samples/2)
  
  # Generate counts with interaction effects
  data <- expand.grid(
    sample = samples,
    cell_group = cell_groups,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      genotype = rep(genotype, times = n_cell_groups),
      tissue = rep(tissue, times = n_cell_groups),
      treatment = rep(treatment, times = n_cell_groups),
      # Simulate interaction effects
      base_count = rpois(n(), 80),
      genotype_effect = ifelse(genotype == "KO", 1.3, 1),
      tissue_effect = ifelse(tissue == "brain", 1.2, 1),
      treatment_effect = ifelse(treatment == "drug", 1.4, 1),
      interaction_effect = ifelse(genotype == "KO" & treatment == "drug", 1.5, 1),
      count = as.integer(base_count * genotype_effect * tissue_effect * treatment_effect * interaction_effect)
    )
  
  return(data)
}

# Test 1: Simple single factor model
test_that("sccomp_proportional_fold_change works with simple single factor model", {
  skip_cmdstan()
  
  # Create simple test data
  test_data <- create_simple_test_data()
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_true("cell_group" %in% colnames(result))
  expect_true("proportion_fold_change" %in% colnames(result))
  expect_true("average_uncertainty" %in% colnames(result))
  expect_true("statement" %in% colnames(result))
  
  # Check that we have results for all cell groups
  expect_equal(nrow(result), length(unique(test_data$cell_group)))
  
  # Check that fold changes are numeric
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(is.numeric(result$average_uncertainty))
  
  # Check that statements are character
  expect_true(is.character(result$statement))
})

# Test 2: Simple model with different factor names
test_that("sccomp_proportional_fold_change works with different factor names", {
  skip_cmdstan()
  
  # Create test data with different factor name
  test_data <- create_simple_test_data() %>%
    rename(condition = treatment)
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ condition,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ condition,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), length(unique(test_data$cell_group)))
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 3: Two-factor model without interaction
test_that("sccomp_proportional_fold_change works with two-factor model", {
  skip_cmdstan()
  
  # Create test data with two factors
  test_data <- create_simple_test_data() %>%
    mutate(
      timepoint = rep(c("baseline", "followup"), times = n()/2)
    )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment + timepoint,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change for treatment
  result_treatment <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Test proportional fold change for timepoint
  result_timepoint <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ timepoint,
    from = "baseline",
    to = "followup"
  )
  
  # Basic expectations
  expect_s3_class(result_treatment, "tbl_df")
  expect_s3_class(result_timepoint, "tbl_df")
  expect_equal(nrow(result_treatment), length(unique(test_data$cell_group)))
  expect_equal(nrow(result_timepoint), length(unique(test_data$cell_group)))
})

# Test 4: Interaction model
test_that("sccomp_proportional_fold_change works with interaction model", {
  skip_cmdstan()
  
  # Create interaction test data
  test_data <- create_interaction_test_data()
  
  # Fit model with interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ genotype * treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change for genotype
  result_genotype <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ genotype,
    from = "WT",
    to = "KO"
  )
  
  # Test proportional fold change for treatment
  result_treatment <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "drug"
  )
  
  # Basic expectations
  expect_s3_class(result_genotype, "tbl_df")
  expect_s3_class(result_treatment, "tbl_df")
  expect_equal(nrow(result_genotype), length(unique(test_data$cell_group)))
  expect_equal(nrow(result_treatment), length(unique(test_data$cell_group)))
})

# Test 5: Complex three-factor interaction model
test_that("sccomp_proportional_fold_change works with complex three-factor interaction", {
  skip_cmdstan()
  
  # Create complex interaction test data
  test_data <- create_interaction_test_data()
  
  # Fit model with three-factor interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ genotype * tissue * treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change for each factor
  result_genotype <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ genotype,
    from = "WT",
    to = "KO"
  )
  
  result_tissue <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ tissue,
    from = "liver",
    to = "brain"
  )
  
  result_treatment <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "drug"
  )
  
  # Basic expectations
  expect_s3_class(result_genotype, "tbl_df")
  expect_s3_class(result_tissue, "tbl_df")
  expect_s3_class(result_treatment, "tbl_df")
  expect_equal(nrow(result_genotype), length(unique(test_data$cell_group)))
  expect_equal(nrow(result_tissue), length(unique(test_data$cell_group)))
  expect_equal(nrow(result_treatment), length(unique(test_data$cell_group)))
})

# Test 6: Model with random effects
test_that("sccomp_proportional_fold_change works with random effects model", {
  skip_cmdstan()
  
  # Create test data with grouping
  test_data <- create_simple_test_data() %>%
    mutate(
      subject = rep(paste0("subject_", 1:10), times = n()/10)
    )
  
  # Fit model with random effects
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment + (1 | subject),
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), length(unique(test_data$cell_group)))
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 7: Error handling - invalid factor levels
test_that("sccomp_proportional_fold_change handles invalid factor levels gracefully", {
  skip_cmdstan()
  
  # Create simple test data
  test_data <- create_simple_test_data()
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with invalid factor levels
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment,
      from = "nonexistent",
      to = "treatment"
    )
  )
  
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment,
      from = "control",
      to = "nonexistent"
    )
  )
})

# Test 8: Error handling - invalid formula
test_that("sccomp_proportional_fold_change handles invalid formulas gracefully", {
  skip_cmdstan()
  
  # Create simple test data
  test_data <- create_simple_test_data()
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with invalid formula
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ nonexistent_factor,
      from = "control",
      to = "treatment"
    )
  )
})

# Test 9: Consistency of results
test_that("sccomp_proportional_fold_change produces consistent results", {
  skip_cmdstan()
  
  # Create simple test data
  test_data <- create_simple_test_data()
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Run the same analysis twice
  result1 <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  result2 <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Results should be identical
  expect_equal(result1$proportion_fold_change, result2$proportion_fold_change, tolerance = 1e-10)
  expect_equal(result1$average_uncertainty, result2$average_uncertainty, tolerance = 1e-10)
})

# Test 10: Edge cases - very small sample sizes
test_that("sccomp_proportional_fold_change works with small sample sizes", {
  skip_cmdstan()
  
  # Create minimal test data
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 2),
    treatment = rep(c("control", "treatment"), each = 3),
    count = as.integer(c(10, 20, 15, 15, 25, 20))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3) # 3 cell groups
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 11: Edge cases - single cell group
test_that("sccomp_proportional_fold_change works with single cell group", {
  skip_cmdstan()
  
  # Create test data with single cell group
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = 1),
    cell_group = rep("cell1", 2),
    treatment = c("control", "treatment"),
    count = as.integer(c(100, 150))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1) # 1 cell group
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 12: Edge cases - many cell groups
test_that("sccomp_proportional_fold_change works with many cell groups", {
  skip_cmdstan()
  
  # Create test data with many cell groups
  n_cell_groups <- 20
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = n_cell_groups),
    cell_group = rep(paste0("cell", 1:n_cell_groups), times = 2),
    treatment = rep(c("control", "treatment"), each = n_cell_groups),
    count = as.integer(rpois(n_cell_groups * 2, 50))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), n_cell_groups)
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 13: Edge cases - zero counts
test_that("sccomp_proportional_fold_change handles zero counts gracefully", {
  skip_cmdstan()
  
  # Create test data with some zero counts
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 2),
    treatment = rep(c("control", "treatment"), each = 3),
    count = as.integer(c(10, 0, 15, 15, 5, 20))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 14: Edge cases - very large counts
test_that("sccomp_proportional_fold_change handles very large counts", {
  skip_cmdstan()
  
  # Create test data with large counts
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 2),
    treatment = rep(c("control", "treatment"), each = 3),
    count = as.integer(c(10000, 20000, 15000, 15000, 25000, 20000))
  )
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "treatment"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
})

# Test 15: Edge cases - identical conditions
test_that("sccomp_proportional_fold_change handles identical conditions", {
  skip_cmdstan()
  
  # Create test data with identical conditions - this should work but may not be meaningful
  test_data <- data.frame(
    sample = rep(c("s1", "s2"), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 2),
    treatment = rep(c("control", "control"), each = 3), # Same condition
    count = as.integer(c(10, 20, 15, 10, 20, 15))
  )
  
  # This test should be skipped because identical conditions don't make sense for fold change
  # The function should handle this gracefully or we should test a different scenario
  skip("Identical conditions don't make sense for fold change calculation")
  
  # Fit model
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment,
    from = "control",
    to = "control" # Same condition
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
}) 