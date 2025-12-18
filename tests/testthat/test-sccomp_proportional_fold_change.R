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
  # Generated quantities may show tiny numerical differences; require close agreement.
  expect_equal(result1$proportion_fold_change, result2$proportion_fold_change, tolerance = 1e-2)
  expect_equal(result1$average_uncertainty, result2$average_uncertainty, tolerance = 1e-2)
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

# Test 16: Interaction in from/to categories - two-factor interaction
test_that("sccomp_proportional_fold_change works with interaction categories in from/to", {
  skip_cmdstan()
  
  # Create test data with two factors that will have interactions
  test_data <- data.frame(
    sample = rep(paste0("s", 1:8), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 8),
    treatment = rep(c("control", "treatment"), each = 12),
    timepoint = rep(c("baseline", "followup"), each = 6, times = 2),
    count = as.integer(c(
      # control-baseline
      10, 20, 15,
      12, 22, 17,
      # control-followup  
      15, 25, 20,
      18, 28, 23,
      # treatment-baseline
      20, 30, 25,
      22, 32, 27,
      # treatment-followup
      30, 40, 35,
      35, 45, 40
    ))
  )
  
  # Fit model with interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change with interaction categories
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint,
    from = "control:baseline",
    to = "treatment:followup"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(all(result$proportion_fold_change > 0)) # Should be positive for treatment effect
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
})

# Test 17: Interaction in from/to categories - three-factor interaction
test_that("sccomp_proportional_fold_change works with three-factor interaction categories", {
  skip_cmdstan()
  
  # Create test data with three factors
  test_data <- data.frame(
    sample = rep(paste0("s", 1:16), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 16),
    treatment = rep(c("control", "treatment"), each = 24),
    timepoint = rep(c("baseline", "followup"), each = 12, times = 2),
    cohort = rep(c("A", "B"), each = 6, times = 4),
    count = as.integer(c(
      # control-baseline-A
      10, 20, 15,
      12, 22, 17,
      # control-baseline-B
      11, 21, 16,
      13, 23, 18,
      # control-followup-A
      15, 25, 20,
      18, 28, 23,
      # control-followup-B
      16, 26, 21,
      19, 29, 24,
      # treatment-baseline-A
      20, 30, 25,
      22, 32, 27,
      # treatment-baseline-B
      21, 31, 26,
      23, 33, 28,
      # treatment-followup-A
      30, 40, 35,
      35, 45, 40,
      # treatment-followup-B
      31, 41, 36,
      36, 46, 41
    ))
  )
  
  # Fit model with three-factor interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint * cohort,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change with three-factor interaction categories
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint * cohort,
    from = "control:baseline:A",
    to = "treatment:followup:B"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(all(result$proportion_fold_change > 0)) # Should be positive for treatment effect
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
})

# Test 18: Interaction categories with different factor orders
test_that("sccomp_proportional_fold_change works with different factor orders in interaction", {
  skip_cmdstan()
  
  # Create test data
  test_data <- data.frame(
    sample = rep(paste0("s", 1:8), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 8),
    treatment = rep(c("control", "treatment"), each = 12),
    timepoint = rep(c("baseline", "followup"), each = 6, times = 2),
    count = as.integer(c(
      # control-baseline
      10, 20, 15,
      12, 22, 17,
      # control-followup  
      15, 25, 20,
      18, 28, 23,
      # treatment-baseline
      20, 30, 25,
      22, 32, 27,
      # treatment-followup
      30, 40, 35,
      35, 45, 40
    ))
  )
  
  # Fit model with interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with different factor orders in the interaction
  result1 <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint,
    from = "control:baseline",
    to = "treatment:followup"
  )
  
  result2 <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint,
    from = "baseline:control",
    to = "followup:treatment"
  )
  
  # Both should work and produce similar results
  expect_s3_class(result1, "tbl_df")
  expect_s3_class(result2, "tbl_df")
  expect_equal(nrow(result1), 3)
  expect_equal(nrow(result2), 3)
  expect_true(is.numeric(result1$proportion_fold_change))
  expect_true(is.numeric(result2$proportion_fold_change))
})

# Test 19: Complex interaction with nested factors
test_that("sccomp_proportional_fold_change works with complex nested interactions", {
  skip_cmdstan()
  
  # Create test data with nested structure
  test_data <- data.frame(
    sample = rep(paste0("s", 1:12), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 12),
    treatment = rep(c("control", "treatment"), each = 18),
    timepoint = rep(c("baseline", "followup"), each = 9, times = 2),
    batch = rep(c("batch1", "batch2"), each = 3, times = 6),
    count = as.integer(c(
      # control-baseline-batch1
      10, 20, 15,
      12, 22, 17,
      11, 21, 16,
      # control-baseline-batch2
      13, 23, 18,
      14, 24, 19,
      15, 25, 20,
      # control-followup-batch1
      18, 28, 23,
      20, 30, 25,
      19, 29, 24,
      # control-followup-batch2
      21, 31, 26,
      22, 32, 27,
      23, 33, 28,
      # treatment-baseline-batch1
      25, 35, 30,
      27, 37, 32,
      26, 36, 31,
      # treatment-baseline-batch2
      28, 38, 33,
      29, 39, 34,
      30, 40, 35,
      # treatment-followup-batch1
      35, 45, 40,
      37, 47, 42,
      36, 46, 41,
      # treatment-followup-batch2
      38, 48, 43,
      39, 49, 44,
      40, 50, 45
    ))
  )
  
  # Fit model with complex interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint * batch,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with complex interaction categories
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint * batch,
    from = "control:baseline:batch1",
    to = "treatment:followup:batch2"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(all(result$proportion_fold_change > 0)) # Should be positive for treatment effect
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
})

# Test 20: Error handling for invalid interaction categories
test_that("sccomp_proportional_fold_change handles invalid interaction categories gracefully", {
  skip_cmdstan()
  
  # Create test data
  test_data <- data.frame(
    sample = rep(paste0("s", 1:4), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 4),
    treatment = rep(c("control", "treatment"), each = 6),
    timepoint = rep(c("baseline", "followup"), each = 3, times = 2),
    count = as.integer(c(10, 20, 15, 15, 25, 20, 12, 22, 17, 18, 28, 23))
  )
  
  # Fit model with interaction
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with invalid interaction category
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment * timepoint,
      from = "invalid:category",
      to = "treatment:followup"
    ),
    regexp = "Error in.*sccomp_predict"
  )
  
  # Test with missing factor in interaction
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment * timepoint,
      from = "control",
      to = "treatment:followup"
    ),
    regexp = "Error in.*sccomp_predict"
  )
})

# Test 21: Random effects in composition model
test_that("sccomp_proportional_fold_change works with random effects in composition", {
  skip_cmdstan()
  
  # Create test data with random effects
  test_data <- data.frame(
    sample = rep(paste0("s", 1:8), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 8),
    treatment = rep(c("control", "treatment"), each = 12),
    batch = rep(c("batch1", "batch2"), each = 6, times = 2),
    count = as.integer(c(
      # control-batch1
      10, 20, 15,
      10, 20, 15,
      # control-batch2
      12, 22, 17,
      12, 22, 17,
      # treatment-batch1
      20, 30, 25,
      20, 30, 25,
      # treatment-batch2
      22, 32, 27,
      22, 32, 27
    ))
  )
  
  # Fit model with random effects in composition
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment + (1|batch),
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change with random effects
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment + (1|batch),
    from = "control:batch1",
    to = "treatment:batch2"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
})

# Test 22: Random effects in variability model (current limitation)
test_that("sccomp_proportional_fold_change handles random effects in variability model", {
  skip_cmdstan()
  
  # Create test data with random effects in variability
  test_data <- data.frame(
    sample = rep(paste0("s", 1:8), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 8),
    treatment = rep(c("control", "treatment"), each = 12),
    batch = rep(c("batch1", "batch2"), each = 6, times = 2),
    count = as.integer(c(
      # control-batch1
      10, 20, 15,
      10, 20, 15,
      # control-batch2
      12, 22, 17,
      12, 22, 17,
      # treatment-batch1
      20, 30, 25,
      20, 30, 25,
      # treatment-batch2
      22, 32, 27,
      22, 32, 27
    ))
  )
  
  # Fit model with random effects in variability
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment,
    formula_variability = ~ (1|batch),
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change
  # Note: This should work but the variability random effects are ignored
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
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
  
  # Note: The variability random effects are ignored due to sccomp_predict limitation
  # This is a known limitation of the current implementation
})

# Test 23: Complex random effects with multiple grouping variables
test_that("sccomp_proportional_fold_change works with complex random effects", {
  skip_cmdstan()
  
  # Create test data with multiple random effects
  test_data <- data.frame(
    sample = rep(paste0("s", 1:12), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 12),
    treatment = rep(c("control", "treatment"), each = 18),
    batch = rep(c("batch1", "batch2", "batch3"), each = 12),
    timepoint = rep(c("baseline", "followup"), each = 6, times = 6),
    count = as.integer(c(
      # control-batch1-baseline
      10, 20, 15,
      10, 20, 15,
      # control-batch1-followup
      12, 22, 17,
      12, 22, 17,
      # control-batch2-baseline
      11, 21, 16,
      11, 21, 16,
      # control-batch2-followup
      13, 23, 18,
      13, 23, 18,
      # control-batch3-baseline
      12, 22, 17,
      12, 22, 17,
      # control-batch3-followup
      14, 24, 19,
      14, 24, 19,
      # treatment-batch1-baseline
      20, 30, 25,
      20, 30, 25,
      # treatment-batch1-followup
      22, 32, 27,
      22, 32, 27,
      # treatment-batch2-baseline
      21, 31, 26,
      21, 31, 26,
      # treatment-batch2-followup
      23, 33, 28,
      23, 33, 28,
      # treatment-batch3-baseline
      22, 32, 27,
      22, 32, 27,
      # treatment-batch3-followup
      24, 34, 29,
      24, 34, 29
    ))
  )
  
  # Fit model with complex random effects
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment * timepoint + (1|batch),
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test proportional fold change with complex random effects
  result <- sccomp_proportional_fold_change(
    estimate,
    formula_composition = ~ treatment * timepoint + (1|batch),
    from = "control:baseline:batch1",
    to = "treatment:followup:batch2"
  )
  
  # Basic expectations
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 3)
  expect_true(is.numeric(result$proportion_fold_change))
  expect_true(is.character(result$statement))
  expect_true(is.numeric(result$average_uncertainty))
})

# Test 24: Random effects parsing and new_data creation
test_that("sccomp_proportional_fold_change correctly parses random effects formulas", {
  skip_cmdstan()
  
  # Test random effects formula parsing
  formula_with_random <- ~ treatment + (1|batch)
  my_factors_random <- sccomp:::parse_formula(formula_with_random)
  
  # The parse_formula function should handle random effects
  expect_true(is.character(my_factors_random))
  expect_true(length(my_factors_random) > 0)
  
  # Test new_data creation for random effects
  .sample <- quo(sample)
  from_random <- "control:batch1"
  to_random <- "treatment:batch2"
  from_parts <- stringr::str_split(from_random, ":")[[1]]
  to_parts <- stringr::str_split(to_random, ":")[[1]]
  
  # Create new_data for random effects
  new_data_random <- tibble(
    !!quo_name(.sample) := c(to_random, from_random)
  )
  
  # Add factor columns
  if (length(my_factors_random) > 1) {
    for (i in seq_along(my_factors_random)) {
      new_data_random <- new_data_random %>%
        mutate(!!my_factors_random[i] := c(to_parts[i], from_parts[i]))
    }
  }
  
  # Verify the structure is correct
  expect_s3_class(new_data_random, "tbl_df")
  expect_equal(nrow(new_data_random), 2)
  expect_true("sample" %in% colnames(new_data_random))
})

# Test 25: Error handling for invalid random effects
test_that("sccomp_proportional_fold_change handles invalid random effects gracefully", {
  skip_cmdstan()
  
  # Create test data
  test_data <- data.frame(
    sample = rep(paste0("s", 1:4), each = 3),
    cell_group = rep(c("cell1", "cell2", "cell3"), times = 4),
    treatment = rep(c("control", "treatment"), each = 6),
    batch = rep(c("batch1", "batch2"), each = 3, times = 2),
    count = as.integer(c(10, 20, 15, 15, 25, 20, 12, 22, 17, 18, 28, 23))
  )
  
  # Fit model with random effects
  estimate <- sccomp_estimate(
    test_data,
    formula_composition = ~ treatment + (1|batch),
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    cores = 1,
    verbose = FALSE
  )
  
  # Test with invalid random effects category
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment + (1|batch),
      from = "invalid:category",
      to = "treatment:batch1"
    ),
    regexp = "Error in.*sccomp_predict"
  )
  
  # Test with missing random effect level
  expect_error(
    sccomp_proportional_fold_change(
      estimate,
      formula_composition = ~ treatment + (1|batch),
      from = "control",
      to = "treatment:batch1"
    ),
    regexp = "Error in.*sccomp_predict"
  )
})

# Test 26: Minimal dataset with hyphenated factor levels
test_that("sccomp_proportional_fold_change works with hyphenated factor levels - minimal dataset", {
  skip_cmdstan()
  
  # Minimal count data for testing hyphenated factor levels **with >2 levels**
  #
  # IMPORTANT (regression target):
  # The historical bug on master most reliably reproduces when the model is trained
  # with >2 levels, but `new_data` is constructed containing only the two levels
  # being compared. In that situation, the `new_data` design matrix can end up
  # missing columns and predictions for `from` and `to` collapse, yielding
  # fold-changes of ~1 for all cell groups.
  #
  # Dataset:
  # - 8 samples across 4 levels of `main_comparison_type`
  # - includes hyphenated levels "Pre-surgery" and "Post-surgery"
  # - still keeps only 3 cell groups for speed
  minimal_counts <- tibble::tribble(
    ~sample_donor, ~cell_group, ~main_comparison_type, ~count,
    # Pre-surgery (hyphenated)
    "donor1", "B_cell", "Pre-surgery", 50L,
    "donor1", "T_cell", "Pre-surgery", 30L,
    "donor1", "NK_cell", "Pre-surgery", 20L,
    "donor2", "B_cell", "Pre-surgery", 48L,
    "donor2", "T_cell", "Pre-surgery", 32L,
    "donor2", "NK_cell", "Pre-surgery", 20L,
    # Bacterial
    "donor3", "B_cell", "Bacterial", 80L,
    "donor3", "T_cell", "Bacterial", 10L,
    "donor3", "NK_cell", "Bacterial", 10L,
    "donor4", "B_cell", "Bacterial", 78L,
    "donor4", "T_cell", "Bacterial", 12L,
    "donor4", "NK_cell", "Bacterial", 10L,
    # Non-Infectious (additional level to force >2 levels in training)
    "donor5", "B_cell", "Non-Infectious", 20L,
    "donor5", "T_cell", "Non-Infectious", 60L,
    "donor5", "NK_cell", "Non-Infectious", 20L,
    "donor6", "B_cell", "Non-Infectious", 22L,
    "donor6", "T_cell", "Non-Infectious", 58L,
    "donor6", "NK_cell", "Non-Infectious", 20L,
    # Post-surgery (hyphenated additional level)
    "donor7", "B_cell", "Post-surgery", 40L,
    "donor7", "T_cell", "Post-surgery", 40L,
    "donor7", "NK_cell", "Post-surgery", 20L,
    "donor8", "B_cell", "Post-surgery", 38L,
    "donor8", "T_cell", "Post-surgery", 42L,
    "donor8", "NK_cell", "Post-surgery", 20L
  ) |>
    # CRITICAL: Convert to factor to match real-world data (like sce_dummy)
    # Without this, the bug may not manifest on buggy master
    mutate(main_comparison_type = factor(
      main_comparison_type,
      levels = c("Pre-surgery", "Post-surgery", "Non-Infectious", "Bacterial")
    ))
  
  # Run sccomp pipeline
  sccomp_result <- minimal_counts |>
    sccomp_estimate(
      formula_composition = ~main_comparison_type, 
      sample = "sample_donor",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      verbose = FALSE
    ) |>
    sccomp_test()
  
  # Calculate proportional fold changes with hyphenated factor levels
  fold_changes <- sccomp_proportional_fold_change(
    sccomp_result, 
    formula_composition = ~main_comparison_type,
    from = "Pre-surgery",  # Contains hyphen - this is the key test
    to = "Bacterial"
  )
  
  # Assertions
  expect_s3_class(fold_changes, "tbl_df")
  expect_equal(nrow(fold_changes), 3)
  expect_gte(dplyr::n_distinct(minimal_counts$main_comparison_type), 3)
  
  # Check that all expected cell groups are present
  expect_true(all(c("B_cell", "T_cell", "NK_cell") %in% fold_changes$cell_group))
  
  # Verify fold changes are non-trivial (not all equal to 1)
  # This was the bug we fixed - without proper factor level handling, all fold changes were 1
  expect_false(all(abs(fold_changes$proportion_fold_change == 1) ))
  
  # Check that we have at least some diversity in fold changes
  expect_gt(length(unique(fold_changes$proportion_fold_change)), 1)
  
  # Verify the structure of results
  expect_true("proportion_fold_change" %in% colnames(fold_changes))
  expect_true("average_uncertainty" %in% colnames(fold_changes))
  expect_true("statement" %in% colnames(fold_changes))
  
  # Check data types
  expect_true(is.numeric(fold_changes$proportion_fold_change))
  expect_true(is.numeric(fold_changes$average_uncertainty))
  expect_true(is.character(fold_changes$statement))
  
  # Verify statements contain fold change descriptions
  expect_true(all(grepl("fold", fold_changes$statement, ignore.case = TRUE)))
  expect_true(all(grepl("increase|decrease", fold_changes$statement, ignore.case = TRUE)))
}) 
