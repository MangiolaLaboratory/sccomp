test_that("convergence metrics are available for both fixed effects and random effects", {
  # Skip if cmdstanr is not available
  skip_if_not(instantiate::stan_cmdstan_exists(), "cmdstanr not available")
  
  # Load test data
  data("counts_obj")
  
  # Test 1: Model with only fixed effects (no random effects)
  # This should have convergence metrics for all parameters
  fixed_only_model <- sccomp_estimate(
    counts_obj,
    ~ type, ~ 1, "sample", "cell_group", "count",
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
  
  # Check that fixed effects have convergence metrics
  fixed_effects <- fixed_only_model |> 
    filter(!str_detect(parameter, "___")) # Exclude random effects (which shouldn't exist here)
  
  expect_true(all(c("c_rhat", "c_ess_bulk", "c_ess_tail") %in% colnames(fixed_effects)))
  expect_true(all(!is.na(fixed_effects$c_rhat)))
  expect_true(all(!is.na(fixed_effects$c_ess_bulk)))
  expect_true(all(!is.na(fixed_effects$c_ess_tail)))
  
  # Test 2: Model with random effects using existing data structure
  # Use the existing 'type' variable as a grouping factor for random effects
  random_effects_model <- sccomp_estimate(
    counts_obj,
    ~ type + (1 | type), ~ 1, "sample", "cell_group", "count",
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
  
  # Separate fixed and random effects
  fixed_effects_re <- random_effects_model |> 
    filter(!str_detect(parameter, "___"))
  
  random_effects_re <- random_effects_model |> 
    filter(str_detect(parameter, "___"))
  
  # Check that fixed effects have convergence metrics
  expect_true(all(c("c_rhat", "c_ess_bulk", "c_ess_tail") %in% colnames(fixed_effects_re)))
  expect_true(all(!is.na(fixed_effects_re$c_rhat)))
  expect_true(all(!is.na(fixed_effects_re$c_ess_bulk)))
  expect_true(all(!is.na(fixed_effects_re$c_ess_tail)))
  
  # Check that random effects NOW HAVE convergence metrics
  expect_true(all(c("c_rhat", "c_ess_bulk", "c_ess_tail") %in% colnames(random_effects_re)))
  expect_true(all(!is.na(random_effects_re$c_rhat)))
  expect_true(all(!is.na(random_effects_re$c_ess_bulk)))
  expect_true(all(!is.na(random_effects_re$c_ess_tail)))
  
  # Test 3: Verify this applies to both composition and variability parameters
  # Model with variability formula
  variability_model <- sccomp_estimate(
    counts_obj,
    ~ type + (1 | type), ~ type, "sample", "cell_group", "count",
    cores = 1, verbose = FALSE
  ) |> 
  sccomp_test()
  
  # Check composition fixed effects
  comp_fixed <- variability_model |> 
    filter(!str_detect(parameter, "___") & str_detect(parameter, "^type"))
  
  expect_true(all(!is.na(comp_fixed$c_rhat)))
  
  # Check composition random effects
  comp_random <- variability_model |> 
    filter(str_detect(parameter, "___"))
  
  expect_true(all(!is.na(comp_random$c_rhat)))
  
  # Check variability fixed effects (should have convergence metrics)
  var_fixed <- variability_model |> 
    filter(!str_detect(parameter, "___") & str_detect(parameter, "^type"))
  
  expect_true(all(!is.na(var_fixed$v_rhat)))
  
  # Summary of findings
  cat("\nTest Results Summary:\n")
  cat("- Fixed effects: Convergence metrics available (rhat, ess_bulk, ess_tail)\n")
  cat("- Random effects: Convergence metrics NOW available (rhat, ess_bulk, ess_tail)\n")
  cat("- This applies to both composition (c_) and variability (v_) parameters\n")
}) 