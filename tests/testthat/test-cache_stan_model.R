test_that("cache_stan_model parameter works correctly", {
  # Skip if cmdstanr is not available
  skip_if_not(instantiate::stan_cmdstan_exists())
  
  # Test that cache_stan_model parameter is available in all relevant functions
  expect_true("cache_stan_model" %in% names(formals(sccomp_estimate)))
  expect_true("cache_stan_model" %in% names(formals(sccomp_replicate)))
  expect_true("cache_stan_model" %in% names(formals(simulate_data)))
  expect_true("cache_stan_model" %in% names(formals(sccomp_remove_outliers)))
  
  # Test that cache_stan_model defaults to sccomp_stan_models_cache_dir
  expect_equal(formals(sccomp_estimate)$cache_stan_model, quote(sccomp_stan_models_cache_dir))
  expect_equal(formals(sccomp_replicate)$cache_stan_model, quote(sccomp_stan_models_cache_dir))
  expect_equal(formals(simulate_data)$cache_stan_model, quote(sccomp_stan_models_cache_dir))
  expect_equal(formals(sccomp_remove_outliers)$cache_stan_model, quote(sccomp_stan_models_cache_dir))
  
  # Test that sccomp:::load_model function handles cache_stan_model correctly
  # This is an internal test to verify the cache directory construction
  test_cache_dir <- tempdir()
  sccomp_version <- as.character(packageVersion("sccomp"))
  expected_cache_dir <- file.path(test_cache_dir, sccomp_version)
  
  # Test with default cache_stan_model (should use default)
  expect_no_error({
    model <- sccomp:::load_model("glm_multi_beta_binomial", threads = 1)
  })
  
  # Test with custom cache_stan_model
  expect_no_error({
    model <- sccomp:::load_model("glm_multi_beta_binomial", cache_dir = test_cache_dir, threads = 1)
  })
  
  # Verify that the cache directory was created with the sccomp version
  expect_true(dir.exists(expected_cache_dir))
  
  # Test that default cache directory includes version
  default_cache <- sccomp:::get_sccomp_cache_dir()
  expect_true(grepl(sccomp_version, default_cache))
  
  # Clean up
  unlink(expected_cache_dir, recursive = TRUE)
}) 
