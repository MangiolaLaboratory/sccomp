test_that("model_cache_directory parameter works correctly", {
  # Test 1: Verify that get_sccomp_cache_dir properly versions cache directories
  test_that("get_sccomp_cache_dir appends version to cache directories", {
    version <- as.character(utils::packageVersion("sccomp"))
    
    # Test default cache directory
    default_cache <- get_sccomp_cache_dir()
    expect_true(grepl(paste0("/", version, "$"), default_cache))
    
    # Test custom cache directory
    custom_cache <- get_sccomp_cache_dir("/my/custom/path")
    expect_equal(custom_cache, file.path("/my/custom/path", version))
    
    # Test that version is always appended
    expect_false(identical(default_cache, file.path(path.expand("~"), ".sccomp_models")))
    expect_false(identical(custom_cache, "/my/custom/path"))
  })
  
  # Test 2: Verify that all main functions accept model_cache_directory parameter
  test_that("all main functions accept model_cache_directory parameter", {
    # Check function signatures
    expect_true("model_cache_directory" %in% names(formals(sccomp_estimate)))
    expect_true("model_cache_directory" %in% names(formals(simulate_data)))
    expect_true("model_cache_directory" %in% names(formals(sccomp_remove_outliers)))
    expect_true("model_cache_directory" %in% names(formals(sccomp_replicate)))
    
    # Check default values
    expect_equal(formals(sccomp_estimate)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
    expect_equal(formals(simulate_data)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
    expect_equal(formals(sccomp_remove_outliers)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
    expect_equal(formals(sccomp_replicate)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
  })
  
  # Test 3: Verify that load_model accepts and uses model_cache_directory
  test_that("load_model accepts model_cache_directory parameter", {
    expect_true("model_cache_directory" %in% names(formals(load_model)))
    expect_equal(formals(load_model)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
  })
  
  # Test 4: Verify that clear_stan_model_cache accepts model_cache_directory
  test_that("clear_stan_model_cache accepts model_cache_directory parameter", {
    expect_true("model_cache_directory" %in% names(formals(clear_stan_model_cache)))
    expect_equal(formals(clear_stan_model_cache)$model_cache_directory, quote(sccomp_stan_models_cache_dir))
  })
  
  # Test 5: Test that custom cache directory can be created
  test_that("custom cache directory can be created", {
    temp_cache_dir <- tempfile("sccomp_test_cache")
    on.exit(unlink(temp_cache_dir, recursive = TRUE, force = TRUE))
    
    expect_no_error({
      clear_stan_model_cache(temp_cache_dir)
    })
    # The function should not error even if directory doesn't exist
    expect_true(TRUE) # If we get here, no error was thrown
  })
  
  # Test 6: Verify default cache directory is accessible
  test_that("default cache directory is accessible", {
    default_cache <- sccomp_stan_models_cache_dir
    expect_true(is.character(default_cache))
    expect_true(length(default_cache) == 1)
    expect_true(nchar(default_cache) > 0)
  })
})

test_that("get_sccomp_cache_dir function works correctly", {
  version <- as.character(utils::packageVersion("sccomp"))
  
  # Test 1: Default behavior (NULL input)
  test_that("get_sccomp_cache_dir with NULL input", {
    result <- get_sccomp_cache_dir()
    expected <- file.path(path.expand("~"), ".sccomp_models", version)
    expect_equal(result, expected)
  })
  
  # Test 2: Custom path input
  test_that("get_sccomp_cache_dir with custom path", {
    result <- get_sccomp_cache_dir("/my/custom/path")
    expected <- file.path("/my/custom/path", version)
    expect_equal(result, expected)
  })
  
  # Test 3: Path with spaces
  test_that("get_sccomp_cache_dir with path containing spaces", {
    result <- get_sccomp_cache_dir("/my custom path")
    expected <- file.path("/my custom path", version)
    expect_equal(result, expected)
  })
  
  # Test 4: Relative path
  test_that("get_sccomp_cache_dir with relative path", {
    result <- get_sccomp_cache_dir("./relative/path")
    expected <- file.path("./relative/path", version)
    expect_equal(result, expected)
  })
  
  # Test 5: Empty string (edge case)
  test_that("get_sccomp_cache_dir with empty string", {
    result <- get_sccomp_cache_dir("")
    expected <- file.path("", version)
    expect_equal(result, expected)
  })
  
  # Test 6: Verify version is always appended
  test_that("get_sccomp_cache_dir always appends version", {
    # Test with various inputs
    inputs <- c(
      NULL,
      "/simple/path",
      "/path/with/spaces",
      "./relative/path",
      "",
      "/very/long/path/with/many/directories"
    )
    
    for (input in inputs) {
      result <- get_sccomp_cache_dir(input)
      # Should end with version
      expect_true(endsWith(result, version))
      # Should not be the same as input (unless input is NULL)
      if (!is.null(input)) {
        expect_false(identical(result, input))
      }
    }
  })
  
  # Test 7: Verify sccomp_stan_models_cache_dir is not versioned
  test_that("sccomp_stan_models_cache_dir is not versioned", {
    base_cache <- sccomp_stan_models_cache_dir
    versioned_cache <- get_sccomp_cache_dir(base_cache)
    
    # The base cache should not contain version
    expect_false(grepl(version, base_cache))
    # The versioned cache should contain version
    expect_true(grepl(version, versioned_cache))
    # They should be different
    expect_false(identical(base_cache, versioned_cache))
  })
}) 