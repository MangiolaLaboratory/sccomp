# Test script for match.arg validation
library(sccomp)

cat("Testing match.arg for inference_method parameter...\n")

# Create a simple test data frame with integer counts
test_data <- data.frame(
  sample = c("s1", "s2", "s3"),
  cell_group = c("A", "B", "C"),
  count = as.integer(c(10, 20, 30)),
  type = c("control", "treatment", "control")
)

# Test 1: Valid values should work
cat("Test 1: Valid values\n")
valid_methods <- c("pathfinder", "hmc", "variational")

for (method in valid_methods) {
  cat("Testing '", method, "'...\n", sep = "")
  tryCatch({
    # This should work with valid inference_method
    # We'll just test that the function accepts the parameter
    result <- sccomp_estimate(
      test_data,
      ~ type,
      ~ 1,
      "sample",
      "cell_group", 
      "count",
      inference_method = method,
      cores = 1
    )
    cat("  ✓ Success with '", method, "'\n", sep = "")
  }, error = function(e) {
    cat("  ✗ Error with '", method, "': ", e$message, "\n", sep = "")
  })
}

# Test 2: Invalid value should give clear error
cat("\nTest 2: Invalid value\n")
tryCatch({
  result <- sccomp_estimate(
    test_data,
    ~ type,
    ~ 1,
    "sample",
    "cell_group",
    "count", 
    inference_method = "invalid_method",
    cores = 1
  )
  cat("  ✗ Should have failed with invalid method\n")
}, error = function(e) {
  cat("  ✓ Correctly caught error: ", e$message, "\n", sep = "")
})

cat("\nTest completed.\n") 