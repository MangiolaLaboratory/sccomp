test_that("sccomp_estimate has a DuckDB table method", {
  expect_true(is.function(getS3method("sccomp_estimate", "tbl_duckdb_connection")))
})
