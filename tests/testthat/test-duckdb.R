test_that("sccomp_estimate has a DuckDB table method", {
  expect_true(is.function(getS3method("sccomp_estimate", "tbl_duckdb_connection")))
})

test_that("sccomp_estimate gives the same results for DuckDB tables", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("dbplyr")
  skip_cmdstan()

  data("counts_obj", package = "sccomp", envir = environment())
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbWriteTable(con, "counts_obj", counts_obj)

  estimate_args <- list(
    formula_composition = ~type,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    inference_method = "pathfinder",
    cores = 1,
    mcmc_seed = 12345,
    max_sampling_iterations = 1000,
    verbose = FALSE
  )

  data_frame_estimate <- do.call(
    sccomp_estimate,
    c(list(counts_obj), estimate_args)
  )
  duckdb_estimate <- do.call(
    sccomp_estimate,
    c(list(dplyr::tbl(con, "counts_obj")), estimate_args)
  )

  expect_equal(
    data_frame_estimate |>
      dplyr::arrange(cell_group, parameter) |>
      dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper),
    duckdb_estimate |>
      dplyr::arrange(cell_group, parameter) |>
      dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper),
      tolerance = 1e-8
  )
})
