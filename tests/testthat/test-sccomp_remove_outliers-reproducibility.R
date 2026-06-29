library(sccomp)
data("counts_obj")
n_iterations <- 1000

test_that("sccomp_remove_outliers is reproducible with fixed mcmc_seed", {
  skip_cmdstan()
  
  fit0 =
    counts_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      inference_method = "pathfinder",
      cores = 1,
      mcmc_seed = 12345,
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )
  
  out1 =
    fit0 |>
    sccomp_remove_outliers(
      inference_method = "pathfinder",
      cores = 1,
      mcmc_seed = 67890,
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )
  
  out2 =
    fit0 |>
    sccomp_remove_outliers(
      inference_method = "pathfinder",
      cores = 1,
      mcmc_seed = 67890,
      max_sampling_iterations = n_iterations,
      verbose = FALSE
    )
  
  expect_identical(
    attr(out1, "model_input")$truncation_not_idx,
    attr(out2, "model_input")$truncation_not_idx
  )
  
  expect_identical(
    attr(out1, "outliers") |>
      dplyr::arrange(sample, cell_group) |>
      dplyr::pull(outlier),
    attr(out2, "outliers") |>
      dplyr::arrange(sample, cell_group) |>
      dplyr::pull(outlier)
  )
})
