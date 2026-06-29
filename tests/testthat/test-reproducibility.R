library(sccomp)
data("counts_obj")
n_iterations <- 1000
methods <- c("pathfinder", "hmc")

test_that("sccomp_estimate is reproducible with fixed mcmc_seed", {
  skip_cmdstan()
  
  for (method in methods) {
    fit1 <-
      counts_obj |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample = "sample",
        cell_group = "cell_group",
        abundance = "count",
        inference_method = method,
        cores = 1,
        mcmc_seed = 12345,
        max_sampling_iterations = n_iterations,
        verbose = FALSE
      )
    
    fit2 <-
      counts_obj |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample = "sample",
        cell_group = "cell_group",
        abundance = "count",
        inference_method = method,
        cores = 1,
        mcmc_seed = 12345,
        max_sampling_iterations = n_iterations,
        verbose = FALSE
      )
    
    expect_identical(
      attr(fit1, "fit")$draws(format = "matrix"),
      attr(fit2, "fit")$draws(format = "matrix"),
      info = method
    )
    
    expect_identical(
      fit1 |>
        dplyr::arrange(cell_group, parameter) |>
        dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper),
      fit2 |>
        dplyr::arrange(cell_group, parameter) |>
        dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper),
      info = method
    )
  }
})

test_that("sccomp_remove_outliers is reproducible with fixed mcmc_seed", {
  skip_cmdstan()
  
  for (method in methods) {
    fit0 <-
      counts_obj |>
      sccomp_estimate(
        formula_composition = ~ type,
        sample = "sample",
        cell_group = "cell_group",
        abundance = "count",
        inference_method = method,
        cores = 1,
        mcmc_seed = 12345,
        max_sampling_iterations = n_iterations,
        verbose = FALSE
      )
    
    out1 <-
      fit0 |>
      sccomp_remove_outliers(
        inference_method = method,
        cores = 1,
        mcmc_seed = 67890,
        max_sampling_iterations = n_iterations,
        verbose = FALSE
      )
    
    out2 <-
      fit0 |>
      sccomp_remove_outliers(
        inference_method = method,
        cores = 1,
        mcmc_seed = 67890,
        max_sampling_iterations = n_iterations,
        verbose = FALSE
      )
    
    expect_identical(
      attr(out1, "model_input")$truncation_not_idx,
      attr(out2, "model_input")$truncation_not_idx,
      info = method
    )
    
    expect_identical(
      attr(out1, "outliers") |>
        dplyr::arrange(sample, cell_group) |>
        dplyr::pull(outlier),
      attr(out2, "outliers") |>
        dplyr::arrange(sample, cell_group) |>
        dplyr::pull(outlier),
      info = method
    )
  }
})
