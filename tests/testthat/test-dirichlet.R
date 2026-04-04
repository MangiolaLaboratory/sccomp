library(testthat)
library(dplyr)
library(sccomp)

test_that("sccomp_estimate with noise_model dirichlet_multinomial runs and returns expected structure", {
  skip_cmdstan()

  data(counts_obj)
  cache_dir <- file.path(tempdir(), "sccomp_test_cache_dirichlet")

  est <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    noise_model = "dirichlet_multinomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  expect_equal(attr(est, "noise_model"), "dirichlet_multinomial")
  expect_s3_class(est, "sccomp_tbl")
  expect_true(all(c("cell_group", "parameter", "factor", "c_effect", "v_effect") %in% names(est)))
  expect_gt(nrow(est), 0)
})

test_that("sccomp_estimate dirichlet_multinomial works with counts", {
  skip_cmdstan()

  data(counts_obj)
  cache_dir <- file.path(tempdir(), "sccomp_test_cache_dirichlet")

  est <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    noise_model = "dirichlet_multinomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  expect_equal(attr(est, "noise_model"), "dirichlet_multinomial")
  expect_s3_class(est, "sccomp_tbl")
  expect_gt(nrow(est), 0)
})

test_that("sccomp_test works with dirichlet_multinomial estimate", {
  skip_cmdstan()

  data(counts_obj)
  cache_dir <- file.path(tempdir(), "sccomp_test_cache_dirichlet")

  est <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    noise_model = "dirichlet_multinomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  tested <- est |> sccomp_test()
  expect_true(all(c("c_FDR", "v_FDR") %in% names(tested)))
})

test_that("composition means are highly correlated between beta-binomial and dirichlet-multinomial", {
  skip_cmdstan()

  data(counts_obj)
  cache_dir <- file.path(tempdir(), "sccomp_test_cache_dirichlet_correlation")

  est_bb <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    noise_model = "multi_beta_binomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  est_dm <- sccomp_estimate(
    counts_obj,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "count",
    noise_model = "dirichlet_multinomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  comparison <- est_bb |>
    select(cell_group, parameter, c_effect_bb = c_effect) |>
    inner_join(
      est_dm |> select(cell_group, parameter, c_effect_dir = c_effect),
      by = c("cell_group", "parameter")
    )

  r2 <- with(comparison, cor(c_effect_bb, c_effect_dir, use = "complete.obs")^2)
  expect_gt(r2, 0.8)
})

test_that("dirichlet_multinomial accepts proportion input", {
  skip_cmdstan()
  data(counts_obj)
  cache_dir <- file.path(tempdir(), "sccomp_test_cache_dirichlet_proportion")

  counts_prop <- counts_obj |>
    group_by(sample) |>
    mutate(proportion = count / sum(count)) |>
    ungroup()

  est <- sccomp_estimate(
    counts_prop,
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample = "sample",
    cell_group = "cell_group",
    abundance = "proportion",
    noise_model = "dirichlet_multinomial",
    cores = 1,
    verbose = FALSE,
    inference_method = "hmc",
    max_sampling_iterations = 500,
    cache_stan_model = cache_dir
  )

  expect_equal(attr(est, "noise_model"), "dirichlet_multinomial")
  expect_gt(nrow(est), 0)
})
