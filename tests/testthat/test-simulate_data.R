library(testthat)
library(dplyr)
library(sccomp)

# Helpers that need no fit

test_that("expand_variability_term_values recycles a scalar and matches names", {
  terms = c("(Intercept)", "typecancer")

  none = sccomp:::expand_variability_term_values(NULL, terms, "x")
  expect_equal(as.vector(none$provided), c(0L, 0L))

  scalar = sccomp:::expand_variability_term_values(-0.7, terms, "x")
  expect_equal(as.vector(scalar$values), c(-0.7, -0.7))
  expect_equal(as.vector(scalar$provided), c(1L, 1L))

  named = sccomp:::expand_variability_term_values(c(typecancer = 2), terms, "x")
  expect_equal(as.vector(named$values), c(0, 2))
  expect_equal(as.vector(named$provided), c(0L, 1L))

  expect_error(sccomp:::expand_variability_term_values(c(foo = 1), terms, "x"), "not in the variability design: foo")
  expect_error(sccomp:::expand_variability_term_values(c(1, 2), terms, "x"), "single number, or a vector named")
  expect_error(sccomp:::expand_variability_term_values("a", terms, "x"), "must be numeric")
})

test_that("coefficients_to_beta aligns cell groups and design columns by name", {
  coefficients = tibble(
    cell_group = c("b", "a", "c"),
    typecancer = c(1, -2, 1),
    unused = 9
  )

  beta = sccomp:::coefficients_to_beta(coefficients, "cell_group", c("(Intercept)", "typecancer"), c("a", "b", "c"))

  expect_equal(dimnames(beta), list(c("(Intercept)", "typecancer"), c("a", "b", "c")))
  expect_equal(unname(beta["(Intercept)", ]), c(0, 0, 0))
  expect_equal(unname(beta["typecancer", ]), c(-2, 1, 1))

  expect_error(
    sccomp:::coefficients_to_beta(coefficients, "cell_group", "typecancer", c("a", "b", "c", "d")),
    "must cover every fitted cell group.*Missing: d"
  )
  expect_error(
    sccomp:::coefficients_to_beta(coefficients, "cell_group", "typecancer", c("a", "b")),
    "not present in the fitted model.*: c"
  )
  expect_error(
    sccomp:::coefficients_to_beta(coefficients |> select(-cell_group), "cell_group", "typecancer", c("a", "b", "c")),
    "must contain the cell group column"
  )
})

test_that("coefficient_columns_to_beta maps by name, then by position", {
  data = tibble(
    sample = rep(c("s1", "s2"), each = 2),
    cell_group = rep(c("a", "b"), 2),
    b_0 = rep(c(1, -1), 2),
    typecancer = rep(c(-3, 3), 2)
  )

  beta = sccomp:::coefficient_columns_to_beta(data, c("b_0", "typecancer"), "cell_group", c("(Intercept)", "typecancer"), c("a", "b"))

  expect_equal(unname(beta["(Intercept)", ]), c(1, -1))
  expect_equal(unname(beta["typecancer", ]), c(-3, 3))
})

test_that("check_coefficients_sum_to_zero", {
  ok = rbind(`(Intercept)` = c(1, -0.5, -0.5), typecancer = c(0.1, 0.2, -0.3))
  expect_silent(sccomp:::check_coefficients_sum_to_zero(ok))

  bad = rbind(`(Intercept)` = c(1, -0.5, -0.5), typecancer = c(1, 0, 0))
  expect_error(sccomp:::check_coefficients_sum_to_zero(bad), "design column 'typecancer' do not sum to zero")
})

# Fitted models

if (instantiate::stan_cmdstan_exists()) {
  data("counts_obj")
  data("seurat_obj")
  set.seed(42)

  simulate_fit_v1 = sccomp_estimate(
    counts_obj, ~ type, ~ 1, "sample", "cell_group", "count",
    cores = 1, inference_method = "pathfinder", verbose = FALSE
  )

  simulate_fit_vt = sccomp_estimate(
    counts_obj, ~ type, ~ type, "sample", "cell_group", "count",
    cores = 1, inference_method = "pathfinder", verbose = FALSE
  )

  simulate_fit_exclude = sccomp_estimate(
    counts_obj, ~ type, ~ type, "sample", "cell_group", "count",
    cores = 1, inference_method = "pathfinder", verbose = FALSE,
    exclude_mean_variability_association = TRUE
  )

  simulate_fit_random_effect = sccomp_estimate(
    seurat_obj, ~ type + (1 | group__), ~ 1, "sample", "cell_group",
    cores = 1, inference_method = "pathfinder", verbose = FALSE
  )

  simulate_cell_groups = colnames(attr(simulate_fit_v1, "model_input")$y)
  simulate_M = length(simulate_cell_groups)

  # Distinct, sum-to-zero coefficients, so that misalignment would show
  simulate_coefficients =
    tibble(cell_group = simulate_cell_groups) |>
    mutate(
      `(Intercept)` = seq(-1, 1, length.out = simulate_M),
      typecancer = rev(seq(-0.5, 0.5, length.out = simulate_M))
    ) |>
    mutate(across(-cell_group, ~ .x - mean(.x)))

  simulate_samples = counts_obj |> distinct(sample, type)
}

alpha_simulated_draws = function(rng, term) {
  draws = rng$fit$draws("alpha_simulated", format = "matrix")
  draws = draws[, grep(sprintf("^alpha_simulated\\[%d,", term), colnames(draws)), drop = FALSE]
  matrix(as.numeric(draws), nrow(draws))
}

test_that("simulated variability is beta * slope + intercept of the provided values", {
  skip_cmdstan()

  rng = sccomp:::simulate_data_rng(
    simulate_fit_v1, ~ type, ~ 1,
    new_data = simulate_samples, coefficients = simulate_coefficients,
    mean_dispersion_slope = -0.7, mean_dispersion_intercept = 5,
    .sample = sample, .cell_group = cell_group, number_of_draws = 3, mcmc_seed = 42
  )

  expected = simulate_coefficients$`(Intercept)` * -0.7 + 5
  alpha = alpha_simulated_draws(rng, 1)
  expect_equal(alpha, matrix(expected, nrow(alpha), simulate_M, byrow = TRUE), tolerance = 1e-8)
})

test_that("a named slope overrides only its variability term", {
  skip_cmdstan()

  rng = sccomp:::simulate_data_rng(
    simulate_fit_vt, ~ type, ~ type,
    new_data = simulate_samples, coefficients = simulate_coefficients,
    mean_dispersion_slope = c(typecancer = 2), mean_dispersion_intercept = c(typecancer = 0.5),
    .sample = sample, .cell_group = cell_group, number_of_draws = 3, mcmc_seed = 42
  )

  alpha_type = alpha_simulated_draws(rng, 2)
  expect_equal(
    alpha_type,
    matrix(simulate_coefficients$typecancer * 2 + 0.5, nrow(alpha_type), simulate_M, byrow = TRUE),
    tolerance = 1e-8
  )

  # The intercept term keeps its posterior association, which varies by draw
  expect_gt(sd(alpha_simulated_draws(rng, 1)[, 1]), 0)
})

test_that("a model fitted without the mean-variability association simulates with slope 0", {
  skip_cmdstan()

  rng = sccomp:::simulate_data_rng(
    simulate_fit_exclude, ~ type, ~ type,
    new_data = simulate_samples, coefficients = simulate_coefficients,
    .sample = sample, .cell_group = cell_group, number_of_draws = 3, mcmc_seed = 42
  )

  alpha = alpha_simulated_draws(rng, 1)
  expect_equal(apply(alpha, 1, sd), rep(0, nrow(alpha)), tolerance = 1e-8)
})

test_that("sccomp_simulate returns one row per sample, cell group and draw", {
  skip_cmdstan()

  result = sccomp_simulate(
    simulate_fit_v1, ~ type, ~ 1,
    new_data = simulate_samples, coefficients = simulate_coefficients,
    number_of_draws = 2, cores = 1, mcmc_seed = 42
  )

  expect_s3_class(result, "tbl")
  expect_true(all(c("sample", "type", "cell_group", "generated_proportions", "generated_counts", "replicate") %in% colnames(result)))
  expect_equal(nrow(result), nrow(simulate_samples) * simulate_M * 2)
  expect_false(anyNA(result$generated_counts))
  expect_setequal(result$cell_group, simulate_coefficients$cell_group)

  proportion_sums =
    result |>
    group_by(sample, replicate) |>
    summarise(total = sum(generated_proportions), .groups = "drop")
  expect_equal(proportion_sums$total, rep(1, nrow(proportion_sums)), tolerance = 1e-8)
})

test_that("simulated proportions follow the provided coefficients", {
  skip_cmdstan()

  # A high, flat precision so the average lands close to softmax(X * beta)
  result = sccomp_simulate(
    simulate_fit_v1, ~ type, ~ 1,
    new_data = simulate_samples, coefficients = simulate_coefficients,
    mean_dispersion_slope = 0, mean_dispersion_intercept = 8,
    number_of_draws = 50, cores = 1, mcmc_seed = 1
  )

  expected =
    simulate_coefficients |>
    tidyr::crossing(type = c("benign", "cancer")) |>
    mutate(eta = `(Intercept)` + typecancer * (type == "cancer")) |>
    group_by(type) |>
    mutate(expected = exp(eta) / sum(exp(eta))) |>
    ungroup() |>
    select(cell_group, type, expected)

  observed =
    result |>
    mutate(type = as.character(type)) |>
    group_by(cell_group, type) |>
    summarise(observed = mean(generated_proportions), .groups = "drop") |>
    left_join(expected, by = c("cell_group", "type"))

  expect_lt(max(abs(observed$observed - observed$expected)), 0.005)
})

test_that("sccomp_simulate is reproducible with the same seeds", {
  skip_cmdstan()

  simulate_once = function() {
    set.seed(7)
    sccomp_simulate(
      simulate_fit_v1, ~ type, ~ 1,
      new_data = simulate_samples, coefficients = simulate_coefficients,
      cores = 1, mcmc_seed = 9
    )
  }

  expect_identical(simulate_once()$generated_counts, simulate_once()$generated_counts)
})

test_that("sccomp_simulate uses the posterior and the fitted data by default", {
  skip_cmdstan()

  result = sccomp_simulate(simulate_fit_v1, ~ type, ~ 1, cores = 1)

  expect_equal(nrow(result), nrow(attr(simulate_fit_v1, "count_data")))
  expect_false(anyNA(result$generated_counts))
})

test_that("sccomp_simulate takes the fitted variability formula when it is NULL", {
  skip_cmdstan()

  result = sccomp_simulate(
    simulate_fit_vt, ~ type, NULL,
    new_data = simulate_samples, coefficients = simulate_coefficients, cores = 1
  )

  expect_equal(nrow(result), nrow(simulate_samples) * simulate_M)
})

test_that("sccomp_simulate works on a subset of samples", {
  skip_cmdstan()

  subset_samples = simulate_samples |> slice(1:3)

  result = sccomp_simulate(
    simulate_fit_v1, ~ type, ~ 1,
    new_data = subset_samples, coefficients = simulate_coefficients, cores = 1
  )

  expect_setequal(result$sample, subset_samples$sample)
  expect_equal(nrow(result), 3 * simulate_M)
})

test_that("sccomp_simulate works with random effects", {
  skip_cmdstan()

  result = sccomp_simulate(simulate_fit_random_effect, ~ type + (1 | group__), ~ 1, cores = 1)

  expect_equal(nrow(result), nrow(attr(simulate_fit_random_effect, "count_data")))
  expect_false(anyNA(result$generated_counts))
})

test_that("sccomp_simulate validates the coefficients", {
  skip_cmdstan()

  expect_error(
    sccomp_simulate(
      simulate_fit_v1, ~ type, ~ 1, new_data = simulate_samples,
      coefficients = simulate_coefficients |> mutate(typecancer = typecancer + 1), cores = 1
    ),
    "do not sum to zero"
  )

  expect_error(
    sccomp_simulate(
      simulate_fit_v1, ~ type, ~ 1, new_data = simulate_samples,
      coefficients = bind_rows(simulate_coefficients, simulate_coefficients |> mutate(cell_group = paste0(cell_group, "_new"))),
      cores = 1
    ),
    "not present in the fitted model"
  )
})

test_that("simulate_data is deprecated and still accepts coefficient columns", {
  skip_cmdstan()

  lifecycle::expect_deprecated(
    result <- simulate_data(
      counts_obj |> mutate(b_0 = 0, b_1 = 0),
      simulate_fit_v1, ~ type, ~ 1, sample, cell_group, c(b_0, b_1),
      cores = 1
    )
  )

  expect_equal(nrow(result), nrow(counts_obj))
  expect_false(anyNA(result$generated_counts))
})
