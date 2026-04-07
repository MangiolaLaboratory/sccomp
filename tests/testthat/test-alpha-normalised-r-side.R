library(dplyr)
library(sccomp)

estimate_for_alpha_normalisation_tests <- function(output_directory, bimodal = FALSE) {
  data("counts_obj", package = "sccomp", envir = environment())
  cache_dir <- file.path(output_directory, "stan_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  counts_obj |>
    sccomp_estimate(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      inference_method = "pathfinder",
      max_sampling_iterations = 500,
      output_samples = 200,
      verbose = FALSE,
      output_directory = output_directory,
      portable = FALSE,
      cache_stan_model = cache_dir,
      bimodal_mean_variability_association = bimodal
    )
}

test_that("alpha_normalised is computed in R for unimodal fits", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_alpha_norm_unimodal_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_alpha_normalisation_tests(test_output_dir, bimodal = FALSE)
  fit <- attr(result, "fit")
  model_input <- attr(result, "model_input")

  # Use one variability coefficient across all cell groups to keep the test fast.
  n_m <- ncol(model_input$y)
  alpha_subset <- sprintf("alpha[%d,%d]", 1L, seq_len(n_m))

  calc <- sccomp:::compute_alpha_normalised_draws(
    fit = fit,
    model_input = model_input,
    alpha_variable_subset = alpha_subset
  ) |>
    arrange(C, M, .chain, .iteration) |>
    select(C, M, .chain, .iteration, calc = .value)

  alpha <- sccomp:::draws_to_tibble_x_y(fit, alpha_subset, "C", "M") |>
    rename(alpha = .value)
  beta <- sccomp:::draws_to_tibble_x_y(fit, sprintf("beta[%d,%d]", 1L, seq_len(n_m)), "C", "M") |>
    rename(C_comp = C, beta = .value)
  slope_1 <- sccomp:::draws_to_tibble_x(fit, "prec_slope_1", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_slope_1 = .value)

  expected <- alpha |>
    mutate(C_comp = model_input$variability_to_composition_map[C]) |>
    left_join(beta, by = c("C_comp", "M", ".chain", ".iteration")) |>
    left_join(slope_1, by = c("C", ".chain", ".iteration")) |>
    mutate(expected = alpha - (beta * prec_slope_1)) |>
    arrange(C, M, .chain, .iteration) |>
    select(C, M, .chain, .iteration, expected)

  joined <- left_join(calc, expected, by = c("C", "M", ".chain", ".iteration"))

  expect_false(any(is.na(joined$calc)))
  expect_false(any(is.na(joined$expected)))
  expect_equal(joined$calc, joined$expected, tolerance = 1e-8)
})

test_that("alpha_normalised soft correction matches manual bimodal computation", {
  skip_cmdstan()

  test_output_dir <- tempfile("sccomp_test_alpha_norm_bimodal_")
  dir.create(test_output_dir)
  on.exit(unlink(test_output_dir, recursive = TRUE), add = TRUE)

  result <- estimate_for_alpha_normalisation_tests(test_output_dir, bimodal = TRUE)
  fit <- attr(result, "fit")
  model_input <- attr(result, "model_input")

  # Use one variability coefficient across all cell groups to keep the test fast.
  n_m <- ncol(model_input$y)
  alpha_subset <- sprintf("alpha[%d,%d]", 1L, seq_len(n_m))

  calc <- sccomp:::compute_alpha_normalised_draws(
    fit = fit,
    model_input = model_input,
    alpha_variable_subset = alpha_subset
  ) |>
    arrange(C, M, .chain, .iteration) |>
    select(C, M, .chain, .iteration, calc = .value)

  alpha <- sccomp:::draws_to_tibble_x_y(fit, alpha_subset, "C", "M") |>
    rename(alpha = .value)
  beta <- sccomp:::draws_to_tibble_x_y(fit, sprintf("beta[%d,%d]", 1L, seq_len(n_m)), "C", "M") |>
    rename(C_comp = C, beta = .value)
  slope_1 <- sccomp:::draws_to_tibble_x(fit, "prec_slope_1", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_slope_1 = .value)
  slope_2 <- sccomp:::draws_to_tibble_x(fit, "prec_slope_2", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_slope_2 = .value)
  intercept_1 <- sccomp:::draws_to_tibble_x(fit, "prec_intercept_1", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_intercept_1 = .value)
  intercept_2 <- sccomp:::draws_to_tibble_x(fit, "prec_intercept_2", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_intercept_2 = .value)
  prec_sd <- sccomp:::draws_to_tibble_x(fit, "prec_sd", "C") |>
    filter(C == 1L) |>
    transmute(C, .chain, .iteration, prec_sd = .value)
  mix_p <- fit$draws(variables = "mix_p", format = "draws_df") |>
    transmute(
      .chain = as.integer(.chain),
      .iteration = as.integer(.iteration),
      mix_p = mix_p
    )

  expected <- alpha |>
    mutate(C_comp = model_input$variability_to_composition_map[C]) |>
    left_join(beta, by = c("C_comp", "M", ".chain", ".iteration")) |>
    left_join(slope_1, by = c("C", ".chain", ".iteration")) |>
    left_join(slope_2, by = c("C", ".chain", ".iteration")) |>
    left_join(intercept_1, by = c("C", ".chain", ".iteration")) |>
    left_join(intercept_2, by = c("C", ".chain", ".iteration")) |>
    left_join(prec_sd, by = c("C", ".chain", ".iteration")) |>
    left_join(mix_p, by = c(".chain", ".iteration")) |>
    mutate(
      log_1 = log(mix_p) +
        stats::dt((alpha - (beta * prec_slope_1 + prec_intercept_1)) / prec_sd, df = 3, log = TRUE) -
        log(prec_sd),
      log_2 = log1p(-mix_p) +
        stats::dt((alpha - (beta * prec_slope_2 + prec_intercept_2)) / prec_sd, df = 3, log = TRUE) -
        log(prec_sd),
      max_log = pmax(log_1, log_2),
      weight_1 = exp(log_1 - max_log) / (exp(log_1 - max_log) + exp(log_2 - max_log)),
      slope_effective = weight_1 * prec_slope_1 + (1 - weight_1) * prec_slope_2,
      expected = alpha - (beta * slope_effective)
    ) |>
    arrange(C, M, .chain, .iteration) |>
    select(C, M, .chain, .iteration, expected)

  joined <- left_join(calc, expected, by = c("C", "M", ".chain", ".iteration"))

  expect_false(any(is.na(joined$calc)))
  expect_false(any(is.na(joined$expected)))
  expect_equal(joined$calc, joined$expected, tolerance = 1e-8)
})

