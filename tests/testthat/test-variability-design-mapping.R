library(testthat)
library(dplyr)
library(sccomp)

test_that("get_variability_to_composition_map matches by column name", {
  X <- matrix(0, nrow = 2, ncol = 4)
  colnames(X) <- c("(Intercept)", "typehealthy", "phenotypep2", "phenotypep3")

  Xa <- matrix(0, nrow = 2, ncol = 4)
  colnames(Xa) <- c("(Intercept)", "phenotypep2", "phenotypep3", "typehealthy")

  expect_equal(
    sccomp:::get_variability_to_composition_map(X, Xa),
    c(1L, 3L, 4L, 2L)
  )
})

test_that("get_variability_to_composition_map errors on missing terms", {
  X <- matrix(0, nrow = 2, ncol = 2)
  colnames(X) <- c("(Intercept)", "typehealthy")

  Xa <- matrix(0, nrow = 2, ncol = 3)
  colnames(Xa) <- c("(Intercept)", "typehealthy", "phenotypep2")

  expect_error(
    sccomp:::get_variability_to_composition_map(X, Xa),
    "every variability design term must also be present in the composition design matrix"
  )

  expect_error(
    sccomp:::get_variability_to_composition_map(X, Xa),
    "Missing terms: phenotypep2"
  )
})

test_that("variability design maps to composition design by name", {
  test_counts <- tibble::tibble(
    sample = rep(c("s1", "s2", "s3", "s4"), each = 2),
    type = rep(c("healthy", "healthy", "cancer", "cancer"), each = 2),
    phenotype = rep(c("p1", "p2", "p1", "p2"), each = 2),
    cell_group = rep(c("cg1", "cg2"), times = 4),
    count = c(10L, 20L, 12L, 18L, 22L, 9L, 19L, 13L)
  )

  formula_composition <- ~ type + phenotype
  formula_variability <- ~ phenotype + type

  model_input <-
    test_counts |>
    mutate(random_effect = "1") |>
    sccomp:::data_to_spread(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_group = !!quo(cell_group),
      .count = !!quo(count),
      .grouping_for_random_effect = "random_effect"
    ) |>
    sccomp:::data_spread_to_model_input(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_group = !!quo(cell_group),
      .count = !!quo(count),
      truncation_ajustment = 1.1,
      approximate_posterior_inference = FALSE,
      formula_variability = formula_variability,
      contrasts = NULL,
      bimodal_mean_variability_association = FALSE,
      use_data = TRUE,
      random_effect_elements = tibble(factor = character(), grouping = character())
    )

  expected_map <- match(colnames(model_input$Xa), colnames(model_input$X))
  expect_equal(model_input$variability_to_composition_map, expected_map)
})

test_that("variability design terms must exist in composition design", {
  test_counts <- tibble::tibble(
    sample = rep(c("s1", "s2", "s3", "s4"), each = 2),
    type = rep(c("healthy", "healthy", "cancer", "cancer"), each = 2),
    phenotype = rep(c("p1", "p2", "p1", "p2"), each = 2),
    cell_group = rep(c("cg1", "cg2"), times = 4),
    count = c(10L, 20L, 12L, 18L, 22L, 9L, 19L, 13L)
  )

  formula_composition <- ~ type + phenotype
  formula_variability <- ~ type * phenotype

  expect_error(
    test_counts |>
      mutate(random_effect = "1") |>
      sccomp:::data_to_spread(
        formula = formula_composition,
        .sample = !!quo(sample),
        .cell_group = !!quo(cell_group),
        .count = !!quo(count),
        .grouping_for_random_effect = "random_effect"
      ) |>
      sccomp:::data_spread_to_model_input(
        formula = formula_composition,
        .sample = !!quo(sample),
        .cell_group = !!quo(cell_group),
        .count = !!quo(count),
        truncation_ajustment = 1.1,
        approximate_posterior_inference = FALSE,
        formula_variability = formula_variability,
        contrasts = NULL,
        bimodal_mean_variability_association = FALSE,
        use_data = TRUE,
        random_effect_elements = tibble(factor = character(), grouping = character())
      ),
    "every variability design term must also be present in the composition design matrix"
  )
})
