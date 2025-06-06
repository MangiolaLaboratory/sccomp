library(testthat)
library(dplyr)
library(tidyr)
library(sccomp)

test_that("replicate_data works correctly", {
  # Load test data
  data("counts_obj")
  
  # Create formulas
  formula_composition = ~ type
  formula_variability = ~ type
  
  # Create .data object
  model_input = 
    counts_obj |> 
      mutate(random_effect = "1") |> 
       sccomp:::data_to_spread(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      .grouping_for_random_effect = "random_effect"
    ) |>
    sccomp:::data_spread_to_model_input(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      truncation_ajustment = 1.1,
      approximate_posterior_inference = FALSE,
      formula_variability = formula_variability,
      contrasts = NULL,
      bimodal_mean_variability_association = FALSE,
      use_data = TRUE,
      random_effect_elements = tibble(factor = character(), grouping = character())
    )
 

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = rlang::quo(sample),
    .cell_group = rlang::quo(cell_group),
    .count = rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = counts_obj
  )

  # Check structure
  expect_type(result, "list")
  expect_true("model_input" %in% names(result))
  expect_true("X_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check dimensions
  expect_equal(nrow(result$model_input$X), nrow(model_input$X))
  expect_equal(ncol(result$model_input$X), ncol(model_input$X))
})

test_that("replicate_data works with random intercept model", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random intercept
  formula_composition = ~ type + (1 | group__)
  formula_variability = ~ type
  
  # Create test data with group information
  test_data = 
    counts_obj |> 
    group_by(sample) |>
    mutate(group__ = sample(c("GROUP1", "GROUP2"), 1)) |>
    ungroup()
  
  # Create .data object
  model_input = 
    test_data |> 
    sccomp:::data_to_spread(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      .grouping_for_random_effect = "group__"
    ) |>
    sccomp:::data_spread_to_model_input(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      truncation_ajustment = 1.1,
      approximate_posterior_inference = FALSE,
      formula_variability = formula_variability,
      contrasts = NULL,
      bimodal_mean_variability_association = FALSE,
      use_data = TRUE,
      random_effect_elements = tibble(factor = "(Intercept)", grouping = "group__")
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = rlang::quo(sample),
    .cell_group = rlang::quo(cell_group),
    .count = rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data
  )

  # Check structure
  expect_type(result, "list")
  expect_true("model_input" %in% names(result))
  expect_true("X_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  expect_true(!is.null(result$model_input$X_random_effect))
  expect_true(ncol(result$model_input$X_random_effect) > 0)
})

test_that("replicate_data works with random slope model", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random slope
  formula_composition = ~ type + (type | group__)
  formula_variability = ~ type
  
  # Create test data with group information
  test_data = 
    counts_obj |> 
    group_by(sample) |>
    mutate(
      group__ = sample(c("GROUP1", "GROUP2"), 1),
      group2__ = paste0(group__, "_", sample(1:2, 1))
    ) |>
    ungroup()
  
  # Create .data object
  model_input = 
    test_data |> 
    sccomp:::data_to_spread(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      .grouping_for_random_effect = "group__"
    ) |>
    sccomp:::data_spread_to_model_input(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      truncation_ajustment = 1.1,
      approximate_posterior_inference = FALSE,
      formula_variability = formula_variability,
      contrasts = NULL,
      bimodal_mean_variability_association = FALSE,
      use_data = TRUE,
      random_effect_elements = tibble(
        factor = c("(Intercept)", "type"),
        grouping = c("group__", "group__")
      )
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = rlang::quo(sample),
    .cell_group = rlang::quo(cell_group),
    .count = rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data
  )

  # Check structure
  expect_type(result, "list")
  expect_true("model_input" %in% names(result))
  expect_true("X_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  expect_true(!is.null(result$model_input$X_random_effect))
  expect_true(ncol(result$model_input$X_random_effect) > 0)
})

test_that("replicate_data works with nested random effects", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with nested random effects
  formula_composition = ~ type + (type | group__) + (1 | group2__)
  formula_variability = ~ type
  
  # Create test data with nested group information
  test_data = 
    counts_obj |> 
    group_by(sample) |>
    mutate(
      group__ = sample(c("GROUP1", "GROUP2"), 1),
      group2__ = paste0(group__, "_", sample(1:2, 1))
    ) |>
    ungroup()
  
  # Create .data object
  model_input = 
    test_data |> 
    sccomp:::data_to_spread(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      .grouping_for_random_effect = c("group__", "group2__")
    ) |>
    sccomp:::data_spread_to_model_input(
      formula = formula_composition,
      .sample = !!quo(sample),
      .cell_type = !!quo(cell_group),
      .count = !!quo(count),
      truncation_ajustment = 1.1,
      approximate_posterior_inference = FALSE,
      formula_variability = formula_variability,
      contrasts = NULL,
      bimodal_mean_variability_association = FALSE,
      use_data = TRUE,
      random_effect_elements = tibble(
        factor = c("(Intercept)", "type", "(Intercept)"),
        grouping = c("group__", "group__", "group2__")
      )
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = rlang::quo(sample),
    .cell_group = rlang::quo(cell_group),
    .count = rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data
  )

  # Check structure
  expect_type(result, "list")
  expect_true("model_input" %in% names(result))
  expect_true("X_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrices
  expect_true(!is.null(result$model_input$X_random_effect))
  expect_true(!is.null(result$model_input$X_random_effect_2))
  expect_true(ncol(result$model_input$X_random_effect) > 0)
  expect_true(ncol(result$model_input$X_random_effect_2) > 0)
}) 
