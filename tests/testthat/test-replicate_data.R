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
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = counts_obj |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check dimensions
  print('nrow(result$X):'); print(nrow(result$X))
  print('nrow(model_input$X):'); print(nrow(model_input$X))
  expect_equal(nrow(result$X), nrow(model_input$X))
  print('ncol(result$X):'); print(ncol(result$X))
  print('ncol(model_input$X):'); print(ncol(model_input$X))
  expect_equal(ncol(result$X), ncol(model_input$X))

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
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
      random_effect_elements = tibble(factor = "(Intercept)", grouping = "group__"),
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"model_input" %in% names(result):'); print("model_input" %in% names(result))
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))
  print('ncol(result$X_random_effect):'); print(ncol(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
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
      ),
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"model_input" %in% names(result):'); print("model_input" %in% names(result))
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))
  print('ncol(result$X_random_effect):'); print(ncol(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
})

test_that("replicate_data works with random slope model and type NA", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random slope
  formula_composition = ~ type + (type | group__)
  formula_variability = ~ type

  # test data 
  test_data = 
    counts_obj |> 
    left_join(
      tibble(
        sample = unique(counts_obj$sample),
        group__ = case_when(
          sample %in% unique(sample)[1:10] ~ "GROUP1",
          TRUE ~ "GROUP2"
        )
      )
    )


  # Create test data with group information and type NA
  

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
      ),
      accept_NA_as_average_effect = TRUE
    )



new_data = 
    test_data |> 
    sccomp:::.subset(!!quo(sample)) |> 
    mutate(type = NA_character_)



  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = new_data,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"model_input" %in% names(result):'); print("model_input" %in% names(result))
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))

  # Find a column with 'type' in its name
  col_type <- grep("type", colnames(result$X_random_effect), value = TRUE)[1]
  col_intercept <- grep("Intercept", colnames(result$X_random_effect), value = TRUE)[1]
  print(paste('result$X_random_effect[1,', col_type, ']:')); print(result$X_random_effect[1, col_type])
  print(paste('result$X_random_effect[1,', col_intercept, ']:')); print(result$X_random_effect[1, col_intercept])
  expect_equal(result$X_random_effect[1, col_type], 0.5)
  expect_equal(result$X_random_effect[1, col_intercept], 1)

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
      ),
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"model_input" %in% names(result):'); print("model_input" %in% names(result))
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrices
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))
  print('is.null(result$X_random_effect_2):'); print(is.null(result$X_random_effect_2))
  expect_true(!is.null(result$X_random_effect_2))
  print('ncol(result$X_random_effect):'); print(ncol(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)
  print('ncol(result$X_random_effect_2):'); print(ncol(result$X_random_effect_2))
  expect_true(ncol(result$X_random_effect_2) > 0)

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
})

test_that("replicate_data works with NA values in grouping", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random intercept
  formula_composition = ~ type + (1 | group__)
  formula_variability = ~ type
  
  # Pick one sample to have NA group__
  unique_samples <- unique(counts_obj$sample)
  sample_with_na <- unique_samples[1]
  
  # Create test data with group information including NA values
  test_data = 
    counts_obj |> 
    mutate(
      group__ = if_else(
        sample == sample_with_na,
        NA_character_,
        sample(c("GROUP1", "GROUP2"), 1)
      ),
      count = if_else(is.na(group__), count, count)
    )

  # Assert only one unique sample has NA in group__
  expect_equal(
    test_data |> filter(is.na(group__)) |> distinct(sample) |> nrow(),
    1
  )
  
  # Use parse_formula_random_effect for random_effect_elements
  random_effect_elements = sccomp:::parse_formula_random_effect(formula_composition)

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
      random_effect_elements = random_effect_elements,
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"model_input" %in% names(result):'); print("model_input" %in% names(result))
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))
  print('ncol(result$X_random_effect):'); print(ncol(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)
  
  # Check that NA values are handled correctly
  print('all(!is.na(result$X_random_effect))'); print(all(!is.na(result$X_random_effect)))
  expect_true(all(!is.na(result$X_random_effect)))
  
  # Check that the design matrix dimensions are correct
  print('nrow(result$X_random_effect):'); print(nrow(result$X_random_effect))
  print('nrow(distinct(test_data, sample))'); print(nrow(test_data |> distinct(sample)))
  expect_equal(nrow(result$X_random_effect), nrow(test_data |> distinct(sample)))

  # Check that X_random_effect has only one 0 for the expected sample
  group1_samples <- test_data |> filter(group__ == "GROUP1") |> distinct(sample) |> pull(sample)
  if ("(Intercept)___GROUP1" %in% colnames(result$X_random_effect)) {
    print('result$X_random_effect[group1_samples, "(Intercept)___GROUP1"]:')
    print(result$X_random_effect[group1_samples, "(Intercept)___GROUP1"])
    expect_true(all(result$X_random_effect[group1_samples, "(Intercept)___GROUP1"] == 1))
  }
  
  # Check that X_random_effect_unseen has the expected number of 1s
  print('result$X_random_effect_unseen[rownames(result$X_random_effect_unseen) == sample_with_na, "(Intercept)___NA"]:')
  print(result$X_random_effect_unseen[rownames(result$X_random_effect_unseen) == sample_with_na, "(Intercept)___NA"])
  expect_equal(result$X_random_effect_unseen[rownames(result$X_random_effect_unseen) == sample_with_na, "(Intercept)___NA"],
    1
  ) 

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
})

test_that("replicate_data works with NA values in grouping and random effects", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random intercept and slope
  formula_composition = ~ type + (type | group__)
  formula_variability = ~ type
  
  # Pick one sample to have NA group__
  unique_samples <- unique(counts_obj$sample)
  sample_with_na <- unique_samples[c(1, 9)]
  
  # Create test data with group information including NA values and random effects
  test_data = 
    counts_obj |> 
    group_by(sample) |>
    mutate(
      group__ = case_when(
        sample == sample_with_na[1] ~ NA_character_,
        sample == sample_with_na[2] ~ NA_character_,
        sample %in% unique_samples[3:10] ~ "GROUP1",
        TRUE ~ "GROUP2"
      ),
      count = if_else(is.na(group__), count, count)
    ) |>
    ungroup()
  
  # Verify we have all three groups
  expect_equal(
    test_data |> distinct(group__) |> pull(group__) |> sort(),
    c("GROUP1", "GROUP2", NA_character_) |> sort()
  )
  
  # Use parse_formula_random_effect for random_effect_elements
  random_effect_elements = sccomp:::parse_formula_random_effect(formula_composition)

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
      random_effect_elements = random_effect_elements,
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = NULL,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  print('"X_which" %in% names(result):'); print("X_which" %in% names(result))
  expect_true("X_which" %in% names(result))
  print('"XA_which" %in% names(result):'); print("XA_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  print('"create_intercept" %in% names(result):'); print("create_intercept" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  print('is.null(result$X_random_effect):'); print(is.null(result$X_random_effect))
  expect_true(!is.null(result$X_random_effect))
  print('ncol(result$X_random_effect):'); print(ncol(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)
  
  # Check that NA values are handled correctly in random effects
  print('all(!is.na(result$X_random_effect))'); print(all(!is.na(result$X_random_effect)))
  expect_true(all(!is.na(result$X_random_effect)))
  
  # Check that the design matrix dimensions are correct
  print('nrow(result$X_random_effect):'); print(nrow(result$X_random_effect))
  print('nrow(distinct(test_data, sample))'); print(nrow(test_data |> distinct(sample)))
  expect_equal(nrow(result$X_random_effect), nrow(test_data |> distinct(sample)))
  
  # Check that the random effect design matrix has the correct structure
  print('any(grepl("type", colnames(result$X_random_effect)))'); print(any(grepl("type", colnames(result$X_random_effect))))
  expect_true(any(grepl("type", colnames(result$X_random_effect))))
  print('any(grepl("Intercept", colnames(result$X_random_effect)))'); print(any(grepl("Intercept", colnames(result$X_random_effect))))
  expect_true(any(grepl("Intercept", colnames(result$X_random_effect))))
  
  # Test the properties of X_random_effect_unseen
  print('all(!is.na(result$X_random_effect_unseen))'); print(all(!is.na(result$X_random_effect_unseen)))
  expect_true(all(!is.na(result$X_random_effect_unseen)))
  print('nrow(result$X_random_effect_unseen):'); print(nrow(result$X_random_effect_unseen))
  print('nrow(distinct(test_data, sample)):'); print(nrow(test_data |> distinct(sample)))
  expect_equal(nrow(result$X_random_effect_unseen), nrow(test_data |> distinct(sample)))
  print('ncol(result$X_random_effect_unseen):'); print(ncol(result$X_random_effect_unseen))
  expect_equal(ncol(result$X_random_effect_unseen), 2)  # Should have (Intercept)___NA and typecancer___NA
  print('rownames(result$X_random_effect_unseen):'); print(rownames(result$X_random_effect_unseen))
  print('rownames(result$X_random_effect):'); print(rownames(result$X_random_effect))
  expect_equal(rownames(result$X_random_effect_unseen), rownames(result$X_random_effect))
  print('all(grepl("___NA$", colnames(result$X_random_effect_unseen)))'); print(all(grepl("___NA$", colnames(result$X_random_effect_unseen))))
  expect_true(all(grepl("___NA$", colnames(result$X_random_effect_unseen))))
  print('all(grepl("___GROUP", colnames(result$X_random_effect)))'); print(all(grepl("___GROUP", colnames(result$X_random_effect))))
  expect_true(all(grepl("___GROUP", colnames(result$X_random_effect))))
  
  # Check that the sums of the NA columns are correct
  print('sum(result$X_random_effect_unseen[, "typecancer___NA"]):'); print(sum(as.vector(result$X_random_effect_unseen[, "typecancer___NA"])))
  expect_equal(
    as.vector(result$X_random_effect_unseen[, "typecancer___NA"]) |> sum(),
    1
  )
  
  print('sum(result$X_random_effect_unseen[, "(Intercept)___NA"]):'); print(sum(as.vector(result$X_random_effect_unseen[, "(Intercept)___NA"])))
  expect_equal(
    as.vector(result$X_random_effect_unseen[, "(Intercept)___NA"]) |> sum(),
    2
  )

  # Robust per-sample checks for random effect matrices
  for (s in rownames(result$X_random_effect)) {
    group_val <- test_data |> filter(sample == s) |> distinct(group__) |> pull(group__)
    if (!is.na(group_val)) {
      colname <- paste0("(Intercept)___", group_val)
      if (colname %in% colnames(result$X_random_effect)) {
        expect_equal(result$X_random_effect[s, colname], 1)
      }
    }
  }
  if (!is.null(result$X_random_effect_unseen)) {
    for (s in rownames(result$X_random_effect_unseen)) {
      group_val <- test_data |> filter(sample == s) |> distinct(group__) |> pull(group__)
      if (is.na(group_val)) {
        expect_equal(result$X_random_effect_unseen[s, "(Intercept)___NA"], 1)
      }
    }
  }

  # Print key objects for debugging
  print("colnames(result$X_random_effect):")
  print(colnames(result$X_random_effect))
  print("rownames(result$X_random_effect):")
  print(rownames(result$X_random_effect))
  print("result$X_random_effect:")
  print(result$X_random_effect)
  print("colnames(result$X_random_effect_unseen):")
  print(colnames(result$X_random_effect_unseen))
  print("result$X_random_effect_unseen:")
  print(result$X_random_effect_unseen)
}) 

test_that("replicate_data works with type NA and group__ NA", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random intercept and slope
  formula_composition = ~ type + (type | group__)
  formula_variability = ~ type

  # Create original data with group information and type NA
  test_data = 
    counts_obj |> 
    left_join(
      tibble(
        sample = unique(counts_obj$sample),
        group__ = case_when(
          sample %in% unique(sample)[1:10] ~ "GROUP1",
          TRUE ~ "GROUP2"
        )
      )
    )   

  # New data is the same as test data but with type = NA and group__ = NA
  new_data = 
    test_data |>
    sccomp:::.subset(!!quo(sample)) |>
    mutate(
      type = NA_character_,
      group__ = NA_character_,
      sample = paste0("new_", sample)  # Give new samples unique names
    )

  # Use parse_formula_random_effect for random_effect_elements
  random_effect_elements = sccomp:::parse_formula_random_effect(formula_composition)

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
      random_effect_elements = random_effect_elements,
      accept_NA_as_average_effect = TRUE
    )

  # Prepare replicate data
  result = sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = new_data,
    original_count_data = test_data |> sccomp:::.subset(!!quo(sample))
  )

  # Check structure
  print("names(result):"); print(names(result))
  expect_type(result, "list")
  expect_true("X_which" %in% names(result))
  expect_true("XA_which" %in% names(result))
  expect_true("create_intercept" %in% names(result))
  
  # Check random effect design matrix
  expect_true(!is.null(result$X_random_effect))
  expect_true(ncol(result$X_random_effect) > 0)
  
  # Check that NA values are handled correctly
  expect_true(all(!is.na(result$X_random_effect)))
  
  # Check that the design matrix dimensions are correct
  expect_equal(nrow(result$X_random_effect), nrow(new_data |> distinct(sample)))
  
  # Check that X_random_effect_unseen has the expected structure
  expect_true(!is.null(result$X_random_effect_unseen))
  expect_true(all(!is.na(result$X_random_effect_unseen)))
  expect_equal(nrow(result$X_random_effect_unseen), nrow(new_data |> distinct(sample)))
  expect_true(all(grepl("___NA$", colnames(result$X_random_effect_unseen))))
}) 

test_that("replicate_data works with new data containing only NA groups", {
  # Load test data
  data("counts_obj")
  
  # Create formulas with random intercept and slope
  formula_composition = ~ type + (type | group__)
  formula_variability = ~ type
  
  # Create original data with group information
  original_data = 
    counts_obj |> 
    left_join(
      tibble(
        sample = unique(counts_obj$sample),
        group__ = case_when(
          sample %in% unique(sample)[1:10] ~ "GROUP1",
          TRUE ~ "GROUP2"
        )
      )
    ) 

  # Create sccomp estimate object
  estimate = sccomp_estimate(
    original_data,
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    .sample = sample,
    .cell_group = cell_group,
    .abundance = count,
    cores = 1,
    max_sampling_iterations = 1000,
    verbose = FALSE
  )

  # Create new data with only NA groups
  new_data = 
    counts_obj |>
    sccomp:::.subset(!!quo(sample)) |> 
    slice(1, 10) |> 
    mutate(
      group__ = NA_character_,  # Set all groups to NA
      sample = paste0("new_", sample)  # Give them new sample names
    )
  
  
  # Test
  result = sccomp:::replicate_data(estimate,
                formula_composition = formula_composition,
                formula_variability = formula_variability,
                new_data = new_data,
                number_of_draws = 1,
                mcmc_seed = sample(1e5, 1),
                cores = 1)  

  print(result)
  print(result$summary())
  expect_s3_class(result, "CmdStanGQ")
  expect_true(any(grepl("counts_uncorrected", result$summary()$variable)))
}) 

test_that("prepare_replicate_data handles complex design with NAs and prints new X and random effect design matrices", {
  # Create test data with complex factors and NAs
  age_raw <- c(25, 20, 35, 40, 45)
  age_scaled <- as.numeric(scale(age_raw, center = TRUE, scale = TRUE))

  # Test data without NAs
  test_data <- tibble(
    sample = c("sample1", "sample1", "sample2", "sample2", "sample3", "sample3", "sample4", "sample4", "sample5", "sample5"),
    type = c("healthy", "healthy", "healthy", "healthy", "cancer", "cancer", "cancer", "cancer", "healthy", "healthy"),
    age = rep(age_scaled, each = 2),
    sex = c("M", "M", "F", "F", "M", "M", "F", "F", "M", "M"),
    group__ = c("GROUP1", "GROUP1", "GROUP2", "GROUP2", "GROUP1", "GROUP1", "GROUP2", "GROUP2", "GROUP1", "GROUP1"),
    count = c(10, 12, 15, 17, 20, 22, 25, 27, 30, 32),
    cell_group = c("A", "B", "A", "B", "A", "B", "A", "B", "A", "B")
  )

  # New data - only sample-related columns
  new_data = tibble(
    sample = c("sample1", "sample2", "sample3", "sample4", "sample5"),
    type = c("healthy", "healthy", "cancer", NA, "healthy"),
    age = c(age_scaled[1], age_scaled[2], NA, age_scaled[4], age_scaled[5]),
    sex = c("M", "F", NA, "M", "F"),
    group__ = c("GROUP1", "GROUP2", "GROUP1", "GROUP2", "GROUP1")
  )

  # Define formulas
  formula_composition <- ~ type * age * sex + (1 + age | group__)
  formula_variability <- ~ type

  # Use parse_formula_random_effect for random_effect_elements
  random_effect_elements <- sccomp:::parse_formula_random_effect(formula_composition)
  
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
      random_effect_elements = random_effect_elements
    )
  
  # Create original_count_data using .subset()
  original_count_data = test_data |>
    sccomp:::.subset(!!quo(sample)) |>
    distinct(sample, type, age, sex, group__)
  
  # Prepare replicate data
  result <- sccomp:::prepare_replicate_data(
    X = model_input$X,
    Xa = model_input$Xa,
    N = model_input$N,
    intercept_in_design = model_input$intercept_in_design,
    X_random_effect = model_input$X_random_effect,
    X_random_effect_2 = model_input$X_random_effect_2,
    .sample = !!rlang::quo(sample),
    .cell_group = !!rlang::quo(cell_group),
    .count = !!rlang::quo(count),
    formula_composition = formula_composition,
    formula_variability = formula_variability,
    new_data = new_data,
    original_count_data = original_count_data
  )

  # Print the new fixed effect design matrix (X)
  cat("\nNew fixed effect design matrix (X):\n")
  print(result$X)

  # Print the new random effect design matrix (X_random_effect)
  cat("\nNew random effect design matrix (X_random_effect):\n")
  print(result$X_random_effect)

  # Print the new random effect design matrix 2 (if present)
  if (!is.null(result$X_random_effect_2)) {
    cat("\nNew random effect design matrix 2 (X_random_effect_2):\n")
    print(result$X_random_effect_2)
  }

  # Evaluate structure
  expect_true(is.matrix(result$X) || is.data.frame(result$X))
  expect_true(is.matrix(result$X_random_effect) || is.data.frame(result$X_random_effect))
  expect_equal(nrow(result$X), nrow(original_count_data))
  expect_equal(nrow(result$X_random_effect), nrow(original_count_data))
  # Check that NAs are handled (no NA in design matrices)
  expect_true(all(!is.na(result$X)))
  expect_true(all(!is.na(result$X_random_effect)))
})

test_that("prepare_replicate_data throws error for duplicate sample names", {
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

  # Test case 1: Duplicate samples in original_count_data
  original_count_data_duplicates = 
    counts_obj |>
    sccomp:::.subset(!!quo(sample)) 
  
  original_count_data_duplicates = 
    original_count_data_duplicates |>
    bind_rows(
      original_count_data_duplicates |> filter(sample == first(sample))  # Duplicate first sample
    ) 

  expect_error(
    sccomp:::prepare_replicate_data(
      X = model_input$X,
      Xa = model_input$Xa,
      N = model_input$N,
      intercept_in_design = model_input$intercept_in_design,
      X_random_effect = model_input$X_random_effect,
      X_random_effect_2 = model_input$X_random_effect_2,
      .sample = !!rlang::quo(sample),
      .cell_group = !!rlang::quo(cell_group),
      .count = !!rlang::quo(count),
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      new_data = NULL,
      original_count_data = original_count_data_duplicates
    ),
    "sccomp says: your original_count_data has duplicated sample names"
  )

  # Test case 2: Duplicate samples in new_data
  new_data_duplicates = 
    counts_obj |>
    sccomp:::.subset(!!quo(sample)) |> 
    bind_rows(
      counts_obj |> sccomp:::.subset(!!quo(sample)) |> filter(sample == first(sample))  # Duplicate first sample
    ) 

  expect_error(
    sccomp:::prepare_replicate_data(
      X = model_input$X,
      Xa = model_input$Xa,
      N = model_input$N,
      intercept_in_design = model_input$intercept_in_design,
      X_random_effect = model_input$X_random_effect,
      X_random_effect_2 = model_input$X_random_effect_2,
      .sample = !!rlang::quo(sample),
      .cell_group = !!rlang::quo(cell_group),
      .count = !!rlang::quo(count),
      formula_composition = formula_composition,
      formula_variability = formula_variability,
      new_data = new_data_duplicates,
      original_count_data = counts_obj |> sccomp:::.subset(!!quo(sample))
    ),
    "sccomp says: your new_data has duplicated sample names"
  )
})
