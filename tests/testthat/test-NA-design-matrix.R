test_that("test NA factors with interactions",{
  
  
  library(dplyr)
 
  n_iterations = 1000
  
  formula = ~ group2__ * type
  print("Formula:")
  print(formula)
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = n_iterations,
      inference_method = "hmc", verbose=FALSE
    )
  
  print("\nOriginal count data structure:")
  print(str(attr(res, "count_data")))
  
  # Create new data using one sample and setting continuous covariate to NA
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(type = NA)
  
  print("\nNew data with NA type:")
  print(new_data)
  
  print("\nCombined data structure:")
  combined_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data)
  print(str(combined_data))
  
  print("\nDesign matrix with NA handling:")
  design_matrix = 
    combined_data |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  print(head(design_matrix))
  print("\nLast row of design matrix:")
  print(tail(design_matrix, 1))
  
  res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample) |>
    tail(1) |>
    as.numeric() |>
    
    # Expect
    #  (Intercept) group2__GROUP22 typehealthy group2__GROUP22:typehealthy
    #          1               0         0.5                           0
    expect_equal(
      c(1, 0, 0.5, 0)
    )
})

test_that("test NA in both group2__ and type", {
  
  library(dplyr)
  
  
  
  formula = ~ group2__ * type
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in both group2__ and type
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(
      type = NA,
      group2__ = NA
    )
  
  print("\nNew data with NA in both group2__ and type:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in both variables:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NAs)
  tail(design_matrix, 1) |>
    as.numeric() |>
    expect_equal(
      c(1, 0.5, 0.5, 0.25)  # Intercept, group2__ effect split, type effect split, interaction split
    )
})

test_that("test NA in three-way interaction", {
  
  library(dplyr)
  
  
  
  # Create a more complex formula with three-way interaction
  formula = ~ group2__ * type * continuous_covariate
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in one variable
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(type = NA)
  
  print("\nNew data with NA in type for three-way interaction:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in three-way interaction:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NA)
  last_row = tail(design_matrix, 1) |> as.numeric()
  
  # Check intercept
  expect_equal(last_row[1], 1)
  
  # Check group2__ effect (should be 0 since it's GROUP21)
  expect_equal(last_row[2], 0)
  
  # Check type effect (should be 0.5 since it's NA)
  expect_equal(last_row[3], 0.5)
  
  # Check continuous_covariate effect (should be preserved)
  expect_equal(last_row[4], new_data$continuous_covariate[1])
  
  # Check two-way interactions
  expect_equal(last_row[5], 0)  # group2__:type
  expect_equal(last_row[6], 0)  # group2__:continuous_covariate
  expect_equal(last_row[7], 0.5 * last_row[4])  # type:continuous_covariate
  
  # Check three-way interaction
  expect_equal(last_row[8], 0)  # group2__:type:continuous_covariate
})

test_that("test NA in nested effects", {
  
  library(dplyr)
  
  
  
  # Create formula with nested effects
  formula = ~ type + group2__/type
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in nested effect
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SCP424_pbmc2") |>
    mutate(type = NA)
  
  print("\nNew data with NA in nested effect:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in nested effect:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NA)
  last_row = tail(design_matrix, 1) |> as.numeric()
  
  # Check intercept
  expect_equal(last_row[1], 1)
  
  # Check main effect
  expect_equal(last_row[2], 0.5)  # type effect split
  
  # Check nested effects
  expect_equal(last_row[3:length(last_row)], c(1, 0.5))  # nested effects split
})

test_that("test NA in three-way interaction with one factor NA", {
  
  library(dplyr)
  
  
  
  # Create a more complex formula with three-way interaction
  formula = ~ group2__ * type * continuous_covariate
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in one variable
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(type = NA)
  
  print("\nNew data with NA in type for three-way interaction:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in three-way interaction:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NA)
  last_row = tail(design_matrix, 1) |> as.numeric()
  
  # Check intercept
  expect_equal(last_row[1], 1)
  
  # Check group2__ effect (should be 0 since it's GROUP21)
  expect_equal(last_row[2], 0)
  
  # Check type effect (should be 0.5 since it's NA)
  expect_equal(last_row[3], 0.5)
  
  # Check continuous_covariate effect (should be preserved)
  expect_equal(last_row[4], new_data$continuous_covariate[1])
  
  # Check two-way interactions
  expect_equal(last_row[5], 0)  # group2__:type
  expect_equal(last_row[6], 0)  # group2__:continuous_covariate
  expect_equal(last_row[7], 0.5 * last_row[4])  # type:continuous_covariate
  
  # Check three-way interaction
  expect_equal(last_row[8], 0)  # group2__:type:continuous_covariate
})

test_that("test NA in three-way interaction with two factors NA", {
  
  library(dplyr)
  
  
  
  # Create a more complex formula with three-way interaction
  formula = ~ group2__ * type * continuous_covariate
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in two variables
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(type = NA, group2__ = NA)
  
  print("\nNew data with NA in type and group2__ for three-way interaction:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in three-way interaction:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NA)
  last_row = tail(design_matrix, 1) |> as.numeric()
  
  # Check intercept
  expect_equal(last_row[1], 1)
  
  # Check group2__ effect (should be 0.5 since it's NA)
  expect_equal(last_row[2], 0.5)
  
  # Check type effect (should be 0.5 since it's NA)
  expect_equal(last_row[3], 0.5)
  
  # Check continuous_covariate effect (should be preserved)
  expect_equal(last_row[4], new_data$continuous_covariate[1])
  
  # Check two-way interactions
  expect_equal(last_row[5], 0.25)  # group2__:type
  expect_equal(last_row[6], 0.5 * last_row[4])  # group2__:continuous_covariate
  expect_equal(last_row[7], 0.5 * last_row[4])  # type:continuous_covariate
  
  # Check three-way interaction
  expect_equal(last_row[8], 0.25 * last_row[4])  # group2__:type:continuous_covariate
})

test_that("test NA in three-way interaction with three factors NA", {
  
  library(dplyr)
  
  # Create a more complex formula with three-way interaction
  formula = ~ group2__ * type * continuous_covariate
  
  res =
    seurat_obj |>
    sccomp_estimate(
      formula_composition = formula,
      formula_variability = ~ 1,
      "sample", "cell_group",
      cores = 1,
      mcmc_seed = 42,   
      max_sampling_iterations = 1000,
      inference_method = "hmc", verbose=FALSE
    )
  
  # Create new data with NA in all three variables
  new_data = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample)  |>
    filter(sample == "SI-GA-H1") |>
    mutate(type = NA, group2__ = NA, continuous_covariate = NA)
  
  print("\nNew data with NA in type, group2__, and continuous_covariate for three-way interaction:")
  print(new_data)
  
  design_matrix = 
    res |> 
    attr("count_data") |> 
    sccomp:::.subset(sample) |> 
    bind_rows(new_data) |> 
    sccomp:::get_design_matrix_with_na_handling(formula, sample)
  
  print("\nDesign matrix with NA in three-way interaction:")
  print(tail(design_matrix, 1))
  
  # Test the last row (with NA)
  last_row = tail(design_matrix, 1) |> as.numeric()
  
  # Check intercept
  expect_equal(last_row[1], 1)
  
  # Check group2__ effect (should be 0.5 since it's NA)
  expect_equal(last_row[2], 0.5)
  
  # Check type effect (should be 0.5 since it's NA)
  expect_equal(last_row[3], 0.5)
  
  # Check continuous_covariate effect (should be the mean of the continuous covariate)
  expect_equal(last_row[4], mean(res |> 
                                   attr("count_data") |> 
                                   sccomp:::.subset(sample) |> 
                                   pull(continuous_covariate), na.rm = TRUE))
  
  # Check two-way interactions
  expect_equal(last_row[5],  0.25)  # group2__:type
  expect_equal(last_row[6], mean(res |> 
                                   attr("count_data") |> 
                                   sccomp:::.subset(sample) |> 
                                   pull(continuous_covariate), na.rm = TRUE) * 0.5)  # group2__:continuous_covariate
  expect_equal(last_row[7], mean(res |> 
                                   attr("count_data") |> 
                                   sccomp:::.subset(sample) |> 
                                   pull(continuous_covariate), na.rm = TRUE) * 0.5)  # type:continuous_covariate
  
  # Check three-way interaction
  expect_equal(last_row[8], mean(res |> 
                                   attr("count_data") |> 
                                   sccomp:::.subset(sample) |> 
                                   pull(continuous_covariate), na.rm = TRUE) * 0.25)  # group2__:type:continuous_covariate
})
