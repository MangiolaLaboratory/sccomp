library(testthat)
library(sccomp)

test_that("covariates are z-scored against the reference samples", {
  training <- tibble::tibble(
    sample = paste0("S", 1:5),
    age = c(10, 20, 30, 40, 50),
    condition = c("a", "b", "a", "b", "a")
  )
  grid <- tibble::tibble(sample = "g1", age = 50, condition = "a")
  
  scaled <- grid |>
    sccomp:::scale_numeric_covariates(
      c("age", "condition"),
      reference = training
    )
  
  expect_equal(scaled$age, (50 - 30) / sd(c(10, 20, 30, 40, 50)))
  expect_equal(scaled$condition, "a")
})

test_that("a covariate constant across samples is not divided by zero", {
  training <- tibble::tibble(sample = paste0("S", 1:3), age = c(40, 40, 40))
  
  scaled <- training |> sccomp:::scale_numeric_covariates("age")
  
  expect_equal(scaled$age, c(0, 0, 0))
})

test_that("design matrix columns do not depend on the prediction grid", {
  training <- tibble::tibble(
    sample = paste0("S", 1:5),
    age = c(10, 20, 30, 40, 50)
  )
  
  # Two grids that share the age of interest but differ in range and density
  narrow_grid <- tibble::tibble(sample = c("N1", "N2"), age = c(35, 45))
  wide_grid <- tibble::tibble(
    sample = paste0("W", 1:4),
    age = c(0, 35, 45, 200)
  )
  
  column_at_35 <- function(grid) {
    design <- training |>
      dplyr::bind_rows(grid) |>
      sccomp:::get_design_matrix(~ age, sample, scaling_reference = training)
    design[grid$sample[grid$age == 35], "age"]
  }
  
  expect_equal(column_at_35(narrow_grid), column_at_35(wide_grid))
  expect_equal(
    unname(column_at_35(narrow_grid)),
    (35 - 30) / sd(c(10, 20, 30, 40, 50))
  )
})

test_that("design matrix scaling is unchanged when no reference is supplied", {
  training <- tibble::tibble(
    sample = paste0("S", 1:5),
    age = c(10, 20, 30, 40, 50)
  )
  
  design <- training |>
    sccomp:::get_design_matrix(~ age, sample)
  
  expect_equal(
    unname(design[, "age"]),
    as.vector(scale(c(10, 20, 30, 40, 50)))
  )
})
