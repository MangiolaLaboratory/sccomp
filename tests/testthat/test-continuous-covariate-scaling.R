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
    scaled <- grid |>
      sccomp:::scale_numeric_covariates("age", reference = training)
    design <- scaled |> sccomp:::get_design_matrix(~ age, sample)
    design[grid$sample[grid$age == 35], "age"]
  }
  
  expect_equal(column_at_35(narrow_grid), column_at_35(wide_grid))
  expect_equal(
    unname(column_at_35(narrow_grid)),
    (35 - 30) / sd(c(10, 20, 30, 40, 50))
  )
})

test_that("get_design_matrix does not scale a second time", {
  training <- tibble::tibble(
    sample = paste0("S", 1:5),
    age = c(10, 20, 30, 40, 50)
  )
  scaled <- training |> sccomp:::scale_numeric_covariates("age")
  design <- scaled |> sccomp:::get_design_matrix(~ age, sample)
  
  expect_equal(unname(design[, "age"]), scaled$age)
})

test_that("an NA covariate resolves to the fitted average, not the grid average", {
  training <- tibble::tibble(
    sample = paste0("S", 1:6),
    age = c(10, 20, 30, 40, 50, 60),
    condition = c("a", "b", "c", "a", "b", "c")
  )
  
  # The grid sits far off the fitted centre and never mentions level "c", so a
  # design read off the grid alone would centre the NA age on 150 and split the
  # NA condition two ways instead of three.
  grid <- tibble::tibble(
    sample = c("g1", "g2", "g3"),
    age = c(100, NA, 200),
    condition = c("a", "b", NA)
  )
  
  design <-
    grid |>
    sccomp:::declare_fitted_levels(training, exclude = "sample") |>
    sccomp:::scale_numeric_covariates(c("age", "condition"), reference = training) |>
    sccomp:::get_design_matrix(~ age + condition, sample,
                               accept_NA_as_average_effect = TRUE)
  
  # The level absent from the grid still gets a column, so the design lines up
  # with the fitted parameters
  expect_equal(colnames(design), c("(Intercept)", "age", "conditionb", "conditionc"))
  
  # An unknown age is the fitted mean, which z-scoring puts at 0
  expect_equal(unname(design["g2", "age"]), 0)
  
  # An unknown condition spreads evenly over the three fitted levels
  expect_equal(unname(design["g3", c("conditionb", "conditionc")]), c(1/3, 1/3))
})
