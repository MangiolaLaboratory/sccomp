library(testthat)
library(sccomp)

test_that("sccomp_identify_covariate_contrasts identifies direct covariates", {
  model_input <- list(
    X = matrix(0, nrow = 1, ncol = 2, dimnames = list(NULL, c("(Intercept)", "treatment"))),
    XA = matrix(0, nrow = 1, ncol = 2, dimnames = list(NULL, c("(Intercept)", "treatment")))
  )

  identified <- sccomp:::sccomp_identify_covariate_contrasts(c("treatment"), model_input)

  expect_false(is.null(identified))
  expect_equal(identified$par_name, "beta")
  expect_equal(identified$full_design, colnames(model_input$X))
  expect_true(isTRUE(identified$variab_ok))
  expect_equal(identified$contrast_mapping$contrast, "treatment")
  expect_equal(unname(identified$contrast_mapping$design_param), "treatment")
})

test_that("sccomp_identify_covariate_contrasts returns NULL for non-atomic or missing terms", {
  model_input <- list(
    X = matrix(0, nrow = 1, ncol = 2, dimnames = list(NULL, c("(Intercept)", "treatment"))),
    XA = matrix(0, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  )

  expect_null(sccomp:::sccomp_identify_covariate_contrasts(NULL, model_input))
  expect_null(sccomp:::sccomp_identify_covariate_contrasts(c("treatment - (Intercept)"), model_input))
  expect_null(sccomp:::sccomp_identify_covariate_contrasts(c("does_not_exist"), model_input))

  identified <- sccomp:::sccomp_identify_covariate_contrasts(c("treatment"), model_input)
  expect_false(is.null(identified))
  expect_false(identified$variab_ok)
})

