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

test_that("build_stan_parameter_subset names terms that sit late in the design", {
  # The matched terms are at design positions 4 and 6, well past `length(matched)`.
  # Subsetting `matched` by those positions instead of by position within
  # `matched` returns NA, and the random-effect path then filters its draws by
  # these names, dropping the columns the contrast needs.
  design_columns <- c(
    "(Intercept)___a", "x1___a", "x2___a", "x3___a", "x1___b", "x3___b"
  )
  model_input <- list(y = matrix(0, nrow = 1, ncol = 2))

  subset <- sccomp:::build_stan_parameter_subset(
    contrasts = c(late = "x3___a + 0.5 * x3___b"),
    design_columns = design_columns,
    stan_parameter = "random_effect_2",
    model_input = model_input
  )

  expect_false(anyNA(subset$parameter))
  expect_setequal(unique(subset$parameter), c("x3___a", "x3___b"))

  # Each name must still carry the Stan index of its own design column
  expect_setequal(
    subset$variable[subset$parameter == "x3___a"],
    c("random_effect_2[4,1]", "random_effect_2[4,2]")
  )
  expect_setequal(
    subset$variable[subset$parameter == "x3___b"],
    c("random_effect_2[6,1]", "random_effect_2[6,2]")
  )
})
