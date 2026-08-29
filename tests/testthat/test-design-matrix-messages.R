library(testthat)
library(sccomp)


test_that("design matrix parameter previews are compact", {
  x <- matrix(0, nrow = 2, ncol = 12)
  colnames(x) <- paste0("parameter_", seq_len(ncol(x)))
  
  expect_equal(
    sccomp:::format_design_matrix_preview(x),
    paste0(
      "12 parameters: ",
      paste0("parameter_", seq_len(10), collapse = ", "),
      ", ... (+2 more)"
    )
  )
})


test_that("design matrix messages identify every modelled matrix", {
  model_input <- list(
    X = matrix(
      0,
      nrow = 2,
      ncol = 2,
      dimnames = list(NULL, c("(Intercept)", "conditiontreated"))
    ),
    Xa = matrix(
      0,
      nrow = 2,
      ncol = 1,
      dimnames = list(NULL, "(Intercept)")
    ),
    ncol_X_random_eff = c(4L, 3L, 0L, 0L),
    X_random_effect_1 = matrix(
      0,
      nrow = 2,
      ncol = 4,
      dimnames = list(
        NULL,
        c(
          "(Intercept)___A", "(Intercept)___B",
          "age___A", "age___B"
        )
      )
    ),
    X_random_effect_2 = matrix(
      0,
      nrow = 2,
      ncol = 3,
      dimnames = list(
        NULL,
        c(
          "s(age, k = 5)___basis01",
          "s(age, k = 5)___basis02",
          "s(age, k = 5)___basis03"
        )
      )
    )
  )
  model_input$random_effect_design_terms <- tibble::tibble(
    slot = c(1L, 2L),
    term = c("(1 + age | donor)", "smooth s(age, k = 5)")
  )
  
  messages <- character(0)
  withCallingHandlers(
    sccomp:::message_design_matrices(model_input),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  
  expect_length(messages, 1L)
  expect_match(messages, "composition X - 2 parameters", fixed = TRUE)
  expect_match(messages, "variability Xa - 1 parameter", fixed = TRUE)
  expect_match(
    messages,
    "random effect X1 [(1 + age | donor)] - 4 parameters",
    fixed = TRUE
  )
  expect_match(
    messages,
    "random effect X2 [smooth s(age, k = 5)] - 3 parameters",
    fixed = TRUE
  )
  expect_false(grepl("random effect X3", messages, fixed = TRUE))
})


test_that("design matrix messages survive fits without recorded terms", {
  model_input <- list(
    X = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "(Intercept)")),
    Xa = matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, "(Intercept)")),
    ncol_X_random_eff = c(1L, 0L, 0L, 0L),
    X_random_effect_1 = matrix(
      0,
      nrow = 2,
      ncol = 1,
      dimnames = list(NULL, "(Intercept)___A")
    )
  )
  
  expect_message(
    sccomp:::message_design_matrices(model_input),
    "random effect X1 [slot 1] - 1 parameter",
    fixed = TRUE
  )
})
