stan_test <- function(label, code) {
  skip_if_not_installed("withr")
  expr <- substitute(
    testthat::test_that(label, code),
    env = list(label = label, code = substitute(code))
  )
  temp <- tempfile()
  dir.create(temp)
  on.exit(unlink(temp))
  withr::with_dir(
    temp,
    suppressMessages(eval(expr, envir = parent.frame()))
  )
}

skip_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  if (!stan_cmdstan_exists()) {
    skip("CmdStan not found.")
  }
}
