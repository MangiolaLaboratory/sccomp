
skip_cmdstan <- function() {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  if (!stan_cmdstan_exists()) {
    skip("CmdStan not found.")
  }
}
