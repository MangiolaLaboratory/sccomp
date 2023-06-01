library(here)
library(readr)
glm_multi_beta_binomial = read_file(here("dev/glm_multi_beta_binomial.stan"))
save(glm_multi_beta_binomial, file = here("data/glm_multi_beta_binomial.rda"), compress = "xz")

glm_multi_beta_binomial_generate = read_file(here("dev/glm_multi_beta_binomial_generate.stan"))
save(glm_multi_beta_binomial_generate, file = here("data/glm_multi_beta_binomial_generate.rda"), compress = "xz")

glm_multi_beta_binomial_simulate = read_file(here("dev/glm_multi_beta_binomial_simulate.stan"))
save(glm_multi_beta_binomial_simulate, file = here("data/glm_multi_beta_binomial_simulate.rda"), compress = "xz")