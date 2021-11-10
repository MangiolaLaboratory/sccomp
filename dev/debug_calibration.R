library(tidyverse)
library(sccomp)
design_matrix_and_coefficients_to_dir_mult_simulation(
  c(rep(1, 100), rep(0, 100)) %>% sample,
  matrix(c(0, 0, 0, 0), ncol = 2),
  seed = seed
)  %>%
sccomp_glm(
  ~ covariate_1,
  .sample = sample,
  .cell_group = cell_type,
  .count = generated_counts,
  check_outliers = F,
  seed = seed,
  prior_mean_variable_association = list(
    intercept = c(0, 5),
    slope = c(0,  5),
    standard_deviation = c(0, 2)
  )
)



data_for_model_simple = data_for_model
data_for_model_simple$X[,2] = data_for_model_simple$X[,2] + 0.5
fit_simple =
  sampling(
    readRDS("dev/simple_model_for_testing.rds"),
    data = data_for_model_simple,
    cores = 4
  )

fit_simple %>% tidybayes::gather_draws(beta[C, M]) %>% filter(C==2 & M==1) %>% filter(.value>0.5) %>% nrow %>% `/` (4000) %>%  `*` (2)

fit_simple %>% traceplot("beta")


f2 %>% tidybayes::gather_draws(beta[C, M]) %>% filter(C==2 & M==1) %>% filter(.value>0) %>% nrow %>% `/` (4000) %>%  `*` (2)


mcmc_parcoord(f, np = nuts_params(f), regex_pars = "beta|alpha|prec_sd|prec_coeff" ) + theme(axis.text.x = element_text(angle=90))
