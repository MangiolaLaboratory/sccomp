library(dplyr)
library(sccomp)
library(tidyverse)
data("seurat_obj")
data("sce_obj")
data("counts_obj")


my_estimate_0_1_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42
  )  |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_1_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_1_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_3_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_3_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_3_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_33_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_33_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_0_33_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group,  prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39
  ) |> sccomp_test(contrasts = "typehealthy")

my_estimate_no_intercept_0_1_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42
  )  |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_1_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_1_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 1), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))


my_estimate_no_intercept_0_3_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_3_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_3_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,1)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_33_a = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 42,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_33_b = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 40,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

my_estimate_no_intercept_0_33_c = 
  seurat_obj |>
  sccomp_estimate(
    formula_composition = ~ 0 + type + continuous_covariate,
    formula_variability = ~ 1,
    sample, cell_group, prior_overdispersion_mean_association = list(intercept = c(4, 3), slope = c(0, 2),     standard_deviation = c(10, 10)),
    prior_mean = list(intercept = c(0, 3), coefficients = c(0,3)),
    approximate_posterior_inference = FALSE,
    cores = 1,  #exclude_priors = T,
    mcmc_seed = 39,
    max_sampling_iterations = 1000
  ) |> sccomp_test(contrasts = c("typehealthy - typecancer"))

df_of_models = 
  tibble(
    test = c(
      "my_estimate_0_1_a", "my_estimate_0_1_b", "my_estimate_0_1_c", "my_estimate_0_3_a", "my_estimate_0_3_b", "my_estimate_0_3_c",
      "my_estimate_0_33_a", "my_estimate_0_33_b", "my_estimate_0_33_c",
      "my_estimate_no_intercept_0_1_a", "my_estimate_no_intercept_0_1_b", "my_estimate_no_intercept_0_1_c", "my_estimate_no_intercept_0_3_a", "my_estimate_no_intercept_0_3_b", "my_estimate_no_intercept_0_3_c",
      "my_estimate_no_intercept_0_33_a", "my_estimate_no_intercept_0_33_b", "my_estimate_no_intercept_0_33_c"
    ),
    estimate = list(
      my_estimate_0_1_a, my_estimate_0_1_b, my_estimate_0_1_c, my_estimate_0_3_a, my_estimate_0_3_b, my_estimate_0_3_c,
      my_estimate_0_33_a, my_estimate_0_33_b, my_estimate_0_33_c,
      my_estimate_no_intercept_0_1_a, my_estimate_no_intercept_0_1_b, my_estimate_no_intercept_0_1_c, my_estimate_no_intercept_0_3_a, my_estimate_no_intercept_0_3_b, my_estimate_no_intercept_0_3_c,
      my_estimate_no_intercept_0_33_a, my_estimate_no_intercept_0_33_b, my_estimate_no_intercept_0_33_c
    )
  ) 


df_of_models |> 
  mutate(estimate = map(estimate, ~ .x |> select(cell_group, starts_with("c_")))) |> 
  unnest(estimate) |> 
  tidyr::extract(test, c("model", "prior", "run"), "([a-z_]+)_[0-9]_([0-9]+)_([a-z])", remove = F) |> 
  mutate(log_c_pH0 = log(c_pH0)) |> 
  
  ggplot(aes(cell_group, c_effect, color = model, shape = prior)) +
  geom_point(size = 2) +
  theme(axis.text.x = element_text(angle=90))

df_of_models |> 
  mutate(estimate = map(estimate, ~ .x |> select(cell_group, starts_with("c_")))) |> 
  unnest(estimate) |> 
  extract(test, c("model", "prior", "run"), "([a-z_]+)_[0-9]_([0-9]+)_([a-z])", remove = F) |> 
  mutate(log_c_pH0 = log(c_pH0)) |> 
  tidybulk::reduce_dimensions(test, cell_group, c_effect, method = "PCA", scale = FALSE) |> 
  
  ggplot(aes(PC1, PC2, color = model, shape = prior)) +
  
  geom_point(size = 2) +
  theme(axis.text.x = element_text(angle=90))
