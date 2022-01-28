library(job)

model_glm_dirichlet_multinomial = read_file("dev/stan_models/glm_dirichlet_multinomial.stan")
save(model_glm_dirichlet_multinomial, file="data/model_glm_dirichlet_multinomial.rda", compress = "xz")


model_imputation = read_file("dev/stan_models/glm_dirichlet_multinomial_imputation.stan")
save(model_imputation, file="data/model_glm_dirichlet_multinomial_imputation.rda", compress = "xz")


glm_dirichlet_multinomial_generate_quantities = read_file("dev/stan_models/glm_dirichlet_multinomial_generate_quantities.stan")
save(glm_dirichlet_multinomial_generate_quantities, file="data/model_glm_dirichlet_multinomial_generate_quantities.rds", compress = "xz")

job({stan_model("dev/stan_models/glm_dirichlet_multinomial.stan") %>% saveRDS("model_glm_dirichlet_multinomial.rds")})
job({stan_model("dev/stan_models/glm_dirichlet_multinomial_imputation.stan") %>% saveRDS("model_glm_dirichlet_multinomial_imputation.rds")})
job({stan_model("dev/stan_models/glm_dirichlet_multinomial_generate_quantities.stan") %>% saveRDS("model_glm_dirichlet_multinomial_generate_quantities.rds")})
