library(job)

job({
  glm_dirichlet_multinomial =
    read_file("dev/stan_models/glm_dirichlet_multinomial.stan")

  save(glm_dirichlet_multinomial, file="data/glm_dirichlet_multinomial.rda", compress = "xz")
})

job({
  glm_dirichlet_multinomial_generate_quantities =
    read_file("dev/stan_models/glm_dirichlet_multinomial_generate_quantities.stan")

  save(glm_dirichlet_multinomial_generate_quantities, file="data/glm_dirichlet_multinomial_generate_quantities.rda", compress = "xz")
})

job({
  glm_dirichlet_multinomial_imputation =
    read_file("dev/stan_models/glm_dirichlet_multinomial_imputation.stan")

  save(glm_dirichlet_multinomial_imputation, file="data/glm_dirichlet_multinomial_imputation.rda", compress = "xz")
})

job({
  glm_multi_beta_generate_data =
    read_file("dev/stan_models/glm_multi_beta_generate_data.stan")

  save(glm_multi_beta_generate_data, file="data/glm_multi_beta_generate_data.rda", compress = "xz")
})

job({
  glm_multi_beta =
    read_file("dev/stan_models/glm_multi_beta.stan")

  save(glm_multi_beta, file="data/glm_multi_beta.rda", compress = "xz")
})

