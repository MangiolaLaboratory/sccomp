library(dplyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj")


# test_that("dirichlet multinomial outliers",{
#
#   library(dplyr)
#   library(sccomp)
#   library(digest)
#
#   res =
#     sccomp::counts_obj  |>
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type, count
#     )
#
#   expect_equal(
#     res |>
#       distinct(cell_type, significant) |>
#       pull(significant) |>
#       digest(algo="md5"),
#     "f6cef772af43198f586e15c96b2f1239"
#   )
#
# })
# test_that("counts dirichlet multinomial outlier VB",{
#
#   if(interactive()){
#     res =
#       counts_obj  |>
#       sccomp_glm(
#         formula = ~ type,
#         sample, cell_group, count,
#         noise_model = "dirichlet_multinomial",
#         cores = 1
#       )
#   }
#
# })

test_that("Generate data",{


  seurat_obj |>
    sccomp_glm(
      formula_composition = ~ type ,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = FALSE,
      cores = 1,
      mcmc_seed = 42
    ) |>


    sccomp_replicate() |>
    nrow() |>
    expect_equal(600)

  # With grouping
  seurat_obj |>
    sccomp_glm(
      formula_composition = ~ 0 + type + (type | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = FALSE,
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 1,
      mcmc_seed = 42
    ) |>


    sccomp_replicate(~ 0 + type) |>
    nrow() |>
    expect_equal(600)


})

test_that("multilevel multi beta binomial from Seurat",{

  res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_glm(
      formula_composition = ~ type + (1 | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

  # Check order
  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1:3) |>
    pull(cell_group) |>
    expect_equal(c("CD4 ribosome" ,        "CD4 cm high cytokine", "Mono NKG7 2"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

  res =
    seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_glm(
      formula_composition = ~ 0 + type + (type | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 1,
      mcmc_seed = 42
    )

  # res |>
  #   filter(parameter == "typecancer - typehealthy") |>
  #   arrange(desc(abs(c_effect))) |>
  #   slice(1:3) |>
  #   pull(cell_group) |>
  #   sort() |>
  #   expect_equal(c("B mem"  ,  "CD4 cm high cytokine" ,"CD4 ribosome"         ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

test_that("multilevel multi beta binomial from Seurat with intercept and continuous covariate",{


  res =
    seurat_obj |>
    sccomp_glm(
      formula_composition = ~ continuous_covariate + (1 + continuous_covariate | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )


  res |>
    filter(parameter == "continuous_covariate") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("CD8 em 1"))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)


})


test_that("multilevel continuous",{


  # seurat_obj =
  #   seurat_obj |>
  #   mutate(group__ = glue::glue("GROUP{group__}")) |>
  #   nest(data = -c(sample, type)) |>
  #   mutate(group2__ = glue("GROUP2{sample(c(1, 2), replace = T, size = n())}") ) |>
  #   mutate(continuous_covariate = rnorm(n())) |>
  #   unnest(data)

  res =
    seurat_obj |>
    sccomp_glm(
      formula_composition = ~ 0 + type + continuous_covariate + (type | group__) + (continuous_covariate | type),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 1,
      mcmc_seed = 42
    )

  res |>
    filter(parameter == "typecancer - typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

# test_that("wrongly-set groups",{
#
#   # library(tidyseurat)
#   # seurat_obj =
#   #   seurat_obj |>
#   #   nest(data = -c(sample, type)) |>
#   #   mutate(group__wrong = c(1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 2,2,2,2,2) |> as.character()) |>
#   #   unnest(data)
#
#     expect_error(
#       object =
#         seurat_obj |>
#         ## filter(cell_group %in% c("NK cycling", "B immature")) |>
#         sccomp_glm(
#           formula_composition = ~ 0 + type + (type | group__wrong),
#           formula_variability = ~ 1,
#           sample, cell_group,
#           check_outliers = FALSE,
#           approximate_posterior_inference = "all",
#           contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
#           cores = 1,
#           mcmc_seed = 42
#         ) ,
#       regexp = "should not be shared"
#     )
#
# })

test_that("multi beta binomial from Seurat",{

  res =
    seurat_obj |>
    sccomp_glm(
      formula_composition = ~  type,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

test_that("multi beta binomial contrasts from Seurat",{

  res = seurat_obj |>
    sccomp_glm(
      formula_composition = ~ 0 + type,
      formula_variability = ~ 1,
      sample, cell_group,
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

    expect_equal(
      res[1,"c_effect"] |> as.numeric(),
      -res[2,"c_effect"] |> as.numeric()
    )

})

test_that("multi beta binomial contrasts from Seurat",{

  res = seurat_obj |>
    sccomp_glm(
      formula_composition = ~ 0 + type,
      formula_variability = ~ 1,
      sample, cell_group,
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

  expect_equal(
    res[1,"c_effect"] |> as.numeric(),
    -res[2,"c_effect"] |> as.numeric()
  )

})

test_that("remove unwanted variation",{

  library(tidyseurat)

  data =
    seurat_obj |>

    # Add batch
    nest(data = -c(sample, type)) |>
    mutate(batch = rep(c(0,1), 10)) |>
    unnest(data)

  # Estimate
  estimate =
    data |>
    sccomp_glm(
      formula_composition = ~ type + batch,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

  estimate |>
    remove_unwanted_variation(~ type)

})

test_that("multi beta binomial from SCE",{

    res =
      sce_obj |>
    sccomp_glm(
      formula_composition = ~ type,
      formula_variability = ~ 1,
      sample,
      cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )

  res |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)
})

res_composition =
  seurat_obj[[]] |>
  sccomp_glm(
    formula_composition = ~ type,
    formula_variability = ~ 1,
    sample,
    cell_group,
    check_outliers = FALSE,
    approximate_posterior_inference = "all",
    cores = 1,
    mcmc_seed = 42
  )

res_composition_variability =
  seurat_obj[[]] |>
  sccomp_glm(
    formula_composition = ~ type,
    formula_variability = ~ type,
    sample,
    cell_group,
    check_outliers = FALSE,
    approximate_posterior_inference = "all",
    cores = 1,
    mcmc_seed = 42
  )

test_that("multi beta binomial from metadata",{

  res_composition  |>
    filter(parameter == "typehealthy") |>
    arrange(desc(abs(c_effect))) |>
    slice(1) |>
    pull(cell_group) |>
    expect_equal(c("B mem"  ))

  # Check convergence
  res_composition |>
    filter(c_R_k_hat > 4) |>
    nrow() |>
    expect_equal(0)

})

test_that("plot test composition",{

  plot_summary(res_composition)


})

test_that("plot test variability",{

  plot_summary(res_composition_variability)


})

test_that("test constrasts",{


  estimate =
    seurat_obj |>
    sccomp_glm(
      formula_composition = ~ type ,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = FALSE,
      cores = 1,
      mcmc_seed = 42
    )

  new_test =
    estimate |>
    test_contrasts() |>
    test_contrasts()

})
