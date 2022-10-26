context('sccomp')

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
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 20,
      mcmc_seed = 42
    ) |>


    replicate_data() |>
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
      cores = 20,
      mcmc_seed = 42
    ) |>


    replicate_data(~ 0 + type) |>
    nrow() |>
    expect_equal(600)

})

test_that("multilevel multi beta binomial from Seurat",{

  seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_glm(
      formula_composition = ~ type + (1 | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 20,
      mcmc_seed = 42
    ) |>


    filter(parameter == "typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(16)


  seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
    sccomp_glm(
      formula_composition = ~ 0 + type + (type | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 20,
      mcmc_seed = 42
    ) |>


    filter(parameter == "typecancer - typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(13)

})

test_that("multilevel multi beta binomial from Seurat with intercept and continuous covariate",{



  #debugonce(sccomp:::data_spread_to_model_input)
  seurat_obj |>
    sccomp_glm(
      formula_composition = ~ continuous_covariate + (1 + continuous_covariate | group__),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 20,
      mcmc_seed = 42
    ) |>


    filter(parameter == "continuous_covariate") |>
    filter(c_pH0<0.1)

  #|>
  #  nrow() |>
  #  expect_equal(1)

})


test_that("multilevel continuous",{


  # seurat_obj =
  #   seurat_obj |>
  #   mutate(group__ = glue::glue("GROUP{group__}")) |>
  #   nest(data = -c(sample, type)) |>
  #   mutate(group2__ = glue("GROUP2{sample(c(1, 2), replace = T, size = n())}") ) |>
  #   mutate(continuous_covariate = rnorm(n())) |>
  #   unnest(data)

  seurat_obj |>
    sccomp_glm(
      formula_composition = ~ 0 + type + continuous_covariate + (type | group__) + (continuous_covariate | type),
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
      cores = 20,
      mcmc_seed = 42
    ) |>


    filter(parameter == "typecancer - typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(16)

})

test_that("wrongly-set groups",{

  # library(tidyseurat)
  # seurat_obj =
  #   seurat_obj |>
  #   nest(data = -c(sample, type)) |>
  #   mutate(group__wrong = c(1,1,1,1,1, 2,2,2,2,2, 1,1,1,1,1, 2,2,2,2,2) |> as.character()) |>
  #   unnest(data)

    expect_error(
      object =
        seurat_obj |>
        ## filter(cell_group %in% c("NK cycling", "B immature")) |>
        sccomp_glm(
          formula_composition = ~ 0 + type + (type | group__wrong),
          formula_variability = ~ 1,
          sample, cell_group,
          check_outliers = FALSE,
          approximate_posterior_inference = "all",
          contrasts = c("typecancer - typehealthy", "typehealthy - typecancer"),
          cores = 20,
          mcmc_seed = 42
        ) ,
      regexp = "should not be shared"
    )

})

test_that("multi beta binomial from Seurat",{

  seurat_obj |>
    sccomp_glm(
      formula_composition = ~  type,
      formula_variability = ~ 1,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = "all",
      cores = 1,
      mcmc_seed = 42
    )  |>
    filter(parameter == "typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(18)

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

  estimate |> remove_unwanted_variation(~ type)

})

test_that("multi beta binomial from SCE",{

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
    )  |>
    filter(parameter == "typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(17)

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
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(17)

})

test_that("plot test composition",{

  plot_summary(res_composition)


})

test_that("plot test composition",{

  plot_summary(res_composition_variability)


})
