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

test_that("multilevel multi beta binomial from Seurat",{

  library(tidyseurat)
  seurat_obj =
    seurat_obj |>
    nest(data = -c(sample, type)) |>
    mutate(group__ = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4,4, 4, 4) |> as.character()) |>
    unnest(data)

  seurat_obj |>
    ## filter(cell_group %in% c("NK cycling", "B immature")) |>
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


    filter(parameter == "typecancer - typehealthy") |>
    filter(c_pH0<0.1) |>
    nrow() |>
    expect_equal(11)

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
    expect_equal(17)

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
