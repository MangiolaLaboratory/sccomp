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

test_that("counts dirichlet multinomial outlier VB",{

  if(interactive()){
    res =
      counts_obj  |>
      sccomp_glm(
        formula = ~ type,
        sample, cell_group, count,
        noise_model = "dirichlet_multinomial"
      )
  }

})

test_that("counts multi beta binomial outlier VB",{

  res =
    counts_obj  |>
    sccomp_glm(
      formula = ~ type,
      sample, cell_group, count,
      approximate_posterior_inference = T
    )


})

test_that("counts multi beta binomial outlier VB",{

  res =
    counts_obj  |>
    sccomp_glm(
      formula = ~ type,
      sample, cell_group, count,
      approximate_posterior_inference = F,
      check_outliers = F
    )


})

test_that("multi beta binomial from Seurat",{

  res =
    seurat_obj |>
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = T
    )

})

test_that("multi beta binomial from SCE",{

  res =
    sce_obj |>
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = T
    )

})

test_that("multi beta binomial from metadata",{

  res =
    seurat_obj[[]] |>
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE,
      approximate_posterior_inference = T
    )

})

