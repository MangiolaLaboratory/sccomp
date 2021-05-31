context('ppcseq')

library(dplyr)
library(sccomp)
data("seurat_obj")
data("sce_obj")
data("counts_obj ")

# test_that("dirichlet multinomial outliers",{
#
#   library(dplyr)
#   library(sccomp)
#   library(digest)
#
#   res =
#     sccomp::counts_obj  %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type, count
#     )
#
#   expect_equal(
#     res %>%
#       distinct(cell_type, significant) %>%
#       pull(significant) %>%
#       digest(algo="md5"),
#     "f6cef772af43198f586e15c96b2f1239"
#   )
#
# })

test_that("dirichlet multinomial",{

  res =
    sccomp::counts_obj  %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_group, count,
      noise_model = "dirichlet_multinomial"
    )

  # [1]  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE

})


test_that("multi beta binomial from Seurat",{

  res =
    seurat_obj %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE
    )

})

test_that("multi beta binomial from SCE",{

  res =
    sce_obj %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE
    )

})

test_that("multi beta binomial from metadata",{

  res =
    seurat_obj[[]] %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_group,
      check_outliers = FALSE
    )

})


test_that("multi beta binomial from counts",{

res =
  counts_obj  %>%
  sccomp_glm(
    formula = ~ type,
    sample, cell_group, count,
    check_outliers = FALSE
  )

})

test_that("multi beta binomial outliers",{

  res =
    seurat_obj %>%
    sccomp_glm(formula = ~ type, sample, cell_group,  approximate_posterior_inference = F )

  res =
    seurat_obj %>%
    sccomp_glm(~ type,  sample, cell_group )
})
