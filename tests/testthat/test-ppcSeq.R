context('ppcseq')

# test_that("dirichlet multinomial outliers",{
#
#   library(dplyr)
#   library(sccomp)
#   library(digest)
#
#   res =
#     sccomp::cell_counts %>%
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

  library(dplyr)
  library(sccomp)
  library(digest)

  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count,
      noise_model = "dirichlet_multinomial"
    )

  expect_equal(
    res %>%
      distinct(cell_type, significant) %>%
      pull(significant) %>%
      digest(algo="md5"),
    "872b45dd0f77c8a355f767294fb44298"
  )

  # [1]  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE

})


test_that("multi beta",{


  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sccomp)
  library(digest)
  library(rstan)

  # Get proportions
  res2 =
    sccomp::cell_counts %>%
    mutate(count = count + 1L) %>%
    nest(data = -sample) %>%
    mutate(exposure = map_int(data, ~ sum(.x$count))) %>%
    unnest(data) %>%
    mutate(frac = (count/exposure) ) %>%
    select(-c(count, exposure, phenotype)) %>%

    sccomp_glm(
      formula = ~ type,
      sample, cell_type, frac,
      noise_model = "multi_beta"
    )

})

test_that("multi beta binomial",{

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sccomp)
  library(digest)
  library(rstan)

  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count,
      noise_model = "multi_beta_binomial"
    )

})
