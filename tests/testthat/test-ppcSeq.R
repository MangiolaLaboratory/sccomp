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
      sample, cell_type, count
    )

  expect_equal(
    res %>%
      distinct(cell_type, significant) %>%
      pull(significant) %>%
      digest(algo="md5"),
    "872b45dd0f77c8a355f767294fb44298"
  )

})


test_that("multi beta",{


  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sccomp)
  library(digest)
  library(rstan)

  input_df =
    sccomp::cell_counts %>%
    mutate(count = count + 1L) %>%
    nest(data = -sample) %>%
    mutate(sample_tot = map_int(data, ~ sum(.x$count))) %>%
    unnest(data) %>%
    mutate(frac = (count/sample_tot) ) %>%
    select(-c(count, sample_tot, phenotype)) %>%
    spread(cell_type, frac)

  f = input_df %>% glm_multi_beta(~ type, sample)

})

test_that("multi beta binomial",{

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sccomp)
  library(digest)
  library(rstan)

  input_df =
    sccomp::cell_counts %>%
    nest(data = -sample) %>%
    mutate(sample_tot = map_int(data, ~ sum(.x$count))) %>%
    unnest(data) %>%
    select(-c(phenotype)) %>%
    spread(cell_type, count)

  f = input_df %>% glm_multi_beta_binomial()

})
