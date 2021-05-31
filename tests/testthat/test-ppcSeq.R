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

  # [1]  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE

})


# test_that("multi beta",{
#
#
#   library(dplyr)
#   library(tidyr)
#   library(purrr)
#   library(sccomp)
#   library(digest)
#   library(rstan)
#
#   # Get proportions
#   res2 =
#     sccomp::cell_counts %>%
#     mutate(count = count + 1L) %>%
#     nest(data = -sample) %>%
#     mutate(exposure = map_int(data, ~ sum(.x$count))) %>%
#     unnest(data) %>%
#     mutate(frac = (count/exposure) ) %>%
#     select(-c(count, exposure, phenotype)) %>%
#
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type, frac,
#       noise_model = "multi_beta"
#     )
#
# })

test_that("multi beta binomial",{

library(dplyr)
library(tidyr)
library(purrr)
library(sccomp)
library(digest)
library(rstan)

#debugonce(multi_beta_binomial_glm)

res =
  sccomp::cell_counts %>%
  sccomp_glm(
    formula = ~ type,
    sample, cell_type, count,
    noise_model = "multi_beta_binomial",
    check_outliers = FALSE
  )

})

test_that("multi beta binomial outliers",{

  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sccomp)
  library(digest)
  library(rstan)
#debugonce(multi_beta_binomial_glm)
  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count,
      noise_model = "multi_beta_binomial",
      approximate_posterior_inference = F
    )

  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count,
      noise_model = "multi_beta_binomial"
    )
})
