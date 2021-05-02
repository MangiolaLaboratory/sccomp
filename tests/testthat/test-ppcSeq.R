context('ppcseq')

test_that("first test",{

  # library(furrr)
  # plan(multisession, workers=20)
  library(dplyr)
  library(sccomp)
  library(digest)
  #debugonce(sccomp_glm)

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
    "f6cef772af43198f586e15c96b2f1239"
  )

})



test_that("multi beta binomial",{

  library(dplyr)
  library(sccomp)
  library(digest)

  input_df =
    sccomp::cell_counts %>%
    mutate(count = count + 1L) %>%
    nest(data = -sample) %>%
    mutate(sample_tot = map_int(data, ~ sum(.x$count))) %>%
    unnest(data) %>%
    select(-c(phenotype)) %>%
    spread(cell_type, count)

  f4 = stan(file = "inst/stan/glm_multi_beta_binomial.stan",
            data = list(
              N = 20,
              M = 36,
              tot = input_df$sample_tot,
              y = input_df %>% select(-type, -sample_tot) %>% nanny::as_matrix(rownames = sample),
              X = input_df %>% select(sample, type) %>% model.matrix(~ type, data=.)
            ),
            cores = 4
  )

})
