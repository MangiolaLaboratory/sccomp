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




