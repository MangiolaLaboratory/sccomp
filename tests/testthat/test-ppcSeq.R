context('ppcseq')

test_that("first test",{

  library(dplyr)
  library(sccomp)
  debugonce(sccomp_glm)

  res =
    sccomp::cell_counts %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type, count,
      cores=1
    )

  expect_equal(

    as.integer(unlist(res[,4])),
    c(0,1,0)
  )

})




