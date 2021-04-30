

test_check("ppcseq")
library(testthat)
library(sccomp)
res =
  sccomp::cell_counts %>%
  sccomp_glm(
    formula = ~ type,
    sample, cell_type, count
  )
