library(testthat)
test_that("old .sample/.cell_group argument style works with placeholders", {
  skip_cmdstan()
  
  suppressPackageStartupMessages({
    library(dplyr)
  })

  test_data <- tibble::tribble(
    ~orig.ident, ~micro_conf, ~sex, ~age, ~manual_fine, ~count,
    "NS031", "yes", "F", 31, "Memory B cells", 155L,
    "NS031", "yes", "F", 31, "Naïve B cells", 18L,
    "NS031", "yes", "F", 31, "IGHG1+ plasma cells", 45L,
    "NS031", "yes", "F", 31, "IGHM+ memory B cells", 36L,
    "NS031", "yes", "F", 31, "IGKC+ plasma cells", 56L,
    "NS031", "yes", "F", 31, "IGHA+ memory B cells", 41L,
    "NS031", "yes", "F", 31, "IGLC1+ memory B cells", 14L,
    "NS031", "yes", "F", 31, "Plasmablasts", 12L,
    "NS031", "yes", "F", 31, "IGHM+ plasma cells", 2L,
    "NS031", "yes", "F", 31, "IGHD+ memory B cells", 4L,
    "MP075", "no", "M", 25, "Memory B cells", 136L,
    "MP075", "no", "M", 25, "Naïve B cells", 40L,
    "MP075", "no", "M", 25, "IGHG1+ plasma cells", 25L,
    "MP075", "no", "M", 25, "IGHM+ memory B cells", 58L,
    "MP075", "no", "M", 25, "IGKC+ plasma cells", 0L,
    "MP075", "no", "M", 25, "IGHA+ memory B cells", 0L,
    "MP075", "no", "M", 25, "IGLC1+ memory B cells", 0L,
    "MP075", "no", "M", 25, "Plasmablasts", 0L,
    "MP075", "no", "M", 25, "IGHM+ plasma cells", 0L,
    "MP075", "no", "M", 25, "IGHD+ memory B cells", 0L
  )

  res <- test_data |> 
       sccomp_estimate(
             formula_composition = ~ micro_conf,
             .sample = orig.ident,
             .cell_group = manual_fine,
             abundance = "count",
             cores = 1, verbose = FALSE
           ) |>
      sccomp_test()

  expect_s3_class(res, "sccomp_tbl")
})

