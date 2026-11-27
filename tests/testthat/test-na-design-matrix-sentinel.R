library(testthat)
library(dplyr)
library(tibble)

test_that("get_design_matrix preserves columns and weights when names or levels contain NA", {
  # Factor level ending in "NA" (e.g., "RNA") must not be dropped as an NA dummy
  df_level_suffix <- tibble(
    sample = c("S1", "S2", "S3"),
    assay = factor(c("ATAC", "RNA", NA), levels = c("ATAC", "RNA"))
  )
  dm_level_suffix <- sccomp:::get_design_matrix(
    df_level_suffix, ~ assay, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(colnames(dm_level_suffix), c("(Intercept)", "assayRNA"))
  expect_equal(unname(dm_level_suffix[, "assayRNA"]), c(0, 1, 0.5))

  # Factor name containing "NA" (e.g., "RNA_type") must distribute weights to its columns
  df_var_contains <- tibble(
    sample = c("S1", "S2", "S3"),
    RNA_type = factor(c("A", "B", NA), levels = c("A", "B"))
  )
  dm_var_contains <- sccomp:::get_design_matrix(
    df_var_contains, ~ RNA_type, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(colnames(dm_var_contains), c("(Intercept)", "RNA_typeB"))
  expect_equal(unname(dm_var_contains[, "RNA_typeB"]), c(0, 1, 0.5))

  # Numeric covariate ending in "NA" (e.g., "scRNA") must not be dropped when another column has NA
  df_num_suffix <- tibble(
    sample = c("S1", "S2", "S3"),
    type = factor(c("A", "B", NA), levels = c("A", "B")),
    scRNA = c(-1, 0, 1)
  )
  dm_num_suffix <- sccomp:::get_design_matrix(
    df_num_suffix, ~ type + scRNA, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(colnames(dm_num_suffix), c("(Intercept)", "typeB", "scRNA"))
  expect_equal(unname(dm_num_suffix[, "typeB"]), c(0, 1, 0.5))
  expect_equal(unname(dm_num_suffix[, "scRNA"]), c(-1, 0, 1))

  # Interaction with a non-NA factor whose level contains "NA" (e.g., "NAIVE")
  df_interaction <- tibble(
    sample = c("S1", "S2", "S3", "S4", "S5"),
    group = factor(c("G1", "G2", "G1", "G2", NA), levels = c("G1", "G2")),
    cell_state = factor(c("CTRL", "CTRL", "NAIVE", "NAIVE", "NAIVE"), levels = c("CTRL", "NAIVE"))
  )
  dm_interaction <- sccomp:::get_design_matrix(
    df_interaction, ~ group * cell_state, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(
    colnames(dm_interaction),
    c("(Intercept)", "groupG2", "cell_stateNAIVE", "groupG2:cell_stateNAIVE")
  )
  expect_equal(unname(dm_interaction["S5", ]), c(1, 0.5, 1, 0.5))

  # Three-way interaction with NA in both categorical variables and a numeric variable ending in NA
  df_3way <- tibble(
    sample = c("S1", "S2", "S3", "S4", "S5"),
    RNA_type = factor(c("A", "B", "A", "B", NA), levels = c("A", "B")),
    cell_state = factor(c("CTRL", "CTRL", "NAIVE", "NAIVE", NA), levels = c("CTRL", "NAIVE")),
    scRNA = c(-1, 0, 1, 2, 0.8)
  )
  dm_3way <- sccomp:::get_design_matrix(
    df_3way, ~ RNA_type * cell_state * scRNA, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(
    colnames(dm_3way),
    c(
      "(Intercept)", "RNA_typeB", "cell_stateNAIVE", "scRNA",
      "RNA_typeB:cell_stateNAIVE", "RNA_typeB:scRNA",
      "cell_stateNAIVE:scRNA", "RNA_typeB:cell_stateNAIVE:scRNA"
    )
  )
  expect_equal(
    unname(dm_3way["S5", ]),
    c(1, 0.5, 0.5, 0.8, 0.25, 0.4, 0.4, 0.2),
    tolerance = 1e-12
  )

  # Negative continuous covariate in an interaction with an NA categorical factor
  # (NA dummy + NA:numeric sums to < 0, so row detection must check != 0, not > 0)
  df_neg_num <- tibble(
    sample = c("S1", "S2", "S3"),
    RNA_type = factor(c("A", "B", NA), levels = c("A", "B")),
    scRNA = c(1, 2, -2)
  )
  dm_neg_num <- sccomp:::get_design_matrix(
    df_neg_num, ~ RNA_type * scRNA, sample, accept_NA_as_average_effect = TRUE
  )
  expect_equal(colnames(dm_neg_num), c("(Intercept)", "RNA_typeB", "scRNA", "RNA_typeB:scRNA"))
  expect_equal(unname(dm_neg_num["S3", ]), c(1, 0.5, -2, -1), tolerance = 1e-12)
})
