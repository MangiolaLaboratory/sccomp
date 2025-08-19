# Test script for sccomp issue #229
# -------------------------------------------------
# This script loads the toy dataset that was causing
# the reframe/full_join error in issue #229 and runs
# the exact analysis pipeline described in the GitHub
# ticket to confirm that the bug is now resolved.
#
# Usage (from R console):
#   source("dev/test_issue_229.R")
# -------------------------------------------------


library(sccomp)
library(dplyr)


sccomp_result <- 
read.csv(file.path("dev", "cellcounts_sccomp_issue.csv")) %>%
  filter(sample_type != "unmapped") %>%
 sccomp_estimate( 
    formula_composition = ~ sample_type, 
    .sample = "sample", 
    .cell_group = "celltype_kt", 
    verbose = FALSE,
    .abundance = "count"
  ) |> 
  sccomp_test()

print(sccomp_result)

message("✔ sccomp test pipeline completed without errors.")
