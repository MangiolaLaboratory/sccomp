library(tidyseurat)

counts = readRDS("~/PostDoc/oligo_breast/expanded_analyses_with_control/integrated_counts.rds")  %>%
  select(cell, sample, cell_type = curated_cell_type, type)
DefaultAssay(counts) = "RNA"
counts[["SCT"]] = NULL
counts[["integrated"]] = NULL

counts %>% filter(type =="benign") %>%
  mutate(dataset_id = "GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K") %>%
  saveRDS("dev/data_integration/GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")


counts %>% filter(type !="benign") %>%
  mutate(dataset_id = "internal_Pal") %>%
  saveRDS("dev/data_integration/internal_Pal_breast_cancer.rds")


