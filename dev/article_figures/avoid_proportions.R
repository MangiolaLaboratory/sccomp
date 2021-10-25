setwd("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/")

library(rlang)
library(tidyverse)




readRDS("dev/data_integration/UVM_single_cell/counts.rds")[[]] %>%
  class_list_to_counts(sample, cell_type)
  with_groups(sample, ~ mutate(.x, proportion = (count+1)/sum(count+1)) ) %>%
