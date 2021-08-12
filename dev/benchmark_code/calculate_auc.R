library(tidyverse)
library(magrittr)
library(bayestestR)

# Read arguments
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_file = args[2]

readRDS(input_file) %>%
  group_by(name,  slope, n_samples , n_cell_type,  max_cell_counts_per_sample , add_outliers, probs) %>%
  summarise(mean_FP_rate = mean(FP_rate), mean_TP_rate = mean(TP_rate)) %>%
  ungroup() %>%
  arrange(mean_FP_rate) %>%
  nest(data = -c(name,  slope, n_samples , n_cell_type,  max_cell_counts_per_sample , add_outliers)) %>%
  mutate(auc = map_dbl(data, ~ bayestestR::auc(.x$mean_FP_rate , .x$mean_TP_rate))) %>%
  select(-data) %>%
  saveRDS(output_file)


