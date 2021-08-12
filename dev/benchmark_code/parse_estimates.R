library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)


# Read arguments
args = commandArgs(trailingOnly=TRUE)

estimate_file_1 = args[6]
estimate_file_2 = args[7]
estimate_file_3 = args[8]
output_file = args[9]

probs = seq(0, 0.1,length.out = 50) %>% c(seq(0.1, 1,length.out = 50))

# Import files
readRDS(estimate_file_1) %>%
  left_join(readRDS(estimate_file_2)) %>%
  left_join(readRDS(estimate_file_3))  %>%

  # Calculate sgnificance
  mutate(hypothesis_edger = map(results_edger , ~ .x %>% arrange(FDR) %>% mutate(probability = 1-PValue) %>% mutate(estimate = logFC))) %>%
  mutate(hypothesis_voom = map(results_voom , ~.x %>% arrange(adj.P.Val) %>% mutate(probability = 1-P.Value) %>% mutate(estimate = logFC))) %>%
  mutate(hypothesis_speckle = map(results_speckle , ~ .x %>% arrange(FDR) %>% mutate(probability = 1-P.Value) %>% mutate(estimate = -Tstatistic    ))) %>%
  mutate(hypothesis_sccomp = map(
    results_sccomp ,
    ~ .x  %>%
      arrange(false_discovery_rate) %>%
      mutate(probability = 1-false_discovery_rate) %>%
      mutate(estimate = .median_type )
  )) %>%
  mutate(hypothesis_DirichletMultinomial  = map(
    results_DirichletMultinomial ,
    ~ .x  %>%
      arrange(false_discovery_rate) %>%
      mutate(probability = 1-false_discovery_rate) %>%
      mutate(estimate = .median_type )
  )) %>%

  # Clean
  dplyr::select(-contains("results")) %>%

  # Reshape
  pivot_longer(contains("hypothesis"),names_prefix = "hypothesis_" ) %>%

  # Calculate ROC
  mutate(df_for_roc = map2(
    data, value,
    ~ left_join(
      .x %>%
        unnest(coefficients) %>%
        distinct(cell_type, beta_1),

      .y  %>%
        select(cell_type,  probability , estimate),
      by = "cell_type"
    )
  )) %>%
  select(run, name, df_for_roc) %>%

  unnest(df_for_roc) %>%
  nest(df_for_roc = -name) %>%
  mutate(df_for_roc = map(
    df_for_roc,
    ~ .x %>%
      mutate(probability = if_else(beta_1 * estimate < 0, 0, probability)) %>%
      mutate(significant = (beta_1 != 0) %>% as.factor())
  )) %>%


# NOT USED ANYMORE
# # Calculate accuracy
#   pivot_longer(contains("hypothesis"),names_prefix = "hypothesis_" ) %>%
#   mutate(accuracy_df = map2(
#     data, value             ,
#     ~ left_join(
#
#       .x %>% unnest(coefficients) %>% dplyr::select(cell_type, beta_1) %>% distinct %>% mutate(cell_type = as.character(cell_type)),
#       .y %>% dplyr::select(cell_type, positive, trend),
#       by="cell_type"
#
#     )
#
#   )) %>%
#   mutate(TP = map_int(accuracy_df, ~ .x %>% filter(positive & (beta_1 != 0) & (beta_1 * trend)>0) %>% nrow())) %>%
#   mutate(FP = map_int(accuracy_df, ~ .x %>% filter(
#
#     # Positive when not
#     (positive & (beta_1 == 0)) |
#
#       # Positive when yes but wrong direction
#       (positive & (beta_1 != 0) & (beta_1 * trend)<0)
#   ) %>% nrow())) %>%
#   mutate(total_true_positive = 6, total_true_negative = 24-6) %>%
#   mutate(FP_rate = FP/total_true_negative) %>%
#   mutate(TP_rate = TP/total_true_positive) %>%


  mutate(
    slope = as.numeric(args[1]),
    n_samples = as.integer(args[2]),
    n_cell_type = as.integer(args[3]),
    max_cell_counts_per_sample = as.integer(args[4]),
    add_outliers = as.integer(args[5])
  ) %>%

  # Save
  saveRDS(output_file)

