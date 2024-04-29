library(sccomp)
library(tidyverse)
library(magrittr)
library(tidybulk)

probs = seq(0, 0.1,length.out = 50) %>% c(seq(0.1, 1,length.out = 50))

# Read arguments
args = commandArgs(trailingOnly=TRUE)

output_file = args[6]

# Estimated files
estimated_files = args[7:length(args)]

print(estimated_files)

# Import files
estimated_files %>%
  map(~ readRDS(.x) ) %>%
  imap(~.x %>% rename(!!as.symbol(sprintf("results_sccomp%s", .y)) := results_sccomp)) %>%
  reduce(left_join) %>%


  # Calculate sgnificance
  mutate(hypothesis_sccomp1 = map(
    results_sccomp1 ,
    ~ .x  %>%
      filter(covariate == "type") |>
      arrange(c_pH0) %>%
      mutate(probability = 1-c_pH0) %>%
      mutate(estimate = c_effect )
  )) %>%

  # Calculate sgnificance
  mutate(hypothesis_sccomp2 = map(
    results_sccomp2 ,
    ~ .x  %>%
      filter(covariate == "type") |>
      arrange(c_pH0) %>%
      mutate(probability = 1-c_pH0) %>%
      mutate(estimate = c_effect )
  )) %>%


  # Calculate sgnificance
  mutate(hypothesis_sccomp3 = map(
    results_sccomp3 ,
    ~ .x  %>%
      filter(covariate == "type") |>
      arrange(c_pH0) %>%
      mutate(probability = 1-c_pH0) %>%
      mutate(estimate = c_effect )
  )) %>%

  # Calculate sgnificance
  mutate(hypothesis_sccomp4 = map(
    results_sccomp4 ,
    ~ .x  %>%
      filter(covariate == "type") |>
      arrange(c_pH0) %>%
      mutate(probability = 1-c_pH0) %>%
      mutate(estimate = c_effect )
  )) %>%

  # Clean
  dplyr::select(-contains("results")) %>%

  # Reshape
  pivot_longer(contains("hypothesis"),names_prefix = "hypothesis_" ) %>%

  # Drop empty inference because failed
  filter(map_int(value, ~ nrow(.x))>0) %>%

  # Calculate ROC
  mutate(df_for_roc = map2(
    data, value,
    ~ left_join(
      .x %>%
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

