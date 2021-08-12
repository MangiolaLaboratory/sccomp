library(tidyverse)
library(magrittr)
library(bayestestR)
library(yardstick)


# Read arguments
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_file = args[2]

input = readRDS(input_file)

input %>%

  # Add random
  bind_rows(
     input %>%
       slice(1) %>%
       mutate(name = "random") %>%
       mutate(df_for_roc  = map(df_for_roc, ~ .x %>% mutate(probability = runif(n(), 0, 1))))
  ) %>%

  # Calculate AUC
  mutate(auc = map_dbl(
    df_for_roc,
    ~ .x %>%
      roc_curve(event_level = 'second', truth = significant, probability) %>%
      filter(1-specificity<=0.1) %>%
      summarise(.estimate = bayestestR::auc(1-specificity , sensitivity)) %>%
      pull(.estimate)
  )) %>%
  select(-df_for_roc) %>%

  saveRDS(output_file)
