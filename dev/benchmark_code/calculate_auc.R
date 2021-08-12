library(tidyverse)
library(magrittr)
library(bayestestR)
library(yardstick)


# Read arguments
args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_file = args[2]

readRDS(input_file) %>%

  mutate(auc = map_dbl(
    df_for_roc,
    ~ .x %>%
      roc_auc(event_level = 'second', truth = significant, probability) %>%
      pull(.estimate)
  )) %>%
  select(-df_for_roc) %>%
  saveRDS(output_file)


