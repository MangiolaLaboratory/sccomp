library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

args = commandArgs(trailingOnly=TRUE)
filter_blood = args[[1]]
input_file = args[[2]]
output_file_1 = args[[3]]
output_file_2 = args[[4]]

differential_composition_sex_absolute =
  readRDS("~/PostDoc/HCAquery/dev/data_for_immune_proportion_sex.rds") |>

  # filter
  filter(development_stage!="unknown") |>
  filter(sex != "unknown") |>
  filter(age_days >= 19 * 365) |>

  # Keep shared tissues
  nest(data = -c(.sample, sex, tissue_harmonised)) |>
  add_count(sex, tissue_harmonised) |>
  filter(n>2) |>
  select(-n) |>
  nest(data = -c(sex, tissue_harmonised)) |>
  add_count(tissue_harmonised) |>
  filter(n==2) |>
  unnest(data) |>
  unnest(data) %>%

  # Filter tissue only having > 2 datasets
  inner_join(
    (.) |>
      distinct( sex, tissue_harmonised, file_id) |>
      count(sex, tissue_harmonised, name = "count_file_id") |>
      with_groups(tissue_harmonised, ~ .x |> mutate(min_count_file_id = min(count_file_id))) |>
      filter(min_count_file_id>1)
  ) |>

  unite("group", c(tissue_harmonised , file_id, sex), remove = FALSE) |>
  unite("tissue_harmonised_sex", c(tissue_harmonised , sex), remove = FALSE) |>

  # Estimate
  sccomp_glm(
    formula_composition = ~ sex + tissue_harmonised + ethnicity  + age_days +  assay + (1 | group) + (sex | tissue_harmonised_sex),
    formula_variability = ~ sex + tissue_harmonised + ethnicity,
    .sample, is_immune,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 20,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
  )

differential_composition_sex_absolute |>
  saveRDS(output_file_1)



differential_composition_sex_absolute |>
		remove_unwanted_variation(~ sex, ~ sex) |>
		saveRDS(output_file_2)
