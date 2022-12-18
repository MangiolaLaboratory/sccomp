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
input_blood_proportion = args[[3]]
output_file_1 = args[[4]]
output_file_2 = args[[5]]

my_data =
  readRDS(input_file) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE) |>

  # Count
  count(
    .sample, sex, ethnicity,
    age_days, assay , group,
    tissue_harmonised, cell_type_harmonised,
    name="counts_from_tissue"
  )

if(filter_blood=="TRUE"){

  counts_to_subtract =
    my_data |>

    with_groups(.sample, ~ .x |> mutate(exposure = sum(counts_from_tissue))) |>
    with_groups(c(.sample, exposure, is_naive), ~ .x |> summarise(n = sum(counts_from_tissue))) |>
    complete(nesting(.sample, exposure), is_naive, fill = list(n = 0)) |>
    left_join(blood_contamination) |>
    mutate(total_blood_count = exposure * blood_contamination) |>
    left_join(predicted_blood_composition) |>
    mutate(counts_from_blood = floor(proportion_mean * total_blood_count)) |>
    select(.sample, cell_type_harmonised, counts_from_blood)

  my_data =
    my_data |>
    left_join(counts_to_subtract) |>
    mutate(counts_blood_free = if_else(
      !tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone", "thymus"),
      counts_from_tissue - counts_from_blood,
      counts_from_tissue
    )) |>
    mutate(counts_blood_free = pmax(counts_blood_free, 0))


}

res_relative =
  my_data |>

  # Estimate
  sccomp_glm(
    formula_composition = ~ 0 + tissue_harmonised + sex + ethnicity  + age_days + assay + (tissue_harmonised | group),
    formula_variability = ~ 0 + tissue_harmonised + sex + ethnicity,
    .sample, cell_type_harmonised, counts_blood_free,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 10,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
  )

res_relative |> saveRDS(output_file_1)

# Remove unwanted variation
res_relative |>
		remove_unwanted_variation(~ 0 + tissue_harmonised, ~ 0 + tissue_harmonised) |>
		saveRDS(output_file_2)

