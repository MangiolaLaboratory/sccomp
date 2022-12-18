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
output_blood = args[[4]]
output_file_2 = args[[5]]

my_data =
  readRDS(input_file) |>

  # Fix groups
  unite("group", c(file_id, assay), remove = FALSE) |>

  # filter assays that not have many tissues
  nest(data = -assay) |>
  mutate(n_tissues = map_int(data, ~ .x |> distinct(tissue_harmonised) |> nrow())) |>
  filter(n_tissues > 3) |>
  unnest(data) |>

  # Order by data size
  add_count(assay) |>
  mutate(assay = assay |> fct_reorder(desc(n)))

if(filter_blood=="TRUE"){


  res_relative_blood =
    my_data |>

    # nly blood
    filter(tissue_harmonised=="blood") |>

    # Estimate
    sccomp_glm(
      formula_composition = ~ 0 + assay + sex + ethnicity + age_days   + (1 | group),
      formula_variability = ~ 0 + assay + tissue_harmonised + sex + ethnicity,
      .sample, cell_type_harmonised, counts_from_tissue,
      check_outliers = F,
      approximate_posterior_inference = FALSE,
      cores = 10,
      mcmc_seed = 42,
      verbose = T,
      prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
    )

  res_relative_blood |> saveRDS(output_blood)

  predicted_blood_composition =
    res_relative_blood |>
    sccomp_predict(
      formula_composition = ~ 0 + assay + sex + ethnicity + age_days   ,
      number_of_draws = 500,
      new_data =
        my_data |>
        distinct(.sample, sex , ethnicity , age_days, assay )
    )

  predicted_proportion_of_naive =
    predicted_blood_composition |>
    mutate(is_naive = cell_type_harmonised |> str_detect("naive|stem")) |>
    with_groups(c(.sample, is_naive), ~ .x |> summarise(predicted_proportion = sum(proportion_mean))) |>
    filter(is_naive)

  counts_to_subtract =
    my_data |>
    mutate(is_naive = cell_type_harmonised |> str_detect("naive|stem")) |>
    with_groups(.sample, ~ .x |> mutate(exposure = sum(counts_from_tissue))) |>
    with_groups(c(.sample, exposure, is_naive), ~ .x |> summarise(n = sum(counts_from_tissue))) |>
    complete(nesting(.sample, exposure), is_naive, fill = list(n = 0)) |>
    with_groups(c(.sample), ~ .x |> mutate(observed_proportion = n/sum(n))) |>
    filter(is_naive) |>
    left_join(predicted_proportion_of_naive) |>
    mutate(blood_contamination = observed_proportion / predicted_proportion) |>
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
			formula_composition = ~ 0 + assay + tissue_harmonised + sex + ethnicity + age_days + (assay | group) ,
			formula_variability = ~ 0 + assay + tissue_harmonised + sex + ethnicity,
			.sample, cell_type_harmonised,
			check_outliers = F,
			approximate_posterior_inference = FALSE,
			contrasts = c(
        "assay10x_3_v3 - assay10x_3_v2" ,
        "assaysci_RNA_seq - assay10x_3_v2" ,
        "assaymicrowell_seq - assay10x_3_v2" ,
        "assay10x_5_v2 - assay10x_3_v2" ,
        "assay10x_5_v1 - assay10x_3_v2" ,
        "assaySmart_seq2  - assay10x_3_v2"
			),
			cores = 10,
			mcmc_seed = 42,
			verbose = T,
			prior_mean_variable_association = list(intercept = c(3.6539176, 0.5), slope = c(-0.5255242, 0.1), standard_deviation = c(20, 40))
		)


res_relative |>
  saveRDS(output_file_1)

# Remove unwanted variation
res_relative |>
  remove_unwanted_variation(~ 0 + assay, ~ 0 + assay) |>
  saveRDS(output_file_2)
