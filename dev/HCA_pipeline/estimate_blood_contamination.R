library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source(
  "https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel"
)

args = commandArgs(trailingOnly = TRUE)
input_file_relative = args[[1]]
input_file_absolute = args[[2]]
output_blood_fit = args[[3]]
output_blood_proportion = args[[4]]

my_data =
  readRDS(input_file_relative) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE)

# Estimate the blood composition according to demography
res_relative_blood =
  my_data |>

  # only blood
  filter(tissue_harmonised == "blood") |>

  # Estimate
  sccomp_glm(
    formula_composition = ~ sex + ethnicity  + age_days + assay + (1 | group),
    formula_variability = ~ sex + ethnicity,
    .sample,
    cell_type_harmonised,
    check_outliers = F,
    approximate_posterior_inference = FALSE,
    cores = 10,
    mcmc_seed = 42,
    verbose = T,
    prior_mean_variable_association = list(
      intercept = c(3.6539176, 0.5),
      slope = c(-0.5255242, 0.1),
      standard_deviation = c(20, 40)
    )
  )

# Save
res_relative_blood |> saveRDS(output_blood_fit)

# Predict blood composition
predicted_blood_composition =
  res_relative_blood |>
  sccomp_predict(
    formula_composition = ~ sex + ethnicity  + age_days + assay,
    number_of_draws = 500,
    new_data =
      my_data |>
      distinct(.sample, sex , ethnicity , age_days, assay)
  )

 # How many naive are in the blood
predicted_proportion_of_naive =
  predicted_blood_composition |>
  mutate(is_naive = cell_type_harmonised |> str_detect("naive")) |>
  with_groups(c(.sample, is_naive),
              ~ .x |> summarise(predicted_proportion = sum(proportion_mean))) |>
  filter(is_naive)

# Predicted proportion of the immune system
predicted_proportion_of_the_immune_system =
  my_data |>

  # Count
  count(
    .sample,
    sex,
    ethnicity,
    age_days,
    assay ,
    group,
    tissue_harmonised,
    cell_type_harmonised,
    name = "counts_from_tissue"
  ) |>

  mutate(is_naive = cell_type_harmonised |> str_detect("naive")) |>
  with_groups(c(.sample, is_naive), ~ .x |> summarise(n = sum(counts_from_tissue))) |>
  complete(nesting(.sample), is_naive, fill = list(n = 0)) |>
  with_groups(c(.sample), ~ .x |> mutate(observed_proportion = n / sum(n))) |>
  filter(is_naive) |>
  left_join(predicted_proportion_of_naive) |>
  mutate(blood_contamination_of_immune_system = observed_proportion / predicted_proportion)  |>

  select(.sample, blood_contamination_of_immune_system)

# Predicted proportion of the whole tissue
predicted_proportion_of_the_whole_tissue =
  readRDS(input_file_absolute) |>

  # Filter high confidence for  immune
  filter(confidence_class %in% c(1, 2, 3) | cell_type_harmonised == "non_immune") |>

  # Count
  count(
    .sample,
    sex,
    ethnicity,
    age_days,
    assay ,
    tissue_harmonised,
    cell_type_harmonised,
    name = "counts_from_tissue"
  ) |>
  mutate(is_naive = cell_type_harmonised |> str_detect("naive")) |>
  with_groups(c(.sample, is_naive), ~ .x |> summarise(n = sum(counts_from_tissue))) |>
  complete(nesting(.sample), is_naive, fill = list(n = 0)) |>
  with_groups(c(.sample), ~ .x |> mutate(observed_proportion = n / sum(n))) |>
  filter(is_naive) |>
  left_join(predicted_proportion_of_naive) |>
  mutate(blood_contamination_of_whole_tissue = observed_proportion / predicted_proportion)  |>

  select(.sample, blood_contamination_of_whole_tissue)

# Join data
predicted_proportion_of_the_immune_system |>
  left_join(predicted_proportion_of_the_whole_tissue) |>

  # Join blood composition
  left_join(
    predicted_blood_composition |>
      select(.sample, cell_type_harmonised, proportion_blood_cell_type = proportion_mean)
  ) |>

  # Save
  saveRDS(output_blood_proportion)

# # Plots for feedback
# predicted_proportion_of_the_immune_system =
#   readRDS("~/PostDoc/HCAquery/dev/run_dec_2022/blood_contamination.rds") |>
#   distinct(.sample, blood_contamination_of_immune_system, blood_contamination_of_whole_tissue) |>
#   left_join(my_data |> distinct(.sample, tissue_harmonised))
#
# predicted_proportion_of_the_immune_system |>
#   left_join(predicted_proportion_of_the_whole_tissue) |>
#   select(.sample, tissue_harmonised, blood_contamination_of_whole_tissue) |>
#   mutate(blood_contamination_of_whole_tissue = pmin(blood_contamination_of_whole_tissue, 1)) |>
#   with_groups(tissue_harmonised, ~ .x |> mutate(mean_contamination = median(blood_contamination_of_whole_tissue, na.rm = T))) |>
#   filter(!tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>
#   ggplot(aes(fct_reorder(tissue_harmonised, desc(mean_contamination)), blood_contamination_of_whole_tissue)) +
#   geom_boxplot(outlier.shape = NA) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=30)) +
#   ggtitle("% of blood in the whole tissue")
#
# predicted_proportion_of_the_immune_system |>
#   left_join(predicted_proportion_of_the_whole_tissue) |>
#   select(.sample, tissue_harmonised, blood_contamination_of_immune_system) |>
#   mutate(blood_contamination_of_immune_system = pmin(blood_contamination_of_immune_system, 1)) |>
#   with_groups(tissue_harmonised, ~ .x |> mutate(mean_contamination = median(blood_contamination_of_immune_system, na.rm = T))) |>
#   filter(!tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen", "thymus")) |>
#   ggplot(aes(fct_reorder(tissue_harmonised, desc(mean_contamination)), blood_contamination_of_immune_system)) +
#   geom_boxplot(outlier.shape = NA) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=30)) +
#   ggtitle("% of blood in the immune system")
