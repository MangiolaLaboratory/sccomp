

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

# Read arguments
args = commandArgs(trailingOnly=TRUE)
input_1 = args[[1]]
input_2 = args[[2]]
output_common = args[[3]]
output_absolute = args[[4]]
output_relative = args[[5]]

# DATA INPUT
common_data =
  readRDS(input_1) |>

  left_join(
    # get_metadata() |>
    get_metadata(input_2) |>
      dplyr::select(
        .cell,
        cell_type,
        file_id,
        assay,
        age_days,
        development_stage,
        sex,
        ethnicity,
        confidence_class
      ) |>
      as_tibble()
  ) |>

  # Fix hematopoietic misclassification
  mutate(
    cell_type_harmonised = if_else(
      cell_type_harmonised == "non_immune" &
        cell_type |> str_detect("hematopoietic"),
      "stem",
      cell_type_harmonised
    )
  ) |>

  # Filter intestine as it does not fit the small vs large paradigm
  filter(tissue != "intestine") |>
  filter(tissue != "appendix") |>

  # Filter out
  filter(!cell_type |> str_detect("erythrocyte")) |>
  filter(!cell_type |> str_detect("platelet")) |>

  mutate(is_immune = cell_type_harmonised != "non_immune") |>

  # Format covatriates
  mutate(assay = assay |> str_replace_all(" ", "_") |> str_replace_all("-", "_")  |> str_remove_all("'")) |>
  mutate(
    ethnicity = case_when(
      ethnicity |> str_detect("Chinese|Asian") ~ "Chinese",
      ethnicity |> str_detect("African") ~ "African",
      TRUE ~ ethnicity
    )
  ) |>

  # Fix samples with multiple assays
  unite(".sample", c(.sample , assay), remove = FALSE) |>

  # Scale age
  mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric())

# Save
common_data |> saveRDS(output_common)

# - Immune proportion per tissue
common_data |>

  # Filter unrepresented organs
  filter(!tissue_harmonised %in% c("rectum", "thyroid gland", "salival_gland")) |>

  # Filter only whole tissue
  filter(!tissue_harmonised %in% c("blood", "lymph node", "bone", "spleen")) |>


  # Filter Immune enriched dataset
  filter(file_id != "e756c34a-abe7-4822-9a35-55ef12270247") |>
  filter(file_id != "ca4a7d56-739b-4e3c-8ecd-28704914cc14") |>
  filter(
    file_id != "59dfc135-19c1-4380-a9e8-958908273756" |
      !tissue_harmonised != c("intestine small", "intestine large")
  ) |>

  # Filter samples which not include non immune cells
  nest(data = -c(.sample, tissue_harmonised)) |>
  filter(
    map_int(
      data,
      ~ .x |> filter(cell_type_harmonised == "non_immune") |> nrow()
    ) > 0 | tissue_harmonised %in% c("blood", "lymph node", "bone")
  ) |>
  unnest(data)|>
  saveRDS(output_absolute)


# Relative
relative_data =
  common_data |>

  # Filter low confidence
  filter(confidence_class %in% c(1, 2, 3)) |>

  # Filter unrepresented organs
  filter(!tissue_harmonised %in% c("rectum", "thyroid gland", "salival_gland")) |>

  # Filter only immune
  filter(is_immune) |>
  filter(development_stage != "unknown") |>
  saveRDS(output_relative)



