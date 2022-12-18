library(tidyverse)
library(glue)
library(here)
library(stringr)
library(Seurat)
library(tidyseurat)

args = commandArgs(trailingOnly=TRUE)
run_directory = args[[1]]

R_code_directory = glue("~/PostDoc/sccomp/dev/HCA_pipeline")
root = "~/PostDoc/HCAquery/dev"
tab = "\t"

run_directory |> dirname() |> dir.create( showWarnings = FALSE, recursive = TRUE)

commands = c()

# Create input
commands =
  commands |>
  c(
    glue("CATEGORY=create_input\nMEMORY=30024\nCORES=1\nWALL_TIME=10000"),
    glue("{run_directory}/input_common.rds {run_directory}/input_absolute.rds {run_directory}/input_relative.rds:{root}/cell_metadata_with_harmonised_annotation.rds {root}/metadata.sqlite\n{tab}Rscript {R_code_directory}/create_common_data.R {root}/cell_metadata_with_harmonised_annotation.rds {root}/metadata.sqlite {run_directory}/input_common.rds {run_directory}/input_absolute.rds {run_directory}/input_relative.rds")
  )

# Estimate blood contamination
commands =
  commands |>
  c(
    glue("CATEGORY=create_input\nMEMORY=30024\nCORES=10"),
    glue("{run_directory}/blood_fit.rds {run_directory}/blood_contamination.rds:{run_directory}/input_relative.rds {run_directory}/input_absolute.rds\n{tab}Rscript {R_code_directory}/estimate_blood_contamination.R {run_directory}/input_relative.rds {run_directory}/input_absolute.rds {run_directory}/blood_fit.rds {run_directory}/blood_contamination.rds")
  )

# estimate
commands =
  commands |>
  c(
    glue("CATEGORY=estimate\nMEMORY=30024\nCORES=10"),
    tibble(
      factor = c("tissue", "sex", "ethnicity", "assay", "age")
    ) |>
      expand_grid(modality = c("absolute", "relative"), filter_blood = c("TRUE", "FALSE")) |>
      filter(!(factor == "assay" & modality=="absolute")) |>
      mutate(analysis = glue("{factor}_{modality}")) |>
      mutate(r_script = glue("{R_code_directory}/sccomp_{analysis}.R")) |>
      mutate(input = glue("input_{modality}.rds")) |>
      mutate(command = if_else(
        filter_blood == "TRUE",
        glue("{run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds:{run_directory}/{input}\n{tab}Rscript {R_code_directory}/{r_script} {filter_blood} {run_directory}/{input} {run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds"),
        glue("{run_directory}/{analysis}.rds {run_directory}/{analysis}_proportion_adjusted.rds:{run_directory}/{input}\n{tab}Rscript {R_code_directory}/{r_script} {filter_blood} {run_directory}/{input} {run_directory}/{analysis}.rds {run_directory}/{analysis}_blood.rds {run_directory}/{analysis}_proportion_adjusted.rds")
      )
      )
  )


commands |>
  write_lines(glue("{run_directory}/pipeline.makeflow"))
