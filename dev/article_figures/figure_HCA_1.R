library(tidyverse)
library(forcats)
library(HCAquery)
library(dittoSeq)
library(sccomp)
library(magrittr)
library(patchwork)
library(glue)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

## from http://tr.im/hH5A


softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }

  exp(x - logsumexp(x))
}

clean_names = function(x){
  x |>  mutate(
    tissue_harmonised =
      tissue_harmonised |>
      str_remove("tissue_harmonised") |>
      str_replace_all("_", " ") |>
      str_replace("gland", "gld") |>
      str_replace("node", "nd") |>
      str_replace("skeletal", "sk")
  )
}

# General stats
# General stats

# set.seed(42)
# library(ggupset)
# upset_summary_plot =
#   get_metadata() |>
#   left_join(
#     read_csv("~/PostDoc/HCAquery/dev/metadata_cell_type.csv"),
#     by = "cell_type",
#     copy=TRUE
#   ) |>
#   mutate(is_immune = if_else(lineage_1 == "immune" & !is.na(lineage_1), "Immune", "Non immune")) |>
#   mutate(is_healthy = if_else(disease == "normal", "Healthy", "Disease")) |>
#   mutate(is_primary = if_else(is_primary_data.x=="TRUE", "Primary", "Secondary")) |>
#   select(.cell, file_id, is_primary, is_immune, is_healthy) |>
#   as_tibble()  |>
#   mutate(category = pmap(
#     list(is_primary, is_immune, is_healthy),
#     ~ c(..1, ..2, ..3)
#   )) |>
#   sample_n(2e5) |>
#   ggplot(aes(x=category)) +
#   geom_bar() +
#   scale_x_upset(n_intersections = 20) +
#   scale_y_continuous(labels = function(x) format(x * 144.8768, scientific = TRUE, digits=2)) +
#   theme_multipanel +
#   theme(axis.title.x = element_blank())
#
# upset_summary_plot |> saveRDS("~/PostDoc/HCAquery/dev/upset_summary_plot.rds")

upset_summary_plot = readRDS("~/PostDoc/HCAquery/dev/upset_summary_plot.rds")


cell_metadata_with_harmonised_annotation =
  readRDS("~/PostDoc/HCAquery/dev/cell_metadata_with_harmonised_annotation.rds")

data_for_plot_1 =
  cell_metadata_with_harmonised_annotation |>

  left_join(
    get_metadata("~/PostDoc/HCAquery/dev/metadata.sqlite") |>
      dplyr::select(.cell, is_primary_data.y, name, cell_type, file_id, assay) |>
      as_tibble()
  ) |>

  # Count samples
  distinct(.sample, tissue_harmonised, file_id, assay) |>
  add_count(tissue_harmonised, name = "Sample count") |>

  # Add colours
  nest(data = -assay) |>
  arrange(assay) |>
  mutate( color = RColorBrewer::brewer.pal(7,"Blues") |> c(dittoSeq::dittoColors()[-2][1:7]) ) |>
  unnest(data)

# - Number of samples per tissue
plot_sample_dataset =
  data_for_plot_1 |>
  ggplot(aes(fct_reorder(tissue_harmonised, dplyr::desc(`Sample count`)))) +
  geom_bar(aes(fill = assay)) +
  xlab("Tissue") +
  ylab("N. samples (log10)") +
  scale_y_log10() +
  scale_fill_manual(values = data_for_plot_1 |> distinct(assay, color) |> deframe()) +
  theme_multipanel +
  theme(axis.text.x = element_blank(),, axis.title.x = element_blank(), axis.ticks.x = element_blank())


# - Number of datasets per tissue
plot_count_dataset =
  data_for_plot_1 |>
  clean_names() |>
  distinct(file_id, tissue_harmonised, `Sample count`, assay) |>
  ggplot(aes(fct_reorder(tissue_harmonised, dplyr::desc(`Sample count`)))) +
  geom_bar(aes(fill = assay)) +
  xlab("Tissue") +
  ylab("N. datasets") +
  scale_y_reverse() +
  scale_fill_manual(values = data_for_plot_1 |> distinct(assay, color) |> deframe()) +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# - Histogram of cells per sample
plot_cell_dataset =
  cell_metadata_with_harmonised_annotation |>

  left_join(
    get_metadata("~/PostDoc/HCAquery/dev/metadata.sqlite")|>
      dplyr::select(.cell, is_primary_data.y, name, cell_type, file_id, assay) |>
      as_tibble()
  )  |>
  dplyr::count(.sample, assay) |>
  ggplot(aes(n)) +
  geom_histogram(aes(fill=assay), bins = 100) +
  scale_fill_manual(values = data_for_plot_1 |> distinct(assay, color) |> deframe()) +
  xlab("Cells in sample") +
  ylab("Count instances") +
  scale_x_log10(limits = c(10, 1e5)) +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7c78b50dce501fc7ce0b2a8d8efd3aded91134aa/color_cell_types.R")


plot_count_confidence =
  get_metadata("~/PostDoc/HCAquery/dev/metadata.sqlite") |>
  filter(cell_type_harmonised != "non_immune") |>
  mutate(confidence_class = if_else(confidence_class == 4 | is.na(confidence_class), "Low", "High")) |>

  ggplot(aes(fct_relevel(confidence_class, c("Low", "High")), fill = cell_type_harmonised)) +
  geom_bar() +
  scale_fill_manual(values = color_array) +
  coord_flip() +
  xlab("Confidence class") +
  ylab("Cell count") +
  theme_multipanel

plot_proportion_confidence =
  get_metadata("~/PostDoc/HCAquery/dev/metadata.sqlite") |>
  filter(cell_type_harmonised != "non_immune") |>
  select(.cell, tissue_harmonised, cell_type_harmonised, confidence_class) |>
  mutate(confidence_class = if_else(confidence_class == 4 | is.na(confidence_class), "Low", "High")) |>
  as_tibble() |>
  dplyr::count(confidence_class, tissue_harmonised) |>
  with_groups(tissue_harmonised, ~ .x |> mutate(proportion = n/(sum(n)))) |>

  nest(data = -tissue_harmonised) |>
  mutate(proportion_level_1 = map_dbl(data, ~ .x |> filter(confidence_class =="High") |> pull(proportion) )) |>
  unnest(data) |>

  clean_names() |>

  ggplot(aes(fct_reorder(tissue_harmonised,1- proportion_level_1), proportion, fill = as.factor(confidence_class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FDE725FF", "grey")) +
  xlab("Tissue") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# Technology
softmax <- function (x) {
  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }

  exp(x - logsumexp(x))
}


# # Relative
res_relative_contrasts = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_assay.rds")

contrasts = c(
  "assay10x_3_v3 - assay10x_3_v2" ,
  "assaysci_RNA_seq - assay10x_3_v2" ,
  "assaymicrowell_seq - assay10x_3_v2" ,
  "assay10x_5_v2 - assay10x_3_v2" ,
  "assay10x_5_v1 - assay10x_3_v2" ,
  "assaySmart_seq2  - assay10x_3_v2"
)

library(ggforce)
library(ggpubr)


circle_plot_slope = function(res){

  logsumexp <- function (x) {
    y = max(x)
    y + log(sum(exp(x - y)))
  }
  softmax <- function (x) {

    exp(x - logsumexp(x))
  }

  res_relative_for_plot =
    res |>
    filter(!parameter |> str_detect("group___")) |>

    # Cell type abundance
    with_groups(cell_type_harmonised, ~ .x |>  mutate(cell_type_mean_change = mean(abs(c_effect)))) |>

    # Filter for visualisation
    filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>

    # Tissue diversity
    with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>

    # First rank
    with_groups(cell_type_harmonised, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

    # Cap
    mutate(c_effect = c_effect |> pmax(-5) |> pmin(5))

  inter_type_diversity_plot =
    res_relative_for_plot |>
    distinct(inter_type_diversity, tissue_harmonised) |>
    ggplot(aes(inter_type_diversity, fct_reorder(tissue_harmonised, inter_type_diversity))) +
    geom_bar(stat = "identity") +
    scale_x_reverse() +
    xlab("Diversity from baseline") +
    theme_multipanel

  cell_type_mean_abundance_plot =
    res_relative_for_plot |>
    distinct(cell_type_mean_change, cell_type_harmonised) |>
    ggplot(aes(fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_change)), cell_type_mean_change)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(position = "right") +
    theme_multipanel +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
    )


  circle_plot =
    res_relative_for_plot |>
    arrange(rank==1) |>
    ggplot() +
    geom_point(aes(
      fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_change)),
      fct_reorder(tissue_harmonised, inter_type_diversity) ,
      fill = c_effect,
      stroke=(c_lower * c_upper)>0
    ), shape=21, size = 3
    ) +
    scale_fill_distiller(palette="Spectral") +
    theme_multipanel +
    theme(
      axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )


  plot_spacer() +
    cell_type_mean_abundance_plot +
    inter_type_diversity_plot +
    circle_plot +
    plot_layout(guides = 'collect', height = c(1,11), width = c(1, 8)) &
    theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

}

library(tidyHeatmap)
library(ComplexHeatmap)
plot_circle_relative_assay =
  res_relative_contrasts |>
  filter(parameter %in% contrasts) |>
  mutate(parameter = parameter |> str_remove_all("assay")) |>
  mutate(parameter = parameter |> str_remove("- 10x_3_v2")) |>
  mutate(tissue_harmonised = parameter) |>

  # Calculate stats
  filter(!parameter |> str_detect("group___")) |>

  # Cell type abundance
  with_groups(cell_type_harmonised, ~ .x |>  mutate(cell_type_mean_change = mean(abs(c_effect)))) |>

  # Filter for visualisation
  filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>

  # Tissue diversity
  with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>

  # First rank
  with_groups(cell_type_harmonised, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

  # # Cap
  # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>

  mutate(Difference = c_effect) |>
  mutate(inter_type_diversity = -inter_type_diversity) |>
  rename(`Mean diff` = cell_type_mean_change) |>
  rename(Diversity =inter_type_diversity ) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

  # Order
  mutate(parameter = fct_reorder(parameter, Diversity)) |>
  mutate(cell_type_harmonised = fct_reorder(cell_type_harmonised, -`Mean diff`)) |>


  # Heatmap
  heatmap(
    parameter, cell_type_harmonised, Difference,
    palette_value = circlize::colorRamp2(
      seq(3, -3, length.out = 11),
      RColorBrewer::brewer.pal(11, "RdBu")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    show_heatmap_legend = FALSE
  ) |>
  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
  annotation_bar(Diversity, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
  layer_point((c_lower * c_upper)>0)

# circle_plot_slope()

# job::job({ plot_summary_res_relative = res_relative |> plot_summary() })

# # PCA
# library(tidybulk)
#
# data_for_assay_pca =
#   res_relative |>
#   sccomp:::get_abundance_contrast_draws(
#     contrasts = c(
#       "assay10x_3_v2",
#       "assay10x_3_v3" ,
#       "assaysci_RNA_seq" ,
#       "assaymicrowell_seq" ,
#       "assay10x_5_v2" ,
#       "assay10x_5_v1" ,
#       "assaySmart_seq2"
#     )
#   )|>
#   unite("sample", c(parameter, .draw), remove = FALSE) |>
#   reduce_dimensions(
#     sample, cell_type_harmonised, .value,
#     method="PCA", action="get", scale=FALSE,
#     transform = identity
#   )
#
# data_for_assay_pca |> saveRDS("~/PostDoc/HCAquery/dev/data_for_assay_pca.rds")

data_for_assay_pca = readRDS("~/PostDoc/HCAquery/dev/data_for_assay_pca.rds")

plot_assay_PCA =
  data_for_assay_pca |>
  mutate(parameter = parameter |> str_remove_all("assay")) |>
  ggplot(aes(PC1, PC2, label = parameter, color = parameter)) +
  stat_density_2d(geom = "polygon", aes(fill = parameter), breaks = c(0.02), alpha = 0.2) +
  geom_point(aes(color = parameter), data = data_for_assay_pca |> with_groups(parameter, ~ .x |> summarise(across(c(PC1, PC2), median) ))) +
  scale_fill_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all(" |-", "_") |> str_remove("'")) |>  deframe() ) +
  scale_color_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all(" |-", "_") |> str_remove("'")) |>  deframe() ) +
  guides(fill="none", color = "none")  +
  theme_multipanel



# Variability
res_relative_for_variability_plot =
  res_relative |>
  test_contrasts(
    contrasts = c(
      "assay10x_3_v2",
      "assay10x_3_v3" ,
      "assaysci_RNA_seq" ,
      "assaymicrowell_seq" ,
      "assay10x_5_v2" ,
      "assay10x_5_v1" ,
      "assaySmart_seq2"
    )
  ) |>
  filter(parameter |> str_detect("^assay")) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  mutate(parameter = parameter |> str_remove("assay")) |>
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>

  # Summarise and rank
  with_groups(parameter, ~ .x |> mutate(gran_mean = mean(v_effect))) |>
  arrange(desc(v_effect)) |>
  mutate(rank = formatC(1:n(), width = 3, format = "d", flag = "0")) |>

  mutate(cell_type_harmonised_label = glue("{rank}__{cell_type_harmonised}"))

source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7c78b50dce501fc7ce0b2a8d8efd3aded91134aa/color_cell_types.R")


plot_variability_error_bar_per_cell =
  res_relative_for_variability_plot |>

  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

  # Select the most variable cell types

  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |> with_groups(parameter, ~ .x |> arrange(desc(v_effect)) |>  slice_head(n=3)) |>

  ggplot(aes(v_effect, cell_type_harmonised_label)) +
  geom_errorbar(aes(xmin=v_lower, xmax=v_upper, color=cell_type_harmonised), ) +
  #geom_point() +
  facet_grid(fct_reorder(parameter, gran_mean) ~ ., scales = "free") +
  scale_color_manual(values = color_array) +
  scale_y_discrete(  labels = function(x) x |> str_remove("^[0-9]+__")   ) +
  xlab("Variability") +
  ylab("Cell group") +
  guides(color = "none") +
  theme_multipanel

library(ggpubr)

plot_variability_density =
  res_relative_for_variability_plot |>

  # Clean
  mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
  mutate(parameter = parameter |> str_replace_all("_", " ")) |>

  mutate(v_sd = (v_upper-v_effect)/qnorm(0.95)) |>
  mutate(distribution = map2(v_effect, v_sd, ~ rnorm(10, mean = .x, sd = .y))) |>
  unnest(distribution) |>
  ggplot(aes(distribution, color = parameter)) +
  geom_density(  alpha = 0.3) +
  stat_central_tendency(aes(color = parameter), type = "median", linetype = 2) +
  scale_color_manual(values = data_for_plot_1 |> distinct(assay, color) |> mutate(assay = assay |> str_replace_all("-", " ") |> str_remove("'")) |>  deframe() ) +
  guides(color = "none") +
  theme_multipanel


first_line =
  (
    plot_spacer() |
      upset_summary_plot |
      plot_spacer()
  ) +
  plot_layout(guides = 'collect', width = c( 119, 37, 23) )

second_line =
  (
    (plot_sample_dataset / plot_count_dataset) |
      plot_cell_dataset |
      (plot_count_confidence / plot_proportion_confidence) |
      plot_spacer()
  )+

  plot_layout(guides = 'collect', width = c( 1, 0.5, 1, 0.5) )  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

third_line =
  (
    plot_assay_PCA |
      wrap_heatmap(plot_circle_relative_assay, padding = unit(c(-30, 0, -3, -30), "points" )) |
      ( ( plot_variability_density / plot_variability_error_bar_per_cell ) +  plot_layout(heights = c(1,4) ) )
  ) +
  plot_layout( guides = 'collect', widths = c(1, 2, 1) ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

p =
  (
    first_line /
    second_line /
    third_line
  ) +
  plot_layout( guides = 'collect', heights = c(1, 0.5, 1) ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

ggsave(
  "~/PostDoc/sccomp_dev/dev/article_figures/figure_HCA_1.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 220 ,
  limitsize = FALSE
)

