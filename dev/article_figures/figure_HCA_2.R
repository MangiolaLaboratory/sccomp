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

## from http://tr.im/hH5A




data_for_immune_proportion = readRDS("~/PostDoc/HCAquery/dev/data_for_immune_proportion.rds")
data_for_immune_proportion_relative = readRDS("~/PostDoc/HCAquery/dev/data_for_immune_proportion_relative.rds")

tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()


source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7c78b50dce501fc7ce0b2a8d8efd3aded91134aa/color_cell_types.R")
cell_type_color = color_array
names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")

# # Save data for third party
# data_for_immune_proportion |>
#
#   # Drop only-immune organs
#   filter(!tissue_harmonised %in% c("blood", "lymph node", "spleen", "bone")) |>
#   mutate(is_immune = as.character(is_immune)) |>
#
#   # Mutate days
#   mutate(age_days = age_days  |> scale(center = FALSE) |> as.numeric()) |>
#   filter(development_stage!="unknown") |>
#   saveRDS("~/PostDoc/sccomp_dev/dev/data_age_absolute_for_third_party.rds")


# Track of immune system in life

differential_composition_age_relative = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_age_relative2.rds")
differential_composition_age = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_age3.rds")

# ggplot() +
#   geom_abline(intercept = -0.281, slope = 0.276)

# generate points
line_age_absolute_mean =
  seq(-3, 3, by = 0.01) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_effect)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_effect)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))

line_age_absolute_lower =
  seq(-3, 3, by = 0.1) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_lower)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_lower)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))

line_age_absolute_upper =
  seq(-3, 3, by = 0.1) |>
  enframe(value = "x") |>
  mutate(
    y = x *
      (
        differential_composition_age |> filter(parameter == "age_days" &
                                                 is_immune == "TRUE") |> pull(c_upper)
      ) +
      (
        differential_composition_age |> filter(parameter == "(Intercept)" &
                                                 is_immune == "TRUE") |> pull(c_upper)
      )
  ) |>
  rowwise() |>
  mutate(proportion = softmax(c(y,-y))[1]) |>
  ungroup() |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))

dropLeadingZero <-
  function(l) {
    stringr::str_replace(l, '0(?=.)', '')
  }
S_sqrt <- function(x) {
  sign(x) * sqrt(abs(x))
}
IS_sqrt <- function(x) {
  x ^ 2 * sign(x)
}
S_sqrt_trans <-
  function()
    scales::trans_new("S_sqrt", S_sqrt, IS_sqrt)

proportions_age_relative = readRDS("~/PostDoc/HCAquery/dev/proportions_age_adjusted_relative.rds")
proportions_age_absolute = readRDS("~/PostDoc/HCAquery/dev/proportions_age_adjusted_absolute.rds")

life_stages = tibble(
  start = c(0, 2, 5, 13, 20, 40, 60),
  end = c(1, 4, 12, 19, 39, 59,100)
) |>
  mutate(stage = c("Infant", "Toddler", "Child", "Teen", "Adult", "Middle age", "Senior"))


# Add life stages
rectangles_age =
  line_age_absolute_mean |>
  mutate(stage = case_when(
  (x_corrected / 365) |> between(0,1) ~ "Infant",
  (x_corrected / 365) |> between(1,4) ~ "Toddler",
  (x_corrected / 365) |> between(4,12) ~ "Child",
  (x_corrected / 365) |> between(12,19) ~ "Teen",
  (x_corrected / 365) |> between(19, 39) ~ "Adult",
  (x_corrected / 365) |> between(39, 59) ~ "Middle age",
  (x_corrected / 365) |> between(59,85) ~ "Senior"
)) |>
  left_join(
    tibble(
      start = c(0, 1, 4, 12, 19, 39, 59),
      end = c(1, 4, 12, 19, 39, 59,85)
    ) |>
      mutate(stage = c("Infant", "Toddler", "Child", "Teen", "Adult", "Middle age", "Senior"))
  ) |>
  with_groups(stage, ~ .x |> mutate(mean_proportion = mean(proportion, na.rm=TRUE))) |>
  distinct(stage, start, end, mean_proportion) |>
  mutate(color =  RColorBrewer::brewer.pal(n = 9, name = "Greys") |> head( n()) )


plot_age_absolute =
  proportions_age_absolute |>
  left_join(data_for_immune_proportion |>
              tidybulk::pivot_sample(.sample)) |>

  filter(development_stage != "unknown") |>

  # Filter
  filter(is_immune == "TRUE") |>
  #filter(tissue_harmonised != "blood") |>
  # Fix samples with multiple assays
  unite(".sample", c(.sample , assay), remove = FALSE) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE) |>

  # Plot
  ggplot() +
  geom_rect(
    aes(xmin = start * 365, xmax = end * 365, ymin = 0, ymax = mean_proportion),
    data = rectangles_age,
    fill = rectangles_age |> select(stage, color) |> deframe(),
    alpha = 0.5
    ) +
  geom_point(
    aes(age_days, adjusted_proportion, fill = tissue_harmonised),
    shape = 21,
    stroke = 0,
    size = 1
  ) +
  geom_line(aes(x_corrected, proportion), data = line_age_absolute_mean) +
  geom_line(aes(x_corrected, proportion),
            data = line_age_absolute_lower,
            color = "grey") +
  geom_line(aes(x_corrected, proportion),
            data = line_age_absolute_upper,
            color = "grey") +
  facet_wrap( ~ is_immune, ncol = 9) +
  scale_fill_manual(values = tissue_color) +
  scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none") +
  theme_multipanel

# Color for human heatmap

# Plot age absolute organ cell type
age_absolute_organ_cell_type =
differential_composition_age |>

  # Find stats of random effect with groups
  sccomp_test(
    contrasts =
      differential_composition_age |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( )
  ) |>
  filter(is_immune == "TRUE")


age_absolute_organ_cell_type |>
  select(parameter, c_effect) |> mutate(color = circlize::colorRamp2(
    seq(1.45,-1.45, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
  )(c_effect)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>
  pull(color) |>
  scales::show_col()

differential_composition_age_relative |>
  sccomp_test(test_composition_above_logit_fold_change = 0.2) |>
  arrange(desc(abs(c_effect))) |> filter(covariate == "age_days") |> filter(c_FDR<0.05)

line_age_relative_mean =

  differential_composition_age_relative |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  nest(data = -cell_type_harmonised) |>
  mutate(x = list(seq(-3, 3, by = 0.1))) |>
  mutate(y = map2(data, x, ~ {
    .y *
      (.x |> filter(parameter == "age_days") |> pull(c_effect)) +
      (.x |> filter(parameter == "(Intercept)") |> pull(c_effect))

  })) |>
  select(-data) |>
  unnest(c(x, y)) |>
  with_groups(x, ~ .x |> mutate(proportion = softmax(y))) |>
  mutate(x_corrected = (x * 9610.807 / 0.6) + 12865.75) |>
  filter(x_corrected |> between(30.0 , 30295.0))



# Plot age relative cell types
plot_age_relative =
  proportions_age_relative |>
  left_join(data_for_immune_proportion_relative |>
              tidybulk::pivot_sample(.sample)) |>

  filter(development_stage != "unknown") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  #filter(tissue_harmonised != "blood") |>
  # Fix samples with multiple assays
  unite(".sample", c(.sample , assay), remove = FALSE) |>

  # Fix groups
  unite("group", c(tissue_harmonised , file_id), remove = FALSE)  |>

  # Filter
  filter(cell_type_harmonised != "immune_unclassified") |>

  # IN FUTURE I WILL TEST BASED ON FOLD CHANGE
  inner_join(
    differential_composition_age_relative |>
      filter(parameter == "age_days") |>
      filter(abs(c_effect)> 0.5) |>
      distinct(cell_type_harmonised)
  ) |>

  ggplot(aes(age_days, adjusted_proportion)) +
  geom_point(
    aes(fill = tissue_harmonised),
    shape = 21,
    stroke = 0,
    size = 0.4
  ) +
  geom_line(
    aes(x_corrected, proportion, color = significant),
    data = line_age_relative_mean |>

      # Join statistics
      inner_join(
        differential_composition_age_relative |>
          sccomp_test(test_composition_above_logit_fold_change = 0.4) |>

          filter(parameter == "age_days") |>
          mutate(significant = c_FDR < 0.05) |>
          filter(cell_type_harmonised != "immune_unclassified") |>
          select(cell_type_harmonised, significant)
      )
  ) +
  facet_wrap( ~ cell_type_harmonised, ncol = 3) +
  scale_y_continuous(trans = S_sqrt_trans(), labels = dropLeadingZero) +
  scale_x_continuous(
    labels = function(x)
      round(x / 356)
  ) +
  scale_fill_manual(values = tissue_color) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  xlab("Years") +
  ylab("Adjusted proportions") +
  guides(fill = "none", color = "none") +
  theme_multipanel


# Plot age relative organ cell type
age_relative_organ_cell_type =

  differential_composition_age_relative |>
  sccomp_test(
    contrasts =
      differential_composition_age_relative |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + {parameter}") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( )
  )



plot_age_relative_by_organ =

  age_relative_organ_cell_type |>

  filter(c_FDR<0.01) |>
  arrange(desc(abs(c_effect))) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  add_count(cell_type_harmonised) |>
  arrange(parameter, desc(n)) |>
  ggplot(aes( fct_reorder(cell_type_harmonised, desc(n)), parameter)) +
  geom_tile(aes(fill = c_effect)) +
  scale_fill_distiller(palette="Spectral") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=20, hjust = 1, vjust = 1))


library(tidyHeatmap)
library(ComplexHeatmap)

df_heatmap_age_relative_organ_cell_type =

  differential_composition_age_relative |>

  # Find stats of random effect with groups
  sccomp_test(
    contrasts =
      differential_composition_age_relative |>
      filter(parameter |> str_detect("___age_days")) |>
      distinct(parameter) |>
      mutate(contrast = glue("age_days + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      deframe( ),
    test_composition_above_logit_fold_change = 0.4
  )  |>

  filter(cell_type_harmonised != "immune_unclassified") |>
  add_count(cell_type_harmonised) |>
  arrange(parameter, desc(n)) |>

  rename(tissue = parameter) |>
  rename(cell_type = cell_type_harmonised) |>

  # Cell type abundance
  with_groups(cell_type, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(cell_type_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # Filter for visualisation
  filter(!cell_type %in% c("non_immune", "immune_unclassified")) |>

  # Tissue diversity
  with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # First rank
  with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

  # # Cap
  # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>

  mutate(Difference = c_effect) |>
  rename(`Mean diff` = cell_type_mean_change) |>
  mutate(`Mean diff tissue` = -tissue_mean_change) |>
  mutate(cell_type = cell_type |> str_replace("macrophage", "macro")) |>
  mutate(tissue = tissue |> str_replace_all("_", " ")) |>

  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>

  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>

  # Order
  mutate(tissue = fct_reorder(tissue, `Mean diff tissue`)) |>
  mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))


plot_heatmap_age_relative_organ_cell_type =

  df_heatmap_age_relative_organ_cell_type |>

  # Heatmap
  heatmap(
    tissue, cell_type, Difference,
    # palette_value = circlize::colorRamp2(
    #   seq(-3, 3, length.out = 11),
    #   RColorBrewer::brewer.pal(11, "Spectral")
    # ),
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

  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>
  annotation_tile(
    tissue, show_legend = FALSE,
    palette =
      df_heatmap_age_relative_organ_cell_type |>
      distinct(tissue, tissue_color) |>
      arrange(tissue) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    cell_type, show_legend = FALSE,
    palette =
      df_heatmap_age_relative_organ_cell_type |>
      distinct(cell_type, cell_type_color)  |>
      arrange(cell_type) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    count_tissue, show_legend = FALSE,
    size = unit(0.2, "cm"),
    palette = circlize::colorRamp2(c(0, 5, 10, 15), viridis::magma(4))
  ) |>
  layer_point((c_lower * c_upper)>0)


plot_heatmap_age_relative_organ_cell_type |>
  save_pdf(
    filename = "~/PostDoc/sccomp_dev/dev/article_figures/plot_heatmap_age_relative_organ_cell_type.pdf",
    width = 78*1.5, height = 78*1.5, units = "mm"
  )

# Color body
df_heatmap_age_relative_organ_cell_type |>
  distinct(tissue, `Mean diff tissue`) |> arrange(`Mean diff tissue`) |>  mutate(color =  circlize::colorRamp2(c(0, 5, 10, 12), viridis::viridis(4))(-`Mean diff tissue`)) |>
  mutate(rgb = map_chr(
    color,
    ~ .x |>
      col2rgb() |>
      paste(collapse = " ")
  )) |>  pull(color) |>
  scales::show_col()

# Significance global statistics
count_significance_age_immune_load =
  differential_composition_age |>
  sccomp_test(test_composition_above_logit_fold_change = 0.1) |>
  filter(covariate=="age_days") |>
  filter(is_immune=="TRUE") |>
  count(c_FDR<0.05)

count_significance_age_immune_load_tissue =
  age_absolute_organ_cell_type1.07
FDR|>
  count(c_FDR<0.05)

count_significance_age_cell_type =
  differential_composition_age_relative |>
  sccomp_test(test_composition_above_logit_fold_change = 0.4) |>
  filter(covariate=="age_days") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)

count_significance_age_cell_type_tissue =
  df_heatmap_age_relative_organ_cell_type |>
  count(c_FDR<0.05)



rm(differential_composition_age_relative , differential_composition_age )
gc()



# Sex
differential_composition_sex_relative = readRDS("~/PostDoc/HCAquery/dev/run_30_nov_2022/immune_non_immune_differential_composition_sex_relative2.rds")
differential_composition_sex_absolute = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_sex_absolute2.rds")


# Summary stats
data_for_immune_proportion_sex |>
  distinct(.sample, sex, tissue_harmonised, file_id) |>

  bind_rows(
    data_for_immune_proportion_sex |>
      distinct(.sample, sex, tissue_harmonised, file_id) |>
      mutate(tissue_harmonised = "overall") %>%
      nest(data = -c(sex, tissue_harmonised)) |>
      mutate(count_file_id = map_int(data, ~ .x |> distinct(file_id) |> nrow())) |>
      unnest(data)
  ) |>
  with_groups(c(sex, tissue_harmonised, count_file_id), ~ .x |> summarise(count = n())) |>
  with_groups(tissue_harmonised, ~ .x |> mutate(count_tissue = sum(count))) |>
  unite("name", c(sex, tissue_harmonised), remove = FALSE) |>
  arrange(desc(count_tissue), tissue_harmonised, sex) %>%
  mutate(name = factor(name, levels = .$name)) |>
  ggplot(aes(name, count)) +
  geom_bar(aes(fill = sex), stat="identity") +
  geom_text(aes(label = count_file_id, y = count)) +
  scale_y_log10() +
  guides(fill="none") +
  ylab("Count (log10)") +
  scale_fill_brewer() +
  theme_multipanel +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    )
  )



proportions_sex_absolute_adjusted = readRDS("~/PostDoc/HCAquery/dev/proportions_sex_absolute_adjusted.rds")

plot_sex_absolute_overall =

  differential_composition_sex_absolute |>
  replicate_data(~ sex, ~ sex, number_of_draws = 20) |>
  left_join(data_for_immune_proportion |> distinct(.sample, sex, age_days)) |>
  mutate(tissue_harmonised = "Main") |>
  rename(proportion = generated_proportions) |>
  filter(is_immune=="TRUE") |>
  ggplot(aes(y = tissue_harmonised, x = proportion)) +
  geom_violin(aes(fill = tissue_harmonised), lwd = 0.2) +
  facet_grid(sex ~.) +
  xlim(NA, 1.0) +
  scale_fill_manual(values = tissue_color) +
  scale_y_discrete(label = function(x){
    x |>
      str_remove("tissue_harmonised") |>
      str_replace_all("_", " ") |>
      str_replace("gland", "gld") |>
      str_replace("node", "nd") |>
      str_replace("skeletal", "sk")
  }) +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()
  )



plot_sex_absolute_tissue =
  proportions_sex_absolute_adjusted |>

  left_join(data_for_immune_proportion |> distinct(.sample, sex, age_days, tissue_harmonised)) |>
  rename(proportion = adjusted_proportion) |>
  filter(is_immune=="TRUE") |>

  ggplot(aes(y = tissue_harmonised, x = proportion)) +
  ggridges::geom_density_ridges(aes(fill = tissue_harmonised), lwd = 0.2) +
  facet_grid(sex ~.) +
  xlim(0, 1.0) +
  scale_fill_manual(values = tissue_color) +
  scale_y_discrete(label = function(x){
    x |>
      str_remove("tissue_harmonised") |>
      str_replace_all("_", " ") |>
      str_replace("gland", "gld") |>
      str_replace("node", "nd") |>
      str_replace("skeletal", "sk")
  }) +
  theme_multipanel


plot_sex_absolute_1D =
  differential_composition_sex_absolute |>
  filter(covariate=="sex") |>
  filter(is_immune == "TRUE") |>
  sccomp:::plot_1d_intervals(is_immune, significance_threshold = 0.025, theme_multipanel) &
  xlim(c(-0.6, 0.6)) &
  geom_vline(xintercept = 0, linetype="dashed") &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank())


# sex_absolute_organ =
#   differential_composition_sex_absolute |>
#
#   # Find stats of random effect with groups
#   sccomp_test(
#     contrasts =
#       differential_composition_sex_absolute |>
#       filter(parameter |> str_detect("___sex")) |>
#       distinct(parameter) |>
#       tidyr::extract( parameter, "sex", "_(.+)___", remove = FALSE) |>
#       mutate(sex = glue("sex{sex}")) |>
#       mutate(contrast = glue("`{sex}`  + `{parameter}`") |> as.character()) |>
#       tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
#       select(tissue_harmonised, contrast) |>
#       filter(!contrast |> str_detect("female")) |>
#       deframe( )
#   )

# # CI Absolute NOTHING TO SEE HERE
# sex_absolute_organ |>
#   filter(is_immune == "TRUE") |>
#   extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>
#   filter(c_FDR<0.05) |>
#
#   # Select significant
#   inner_join(
#     differential_composition_sex_absolute_one_vs_all |>
#       extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>
#       distinct(tissue)
#   ) |>
#   extract(parameter, "tissue", "(.+)_", remove = FALSE) |>
#   ggplot(aes(c_effect, fct_reorder(parameter, dplyr::desc(c_effect)))) +
#   geom_linerange(aes(xmin = c_lower, xmax = c_upper), lwd = 0.2,) +
#   geom_point(size = 0.5, color = "red") +
#   facet_wrap(~tissue, scale="free_y") +
#   xlab("Immune proportion") +
#   ylab("sex") +
#   theme_multipanel


# differential_composition_sex_relative |>
#   sccomp_remove_unwanted_variation(~ sex, ~ sex) |>
#   saveRDS("~/PostDoc/HCAquery/dev/proportions_sex_relative_adjusted.rds")


# proportions_sex_absolute_adjusted |>
#   filter(cell_type_harmonised == "cd14 mono") |>
#   left_join(data_for_immune_proportion_relative |> distinct(.sample, sex, age_days, tissue_harmonised)) |>
#
#   rename(proportion = adjusted_proportion) |>
#
#   ggplot(aes(y = proportion, x = sex)) +
#   geom_violin() +
#   scale_fill_manual(values = tissue_color) +
#   theme_multipanel

sex_relative_organ_cell_type =
  differential_composition_sex_relative |>

  # Find stats of random effect with groups
  sccomp_test(
    contrasts =
      differential_composition_sex_relative |>
      filter(parameter |> str_detect("_male___sex")) |>
      distinct(parameter) |>
      mutate(contrast = glue("sexmale  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
      filter(!contrast |> str_detect("female")) |>
      mutate(tissue_harmonised = tissue_harmonised |> str_remove("_male")) |>
      deframe( )
  )

#
df_heatmap_sex_relative_organ_cell_type =

  sex_relative_organ_cell_type |>

  rename(tissue = parameter) |>
  rename(cell_type = cell_type_harmonised) |>

  # Cell type abundance
  with_groups(cell_type, ~ .x |>  mutate(cell_type_mean_change = mean(abs(c_effect)))) |>

  # Filter for visualisation
  filter(!cell_type %in% c("non_immune", "immune_unclassified")) |>

  # Tissue diversity
  with_groups(tissue, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(tissue_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # First rank
  with_groups(cell_type, ~ .x |> mutate(c_effect_significant = case_when(c_FDR<0.05 ~ c_effect)) |>   mutate(cell_type_mean_change = sum(abs(c_effect_significant), na.rm = TRUE))) |>

  # # Cap
  # mutate(c_effect = c_effect |> pmax(-5) |> pmin(5)) |>

  mutate(Difference = c_effect) |>
  rename(`Mean diff` = cell_type_mean_change) |>
  mutate(`Mean diff tissue` = -tissue_mean_change) |>
  mutate(cell_type = cell_type |> str_replace("macrophage", "macro")) |>
  mutate(tissue = tissue |> str_replace_all("_", " ")) |>

  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  )  |>

  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>

  # Order
  mutate(tissue = fct_reorder(tissue, `Mean diff tissue`)) |>
  mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))
#
#
# library(tidyHeatmap)
# library(ComplexHeatmap)
#
# plot_heatmap_sex_relative_organ_cell_type =
#
#   df_heatmap_sex_relative_organ_cell_type |>
#   # Heatmap
#   heatmap(
#     tissue, cell_type, Difference,
#     # palette_value = circlize::colorRamp2(
#     #   seq(-3, 3, length.out = 11),
#     #   RColorBrewer::brewer.pal(11, "Spectral")
#     # ),
#     palette_value = circlize::colorRamp2(
#       seq(-3, 3, length.out = 11),
#       RColorBrewer::brewer.pal(11, "RdBu")
#     ),
#     cluster_rows = FALSE,
#     cluster_columns = FALSE,
#     row_names_gp = gpar(fontsize = 6),
#     column_names_gp = gpar(fontsize = 6),
#     column_title_gp = gpar(fontsize = 0),
#     row_title_gp = gpar(fontsize = 0),
#     show_heatmap_legend = FALSE
#   ) |>
#
#   annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
#   annotation_bar(`Mean diff tissue`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.8, "cm")) |>
#   annotation_tile(
#     tissue, show_legend = FALSE,
#     palette =
#       df_heatmap_sex_relative_organ_cell_type |>
#       distinct(tissue, tissue_color) |>
#       arrange(tissue) |>
#       deframe(),
#     size = unit(0.3, "cm")
#   ) |>
#   annotation_tile(
#     cell_type, show_legend = FALSE,
#     palette =
#       df_heatmap_sex_relative_organ_cell_type |>
#       distinct(cell_type, cell_type_color)  |>
#       arrange(cell_type) |>
#       deframe(),
#     size = unit(0.3, "cm")
#   ) |>
#   annotation_tile(
#     count_tissue, show_legend = FALSE,
#     size = unit(0.2, "cm"),
#     palette = circlize::colorRamp2(c(0, 5, 10, 15), viridis::viridis(4))
#   ) |>
#   layer_point((c_lower * c_upper)>0)
#


# Significance global statistics
count_significance_sex_immune_load =
  differential_composition_sex_absolute |>
  sccomp_test(test_composition_above_logit_fold_change = 0.1) |>
  filter(covariate=="sex") |>
  filter(is_immune=="TRUE") |>
  count(c_FDR<0.05)

count_significance_sex_immune_load_tissue =
  differential_composition_sex_absolute |>
  sccomp_test(
    contrasts =
      differential_composition_sex_absolute |>
        filter(parameter |> str_detect("_male___sex")) |>
        distinct(parameter) |>
        mutate(contrast = glue("sexmale  + `{parameter}`") |> as.character()) |>
        tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+") |>
        filter(!contrast |> str_detect("female")) |>
        mutate(tissue_harmonised = tissue_harmonised |> str_remove("_male")) |>
        deframe( )
  ) |>
  filter(is_immune == "TRUE") |>
  count(c_FDR<0.05)

count_significance_sex_cell_type =
  differential_composition_sex_relative |>
  sccomp_test(test_composition_above_logit_fold_change = 0.4) |>
  filter(covariate=="sex") |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)

count_significance_sex_cell_type_tissue =
  df_heatmap_sex_relative_organ_cell_type |>
  count(c_FDR<0.05)


rm(differential_composition_sex_relative, differential_composition_sex_absolute)
gc()


# Ethnicicy
circle_plot = function(res) {
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
    with_groups(tissue_harmonised, ~ .x |>  mutate(proportion = softmax(c_effect))) |>
    with_groups(cell_type_harmonised,
                ~ .x |>  mutate(cell_type_mean_abundance = mean(proportion))) |>

    # Filter for visualisation
    filter(!cell_type_harmonised %in% c("non_immune", "immune_unclassified")) |>

    # Tissue diversity
    with_groups(tissue_harmonised, ~ .x |>  mutate(inter_type_diversity = sd(c_effect))) |>

    # First rank
    with_groups(cell_type_harmonised,
                ~ .x |> arrange(desc(c_effect)) |>  mutate(rank = 1:n())) |>

    # Cap
    mutate(c_effect = c_effect |> pmax(-5) |> pmin(5))

  inter_type_diversity_plot =
    res_relative_for_plot |>
    distinct(inter_type_diversity, tissue_harmonised) |>
    ggplot(aes(
      inter_type_diversity,
      fct_reorder(tissue_harmonised, inter_type_diversity)
    )) +
    geom_bar(stat = "identity") +
    scale_x_reverse() +
    xlab("Diversity") +
    ylab("Tissue") +
    theme_multipanel

  cell_type_mean_abundance_plot =
    res_relative_for_plot |>
    distinct(cell_type_mean_abundance, cell_type_harmonised) |>
    ggplot(aes(
      fct_reorder(
        cell_type_harmonised,
        dplyr::desc(cell_type_mean_abundance)
      ),
      cell_type_mean_abundance
    )) +
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
    arrange(rank == 1) |>
    ggplot() +
    geom_point(aes(
      fct_reorder(
        cell_type_harmonised,
        dplyr::desc(cell_type_mean_abundance)
      ),
      fct_reorder(tissue_harmonised, inter_type_diversity) ,
      fill = rank,
      size = c_effect,
      stroke = rank == 1
    ),
    shape = 21) +
    scale_size_continuous(range = c(0.5, 3)) +
    xlab("Cell type") +
    theme_multipanel +
    theme(
      axis.text.x = element_text(
        angle = 30,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank()
    )


  plot_spacer() +
    cell_type_mean_abundance_plot +
    inter_type_diversity_plot +
    circle_plot +
    plot_layout(guides = 'collect',
                height = c(1, 4),
                width = c(1, 8)) &
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          legend.position = "bottom")

}


differential_composition_ethnicity_absolute = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_ethnicity_absolute2.rds")


data_for_ethinicity_absolute_plot =
  differential_composition_ethnicity_absolute |>
  filter(parameter %in% c("ethnicityAfrican",  "ethnicityHispanic or Latin American", "ethnicityEuropean", "ethnicityChinese")) |>

  rowwise() |>
  mutate(
    proportion_mean = softmax(c(c_effect,-c_effect))[1],
    proportion_lower = softmax(c(c_lower,-c_lower))[1],
    proportion_upper = softmax(c(c_upper,-c_upper))[1]
  ) |>
  ungroup() |>
  filter(is_immune == "TRUE")


plot_diff_abundance =
  data_for_ethinicity_absolute_plot |>
  mutate(parameter = parameter |> str_remove(" or Latin American") |> str_remove("ethnicity") ) |>
  ggplot(aes(proportion_mean, fct_reorder(parameter, dplyr::desc(proportion_mean)))) +
  geom_linerange(aes(xmin = proportion_lower, xmax = proportion_upper), lwd = 0.2,) +
  geom_point(size = 0.5, color = "red") +
  xlab("Immune proportion") +
  ylab("Ethnicity") +
  theme_multipanel


plot_diff_variability =
  data_for_ethinicity_absolute_plot |>
  mutate(parameter = parameter |> str_remove(" or Latin American") |> str_remove("ethnicity") ) |>
  ggplot(aes(v_effect, fct_reorder(parameter, desc(proportion_mean)))) +
  geom_linerange(aes(xmin = v_lower, xmax = v_upper), lwd = 0.2,) +
  geom_point(size = 0.5, color = "red") +
  xlab("Variability log-fold change") +
  ylab("Ethnicity") +
  theme_multipanel

ethnicity_absolute_organ =
  differential_composition_ethnicity_absolute |>

  # Find stats of random effect with groups
  sccomp_test(
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)___.+", remove = FALSE) |>
      select(tissue_harmonised, contrast) |>
      deframe( )
  )

# CI Absolute
ethnicity_absolute_organ |>
  filter(is_immune == "TRUE") |>
  extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>

  # Select significant
  inner_join(
    differential_composition_ethnicity_absolute_one_vs_all |>
      extract(parameter, "tissue", "(.+)_.+", remove = FALSE) |>
      distinct(tissue)
  ) |>
  extract(parameter, "tissue", "(.+)_", remove = FALSE) |>
  ggplot(aes(c_effect, fct_reorder(parameter, dplyr::desc(c_effect)))) +
  geom_linerange(aes(xmin = c_lower, xmax = c_upper), lwd = 0.2,) +
  geom_point(size = 0.5, color = "red") +
  facet_wrap(~tissue, scale="free_y") +
  xlab("Immune proportion") +
  ylab("Ethnicity") +
  theme_multipanel


proportions_ethnicity_tissue_absolute_adjusted =
  readRDS("~/PostDoc/HCAquery/dev/proportions_ethnicity_tissue_absolute_adjusted.rds")

plot_ethnicity_absolute_organ_boxoplot_adjusted =
  proportions_ethnicity_tissue_absolute_adjusted |>
  left_join(
    data_for_immune_proportion |>
      distinct(.sample, tissue_harmonised, ethnicity, tissue, file_id)
  ) |>
  inner_join(
    differential_composition_ethnicity_absolute_one_vs_all |>
      tidyr::extract(parameter, "tissue_harmonised", "(.+)_.+", remove = FALSE) |>
      distinct(tissue_harmonised)
  ) |>
  mutate(ethnicity = ethnicity |> str_remove(" or Latin American") ) |>
  mutate(tissue_harmonised = tissue_harmonised |> str_to_sentence()) |>
  filter(is_immune =="TRUE") |>
  ggplot(aes(ethnicity, adjusted_proportion )) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2, fatten = 0.2) +
  geom_jitter(aes(color = file_id), width = 0.1, size=0.1) +
  facet_wrap(~ tissue_harmonised, scale="free_x", nrow = 1) +
  guides(color = "none") +
  ylab("Adjusted proportion") +
  theme_multipanel +
  theme(axis.text.x = element_text(angle=30, hjust = 1, vjust = 1))


# Color bodies
differential_composition_ethnicity_absolute_one_vs_all =
  differential_composition_ethnicity_absolute |>

  # Find stats of random effect with groups
  sccomp_test(
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue", "(.+)_.+___.+", remove = FALSE) |>
      select(-parameter) |>
      with_groups(tissue, ~ .x |> mutate(how_many = n())) |>
      filter(how_many > 1) |>
      mutate(ethnicity = ethnicity |> str_remove("ethnicity") |> str_remove(" or Latin American")) |>
      pivot_wider(names_from = ethnicity, values_from = contrast) |>
      replace_na(list(European = "0", Chinese = "0", African = "0", Hispanic = "0")) |>
      mutate(
        contrast_European = glue("{European} - 1/{how_many-1} * ({Chinese} + {African} + {Hispanic})"),
        contrast_Chinese = glue("{Chinese} - 1/{how_many-1} * ({European} + {African} + {Hispanic})"),
        contrast_African = glue("{African} - 1/{how_many-1} * ({Chinese} + {European} + {Hispanic})"),
        contrast_Hispanic = glue("{Hispanic} - 1/{how_many-1} * ({Chinese} + {African} + {European})")
      ) |>
      select(tissue, contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic) |>
      pivot_longer(
        c(contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic),
        names_to = "ethnicity",
        values_to = "contrast"
      ) |>
      mutate(ethnicity = ethnicity |> str_remove("contrast_")) |>
      unite("name", c(tissue,  ethnicity )) |>
      deframe()
  ) |>
  filter(is_immune=="TRUE") |>
  filter(c_FDR<0.01) |>
  arrange(desc(abs(c_effect))) |>
  mutate(color = circlize::colorRamp2(
    seq(2.5,-2.5, length.out = 11),
    RColorBrewer::brewer.pal(11, "RdBu")
      )(c_effect)) |>
      mutate(rgb = map_chr(
        color,
        ~ .x |>
          col2rgb() |>
          paste(collapse = " ")
      ))

differential_composition_ethnicity_absolute_one_vs_all |>
      pull(color) |>
      scales::show_col()


# # Relative
differential_composition_ethnicity_relative = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_ethnicity_relative2.rds")
#
# plot_ethinicy_bubble =
#   ethnicity_relative_organ_cell_type |>
#   tidyr::extract(parameter, c("tissue_harmonised", "ethnicity"), "(.+)_(.+)", remove = FALSE) |>
#   mutate(cell_type_harmonised = cell_type_harmonised |> str_replace("macrophage", "macro")) |>
#   circle_plot() +
#   scale_fill_viridis_c(direction = -1)
#
#
# significant_cell_types =
#   differential_composition_ethnicity_relative |>
#   sccomp_test(
#     c(
#       African = "`(Intercept)`",
#       #Pacific = "`ethnicityPacific Islander` + ((Intercept))",
#       Hispanic = "`ethnicityHispanic or Latin American` + `(Intercept)`",
#       European = "ethnicityEuropean + `(Intercept)`",
#       Chinese = "ethnicityChinese + `(Intercept)`"
#     )
#   ) |>
#
#   filter(parameter %in% c("African",  "Hispanic", "European", "Chinese")) |>
#   nest(data = -cell_type_harmonised) |>
#   mutate(sd_effect = map_dbl(data, ~ sd(.x$c_effect))) |>
#   filter(cell_type_harmonised != "immune_unclassified") |>
#   arrange(desc(sd_effect)) |>
#   head(6) |>
#   pull(cell_type_harmonised)

# Significance global statistics
count_significance_ethnicity_immune_load =
  differential_composition_ethnicity_absolute |>
  sccomp_test(
    contrasts =
      c(
        European = "ethnicityEuropean - 1/3 * (ethnicityAfrican + ethnicityChinese + `ethnicityHispanic or Latin American`)",
        African = "ethnicityAfrican - 1/3 * (ethnicityEuropean + ethnicityChinese + `ethnicityHispanic or Latin American`)",
        Chinese = "ethnicityChinese - 1/3 * (ethnicityAfrican + ethnicityEuropean + `ethnicityHispanic or Latin American`)",
        Hispanic = "`ethnicityHispanic or Latin American` - 1/3 * (ethnicityAfrican + ethnicityChinese + ethnicityEuropean)"
      ),
    test_composition_above_logit_fold_change = 0.1
  ) |>
  filter(is_immune=="TRUE") |>
  count(c_FDR<0.05)

count_significance_ethnicity_immune_load_tissue =
  differential_composition_ethnicity_absolute |>
  sccomp_test(
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue", "(.+)_.+___.+", remove = FALSE) |>
      select(-parameter) |>
      with_groups(tissue, ~ .x |> mutate(how_many = n())) |>
      filter(how_many > 1) |>
      mutate(ethnicity = ethnicity |> str_remove("ethnicity") |> str_remove(" or Latin American")) |>
      pivot_wider(names_from = ethnicity, values_from = contrast) |>
      replace_na(list(European = "0", Chinese = "0", African = "0", Hispanic = "0")) |>
      mutate(
        contrast_European = glue("{European} - 1/{how_many-1} * ({Chinese} + {African} + {Hispanic})"),
        contrast_Chinese = glue("{Chinese} - 1/{how_many-1} * ({European} + {African} + {Hispanic})"),
        contrast_African = glue("{African} - 1/{how_many-1} * ({Chinese} + {European} + {Hispanic})"),
        contrast_Hispanic = glue("{Hispanic} - 1/{how_many-1} * ({Chinese} + {African} + {European})")
      ) |>
      select(tissue, contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic) |>
      pivot_longer(
        c(contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic),
        names_to = "ethnicity",
        values_to = "contrast"
      ) |>
      mutate(ethnicity = ethnicity |> str_remove("contrast_")) |>
      unite("name", c(tissue,  ethnicity )) |>
      deframe()
  ) |>
  filter(is_immune == "TRUE") |>
  count(c_FDR<0.05)

count_significance_ethnicity_cell_type =
  differential_composition_ethnicity_relative |>
  sccomp_test(
    contrasts =
      c(
        European = "ethnicityEuropean - 1/3 * (ethnicityAfrican + ethnicityChinese + `ethnicityHispanic or Latin American`)",
        African = "ethnicityAfrican - 1/3 * (ethnicityEuropean + ethnicityChinese + `ethnicityHispanic or Latin American`)",
        Chinese = "ethnicityChinese - 1/3 * (ethnicityAfrican + ethnicityEuropean + `ethnicityHispanic or Latin American`)",
        Hispanic = "`ethnicityHispanic or Latin American` - 1/3 * (ethnicityAfrican + ethnicityChinese + ethnicityEuropean)"
      ),
    test_composition_above_logit_fold_change = 0.4
  ) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)

count_significance_ethnicity_cell_type_tissue =
  differential_composition_ethnicity_relative |>
  sccomp_test(
    contrasts =
      differential_composition_ethnicity_absolute |>
      filter(parameter |> str_detect("___ethnicity")) |>
      distinct(parameter) |>
      tidyr::extract( parameter, "ethnicity", "_(.+)___", remove = FALSE) |>
      mutate(ethnicity = glue("ethnicity{ethnicity}")) |>
      mutate(contrast = glue("`{ethnicity}`  + `{parameter}`") |> as.character()) |>
      tidyr::extract(parameter, "tissue", "(.+)_.+___.+", remove = FALSE) |>
      select(-parameter) |>
      with_groups(tissue, ~ .x |> mutate(how_many = n())) |>
      filter(how_many > 1) |>
      mutate(ethnicity = ethnicity |> str_remove("ethnicity") |> str_remove(" or Latin American")) |>
      pivot_wider(names_from = ethnicity, values_from = contrast) |>
      replace_na(list(European = "0", Chinese = "0", African = "0", Hispanic = "0")) |>
      mutate(
        contrast_European = glue("{European} - 1/{how_many-1} * ({Chinese} + {African} + {Hispanic})"),
        contrast_Chinese = glue("{Chinese} - 1/{how_many-1} * ({European} + {African} + {Hispanic})"),
        contrast_African = glue("{African} - 1/{how_many-1} * ({Chinese} + {European} + {Hispanic})"),
        contrast_Hispanic = glue("{Hispanic} - 1/{how_many-1} * ({Chinese} + {African} + {European})")
      ) |>
      select(tissue, contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic) |>
      pivot_longer(
        c(contrast_European, contrast_Chinese, contrast_African, contrast_Hispanic),
        names_to = "ethnicity",
        values_to = "contrast"
      ) |>
      mutate(ethnicity = ethnicity |> str_remove("contrast_")) |>
      unite("name", c(tissue,  ethnicity )) |>
      deframe()
  ) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  count(c_FDR<0.05)




rm(differential_composition_ethnicity_relative , differential_composition_ethnicity_absolute )
gc()

plot_significance_overall =
  count_significance_age_immune_load |>
  mutate(name = c(
    "count_significance_age_immune_load"
  )) |>
  bind_rows(
    count_significance_age_immune_load_tissue |>
      mutate(name = c(
        "count_significance_age_immune_load_tissue"
      ))
    ) |>
  bind_rows(count_significance_age_cell_type |>
              mutate(name = c(

                "count_significance_age_cell_type"
              ))) |>
  bind_rows(count_significance_age_cell_type_tissue |>
              mutate(name = c(

                "count_significance_age_cell_type_tissue"
              ))) |>

  bind_rows(count_significance_sex_immune_load |>
              mutate(name = c(

                "count_significance_sex_immune_load"
              ))) |>
  bind_rows(count_significance_sex_immune_load_tissue |>
              mutate(name = c(

                "count_significance_sex_immune_load_tissue"
              ))) |>
  bind_rows(count_significance_sex_cell_type |>
              mutate(name = c(

                "count_significance_sex_cell_type"
              ))) |>
  bind_rows(count_significance_sex_cell_type_tissue |>
              mutate(name = c(

                "count_significance_sex_cell_type_tissue"
              ))) |>

  bind_rows(count_significance_ethnicity_immune_load |>
              mutate(name = c(

                "count_significance_ethnicity_immune_load"
              ))) |>
  bind_rows(count_significance_ethnicity_immune_load_tissue |>
              mutate(name = c(

                "count_significance_ethnicity_immune_load_tissue"
              ))) |>
  bind_rows(count_significance_ethnicity_cell_type |>
              mutate(name = c(

                "count_significance_ethnicity_cell_type"
              ))) |>
  bind_rows(count_significance_ethnicity_cell_type_tissue |>
              mutate(name = c(

                "count_significance_ethnicity_cell_type_tissue"
              ))) |>
    tidyr::extract(name, c("factor", "variable", "resolution"), "count_significance_([a-zA-Z]+)_([a-zA-Z]+_[a-zA-Z]+)_?(.*)", remove = FALSE) |>
    mutate(resolution = if_else(resolution == "", "overall", resolution)) |>
    unite("xlab", c(variable, resolution), remove = FALSE) |>
    mutate(xlab = xlab |> fct_relevel(c("immune_load_overall", "immune_load_tissue", "cell_type_overall", "cell_type_tissue"))) |>
    with_groups(name, ~ .x |> mutate(sum_n = sum(n))) |>
    mutate(proportion = n/sum_n) |>
    mutate(factor = factor |> str_to_sentence()) |>
    ggplot(aes(xlab, proportion, fill=`c_FDR < 0.05`)) +
    geom_bar(stat = "identity")+
    geom_text(aes(y = 0.5, label = sum_n), size = 2.5, angle=90) +
    facet_wrap( ~ factor,  nrow=1) +
    scale_fill_manual(values = c("FALSE"="grey", "TRUE"="#D5C711")) +
    ylab("Proportion of significant tests") +
    xlab("Hypotheses") +
    theme_multipanel +
    theme(axis.text.x = element_text(angle=20, hjust = 1, vjust = 1))



# job::job({ readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_4.rds") |> sccomp_remove_unwanted_variation(~ age_days) })


plot_first_line =
  (
    (
      (
        (
          ( plot_significance_overall | plot_age_absolute ) + plot_layout(widths = c(46,96))
        ) /
          plot_spacer()
      ) + plot_layout(heights = c(44,15))
    ) |
      plot_spacer()
  ) +
  plot_layout(width = c(144,36), guides = 'collect')

plot_second_line =
  (
    plot_age_relative  |
    plot_spacer() |
    wrap_heatmap(plot_heatmap_age_relative_organ_cell_type, padding = unit(c(-40, 0, -10, -30), "points" ))
    ) +
  plot_layout( widths = c(66,  34, 78), guides = 'collect')

plot_third_line =
    (
      (
        (plot_sex_absolute_overall / plot_sex_absolute_tissue / plot_diff_abundance / plot_diff_variability ) +
          plot_layout(heights  = c(10, 42, 13, 13))
      ) |
        (
          ( plot_spacer() / plot_ethnicity_absolute_organ_boxoplot_adjusted ) +
            plot_layout(heights  = c(38, 42))
        )
    ) +
      plot_layout(widths = c(80, 100), guides = 'collect')



# plot_sex_ethnicity |> saveRDS("~/PostDoc/sccomp_dev/dev/plot_sex_ethnicity.rds")
plot_sex_ethnicity = readRDS("~/PostDoc/sccomp_dev/dev/plot_sex_ethnicity.rds")

# Plotting
p =
  (
    plot_first_line /
      plot_second_line /
      plot_third_line
  ) +
  plot_layout(guides = 'collect', heights = c(57, 57, 92))  &
  theme(
    plot.margin = margin(0, 0, 0, 0, "pt"),
    legend.key.size = unit(0.2, 'cm'),
    legend.position = "bottom"
  )


ggsave(
  "~/PostDoc/sccomp_dev/dev/article_figures/HCA_demography.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 200 ,
  limitsize = FALSE
)

