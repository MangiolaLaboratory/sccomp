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


source("https://gist.githubusercontent.com/stemangiola/cfa08c45c28fdf223d4996a6c1256a39/raw/7c78b50dce501fc7ce0b2a8d8efd3aded91134aa/color_cell_types.R")
cell_type_color = color_array
names(cell_type_color) = names(cell_type_color) |>  str_replace("macrophage", "macro")



tissue_color =
  data_for_immune_proportion_relative |>
  distinct(tissue_harmonised ) |>
  arrange(tissue_harmonised) |>
  mutate(color = dittoSeq::dittoColors()[1:n()]) |>
  deframe()


# LOADING RESULTS
result_directory = "/home/users/allstaff/mangiola.s/PostDoc/HCAquery/dev/run_5_dec_2022/"

res_absolute = readRDS(glue("{result_directory}/immune_non_immune_differential_composition.rds"))

res_absolute = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition.rds")

lymphoid_organs = c("blood", "spleen", "bone", "thymus", "lymph node")

res_generated_proportions =
  res_absolute |>
  replicate_data(number_of_draws = 20) |>
  filter(is_immune=="TRUE") |>
  left_join(
    data_for_immune_proportion |>
      select(.sample, tissue_harmonised)
  ) |>
  with_groups(tissue_harmonised, ~ .x |> sample_n(30, replace = T))

# -- Rank immune
dropLeadingZero <- function(l){  stringr::str_replace(l, '0(?=.)', '') }
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() scales::trans_new("S_sqrt",S_sqrt,IS_sqrt)

library(ggforestplot)

plot_immune_proportion_dataset =
  data_for_immune_proportion |>

  # Stats
  dplyr::count(.sample, tissue_harmonised, is_immune, file_id) |>
  with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n), sum = sum(n))) |>
  filter(is_immune) |>
  with_groups(tissue_harmonised, ~ .x |> mutate( median_proportion = mean(proportion))) |>

  # Add multilevel proportion medians
  left_join(
    res_generated_proportions |>
      with_groups(tissue_harmonised, ~ .x |> summarise(median_generated = median(generated_proportions, na.rm = TRUE)))
  ) |>

  # Label lymphoid organs
  mutate(is_lymphoid = tissue_harmonised %in% c("blood", "spleen", "bone", "thymus", "lymph node")) |>

  clean_names() |>

  # Arrange
  arrange(is_lymphoid, median_generated) %>%
  mutate(tissue_harmonised = factor(tissue_harmonised, levels = unique(.$tissue_harmonised))) |>

  # Cap
  mutate(sum = sum |> pmax(500)) |>
  mutate(sum = sum |> pmin(100000)) |>

  # Plot
  ggplot(aes( proportion, tissue_harmonised)) +
  ggforestplot::geom_stripes(odd = "#33333333", even = "#00000000") +
  geom_jitter(aes(size = sum, color=file_id), width = 0) +
  geom_boxplot(aes(generated_proportions, fct_reorder(tissue_harmonised, median_generated)), color="black", data =
                 res_generated_proportions |>
                 with_groups(tissue_harmonised, ~ .x |> mutate(median_generated = median(generated_proportions, na.rm = TRUE))) |>
                 clean_names(),
               fill = NA, outlier.shape = NA, lwd=0.2
  ) +
  guides(color="none") +
  scale_size(trans = "sqrt", range = c(0.1, 1.5), limits = c(500, 100000)) +
  scale_color_manual(values = dittoSeq::dittoColors()) +
  scale_x_continuous(trans=S_sqrt_trans(), labels = dropLeadingZero) +
  xlab("Immune proportion (sqrt)") +
  ylab("Tissue") +
  theme_multipanel

coefficients_regression =
  res_absolute |>
  filter(covariate == "tissue_harmonised") |>
  filter(c_effect<1) |>
  filter(is_immune == "TRUE") |>
  filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) %>%
  lm( v_effect ~ c_effect, data = .) %$%
  coefficients

library(ggrepel)

median_composition =
  res_absolute |>
  filter(covariate == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) |>
  pull(c_effect) |>
  median()

median_variability =
  res_absolute |>
  filter(covariate == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  filter(!parameter %in% c("spleen", "bone", "blood", "lymph node") ) |>
  pull(v_effect) |>
  median()

# - scatter plot of abundance vs variability per tissue
res_for_plot =
  res_absolute |>
  filter(covariate == "tissue_harmonised") |>
  filter(is_immune == "TRUE") |>
  mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
  mutate(intercept = coefficients_regression[1], slope = coefficients_regression[2]) |>

  # # Normalise effects
  # mutate(
  # 	v_effect = v_effect - (c_effect * slope + intercept),
  # 	v_lower = v_lower - (c_effect * slope + intercept),
  # 	v_upper = v_upper - (c_effect * slope + intercept)
  # ) |>
  # mutate(
  # 	c_effect = c_effect - median_composition,
  # 	c_lower = c_lower - median_composition,
  # 	c_upper = c_upper - median_composition
# ) |>

# Define significance
mutate(tissue_harmonised = parameter |> str_remove("tissue_harmonised")) |>
  arrange(desc(c_effect)) |>
  mutate(	c_significant = ( row_number() <=7 | row_number() >= n() -3 ) & !tissue_harmonised %in%	c("blood", "lymph node", "spleen", "bone")	) |>
  arrange(desc(v_effect)) |>
  mutate(	v_significant = ( row_number() <=3 | row_number() >= n() -3) & !tissue_harmonised %in%	c("blood", "lymph node", "spleen", "bone")) |>

  # Define quadrants
  mutate(quadrant = case_when(
    c_effect > median_composition & v_effect > median_variability & (c_significant | v_significant) ~ "Hot and variable",
    c_effect > median_composition & v_effect < median_variability & (c_significant | v_significant) ~ "Hot and consistent",
    c_effect < median_composition & v_effect > median_variability & (c_significant | v_significant) ~ "Cold and variable",
    c_effect < median_composition & v_effect < median_variability & (c_significant | v_significant) ~ "Cold and consistent"

  )) |>

  # Limit the values
  mutate(
    v_effect = pmax(v_effect, -5),
    c_effect = pmin(c_effect, 1.5),
    v_lower = pmax(v_lower, -5),
    v_upper = pmax(v_upper, -5),
    c_lower = pmin(c_lower, 1.5),
    c_upper = pmin(c_upper, 1.5),
  )

# Plot abundance variability
plot_abundance_variability =
  res_for_plot |>
  ggplot(aes(c_effect, v_effect, label = parameter)) +
  geom_vline(xintercept = median_composition, linetype = "dashed",  alpha = 0.5, lwd = 0.2) +
  geom_hline(yintercept = median_variability, linetype = "dashed",  alpha = 0.5, lwd = 0.2) +
  geom_errorbar(aes(ymin = v_lower, ymax = v_upper, color = v_significant), alpha = 0.5, lwd = 0.2) +
  geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_significant), alpha = 0.5, lwd = 0.2) +
  geom_point(aes(fill = quadrant), shape = 21, stroke = 0.2) +
  #geom_text_repel(data =  res_for_plot |> filter(!v_significant & !c_significant)) +
  geom_text_repel(
    data = res_for_plot |>
      filter(v_significant | c_significant),
    size = 2
  ) +
  xlab("Composition") +
  ylab("Variability") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "black")) +
  scale_fill_brewer(palette = "Set1") +
  theme_multipanel

# cell_metadata_with_harmonised_annotation =
#   readRDS("~/PostDoc/HCAquery/dev/cell_metadata_with_harmonised_annotation.rds")
#
# # Stats
# cell_metadata_with_harmonised_annotation |>
#
#   mutate(is_immune = cell_type_harmonised!="non_immune") |>
#
#   # Stats
#   dplyr::count(.sample, tissue_harmonised, is_immune) |>
#   with_groups(.sample, ~ .x |> mutate(proportion = n/sum(n))) |>
#   filter(tissue_harmonised=="heart" & proportion > 0.75)



proportions_tissue_adjusted = readRDS("~/PostDoc/HCAquery/dev/proportions_tissue_adjusted_5.rds")

# proportions_tissue_replicate = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_5.rds") |> replicate_data (~ 0 + tissue_harmonised)
# proportions_tissue_replicate |> saveRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_5_replicate.rds")

proportions_tissue_replicate = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_relative_5_replicate.rds")

res_absolute_cell_types = readRDS("dev/immune_non_immune_differential_composition_relative_5.rds")

proportions_tissue_adjusted = readRDS("~/PostDoc/HCAquery/dev/proportions_tissue_adjusted_5.rds")
immune_non_immune_differential_composition_cell_types_adjusted = readRDS("~/PostDoc/HCAquery/dev/immune_non_immune_differential_composition_cell_types_adjusted.rds")


library(tidybulk)
observed_proportion_PCA_df =
  data_for_immune_proportion |>

  # Mutate days
  filter(development_stage!="unknown") |>

	add_count(tissue_harmonised) |>
	filter(n > 5) |>
	select(-n) |>
	dplyr::count(.sample, cell_type_harmonised, tissue_harmonised, assay, sex, file_id) |>
	with_groups(.sample, ~ .x |> mutate(observed_proportion = n/sum(n))) |>
	tidyr::complete(nesting(.sample, tissue_harmonised, assay, sex, file_id), cell_type_harmonised, fill = list(observed_proportion = 0)) |>
	reduce_dimensions(.sample, cell_type_harmonised, observed_proportion, method="tSNE", action="get")

observed_proportion_PCA_tissue =
  observed_proportion_PCA_df |>
	ggplot(aes(tSNE1, tSNE2)) +
	geom_point(aes(fill = tissue_harmonised), shape=21, stroke = NA, size=0.2) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none") +
  ggtitle("Tissue") +
	theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

observed_proportion_PCA_batch =
  observed_proportion_PCA_df |>
  ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = file_id), shape=21, stroke = NA, size=0.2) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none")  +
  ggtitle("Dataset") +
  theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

adjusted_proportion_PCA =
  immune_non_immune_differential_composition_cell_types_adjusted |>
	left_join(
		data_for_immune_proportion_relative |>
			distinct(.sample, tissue_harmonised, assay, file_id, sex, ethnicity)
	) |>
	add_count(tissue_harmonised) |>
	filter(n > 5) |>
	reduce_dimensions(.sample , cell_type_harmonised, adjusted_proportion, method="tSNE", action="get") |>
	ggplot(aes(tSNE1, tSNE2)) +
  geom_point(aes(fill = tissue_harmonised), shape=21, stroke = NA, size=0.2) +
  #ggdensity::geom_hdr_lines(aes(color = tissue_harmonised)) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  ggtitle("Adjusted tissue") +
	theme_multipanel +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

res_absolute_for_PCA =
  res_absolute_cell_types |>
  filter(covariate == "tissue_harmonised") |>
  mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
  filter(cell_type_harmonised != "immune_unclassified") |>
  rename(feature=cell_type_harmonised) |>
  bind_rows(
    res_absolute |>
      filter(covariate=="tissue_harmonised") |>
      mutate(parameter = parameter |> str_remove("tissue_harmonised")) |>
      filter(is_immune=="TRUE") |>
      mutate(is_immune = "immune") |>
      rename(feature=is_immune)
  ) |>
  filter(parameter != "skeletal_muscle") |>
  reduce_dimensions(
    parameter, feature, c_effect,
    method="PCA", action="get", scale=FALSE,
    transform = identity
  )

plot_tissue_PCA =
  res_absolute_for_PCA |>
  ggplot(aes(PC1, PC2, label = parameter)) +
  geom_point(aes(fill = parameter), shape=21, stroke = 0.2) +
  ggrepel::geom_text_repel(size=2, max.overlaps = 10, min.segment.length = unit(10, "pt")) +
  scale_fill_manual(values = dittoSeq::dittoColors()) +
  guides(fill="none")  +
  theme_multipanel

# res_absolute_for_PCA |> attr("internals") %$% PCA |>
#   ggplot2::autoplot( data = res_absolute_for_PCA, colour = 'parameter',
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 3)

library(ggforce)
library(ggpubr)

circle_plot = function(res){

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
		with_groups(cell_type_harmonised, ~ .x |>  mutate(cell_type_mean_abundance = mean(proportion))) |>

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
	  xlab("Diversity") +
	  ylab("Tissue") +
		theme_multipanel

	cell_type_mean_abundance_plot =
		res_relative_for_plot |>
		distinct(cell_type_mean_abundance, cell_type_harmonised) |>
		ggplot(aes(fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_abundance)), cell_type_mean_abundance)) +
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
			fct_reorder(cell_type_harmonised, dplyr::desc(cell_type_mean_abundance)),
			fct_reorder(tissue_harmonised, inter_type_diversity) ,
			fill = rank, size = c_effect, stroke=rank==1), shape=21
		) +
	  scale_size_continuous(range = c(0.5, 3)) +
	  xlab("Cell type") +
		theme_multipanel +
		theme(
			axis.text.x = element_text(angle=30, hjust = 1, vjust = 1),
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


df_heatmap_relative_organ_cell_type =
  res_relative |>
  filter(covariate == "tissue_harmonised") |>
  mutate(tissue_harmonised =   parameter ) |>
  clean_names() |>
  mutate(tissue =   tissue_harmonised ) |>

  mutate(cell_type = cell_type_harmonised |> str_replace("macrophage", "macro")) |>







  filter(cell_type !="immune_unclassified") |>

  # Color
  left_join(tissue_color |> enframe(name = "tissue", value = "tissue_color")  ) |>
  left_join(cell_type_color |> enframe(name = "cell_type", value = "cell_type_color")  ) |>

  # Counts
  left_join(
    data_for_immune_proportion_relative |>
      count(tissue_harmonised, name = "count_tissue") |>
      rename(tissue = tissue_harmonised) |>
      mutate(count_tissue = log(count_tissue))
  ) |>

  with_groups(cell_type, ~ .x |> arrange(desc(c_effect)) |> mutate(top_3 = row_number() |> between(1, 3))) |>

  # Cell type abundance
  with_groups(tissue, ~ .x |> mutate(proportion = softmax(c_effect))) |>
  with_groups(cell_type, ~ .x |>  mutate(`Mean diff` = mean(proportion, na.rm = TRUE))) |>
  mutate(cell_type = fct_reorder(cell_type, -`Mean diff`))


plot_circle_relative_tissue =
  df_heatmap_relative_organ_cell_type |>

  # Heatmap
  heatmap(
    tissue, cell_type, c_effect,
    # palette_value = circlize::colorRamp2(
    #   seq(-3, 3, length.out = 11),
    #   RColorBrewer::brewer.pal(11, "Spectral")
    # ),
    palette_value = circlize::colorRamp2(c(-5, -2.5, 0, 2.5, 5), viridis::viridis(5)),
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 6),
    column_title_gp = gpar(fontsize = 0),
    row_title_gp = gpar(fontsize = 0),
    show_heatmap_legend = FALSE,
    row_km = 2
  ) |>
  split_rows(6) |>
  annotation_bar(`Mean diff`, annotation_name_gp= gpar(fontsize = 8), size = unit(0.4, "cm")) |>

  annotation_tile(
    tissue, show_legend = FALSE,
    palette =
      df_heatmap_relative_organ_cell_type |>
      distinct(tissue, tissue_color) |>
      arrange(tissue) |>
      deframe(),
    size = unit(0.2, "cm")
  ) |>
  annotation_tile(
    cell_type, show_legend = FALSE,
    palette =
      df_heatmap_relative_organ_cell_type |>
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
  layer_point(top_3)


# 	circle_plot() +
#   scale_fill_viridis_c(direction = -1)


#job::job({ plot_summary_res_relative = res_relative |> plot_summary() })



# Compose plot

second_line_first_column =

  plot_immune_proportion_dataset /
  (
    observed_proportion_PCA_tissue |
      observed_proportion_PCA_batch |
      adjusted_proportion_PCA
  ) /
  plot_tissue_PCA +
  plot_layout( guides = 'collect', height = c(4, 1, 2)  )   &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")


second_line_second_column =
  (
    plot_abundance_variability /
      wrap_heatmap(plot_circle_relative_tissue, padding = unit(c(-67, -10, -0, -30), "points" ))
  ) +
  plot_layout( guides = 'collect', height = c(1,2) )   &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")


p =
  (
    second_line_first_column |
      second_line_second_column
  ) +
  plot_layout(ncol=2,  width = c(1,1.5) )  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")




ggsave(
  "~/PostDoc/sccomp_dev/dev/article_figures/HCA_tissue.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 150 ,
  limitsize = FALSE
)





