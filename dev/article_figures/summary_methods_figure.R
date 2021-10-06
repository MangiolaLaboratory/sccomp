library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(stringr)

library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)
library(forcats)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/305ed9ba2af815fdef3214b9e6171008d5917400/ggplot_theme_multipanel")
theme_UMAP =
  list(
    multipanel_theme,
    theme(
      axis.text.x  = element_blank(),
      axis.ticks.x =  element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y =  element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  )
friendly_cols <- dittoSeq::dittoColors()

set.seed(42)

job({
  counts  =
    dir("dev/data_integration/", pattern = "^UMAP_", full.names = TRUE) %>%
    grep("UMAP_oligo", ., invert = TRUE, value = TRUE) %>%
    map(~ readRDS(.x) %>%
          as_tibble() %>%
          select(cell, sample, cell_type, UMAP_1, UMAP_2) %>%
          mutate(file = .x)
    ) %>%
    reduce(bind_rows)
})

oligo_adjusted_cell_type = readRDS("dev/data_integration/UMAP_oligo.rds")

umap_1 =
  oligo_adjusted_cell_type %>%
  sample_n(10000) %>%
  ggplot(aes(UMAP_1, UMAP_2, fill=cell_type)) +
  geom_point(aes(fill=cell_type, color=cell_type),size=0.2, alpha=1, shape=21, stroke=0) +
  scale_fill_manual(values = friendly_cols) +
  scale_color_manual(values = friendly_cols) +
  guides(color="none", fill="none") +
  multipanel_theme +
  theme_UMAP

# composition_1 =
#   oligo_adjusted_cell_type %>%
#   count(cell_type, sample, file) %>%
#   with_groups(sample, ~ mutate(.x, proportion = n/sum(n))) %>%
#   ggplot(aes(x = sample, y = proportion, fill = cell_type)) +
#   geom_bar(stat = "identity") +
#   scale_fill_manual(values = friendly_cols) +
#   guides(fill="none") +
#   multipanel_theme  +
#   theme_UMAP

composition_1 =
  sccomp::oligo_breast_estimate %>%
  arrange(desc(abs(.median_typecancer))) %>%
  mutate(rank = 1:n()) %>%
  unnest(count_data) %>%
  with_groups(sample, ~ mutate(.x, proportions = (count)/sum(count))) %>%
  filter(cell_group %in% c("M4", "M5", "M8", "B1")) %>%
  left_join(
    sccomp::oligo_breast_estimate %>%
      replicate_data() %>%
      with_groups(sample, ~ mutate(.x, generated_proportions = (generated_counts )/sum(generated_counts ))),
    by = c("cell_group", "sample")
  ) %>%
  pivot_longer(c(proportions,generated_proportions), names_to = "simulation", values_to = "proportion") %>%
  ggplot(aes(x = type, y = proportion)) +
  geom_boxplot(aes( fill = cell_group), outlier.shape = NA, lwd =0.3, fatten=0.3) +
  geom_jitter(size=0.5, shape=21, stroke = 0) +
  facet_wrap(fct_relevel(simulation, c("proportions", "generated_proportions"))~cell_group, nrow=1) +
  scale_fill_manual(values = friendly_cols) +
  scale_y_continuous(trans="logit") +
  guides(fill="none") +
  theme_UMAP

cell_counts =
  sccomp::counts_obj %>%
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count))) %>%
  ggplot(aes(x = sample, y = proportion, fill = cell_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = friendly_cols) +
  guides(fill="none") +
  multipanel_theme +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

umap_2 =
  counts  %>%
  as_tibble() %>%
  with_groups(file, ~ sample_n(.x, min(5000, nrow(.x))) ) %>%
  ggplot(aes(UMAP_1, UMAP_2, fill=cell_type)) +
  geom_point(aes(fill=cell_type, color=cell_type),size=0.2, alpha=1, shape=21, stroke=0) +
  facet_wrap(~file, scale="free") +
  scale_fill_manual(values = friendly_cols) +
  scale_color_manual(values = friendly_cols) +
  guides(color="none", fill="none") +
  theme_UMAP

observed_proportions =
  sccomp::oligo_breast_estimate %>%
  arrange(desc(abs(.median_typecancer))) %>%
  mutate(rank = 1:n()) %>%
  unnest(count_data) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count))) %>%
  filter(cell_group %in% c("M4", "M5", "M8", "B1", "B2", "B3", "BM", "CD4 1", "CD4 2")) %>%
  mutate(cell_group = factor(cell_group, levels = c("B1","M4", "M5", "M8",  "B2", "B3", "BM", "CD4 1", "CD4 2")))


generated_proportions =
  sccomp::oligo_breast_estimate %>%
  replicate_data(number_of_draws = 100) %>%
  with_groups(c(sample, replicate), ~ mutate(.x, proportion = (generated_counts )/sum(generated_counts ))) %>%
  left_join(observed_proportions %>% distinct(sample, type))  %>%
  filter(cell_group %in% c("B1","M4", "M5", "M8",  "B2", "B3", "BM", "CD4 1", "CD4 2")) %>%
  mutate(cell_group = factor(cell_group, levels = c("B1","M4", "M5", "M8", "B2", "B3", "BM", "CD4 1", "CD4 2")))


composition_2 =

  ggplot() +
  geom_jitter(aes(x = type, y = proportion), data =generated_proportions , fill="grey" , outlier.shape = NA , alpha=0.1, size=0.1, shape=21, stroke = 0) +
  geom_boxplot(aes(x = type, y = proportion, fill=cell_group), outlier.shape = NA, data = observed_proportions, alpha=0.5, lwd =0.3, fatten=0.3) +
  geom_jitter(aes(x = type, y = proportion), data = observed_proportions,size=0.5, shape=21, stroke = 0) +
  #geom_jitter(aes(x = type, y = proportion), size=0.2,data = generated_proportions, fill=NA , color="#cc6666" ) +


  facet_wrap(~cell_group, nrow=1) +
  scale_fill_manual(values = friendly_cols) +
  scale_y_continuous(trans="logit") +
  guides(fill="none") +
  ylab("Logit proportions") +
  multipanel_theme  +
  theme_UMAP


p =
  (
    (umap_2   + plot_spacer() +  plot_layout(widths = c(5,1))) /
      ( umap_1 + cell_counts + composition_1 +  plot_layout(widths = c(2, 2, 2))) /
      (composition_2   + plot_spacer() +  plot_layout(widths = c(5,1)))
  ) +

  # Style
  plot_layout(guides = 'collect')  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")


ggsave(
  "dev/article_figures/summary_plot.pdf",
  plot = p,
  units = c("mm"),
  width = 183/2 ,
  height = 130 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/summary_plot.png",
  plot = p,
  units = c("mm"),
  width = 183/2 ,
  height = 130 ,
  limitsize = FALSE
)
