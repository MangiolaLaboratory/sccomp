library(patchwork)
library(bayestestR)
library(yardstick)
library(forcats)
library(glue)
library(purrr)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/7a6df4a80d38b0427a263c21ac2480c4280cfe4b/ggplot_theme_multipanel")

set.seed(42)

result_directory = "dev/benchmark_results_0568ae31db577b7a39976a96327e2ac741247c77"

custom_paired_palette = c("#E31A1C" ,"#A6CEE3" ,"#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99","#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")

cool_palette = c("#b58b4c", "#74a6aa", "#a15259",  "#37666a", "#79477c", "#cb9f93", "#9bd18e", "#eece97", "#8f7b63", "#4c474b", "#415346")
color_palette_link = cool_palette %>% setNames(c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))
scales::show_col(cool_palette)


plot_auc =
  dir(result_directory, pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = c(name, auc)) %>%
  mutate(random_auc = map_dbl(data, ~ .x %>% filter(name=="random") %>% pull(auc))) %>%
  unnest(data) %>%

  filter(name !="speckle") %>%
  filter(n_samples >2) %>%
  mutate(name = if_else(name=="random", "zrandom", name)) %>%
  # Choose one mode
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  mutate(plot = map(
    data,
    ~ .x %>%

      mutate(`Accuracy (AUC)` = auc ) %>%
      mutate(n_cell_type_ = glue("{n_cell_type} cell types"), n_samples_ = glue("{n_samples} samples")) %>%
      rename(Algorithm = name) %>%

      ggplot(aes(slope, `Accuracy (AUC)`, color=Algorithm)) +
      geom_line(size=0.3) +
      facet_grid( fct_reorder(n_cell_type_,n_cell_type)  ~ fct_reorder(n_samples_, n_samples), scale="free_y") +
      geom_hline(yintercept = 0.1, linetype="dashed", color="grey") +
      scale_y_continuous( labels = dropLeadingZero) +
      scale_color_manual(values = cool_palette) +
      theme_bw() +
      ylab("Accuracy (AUC) up to 0.1 false-positive-rate") +
      xlab("Slope") +
      multipanel_theme
  )) %>%
  pull(plot) %>%
  .[[2]]
  # patchwork::wrap_plots(ncol = 1)

plot_ROC =
  readRDS(glue("{result_directory}/parsed_1.5_10_10_1000_1.rds")) %>%
  mutate(roc = map(
    df_for_roc,
    ~ .x %>%
      roc_curve(event_level = 'second', truth = significant, probability)
  )) %>%
  select(-df_for_roc) %>%
  unnest(roc) %>%

  filter( name !="speckle") %>%
  rename(Algorithm = name) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color=Algorithm)) +
  geom_rect( mapping=aes(xmin=0, xmax=0.1, ymin=0, ymax=1), fill="#EBECF0", color=NA, alpha=0.1) +
  geom_path(size=0.3) +
  geom_abline(lty = 3) +
  scale_color_manual(values = cool_palette) +
  scale_y_continuous( labels = dropLeadingZero) +
  guides(color="none") +
  multipanel_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("True positive rate") +
  xlab("False positive rate")

simulated_data =
  readRDS(glue("{result_directory}/input_1.5_15_10_1000_1.rds")) %>%
  pull(data) %>%
  .[[1]] %>%
  mutate(Significant = beta_1!=0) %>%
  # mutate(Outlier = ratio!=1)
  mutate(Outlier = FALSE)

plot_simulated =
  ggplot() +
  geom_boxplot(
    aes(factor(type), generated_proportions, fill=Significant),
    outlier.shape = NA,
    data =  filter(simulated_data, !Outlier),
    lwd = 0.3, fatten = 0.5
  ) +
  geom_jitter(aes(factor(type), generated_proportions, color=Outlier), size = 0.5, data = simulated_data) +
  facet_wrap(~ interaction(cell_type), scale="free_y", ncol = 5) +
  scale_y_continuous(trans="logit", labels = dropLeadingZero) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  guides(color="none") +
  xlab("Biological condition") +
  ylab("Cell-group proportion (decimal)") +
  multipanel_theme +
  theme(strip.background =element_rect(fill="white", color="white"), legend.position = "bottom")

p =
  (
( plot_simulated + plot_ROC ) /
plot_auc
)+

  # Style
  plot_layout(guides = 'collect', heights  = c(1.5,2)) + plot_annotation(tag_levels = c('A')) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

ggsave(
  "dev/article_figures/plot_benchmark.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 150 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/plot_benchmark.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 150 ,
  limitsize = FALSE
)
