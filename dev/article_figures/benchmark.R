library(patchwork)
library(bayestestR)
library(yardstick)
library(forcats)
library(glue)
library(purrr)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/7a6df4a80d38b0427a263c21ac2480c4280cfe4b/ggplot_theme_multipanel")

set.seed(42)

#result_directory = "dev/benchmark_results_0568ae31db577b7a39976a96327e2ac741247c77"
result_directory = "dev/benchmark_results_3573b63473a94fcc9223277b17e6d59bad3890b7_simulation_logitLinear/"

cool_palette = c("#cc374a","#e29a3b", "#3a9997", "#a15259",   "#94379b", "#ff8d73", "#7be05f", "#ffc571", "#724c24", "#4c474b", "#236633")
scales::show_col(cool_palette)

levels_algorithm =
  c("sccomp", "speckle" , "logitLinear" , "ttest"  , "quasiBinomial" , "rlm" , "DirichletMultinomial" , "random"   )

plot_auc =
  dir(result_directory, pattern = "auc", full.names = TRUE) %>%
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = c(name, auc)) %>%
  mutate(random_auc = map_dbl(data, ~ .x %>% filter(name=="random") %>% pull(auc))) %>%
  unnest(data) %>%

  filter( name !="random") %>%
  filter(n_samples >2) %>%
  mutate(name = if_else(name=="random", "zrandom", name)) %>%

  # Filter too small slope
  #filter(slope>0.1) %>%

  # Choose one mode
  nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
  mutate(plot = map(
    data,
    ~ .x %>%

      mutate(`Accuracy (AUC)` = auc ) %>%
      mutate(n_cell_type_ = glue("{n_cell_type} cell types"), n_samples_ = glue("{n_samples} samples")) %>%
      rename(Algorithm = name) %>%
      ggplot(aes(slope, `Accuracy (AUC)`, color=fct_relevel(Algorithm, levels_algorithm))) +
      geom_line(size=0.3) +
      facet_grid( fct_reorder(n_cell_type_,n_cell_type)  ~ fct_reorder(n_samples_, n_samples), scale="free_y") +
      geom_hline(yintercept = 0.1, linetype="dashed", color="grey") +
      scale_color_manual(values = cool_palette) +
      scale_y_continuous( labels = dropLeadingZero) +
      labs(color="Algorithm") +
      theme_bw() +
      ylab("Accuracy (AUC) to 0.1 false-positive rate (decimal)") +
      xlab("Slope") +
      multipanel_theme
  )) %>%
  pull(plot)


plot_ROC =
  readRDS(glue("{result_directory}/parsed_1.4_20_20_1000_1.rds")) %>%
  mutate(roc = map(
    df_for_roc,
    ~ .x %>%
      roc_curve(event_level = 'second', truth = significant, probability)
  )) %>%
  select(-df_for_roc) %>%
  unnest(roc) %>%

  #filter( name !="speckle") %>%
  rename(Algorithm = name) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity, color=fct_relevel(Algorithm, levels_algorithm))) +
  geom_rect( mapping=aes(xmin=0, xmax=0.1, ymin=0, ymax=1), fill="#EBECF0", color="#EBECF0", size=0, alpha=0.01) +
  geom_path(size=0.3) +
  geom_abline(lty = 3) +
  scale_color_manual(values = cool_palette) +
  multipanel_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed()+
  guides(color="none") +
  ylab("True positive rate") +
  xlab("False positive rate")



simulated_data =
  readRDS(glue("{result_directory}/input_1.4_20_20_1000_1.rds")) %>%
  pull(data) %>%
  .[[5]] %>%
  mutate(Significant = beta_1!=0) %>%
  # mutate(Outlier = ratio!=1)
  mutate(Outlier = ratio!=1)

plot_simulated =
  ggplot() +
  geom_boxplot(
    aes(factor(type), generated_proportions, fill=Significant),
    outlier.shape = NA,
    data =  filter(simulated_data, !Outlier),
    lwd = 0.3, fatten = 0.5
  ) +
  geom_jitter(aes(factor(type), generated_proportions, color=Outlier), size = 0.5, data = simulated_data) +
  facet_wrap(~ fct_reorder(cell_type, !Significant),  ncol = 10) +
  scale_y_continuous(trans="logit", labels = dropLeadingZero) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  #guides(color="none") +
  xlab("Factor of interest") +
  ylab("Cell-group proportion (decimal)") +
  multipanel_theme +
  theme(strip.background =element_rect(fill="white", color="white"), legend.position = "bottom")

p =
  (
( plot_simulated + plot_ROC ) /
plot_auc[[2]]
)+

  plot_layout(guides = 'collect', height=c(1, 1.5)) +
  plot_annotation(tag_levels = c('A')) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position="bottom")

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


p = plot_auc[[1]]


ggsave(
  "dev/article_figures/plot_benchmark_SUPPLEMENTARY_NO_OUTLIER.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 150/2.5*1.5 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/plot_benchmark_SUPPLEMENTARY_NO_OUTLIER.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 150/2.5*1.5 ,
  limitsize = FALSE
)
