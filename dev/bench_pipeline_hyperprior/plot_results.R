library(tidyverse)
library(patchwork)
library(bayestestR)
library(yardstick)
library(forcats)
library(glue)
library(ggplot2)
library(purrr)
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/a0853a1a4e8a46baf33bad6268b09001d49faf51/ggplot_theme_multipanel")

#
# dir("dev/benchmark_results_f58cd58973261636d29f0462eb9b406d762f3895 345a2afa986458672dd95f5dff3f5d973bf1d994//", pattern = "auc", full.names = TRUE) %>%
#   map_dfr(~ .x %>% readRDS()) %>%
#   nest(data = -c(add_outliers, max_cell_counts_per_sample)) %>%
#   pull(data) %>% .[[2]] %>%
#
#
#   nest(data = -c(  slope, n_samples, n_cell_type )) %>%
#   mutate(min_auc = map_dbl(data, ~ min(.x$auc))) %>%
#   unnest(data) %>%
#
#   mutate(`Accuracy (AUC) relative to the smallest` = auc / min_auc) %>%
#
#
#   ggplot(aes(slope, `Accuracy (AUC) relative to the smallest`, color=name)) +
#   geom_line() +
#   facet_grid( n_cell_type ~ n_samples, scale="free_y") +
#   theme_bw()

cool_palette = c("#cc374a","#e29a3b", "#3a9997", "#a15259",   "#94379b", "#ff8d73", "#7be05f", "#ffc571", "#724c24", "#4c474b", "#236633")

levels_algorithm =
  c("sccomp", "speckle" , "logitLinear" , "ttest"  , "quasiBinomial" , "rlm" , "DirichletMultinomial" , "random"   )

plot_auc =
  dir("~/PostDoc/sccomp_dev/dev/benchmark_results_hyperprior_19Dec2022/", pattern = "auc", full.names = TRUE) %>%
  c(dir("~/PostDoc/sccomp_dev/dev/benchmark_results_hyperprior/", pattern = "auc", full.names = TRUE) ) |>
  map_dfr(~ .x %>% readRDS()) %>%
  nest(data = c(name, auc)) %>%
  mutate(random_auc = map_dbl(data, ~ .x %>% filter(name=="random") %>% pull(auc) %>% head(1))) %>%
  unnest(data) %>%

  # filter( name !="speckle") %>%
  #filter(n_samples >2) %>%
  mutate(name = if_else(name=="random", "zrandom", name)) %>%
  filter(name != "zrandom") |>
  # Filter too small slope
  filter(slope>0.1) %>%

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
      theme_bw() +
      ylab("Accuracy (AUC) up to 0.1 false-positive-rate") +
      xlab("Slope") +
      theme_multipanel +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  )) %>%
  pull(plot) %>%
  patchwork::wrap_plots(ncol = 1)

ggsave(
  "~/PostDoc/sccomp_dev/dev/plot_auc_hyperprior.pdf",
  plot = plot_auc,
  units = c("mm"),
  width = 183 ,
  height = 140 ,
  limitsize = FALSE
)

plot_ROC =
  readRDS("~/PostDoc/sccomp_dev/dev/benchmark_results/parsed_1.4_20_20_1000_1.rds") %>%
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
  geom_rect( mapping=aes(xmin=0, xmax=0.1, ymin=0, ymax=1), fill="#EBECF0", color="#EBECF0", alpha=0.01) +
  geom_path() +
  geom_abline(lty = 3) +
  scale_color_brewer(palette="Paired") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_fixed()+
  ylab("True positive rate") +
  xlab("False positive rate")

ggsave(
  "dev/plot_ROC_WEHI_seminar.pdf",
  plot = plot_ROC,
  units = c("mm"),
  width = 183 ,
  height = 183/2 ,
  limitsize = FALSE
)

