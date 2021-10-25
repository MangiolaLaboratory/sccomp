# simulation QQ plots

# library(tidyverse)

library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)
library(stringr)

library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/305ed9ba2af815fdef3214b9e6171008d5917400/ggplot_theme_multipanel")

set.seed(42)


data_for_plot =
  bind_rows(

    # Import data multi beta
    tribble(
      ~dataset, ~data,
      "oligo",
      readRDS("dev/data_integration/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds") ,
      "UVM",
      readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma.rds"),
      "renal_cell_carcinoma",
      readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma.rds"),
      "bc_cells",
      readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
        mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type)))),
      "COVID",
      readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds"),
      "melanoma",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma.rds")
    ) %>%
      mutate(method = "Simplex Beta-binomial"),


    # Import data multi beta
    tribble(
      ~dataset, ~data,
      "oligo",
      readRDS("dev/data_integration/estimate_dirichlet_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds"),
      "UVM",
      readRDS("dev/data_integration/estimate_dirichlet_GSE139829_uveal_melanoma.rds"),
      "renal_cell_carcinoma",
      readRDS("dev/data_integration/estimate_dirichlet_SCP1288_renal_cell_carcinoma.rds"),
      "bc_cells",
      readRDS("dev/data_integration/estimate_dirichlet_SCP1039_bc_cells.rds") %>%
        mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type)))),
      "COVID",
      readRDS("dev/data_integration/estimate_dirichlet_s41587-020-0602-4_COVID_19.rds"),
      "melanoma",
      readRDS("dev/data_integration/estimate_dirichlet_GSE120575_melanoma.rds")
    ) %>%
      mutate(method = "Dirichlet-multinomial")
  ) %>%

  mutate(method = factor(method, levels = c("Simplex Beta-binomial", "Dirichlet-multinomial"))) %>%

  # Simulate
  mutate(simulation = map(
    data,
    ~  .x %>%
      replicate_data() %>%
      with_groups(
        cell_type,
        ~.x %>%
          arrange(generated_proportions ) %>%
          mutate(sample = 1:n() )
      )
  )) %>%

  # Calculate proportion
  mutate(data = map(
    data,
    ~ .x %>% distinct() %>%
      unnest(count_data) %>%
      select(cell_type, sample, outlier, count) %>%
      with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) %>%
      with_groups(cell_type, ~.x %>% arrange(proportion) %>%  mutate(sample = 1:n() ))
  )) %>%

  # Join data
  mutate(data = map2(
    data, simulation,
    ~ left_join(.x, .y, by = c("cell_type", "sample") ) %>%
      mutate(difference_proportion = generated_proportions - proportion )
  )) %>%
  select(-simulation) %>%


  # linear model for qq plot
  unnest(data) %>%
  filter(!outlier) %>%
  nest(data = -c(dataset, cell_type, method)) %>%
  mutate(slope = map_dbl(
    data,
    ~ lm(generated_proportions~proportion, data=.x) %>%
      broom::tidy() %>%
      filter(term!="(Intercept)") %>%
      pull(estimate )
  ))  %>%
  mutate(median_proportion = map_dbl( data, ~ median(.x$proportion) ))

# Boxplot posterior predictive check
data_proportion =
  readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds") %>%
  distinct() %>%
  unnest(count_data) %>%
  mutate(significant = composition_prob_H0 < 0.025) %>%
  select(cell_type, sample, outlier, count, is_critical, significant, composition_effect_is_criticalTRUE) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) )

simulated_proportion =
  readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds") %>%
  replicate_data( number_of_draws = 10) %>%
  left_join(data_proportion %>% distinct(is_critical, sample, cell_type, composition_effect_is_criticalTRUE))

data_simulation_process =
  list(
    data_proportion %>% mutate(step = "Data (proportion representation)") %>% mutate(outlier = FALSE) %>% mutate(composition_effect_is_criticalTRUE = NA),
    data_proportion %>% mutate(step = "Outlier identification") %>% mutate(composition_effect_is_criticalTRUE = NA),
    data_proportion %>% mutate(step = "Model fitting"),
    simulated_proportion %>% rename(proportion = generated_proportions) %>% mutate(step = "Data simulation") %>% filter(replicate==1) %>% mutate(outlier = FALSE) %>% mutate(composition_effect_is_criticalTRUE = NA)
  ) %>%
  reduce(bind_rows) %>%
  filter(cell_type == "Neu") %>%
  mutate(step = forcats::fct_relevel(step, unique(.$step)))


# PLOTS


plot_simulation_process =
  ggplot() +

  geom_boxplot(
    aes(is_critical, proportion, fill = factor(composition_effect_is_criticalTRUE)),
    outlier.shape = NA,
    data = data_simulation_process |> filter(!outlier),
    fatten = 0.5, size=0.5,
  ) +
  geom_point(aes(is_critical, proportion, color=outlier, shape=outlier), position = position_jitter(seed = 41), size = 1, data = data_simulation_process) +
  facet_wrap(~step, ncol = 1) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = "#dd6572", na.value = "white") +
  scale_y_continuous(labels = dropLeadingZero, trans="logit") +
  xlab("Biological condition") +
  ylab("Cell-group proportion (decimal)") +
  guides(fill = "none", shape="none", color="none") +
  coord_cartesian( clip = "off") +
  multipanel_theme
# +
#   theme(strip.clip = "off")


plot_boxplot =
  ggplot() +

  geom_boxplot(
    aes(is_critical, proportion, fill=composition_effect_is_criticalTRUE),
    outlier.shape = NA,
    data = data_proportion |> filter(!outlier), fatten = 0.5, size=0.5,
  ) +
  geom_jitter(aes(is_critical, proportion, shape=outlier, color = composition_effect_is_criticalTRUE),  data = data_proportion) +

  geom_boxplot(
    aes(is_critical, generated_proportions),
    outlier.shape = NA, alpha=0.2,
    data = simulated_proportion, fatten = 0.5, size=0.5,
  ) +
  geom_jitter(aes(is_critical, generated_proportions), color="black" ,alpha=0.2, size = 0.2, data = simulated_proportion) +

  facet_wrap(~ forcats::fct_reorder(cell_type, abs(composition_effect_is_criticalTRUE), .desc = TRUE), scales = "free_y", nrow = 4) +
  #scale_color_manual(values = c("black", "#e11f28")) +
  #scale_fill_manual(values = c("white", "#E2D379")) +
  scale_fill_distiller(palette = "Spectral") +
  scale_color_distiller(palette = "Spectral") +
  scale_y_continuous(labels = dropLeadingZero, trans="logit") +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  multipanel_theme +
  theme(axis.title.y = element_blank())



plot_qq =
  data_for_plot %>%
  filter(method=="Simplex Beta-binomial") %>%
  unnest(data) %>%
  ggplot(aes(proportion, generated_proportions, group=cell_type)) +
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(aes(color=dataset), method = "lm", se = FALSE, size=0.2, alpha=0.3) +
  facet_wrap(~dataset) +
  geom_abline(linetype="dashed", color="grey") +
  scale_color_manual(values = friendly_cols) +
  scale_x_continuous(trans="logit", labels = dropLeadingZero) +
  scale_y_continuous(trans="logit", labels = dropLeadingZero) +
  xlab("Observed proportion (decimal)") +
  ylab("Simulated proportion (decimal)") +
  guides(color = "none") +
  multipanel_theme

plot_slopes =
  data_for_plot  %>%
  filter(median_proportion > 0) %>%
  filter(method=="Simplex Beta-binomial") %>%
  ggplot(aes(slope, color=dataset))+
  geom_vline(xintercept  = 1,linetype="dashed", color="grey" ) +
  geom_density() +
  #facet_wrap(~method) +
  scale_color_manual(values = friendly_cols) +
   scale_x_log10() +
  xlab("Slope qq-plot") +
  ylab("Density") +
  guides(color = "none") +
  multipanel_theme

plot_slopes_median_proportion =
  data_for_plot %>%
  filter(median_proportion > 0) %>%
  ggplot(aes(median_proportion, slope)) +
  geom_hline(yintercept  = 1,linetype="dashed", color="grey" ) +
  geom_point(aes(color=dataset), size=0.4) +
  geom_smooth(color="black", aes=0.4, size=0.4, method="lm") +
  facet_wrap(~method) +
  scale_x_continuous(trans="logit", labels = dropLeadingZero) +
  scale_y_log10() +
  ylab("Slope qq-plot") +
  xlab("Median observed proportion (decimal)") +
  scale_color_manual(values = friendly_cols)  +
  multipanel_theme
# + theme(strip.clip = "off")

p =

  (
    # Boxplots
    ( ( plot_simulation_process | plot_boxplot )  +  plot_layout(widths = c(1,5)) ) /


   (
     # Dotplot
     plot_qq +  plot_slopes + plot_slopes_median_proportion  +
      plot_layout(widths = c(3,1,2))
   )
  ) +

  # Style
  plot_layout(guides = 'collect', heights  = c(2, 1)) + plot_annotation(tag_levels = c('A')) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")

ggsave(
  "dev/article_figures/qq_plot.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 130 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/qq_plot.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 130 ,
  limitsize = FALSE
)
