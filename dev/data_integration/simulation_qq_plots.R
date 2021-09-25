# simulation QQ plots

library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)

multipanel_theme =
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(size=0.5),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1),
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 7),
    panel.spacing.x=unit(0.1, "lines"),
    axis.text.x = element_text(size=6),
    axis.text.y = element_text(size=6),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7)
  )

estimate_beta_binomial = 
  
  # Import data multi beta
  tribble(
    ~dataset, ~data,
    "oligo",
    readRDS("dev/data_integration/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds") %>%
      rename(cell_type = cell_group),
    # "UVM",
    # readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma.rds"),
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
  mutate(method = "beta_binomial") 

estimate_dirichlet_multinomial = 
  
  # Import data multi beta
  tribble(
    ~dataset, ~data,
    "oligo",
    readRDS("dev/data_integration/estimate_dirichlet_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds") %>%
      rename(cell_type = cell_group),
    # "UVM",
    # readRDS("dev/data_integration/estimate_dirichlet_GSE139829_uveal_melanoma.rds"),
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
  mutate(method = "dirichlet_multinomial") 

data_for_plot =  
  bind_rows(
    estimate_beta_binomial,
    estimate_dirichlet_multinomial
  ) %>%
  
  # Simulate
  mutate(simulation = map(
    data,
    ~  .x %>%
      simulate_data(sample, cell_type) %>%
      unnest(generated_data) %>%
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
    ~ lm(difference_proportion~proportion, data=.x) %>%
      broom::tidy() %>%
      filter(term!="(Intercept)") %>%
      pull(estimate )
  ))  %>%
  mutate(median_proportion = map_dbl( data, ~ median(.x$proportion) ))

# Import data


plot_qq =
  data_for_plot %>%
  unnest(data) %>%
  ggplot(aes(proportion, generated_proportions, group=cell_type)) +
  geom_point(size=0.2, alpha=0.5) +
  geom_smooth(aes(color=dataset), method = "lm", se = FALSE, size=0.2, alpha=0.5) +
  facet_wrap(method~dataset) +
  geom_abline(linetype="dashed", color="grey") +
  scale_color_brewer(palette="Set1") +
  scale_x_continuous(trans="logit") +
  scale_y_continuous(trans="logit") +
  multipanel_theme

plot_slopes =
  data_for_plot  %>%
  ggplot(aes(slope, color=dataset))+
  geom_density() +
  facet_wrap(~method) +
  scale_color_brewer(palette="Set1") +
  scale_x_log10() +
  multipanel_theme

plot_slopes_median_proportion =
  data_for_plot %>%
  ggplot(aes(median_proportion, slope)) +
  geom_hline(yintercept  = 0,linetype="dashed", color="grey" ) +
  geom_point(aes(color=dataset)) +
  geom_smooth(color="black", aes=0.4, size=0.4, method="lm") +
  facet_wrap(~method) +
  scale_x_continuous(trans="logit") +
  #scale_y_log10() +
  scale_color_brewer(palette="Set1")  +
  multipanel_theme

( plot_qq / plot_slopes ) +
  plot_layout(guides = 'collect', ) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"),  legend.key.size = unit(0.2, 'cm'), legend.position = "bottom")
