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
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/7a6df4a80d38b0427a263c21ac2480c4280cfe4b/ggplot_theme_multipanel")

set.seed(42)

library(magrittr)
library(scales)
library(utils)
logit2_trans <- function(){


  if (find.package("functional", quiet = TRUE) %>% length %>% equals(0)) {
    message("Installing functional needed for analyses")
    install.packages("functional", repos = "https://cloud.r-project.org")
  }

  trans <- qlogis
  inv <- plogis

  trans_new("logit2",
            transform = trans,
            inverse = inv,
            breaks = functional::Compose(trans, extended_breaks(), inv),
            format = label_scientific(digits = 2)
  )
}

cool_palette = c("#b58b4c", "#74a6aa", "#a15259",  "#37666a", "#79477c", "#cb9f93", "#9bd18e", "#eece97", "#8f7b63", "#4c474b", "#415346")

library(glue)
library(forcats)

S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

color_df =
  dir("dev/study_of_association", pattern = "estimate", full.names = T) %>%
  enframe(value = "file") %>%
  mutate(data_type = "RNA") %>%

  bind_rows(
    dir("dev/metagenomics", pattern = "estimate", full.names = T) %>%
      enframe(value = "file") %>%
      mutate(data_type = "metagenomics")
  ) %>%

  bind_rows(
    dir("dev/cytof", pattern = "estimate", full.names = T) %>%
      enframe(value = "file") %>%
      mutate(data_type = "cytof")
  ) %>%
  filter(!grepl(".R$", file)) %>%

  # Filter out hyperprior because not canging
  filter(!grepl("hyperprior", file)) %>%
  filter(!grepl("dirichlet", file)) %>%

  tidyr::extract(file, "dataset", regex = ".*_?estimate_([^_]+)_?.*.rds", remove = F)  %>%
  nest(data = -data_type) %>%
  mutate(color = cool_palette[1:n()]) %>%
  unnest(data) %>%

  # Add colors
  nest(data = -c(data_type, dataset)) %>%
  arrange(fct_relevel(data_type, c("RNA", "cytof", "metagenomics")), dataset) %>%
  mutate(.color = friendly_cols[1:n()]) %>%
  select(-data)



# Ex figure 2, ow in supplementary mean-variability association

# Read input
df_for_plot =
  dir("dev/study_of_association", pattern = "estimate", full.names = T) %>%
  enframe(value = "file") %>%
  mutate(data_type = "RNA") %>%

  bind_rows(
    dir("dev/metagenomics", pattern = "estimate", full.names = T) %>%
      enframe(value = "file") %>%
      mutate(data_type = "metagenomics")
  ) %>%

  bind_rows(
    dir("dev/cytof", pattern = "estimate", full.names = T) %>%
      enframe(value = "file") %>%
      mutate(data_type = "cytof")
  ) %>%
  filter(!grepl(".R$", file)) %>%

  # Filter out hyperprior because not canging
  filter(!grepl("hyperprior", file)) %>%
  filter(!grepl("dirichlet", file)) %>%

  tidyr::extract(file, "dataset", regex = ".*_?estimate_([^_]+)_?.*.rds", remove = F)  %>%

  # Add cell_type for cytof
  mutate(
    data = map( file,  ~ readRDS(.x) )
  ) %>%

  # Process
  mutate(file = (file)) %>%
  mutate(prior = case_when(
    grepl("priorFree", file) ~ "none",
    grepl("hyperprior", file) ~ "hyperprior",
    TRUE ~ "prior"
  ) ) %>%
  mutate(prior = factor(prior, levels = c("none", "prior", "hyperprior"))) %>%

  # Add lines
  mutate(correlation = map(
    data,
    ~ .x %>%
      attr("mean_concentration_association") %>%
      t() %>%
      as.data.frame %>%
      setNames(c("intercept", "slope", "standard_deviation"))
  )) %>%
  unnest(correlation) %>%

  # hide
  mutate( intercept = if_else(prior!="none", intercept, 99) ) %>%

  # Get dataset size
  mutate(dataset_size = map_int(
    data,
    ~ .x %>% pull(count_data) %>% .[[1]] %>% nrow
  )) %>%
  unite("dataset_size", c(dataset, dataset_size), sep=" S=", remove = FALSE) %>%

  # Get data
  mutate(data = map(data,  ~ select( .x, composition_CI,concentration, cell_type )  )) %>%
  unnest(data) %>%
  unnest(c(composition_CI ,  concentration  ))


data_residuals =
  df_for_plot %>%
  select(-intercept, -slope) %>%
  filter(prior=="none") %>%
  nest(data = -c(dataset, data_type)) %>%
  mutate(rlm_results = map(
    data,
    ~ .x %>%
      rename(logit_mean = `.median_(Intercept)`) %>%
      MASS::rlm(mean ~ logit_mean , data = .)
  )) %>%
  mutate(residuals = map2(
    rlm_results, data,
    ~ resid(.x)
  )) %>%
  mutate(weights = map(
    rlm_results,
    ~ .x$w
  )) %>%
  mutate(coefficients = map(
    rlm_results,
    ~ .x$coefficients %>%
      setNames(c("intercept", "slope")) %>%
      enframe() %>%
      spread(name, value)
  )) %>%
  unnest(coefficients)

plot_residuals =
  data_residuals %>%

  # Group
  mutate(data_plot = pmap(
    list(residuals, weights, dataset, data),
    ~ enframe(..1) %>%
      mutate(m = ..4 %>% pull(`.median_(Intercept)`)) %>%
      mutate(weights = ..2) %>%
      mutate(name = as.numeric(name))
  )) %>% select(data_plot, data_type, dataset) %>%
  unnest(data_plot) %>%
  nest(data_plot = -data_type) %>%

  # Plot
  mutate(plot = map(
    data_plot,
    ~ .x %>%
      ggplot(aes(m, -value, group=dataset, color = dataset)) +
      geom_hline(yintercept=0,linetype="dashed") +
      geom_point(size=0.1, alpha=0.6) +
      geom_smooth(se=FALSE, span=1, mapping = aes(weight = weights), size=0.5, linetype = "dashed", color="#333333") +

      ggside::geom_ysidedensity() +
      ggside::scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
      ggside::scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
      scale_color_manual(values = color_df %>% distinct(dataset, .color) %>% deframe()) +
      xlab("Cell groups/taxa") +
      ylab("Residuals") +
      guides(color="none") +
      multipanel_theme +
      theme(axis.title.x = element_blank())
  ))


# Plot no prior
plot_no_prior =
  data_residuals %>%
  unnest(data) %>%
  filter(prior=="none") %>%
  nest(data = -data_type) %>%
  mutate(plot = map(
    data,
    ~ .x %>%
      ggplot(aes(`.median_(Intercept)`, -mean, group=dataset, color=dataset)) +
      geom_point(size=0.1, alpha=0.6) +
      geom_abline(
        aes(intercept =  -intercept, slope = -slope, color=dataset),
        linetype = "dotted",
        alpha=0.5,
        color="#333333"
      ) +
      scale_color_manual(values = color_df %>% distinct(dataset, .color) %>% deframe()) +
      guides(color="none") +
      xlab("Inverse-multinomial-logit mean") +
      ylab("Log-concentration")+
      multipanel_theme +
      theme(axis.title.x = element_blank())
  ))

# Shrinkage
plot_shrinkage =
  df_for_plot %>%
  mutate(
    diff_in_concentration = abs(`97.5%`-`2.5%`) ,
    diff_in_mean = abs(`.upper_(Intercept)`-`.lower_(Intercept)`)
  ) %>%

  # Scale
  with_groups(dataset, ~ .x %>% mutate(
    diff_in_mean= diff_in_mean %>% scale(center = FALSE) %>% as.numeric(),
    diff_in_concentration= diff_in_concentration %>% scale(center = FALSE) %>% as.numeric(),
  )) %>%

  # Add diff of diff
  nest(data = -c( cell_type, dataset)) %>%

  mutate(log_diff_diff_mean = map_dbl(
    data,
    ~ { v = .x %>%
      arrange(prior) %>%
      pull(diff_in_mean)
    log(v[2]/v[1])

    }

  )) %>%

  mutate(log_diff_diff_concentration = map_dbl(
    data,
    ~ { v = .x %>%
      arrange(prior) %>%
      pull(diff_in_concentration)
    log(v[2]/v[1])

    }

  )) %>%
  unnest(data) %>%

  # Plot
  nest(data = -data_type) %>%

  mutate(plot_diff_concentration = map(
    data,
    ~ .x %>%
      filter(abs(log_diff_diff_concentration)>0.5) %>%
      ggplot( aes(prior, diff_in_concentration)) +
      geom_point(alpha=0.5, size=0.3) +
      geom_line(aes(group=cell_type, color=dataset), alpha=0.5, size=0.5) +
      scale_color_manual(values = color_df %>% distinct(dataset, .color) %>% deframe()) +
      scale_size(range = c(0.01, 0.1)) %>%
      guides(color="none") +
      multipanel_theme +
      theme(axis.title.x = element_blank())
  )) %>%
  mutate(plot_diff_mean = map(
    data,
    ~  .x %>%
      filter(abs(log_diff_diff_mean)>0.5) %>%
      ggplot( aes(prior, diff_in_mean)) +
      geom_point(alpha=0.5, size=0.3) +
      geom_line(aes(group=cell_type, color=dataset), alpha=0.5, size=0.5) +
      scale_color_manual(values = color_df %>% distinct(dataset, .color) %>% deframe()) +
      scale_size(range = c(0.01, 0.1)) %>%
      guides(color="none") +
      multipanel_theme +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ))

rm(data_estimates)
rm(data_for_plot)
rm(df_for_plot)
gc()

# Figure 2 Q-Q plots

data_for_plot =
    dir("dev/data_integration", pattern = "estimate", full.names = T) %>%
    enframe(value = "file") %>%
    mutate(data_type = "RNA") %>%

    bind_rows(
      dir("dev/metagenomics", pattern = "estimate", full.names = T) %>%
        enframe(value = "file") %>%
        mutate(data_type = "metagenomics")
    ) %>%

    bind_rows(
      dir("dev/cytof", pattern = "estimate", full.names = T) %>%
        enframe(value = "file") %>%
        mutate(data_type = "cytof")
    ) %>%
    filter(!grepl(".R$", file)) %>%

    # Filter out hyperprior because not canging
    filter(!grepl("hyperprior", file)) %>%
    filter(!grepl("priorFree", file)) %>%
    filter(!grepl("VarianceApprox", file)) %>%

    # Set method
    mutate(method = if_else(grepl("dirichlet", file), "Dm", "scBb")) %>%
    mutate(method = factor(method, levels = c("scBb", "Dm"))) %>%

    # Create dataset
    tidyr::extract(file, c("dummy", "dataset"), regex = ".*_?estimate_(dirichlet_)?([^_]+)_?.*.rds", remove = F)  %>%

    # Color
    nest(data = -dataset) %>%
    mutate(color = cool_palette[1:n()]) %>%
    unnest(data) %>%

    #filter(method=="Dm") %>%

    # Add cell_type for cytof
    mutate(  data = map( file,  ~ readRDS(.x) ) ) %>%

    # Simulate
    mutate(simulation = imap(
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
      ~ .x %>%
        distinct() %>%
        unnest(count_data) %>%
        when("abundance" %in% colnames(.) ~ rename(., count = abundance), ~ (.)) %>%
        when("Absolute.Abundance" %in% colnames(.) ~ rename(., count = Absolute.Abundance), ~ (.)) %>%
        select(cell_type, sample, one_of("outlier"), count) %>%

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
    filter(is.na(outlier) | !outlier) %>%
    nest(data = -c(data_type, dataset, cell_type, method)) %>%
    mutate(slope = map_dbl(
      data,
      ~ lm(generated_proportions~proportion, data=.x) %>%
        broom::tidy() %>%
        filter(term!="(Intercept)") %>%
        pull(estimate )
    ))  %>%
    mutate(median_proportion = map_dbl( data, ~ median(.x$proportion) ))


data_for_plot %>% saveRDS("dev/data_for_plot.rds")

# data_for_plot =
#   bind_rows(
#
#     # Import data multi beta
#     tribble(
#       ~dataset, ~data,
#       "PBMC",
#       readRDS("dev/data_integration/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds") ,
#       "UVM",
#       readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma.rds"),
#       "RCC",
#       readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma.rds"),
#       "BRCA",
#       readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
#         mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type)))),
#       "COVID",
#       readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds"),
#       "SKCM",
#       readRDS("dev/data_integration/estimate_GSE120575_melanoma.rds")
#     ) %>%
#       mutate(method = "Simplex Beta-binomial"),
#
#
#     # Import data multi beta
#     tribble(
#       ~dataset, ~data,
#       "PBMC",
#       readRDS("dev/data_integration/estimate_dirichlet_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds"),
#       "UVM",
#       readRDS("dev/data_integration/estimate_dirichlet_GSE139829_uveal_melanoma.rds"),
#       "RCC",
#       readRDS("dev/data_integration/estimate_dirichlet_SCP1288_renal_cell_carcinoma.rds"),
#       "BRCA",
#       readRDS("dev/data_integration/estimate_dirichlet_SCP1039_bc_cells.rds") %>%
#         mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type)))),
#       "COVID",
#       readRDS("dev/data_integration/estimate_dirichlet_s41587-020-0602-4_COVID_19.rds"),
#       "SKCM",
#       readRDS("dev/data_integration/estimate_dirichlet_GSE120575_melanoma.rds")
#     ) %>%
#       mutate(method = "Dirichlet-multinomial")
#   ) %>%
#
#   mutate(method = factor(method, levels = c("Simplex Beta-binomial", "Dirichlet-multinomial"))) %>%
#
#   # Simulate
#   mutate(simulation = map(
#     data,
#     ~  .x %>%
#       replicate_data() %>%
#       with_groups(
#         cell_type,
#         ~.x %>%
#           arrange(generated_proportions ) %>%
#           mutate(sample = 1:n() )
#       )
#   )) %>%
#
#   # Calculate proportion
#   mutate(data = map(
#     data,
#     ~ .x %>% distinct() %>%
#       unnest(count_data) %>%
#       select(cell_type, sample, outlier, count) %>%
#       with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) %>%
#       with_groups(cell_type, ~.x %>% arrange(proportion) %>%  mutate(sample = 1:n() ))
#   )) %>%
#
#   # Join data
#   mutate(data = map2(
#     data, simulation,
#     ~ left_join(.x, .y, by = c("cell_type", "sample") ) %>%
#       mutate(difference_proportion = generated_proportions - proportion )
#   )) %>%
#   select(-simulation) %>%
#
#
#   # linear model for qq plot
#   unnest(data) %>%
#   filter(!outlier) %>%
#   nest(data = -c(dataset, cell_type, method)) %>%
#   mutate(slope = map_dbl(
#     data,
#     ~ lm(generated_proportions~proportion, data=.x) %>%
#       broom::tidy() %>%
#       filter(term!="(Intercept)") %>%
#       pull(estimate )
#   ))  %>%
#   mutate(median_proportion = map_dbl( data, ~ median(.x$proportion) ))

# Boxplot posterior predictive check
gc()

covid =  readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds")

data_proportion =
    covid %>%
  distinct() %>%
  unnest(count_data) %>%
  mutate(significant = composition_pH0_is_criticalTRUE < 0.025) %>%
  select(cell_type, sample, outlier, count, is_critical, significant, composition_effect_is_criticalTRUE) %>%
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) %>%
  mutate(Condition = if_else(is_critical, "Critical", "Moderate")) %>%
  mutate(Effect =composition_effect_is_criticalTRUE  )


simulated_proportion =
    covid %>%
  replicate_data( number_of_draws = 100) %>%
  left_join(data_proportion %>% distinct(is_critical, sample, cell_type, composition_effect_is_criticalTRUE))  %>%
  mutate(Condition = if_else(is_critical, "Critical", "Moderate")) %>%
  mutate(Effect =composition_effect_is_criticalTRUE  )

rm(covid)
gc()

data_simulation_process =
  list(
    data_proportion %>% mutate(step = "Data") %>% mutate(outlier = FALSE) %>% mutate(Effect = NA),
    data_proportion %>% mutate(step = "Find outlier") %>% mutate(Effect = NA),
    data_proportion %>% mutate(step = "Fit model"),
    simulated_proportion %>% rename(proportion = generated_proportions) %>% mutate(step = "Simulate data") %>% filter(replicate==1) %>% mutate(outlier = FALSE) %>% mutate(Effect = NA)
  ) %>%
  purrr::reduce(bind_rows) %>%
  filter(cell_type == "Neu") %>%
  mutate(step = forcats::fct_relevel(step, unique(.$step))) %>%
  mutate(composition_effect_is_criticalTRUE = case_when(step=="Fit model" ~ composition_effect_is_criticalTRUE))


# PLOTS
plot_simulation_process =
  ggplot() +

  geom_boxplot(
    aes(Condition, proportion, fill = factor(composition_effect_is_criticalTRUE)),
    outlier.shape = NA,
    data = data_simulation_process |> filter(!outlier),
    fatten = 0.5, size=0.5,
  ) +
  geom_point(aes(Condition, proportion, color=outlier, shape=outlier), position = position_jitter(seed = 41), size = 0.5, data = data_simulation_process) +
  facet_wrap(~step, ncol = 1) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = "#dd6572", na.value = "white") +
  scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +
  #scale_y_continuous(labels = dropLeadingZero, trans="logit") +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  guides(fill = "none", shape="none", color="none") +
  coord_cartesian( clip = "off") +
  multipanel_theme +
  theme(axis.text.x =  element_text(angle=20, hjust = 1), axis.title.x = element_blank())
# +
#   theme(strip.clip = "off")

# data_for_boxplot =
#   data_proportion |>
#   mutate(which="observed") |>
#   bind_rows(
#     simulated_proportion |>
#       rename(proportion = generated_proportions) |>
#       mutate(which="simulated") |>
#       mutate(outlier = FALSE) |>
#       mutate(Effect=NA)
#   )

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

plot_boxplot =
  ggplot() +

  stat_summary(
    aes(Condition, (generated_proportions)),
    fun.data = calc_boxplot_stat, geom="boxplot",
    fatten = 0.5, lwd=0.2,
    data =
      simulated_proportion %>%

      # Filter uanitles because of limits
      inner_join( data_proportion %>% distinct(Condition, cell_type, Effect)) ,
    color="blue"

  ) +

  geom_boxplot(
    aes(Condition, proportion, fill=Effect, group=interaction( Condition)),
    outlier.shape = NA,
    data = data_proportion |> filter(!outlier), fatten = 0.5, lwd=0.5,
  ) +
  geom_jitter(
    aes(Condition, proportion, shape=outlier,  group=interaction( Condition)),
    data = data_proportion,
    position=position_jitterdodge(jitter.height = 0, jitter.width = 0.2),
    size = 0.5
  ) +

  # geom_boxplot(
  #   aes(Condition, generated_proportions),
  #   outlier.shape = NA, alpha=0.2,
  #   data = simulated_proportion, fatten = 0.5, size=0.5,
  # ) +
  # geom_jitter(aes(Condition, generated_proportions), color="black" ,alpha=0.2, size = 0.2, data = simulated_proportion) +

  facet_wrap(
    ~ forcats::fct_reorder(cell_type, abs(Effect), .desc = TRUE, na.rm=TRUE),
    scales = "free_y", nrow = 4
  ) +
  #scale_color_manual(values = c("black", "#e11f28")) +
  #scale_fill_manual(values = c("white", "#E2D379")) +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  #scale_color_distiller(palette = "Spectral") +
  scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +
  #scale_y_continuous(labels = dropLeadingZero, trans="logit") +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  guides(color="none", alpha="none", size="none") +
  labs(fill="Compositional difference") +
  multipanel_theme +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x =  element_text(angle=20, hjust = 1))

get_plot_qq = function(.data, colors){
  .data=
    .data  %>%
    filter(method=="scBb") %>%
    unnest(data)

  .data%>%
    ggplot(aes(proportion, generated_proportions, group=cell_type)) +
    geom_point(size=0.5, alpha=0.5, shape=21, stroke=0, fill="black" ) +
    geom_smooth(
      aes(color=fct_reorder( dataset, data_type)),
      method = "lm", se = FALSE, size=0.2, alpha=0.3
    ) +
    facet_wrap(  fct_reorder( dataset, data_type) ~ ., nrow=2) +
    geom_abline(linetype="dashed", color="grey") +
    scale_color_manual(values = colors) +

    # scale_x_continuous(trans="S_sqrt", labels = dropLeadingZero, breaks = c(0, max(.data$proportion))) +
    # scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero, breaks = c(0, max(.data$generated_proportions))) +

    scale_x_continuous(trans="logit2", breaks = c(min(.data$proportion[.data$proportion>0]), max(.data$proportion))) +
    scale_y_continuous(trans="logit2", breaks = c(min(.data$generated_proportions[.data$generated_proportions>0]), max(.data$generated_proportions))) +

    xlab("Observed proportion") +
    ylab("Simulated proportion") +
    guides(color = "none") +
    multipanel_theme +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
}

plot_qq_RNA =
  data_for_plot %>%
  filter(data_type=="RNA") %>%
  filter(dataset!="GSE115189") %>% # Drop the made up dataset
  get_plot_qq(colors = friendly_cols[1:6]) +
  theme(axis.title.x = element_blank()) +
  ggtitle("RNA")

plot_qq_cytof =
  data_for_plot %>%
  filter(data_type=="cytof") %>%
  get_plot_qq(colors = friendly_cols[7:12]) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  ggtitle("CyTOF")

plot_qq_metagenomics =
  data_for_plot %>%
  filter(data_type=="metagenomics") %>%
  get_plot_qq(colors = friendly_cols[13:18]) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  theme(axis.title.x = element_blank()) +
  ggtitle("Microbiome")

plot_slopes =
  data_for_plot  %>%
  filter(median_proportion > 0) %>%
  filter(method=="scBb") %>%
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
  {
    .x = (.)

    ggplot(.x, aes(median_proportion, slope)) +
      geom_hline(yintercept  = 1,linetype="dashed", color="grey" ) +
      geom_point(aes(color=dataset), size=0.2) +
      geom_smooth(color="black", aes=0.4, size=0.4, method="lm") +
      facet_grid(method~data_type, scales = "free_x") +
      scale_x_continuous(trans="logit2", breaks = c(min(.x$median_proportion), max(.x$median_proportion))) +
      scale_y_log10() +
      ylab("Slope qq-plot") +
      xlab("Observed proportion") +
      scale_color_manual(values = friendly_cols)  +
      multipanel_theme  +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
  }



plot_mean_variability_assoc =
(
  plot_no_prior %>%
    arrange(fct_relevel(data_type, c("RNA", "cytof", "metagenomics"))) %>%
    pull(plot) %>%
    wrap_plots(ncol = 1)
) |

(
  plot_residuals %>%
    arrange(fct_relevel(data_type, c("RNA", "cytof", "metagenomics"))) %>%
    pull(plot) %>%
    wrap_plots(ncol = 1)
) |

(
  plot_shrinkage %>%
    arrange(fct_relevel(data_type, c("RNA", "cytof", "metagenomics"))) %>%
    pull(plot_diff_concentration ) %>%
    wrap_plots(ncol = 1)
) |

(
  plot_shrinkage %>%
    arrange(fct_relevel(data_type, c("RNA", "cytof", "metagenomics"))) %>%
    pull(plot_diff_mean) %>%
    wrap_plots(ncol = 1)
)


p =

  (
    # Mean variability assoc
    plot_mean_variability_assoc /

    # Boxplots
    ( ( plot_simulation_process | plot_boxplot )  +  plot_layout(widths = c(1,8)) ) /

  (
    ( (
      (
       # Dotplot
       plot_qq_RNA + plot_qq_cytof + plot_qq_metagenomics
      ) /
       plot_slopes_median_proportion
    ) |
      plot_slopes )+  plot_layout(width  = c(3,1))
  )

  ) +

  # Style
  plot_layout(guides = 'collect', heights  = c(2, 3, 2))  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom",  legend.key.size = unit(0.5, 'cm'))


ggsave(
  "dev/article_figures/qq_plot.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 225 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/qq_plot.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 225 ,
  limitsize = FALSE
)


