library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)
library(forcats)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/2a5e64e29f79a1e30405903fec1280df8662eff4/ggplot_theme_multipanel")

prior_mean_variable_association = list(
  intercept = c(5,10),
  slope = c(0,  5),
  standard_deviation = c(1,1)
)

friendly_cols <- dittoSeq::dittoColors()
cool_palette = c("#b58b4c", "#74a6aa", "#a15259",  "#37666a", "#79477c", "#cb9f93", "#9bd18e", "#eece97", "#8f7b63", "#4c474b", "#415346")

library(shades)
# http://hughjonesd.github.io/tweaking-colours-with-the-shades-package.html
fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")

# job({
#
#   load("data/counts_obj.rda")
#
#   counts_obj  |>
#     mutate(is_benign = type=="benign") |>
#     rename(cell_type = cell_group) %>%
#     sccomp_glm(
#       formula = ~ is_benign,
#       sample, cell_type, count,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#     ) %>%
#     saveRDS("dev/study_of_association/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
# })
#
# job({
#   readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
#     rename(type = `Sample Type`) %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#       ) %>%
#     saveRDS("dev/study_of_association/estimate_GSE139829_uveal_melanoma.rds")
#
# })
#
# job({
#   readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
#     tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
#     sccomp_glm(
#       formula = ~ sex,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#       ) %>%
#     saveRDS("dev/study_of_association/estimate_SCP1288_renal_cell_carcinoma.rds")
# })
#
# job({
#   readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
#     mutate(type = subtype=="TNBC") %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#       ) %>%
#     saveRDS("dev/study_of_association/estimate_SCP1039_bc_cells.rds")
# })
#
# job({
#   readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
#     mutate(is_critical = severity=="critical") %>%
#     sccomp_glm(
#       formula = ~ is_critical,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#       ) %>%
#     saveRDS("dev/study_of_association/estimate_s41587-020-0602-4_COVID_19.rds")
# })
#
# job({
#   readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
#     sccomp_glm(
#       formula = ~ time,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#       ) %>%
#     saveRDS("dev/study_of_association/estimate_GSE120575_melanoma.rds")
# })
#
# job({
#
#   library(tidySingleCellExperiment)
#
#   readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
#     filter(ptime %>% is.na() %>% `!`) %>%
#
#     # Scale ptime
#     mutate(ptime = scales::rescale(ptime)) %>%
#     rename(cell_type = CellTypesFinal) %>%
#     rename(sample = Sample) %>%
#     sccomp_glm(
#       formula = ~ ptime,
#       sample, cell_type ,
#       approximate_posterior_inference = FALSE,
#       variance_association = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/estimate_BRCA1_s41467-021-21783-3.rds")
# })
#
# # Prior free
# job({
#
#   load("data/counts_obj.rda")
#
#
#   counts_obj  |>
#     mutate(is_benign = type=="benign") |>
#     rename(cell_type = cell_group) %>%
#     sccomp_glm(
#       formula = ~ is_benign,
#       sample, cell_type, count,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
# })
#
# job({
#   readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
#     rename(type = `Sample Type`) %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_GSE139829_uveal_melanoma.rds")
#
# })
#
# job({
#   readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
#     tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
#     sccomp_glm(
#       formula = ~ sex,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_SCP1288_renal_cell_carcinoma.rds")
# })
#
# job({
#   readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
#     mutate(type = subtype=="TNBC") %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_SCP1039_bc_cells.rds")
# })
#
# job({
#   readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
#     mutate(is_critical = severity=="critical") %>%
#     sccomp_glm(
#       formula = ~ is_critical,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_s41587-020-0602-4_COVID_19.rds")
# })
#
# job({
#   readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
#     sccomp_glm(
#       formula = ~ time,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_GSE120575_melanoma.rds")
# })
#
# job({
#
#   library(tidySingleCellExperiment)
#
#   readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
#     filter(ptime %>% is.na() %>% `!`) %>%
#
#     # Scale ptime
#     mutate(ptime = scales::rescale(ptime)) %>%
#     rename(cell_type = CellTypesFinal) %>%
#     rename(sample = Sample) %>%
#     sccomp_glm(
#       formula = ~ ptime,
#       sample, cell_type ,
#       approximate_posterior_inference = FALSE,
#       variance_association = FALSE,
#       exclude_priors = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#
#     ) %>%
#     saveRDS("dev/study_of_association/priorFree_estimate_BRCA1_s41467-021-21783-3.rds")
# })

# Calculate hyperpriors
associations =
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
  filter(!grepl("priorFree", file)) %>%
  filter(!grepl("dirichlet", file)) %>%

  tidyr::extract(file, "dataset", regex = ".*_?estimate_([^_]+)_?.*.rds", remove = F)  %>%

  mutate(data = map(file, ~readRDS(.x))) %>%

  # Add lines
  mutate(correlation = map(
    data,
    ~ .x %>%
      attr("mean_concentration_association") %>%
      as.numeric() %>%
      t() %>%
      as.data.frame %>%
      setNames(c("intercept", "slope", "standard_deviation"))
  )) %>%
  unnest(correlation)

associations %>%
with_groups(data_type, ~summarise(.x,
                                  mean_i = mean(intercept),
                                  sd_i = sd(intercept),
                                  min_i = min(intercept),
                                  max_i = max(intercept),
                                  mean_s = mean(slope),
                                  sd_s = sd(slope),
                                  min_s = min(slope),
                                  max_s = max(slope)
                                  ))




# -0.8452375 0.09860897

fitdistrplus::fitdist(associations$standard_deviation, distr = "gamma", method = "mle") %>%  summary()

# 14.14977 21.98834


# With hyperprior
hyperprior = list(
  intercept = c(4.32452, 1.050276),
  slope = c(-0.8428576, 0.1029212),
  standard_deviation = c(5.390190, 8.746909)
)

# job({
#
#   load("data/counts_obj.rda")
#
#   counts_obj  |>
#     mutate(is_benign = type=="benign") |>
#     rename(cell_type = cell_group) %>%
#     sccomp_glm(
#       formula = ~ is_benign,
#       sample, cell_type, count,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
# })
#
# job({
#   readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
#     rename(type = `Sample Type`) %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_GSE139829_uveal_melanoma.rds")
#
# })
#
# job({
#   readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
#     tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
#     sccomp_glm(
#       formula = ~ sex,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_SCP1288_renal_cell_carcinoma.rds")
# })
#
# job({
#   readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
#     mutate(type = subtype=="TNBC") %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_SCP1039_bc_cells.rds")
# })
#
# job({
#   readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
#     mutate(is_critical = severity=="critical") %>%
#     sccomp_glm(
#       formula = ~ is_critical,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_s41587-020-0602-4_COVID_19.rds")
# })
#
# job({
#   readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
#     sccomp_glm(
#       formula = ~ time,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_GSE120575_melanoma.rds")
# })
#
# job({
#
#   library(tidySingleCellExperiment)
#
#   readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
#     filter(ptime %>% is.na() %>% `!`) %>%
#
#     # Scale ptime
#     mutate(ptime = scales::rescale(ptime)) %>%
#     rename(cell_type = CellTypesFinal) %>%
#     rename(sample = Sample) %>%
#     sccomp_glm(
#       formula = ~ ptime,
#       sample, cell_type ,
#       approximate_posterior_inference = FALSE,
#       variance_association = FALSE,
#       prior_mean_variable_association = hyperprior
#
#     ) %>%
#     saveRDS("dev/study_of_association/hyperprior_estimate_BRCA1_s41467-021-21783-3.rds")
# })


# Plot

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
  nest(data = -data_type) %>%
  mutate(color = cool_palette[1:n()]) %>%
  unnest(data) %>%

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

df_for_plot %>% saveRDS("dev/df_for_plot.rds")

plot_shrinkage =
  df_for_plot %>%
  mutate(
    diff_in_concentration = abs(`97.5%`-`2.5%`),
    diff_in_mean = abs(`.upper_(Intercept)`-`.lower_(Intercept)`)
  ) %>%

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
    ~ ggplot(.x, aes(prior, diff_in_concentration)) +
      geom_point(alpha=0.5, size=0.3) +
      geom_line(aes(group=cell_type, color=log_diff_diff_concentration),  size=0.1) +
      facet_wrap(~ dataset, scales="free_y", ncol=1, strip.position="right") +
      scale_colour_gradient2(
        low="#053061",mid= "grey", high="#67001f",  midpoint = 0,
        limits = c(-quantile(abs(.x$log_diff_diff_concentration), 0.9), quantile(abs(.x$log_diff_diff_concentration), 0.9))
      ) +
      guides(color="none") +
      multipanel_theme +
      theme(axis.title.y = element_blank(), strip.text.y = element_text(angle=0, hjust = 0))
  )) %>%
  mutate(plot_diff_mean = map(
    data,
    ~ ggplot(.x, aes(prior, diff_in_mean)) +
      geom_point(alpha=0.5, size=0.3) +
      geom_line(aes(group=cell_type, color=log_diff_diff_mean),  size=0.1) +
      facet_wrap(~ dataset, scales="free_y", ncol=1, strip.position="right") +
      scale_colour_gradient2(
        low="#053061",mid= "grey", high="#67001f",  midpoint = 0,
        limits = c(-quantile(abs(.x$log_diff_diff_concentration), 0.9), quantile(abs(.x$log_diff_diff_concentration), 0.9))
      ) +
      guides(color="none") +
      multipanel_theme +
      theme(strip.background.y  = element_blank(), strip.text.y = element_blank())
  ))


data_residuals =
  df_for_plot %>%
  select(-intercept, -slope) %>%
  filter(prior=="none") %>%
  nest(data = -c(dataset, data_type, color)) %>%
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
  unnest(coefficients) %>%
  # # library("metRology")
  # mutate(t =map(
  #   residuals,
  #   ~ fitdistrplus::fitdist(.x, distr = "t.scaled", start=list(df=3,mean=mean(.x),sd=sd(.x))) %>%  summary()
  # )) %>%

  mutate(plot = pmap(
    list(residuals, weights, dataset, data),
    ~ enframe(..1) %>%
      mutate(m = ..4 %>% pull(`.median_(Intercept)`)) %>%
      mutate(weights = ..2) %>%
      mutate(name = as.numeric(name)) %>%
      ggplot(aes(m, value)) +
      geom_hline(yintercept=0,linetype="dashed") +
      geom_point(size = 0.5) +
      geom_smooth(se=FALSE, span=1, mapping = aes(weight = weights), size=0.5) +
      # stat_smooth(method=function(formula,data,weights=weight) {MASS::rlm(formula,
      #                                                              data,
      #                                                              weights=weight,
      #                                                              method="MM")},
      #             fullrange=TRUE) +
      ggside::geom_ysidedensity() +
      ggside::scale_xsidey_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
      ggside::scale_ysidex_continuous(breaks = NULL, labels = "", expand = expansion(c(0,.1))) +
      #ggtitle(..3) +
      xlab("Cell groups/taxa") +
      ylab("Residuals") +
      multipanel_theme +
      theme(axis.title.y  = element_blank())
  ))





plot_residuals =
  data_residuals %>%
  nest(data = -data_type) %>%
  mutate(plot = map(
    data,
    ~ .x %>%
      pull(plot) %>%
      wrap_plots(nrow=1)
  ))

# Plot no prior
plot_no_prior =
  data_residuals %>%
  nest(data = -data_type) %>%
  mutate(plot = map(
    data,
    ~ .x %>%
      unnest(data) %>%
      filter(prior=="none") %>%
      ggplot(aes(`.median_(Intercept)`, mean)) +
      geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`),  color=cool_palette[3],  alpha = 0.4, size = 0.5 ) +
      geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`),  color=cool_palette[3], alpha = 0.4, size = 0.5) +
      geom_point(size=0.1) +
      geom_abline(
        aes(intercept =  intercept, slope = slope),
        linetype = "dotted",
        alpha=0.5
      ) +
      # geom_abline(
      #   aes(intercept =  1, slope = slope),
      #   linetype = "dashed"
      # ) +
      facet_wrap( ~ dataset_size, scales = "free", nrow=1) +
      #scale_color_brewer(palette="Set1") +
      #scale_color_manual(values = unique(.x$color)) +
      guides(color="none") +
      xlab("Inverse-multinomial-logit mean") +
      ylab("Log-concentration")+
      multipanel_theme
  ))


plot_prior =
  df_for_plot %>%
  nest(data = -data_type) %>%
  mutate(plot = map(
    data,
    ~ .x %>%
      filter(prior!="none") %>%
      ggplot(aes(`.median_(Intercept)`, mean)) +
      geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`),  color=cool_palette[2],  alpha = 0.4, size = 0.5) +
      geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`),  color=cool_palette[2], alpha = 0.4, size = 0.5) +
      geom_point(size=0.1) +
      geom_abline(
        aes(intercept =  intercept, slope = slope),
        linetype = "dashed",
        alpha=0.5
      ) +
      geom_abline(
        aes(intercept =  1, slope = slope),
        linetype = "dashed",
        alpha=0.5
      ) +
      facet_wrap(prior ~ dataset_size, scales = "free", ncol=7) +
      #scale_color_manual(values = unique(.x$color)) +
      guides(color="none") +
      xlab("Inverse-multinomial-logit mean") +
      ylab("Log-concentration")+
      multipanel_theme +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
  ))




 p =
 (

 (
    filter(plot_no_prior, data_type=="RNA")$plot[[1]] /
      filter(plot_residuals, data_type=="RNA")$plot[[1]]  /
      filter(plot_prior, data_type=="RNA")$plot[[1]] /

     filter(plot_no_prior, data_type=="cytof")$plot[[1]]   /
  filter(plot_residuals, data_type=="cytof")$plot[[1]]     /
   filter(plot_prior, data_type=="cytof")$plot[[1]]   /

    filter(plot_no_prior, data_type=="metagenomics")$plot[[1]] /
    filter(plot_residuals, data_type=="metagenomics")$plot[[1]]  /
    filter(plot_prior, data_type=="metagenomics")$plot[[1]]
 ) | (
   filter(plot_shrinkage, data_type=="RNA")$plot_diff_mean[[1]] +
     filter(plot_shrinkage, data_type=="RNA")$plot_diff_concentration[[1]] +
     filter(plot_shrinkage, data_type=="cytof")$plot_diff_mean[[1]] +
     filter(plot_shrinkage, data_type=="cytof")$plot_diff_concentration[[1]] +
     filter(plot_shrinkage, data_type=="metagenomics")$plot_diff_mean[[1]] +
     filter(plot_shrinkage, data_type=="metagenomics")$plot_diff_concentration[[1]] +
     plot_layout(ncol=2 ,byrow = TRUE)
 )
 )+
  plot_layout(guides = "collect", width = c( 4,1) )  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'))


 ggsave(
   "dev/article_figures/mean_concentration_association_plot.pdf",
   plot = p,
   units = c("mm"),
   width = 183 ,
   height = 183 ,
   limitsize = FALSE
 )

 ggsave(
   "dev/article_figures/mean_concentration_association_plot.png",
   plot = p,
   units = c("mm"),
   width = 183 ,
   height = 183 ,
   limitsize = FALSE
 )


# slopes_df =
#   data_mean_variance_association %>%
#   mutate(data = map(
#     data,
#     ~ .x %>%
#       attr("mean_concentration_association") %>%
#       t() %>%
#       as.data.frame %>%
#       setNames(c("intercept", "slope", "standard_deviation"))
#   )) %>%
#   unnest(data)

plot_association_all =
  data_mean_variance_association %>%
  mutate(data = map(data,  ~ .x %>% select(composition_CI,concentration)  )) %>%
  unnest(data) %>%
  unnest(c(composition_CI ,  concentration  )) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  intercept, slope = slope, color = dataset), linetype = "dashed", data = slopes_df) +
  scale_color_brewer(palette="Set1") +
  xlab("Category logit mean") +
  ylab("Category log-concentration") +
  multipanel_theme

plot_association_all =
  list(
  estimate_benign %>% mutate(dataset = "estimate_benign"),
  estimate_oligo %>% mutate(dataset = "estimate_oligo"),
  estimate_UVM %>% mutate(dataset = "estimate_UVM"),
  estimate_renal_cell_carcinoma %>% mutate(dataset = "estimate_renal_cell_carcinoma"),
  estimate_bc_cells %>% mutate(dataset = "estimate_bc_cells"),
  estimate_COVID %>% mutate(dataset = "estimate_COVID"),
  estimate_melanoma %>% mutate(dataset = "estimate_melanoma")
) %>%
  reduce(bind_rows) %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  intercept, slope = slope, color = dataset), linetype = "dashed", data = slopes_df) +
  scale_color_brewer(palette="Set1") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw() +
  theme(legend.position = "none")

plot_association_one =
  list(
  estimate_benign %>% mutate(dataset = "estimate_benign")
) %>%
  reduce(bind_rows) %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  intercept, slope = slope), linetype = "dashed", color="grey", data = slopes_df %>% filter(dataset == "estimate_benign")) +
  scale_color_brewer(palette="Set1") +
  guides(color='none') +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw()+
  theme(legend.position = "bottom")


# CyTOF
estimate_CyTOF =
  readRDS("dev/data_integration/Cytof_raw.rds") %>%
  as.data.frame() %>%
  setNames(c("cell_cluster", "sample", "count")) %>%
  as_tibble() %>%
  tidyr::extract(sample, c("factor_of_interest", "patient", "batch"), "([A-Z]+)([0-9]+)_B(.+)", remove = FALSE) %>%
  filter(batch=="1") %>%
  sccomp_glm(
    formula = ~ 1,
    sample, cell_cluster,count,
    approximate_posterior_inference = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(5,5))
  )
# estimate_CyTOF %>% attr("mean_concentration_association")
# [1]  4.2623500    -1.1306609
# prec_sd  = 0.4845913

estimate_CyTOF %>%
  unnest(concentration, composition_CI) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(intercept =  4.2623500, slope =  -1.1306609, linetype = "dashed") +
  scale_color_brewer(palette="Set1") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw() +
  theme(legend.position = "none")






(
  plot_association_one +  plot_association_all &
    plot_layout(guides = "collect") &
    theme(legend.position = 'bottom', legend.key.size = unit(0.2, 'cm'))
) %>%
  ggsave(
    "dev/data_integration/association_plot.pdf",
    plot = .,
    units = c("mm"),
    width = 183 ,
    height = 110 ,
    limitsize = FALSE
  )
