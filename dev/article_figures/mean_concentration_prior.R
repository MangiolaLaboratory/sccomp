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

job({

  load("data/counts_obj.rda")

  counts_obj  |>
    mutate(is_benign = type=="benign") |>
    rename(cell_type = cell_group) %>%
    sccomp_glm(
      formula = ~ is_benign,
      sample, cell_type, count,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association
    ) %>%
    saveRDS("dev/study_of_association/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

      ) %>%
    saveRDS("dev/study_of_association/estimate_GSE139829_uveal_melanoma.rds")

})

job({
  readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

      ) %>%
    saveRDS("dev/study_of_association/estimate_SCP1288_renal_cell_carcinoma.rds")
})

job({
  readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    mutate(type = subtype=="TNBC") %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

      ) %>%
    saveRDS("dev/study_of_association/estimate_SCP1039_bc_cells.rds")
})

job({
  readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

      ) %>%
    saveRDS("dev/study_of_association/estimate_s41587-020-0602-4_COVID_19.rds")
})

job({
  readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

      ) %>%
    saveRDS("dev/study_of_association/estimate_GSE120575_melanoma.rds")
})

job({

  library(tidySingleCellExperiment)

  readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
    filter(ptime %>% is.na() %>% `!`) %>%

    # Scale ptime
    mutate(ptime = scales::rescale(ptime)) %>%
    rename(cell_type = CellTypesFinal) %>%
    rename(sample = Sample) %>%
    sccomp_glm(
      formula = ~ ptime,
      sample, cell_type ,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/estimate_BRCA1_s41467-021-21783-3.rds")
})

# Prior free
job({

  load("data/counts_obj.rda")


  counts_obj  |>
    mutate(is_benign = type=="benign") |>
    rename(cell_type = cell_group) %>%
    sccomp_glm(
      formula = ~ is_benign,
      sample, cell_type, count,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_GSE139829_uveal_melanoma.rds")

})

job({
  readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_SCP1288_renal_cell_carcinoma.rds")
})

job({
  readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    mutate(type = subtype=="TNBC") %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_SCP1039_bc_cells.rds")
})

job({
  readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_s41587-020-0602-4_COVID_19.rds")
})

job({
  readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_GSE120575_melanoma.rds")
})

job({

  library(tidySingleCellExperiment)

  readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
    filter(ptime %>% is.na() %>% `!`) %>%

    # Scale ptime
    mutate(ptime = scales::rescale(ptime)) %>%
    rename(cell_type = CellTypesFinal) %>%
    rename(sample = Sample) %>%
    sccomp_glm(
      formula = ~ ptime,
      sample, cell_type ,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      exclude_priors = TRUE,
      prior_mean_variable_association = prior_mean_variable_association

    ) %>%
    saveRDS("dev/study_of_association/priorFree_estimate_BRCA1_s41467-021-21783-3.rds")
})

# Calculate hyperpriors
associations =
  dir("dev/study_of_association/", pattern = "estimate", full.names = T) %>%
  enframe(value = "file") %>%
  filter(!grepl("hyperprior", file)) %>%
  filter(!grepl("priorFree", file)) %>%

  mutate(data = map(file, ~readRDS(.x))) %>%

  # Process
  mutate(file = (file)) %>%
  mutate(prior = !grepl("priorFree", file)) %>%
  tidyr::extract(file, "dataset", regex = "estimate_([^_]+)_.*.rds") %>%

  # Add lines
  mutate(correlation = map(
    data,
    ~ .x %>%
      attr("mean_concentration_association") %>%
      t() %>%
      as.data.frame %>%
      setNames(c("intercept", "slope", "standard_deviation"))
  )) %>%
  unnest(correlation)

mean(associations$intercept)
sd(associations$intercept)
# 4.136855 1.318474

mean(associations$slope)
sd(associations$slope)
# -0.8452375 0.09860897

fitdistrplus::fitdist(associations$standard_deviation, distr = "gamma", method = "mle") %>%  summary()

# 14.14977 21.98834


# With hyperprior
hyperprior = list(
  intercept = c(4.32452, 1.050276),
  slope = c(-0.8428576, 0.1029212),
  standard_deviation = c(5.390190, 8.746909)
)

job({

  load("data/counts_obj.rda")

  counts_obj  |>
    mutate(is_benign = type=="benign") |>
    rename(cell_type = cell_group) %>%
    sccomp_glm(
      formula = ~ is_benign,
      sample, cell_type, count,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior
    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
})

job({
  readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_GSE139829_uveal_melanoma.rds")

})

job({
  readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(sex))  |>
    sccomp_glm(
      formula = ~ sex,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_SCP1288_renal_cell_carcinoma.rds")
})

job({
  readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
    mutate(type = subtype=="TNBC") %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_SCP1039_bc_cells.rds")
})

job({
  readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
    mutate(is_critical = severity=="critical") %>%
    sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_s41587-020-0602-4_COVID_19.rds")
})

job({
  readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    sccomp_glm(
      formula = ~ time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_GSE120575_melanoma.rds")
})

job({

  library(tidySingleCellExperiment)

  readRDS("/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/sccomp/dev/data_integration/BRCA1_s41467-021-21783-3.rds") %>%
    filter(ptime %>% is.na() %>% `!`) %>%

    # Scale ptime
    mutate(ptime = scales::rescale(ptime)) %>%
    rename(cell_type = CellTypesFinal) %>%
    rename(sample = Sample) %>%
    sccomp_glm(
      formula = ~ ptime,
      sample, cell_type ,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      prior_mean_variable_association = hyperprior

    ) %>%
    saveRDS("dev/study_of_association/hyperprior_estimate_BRCA1_s41467-021-21783-3.rds")
})


# Plot

# Read input
df_for_plot =
  dir("dev/study_of_association/", pattern = "estimate", full.names = T) %>%
  enframe(value = "file") %>%
  mutate(data = imap(file, ~{ print(.y); readRDS(.x)} )) %>%

  # Process
  mutate(file = (file)) %>%
  mutate(prior = case_when(
    grepl("priorFree", file) ~ "none",
    grepl("hyperprior", file) ~ "hyperprior",
    TRUE ~ "prior"
  ) ) %>%
  mutate(prior = factor(prior, levels = c("none", "prior", "hyperprior"))) %>%

  tidyr::extract(file, "dataset", regex = "estimate_([^_]+)_.*.rds") %>%

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
  mutate(datase_size = map_int(
    data,
    ~ .x %>% pull(count_data) %>% .[[1]] %>% nrow
  )) %>%
  unite("dataset", c(dataset, datase_size), sep=" S=") %>%

  # Get data
  mutate(data = map(data,  ~ .x %>% select(composition_CI,concentration)  )) %>%
  unnest(data) %>%
  unnest(c(composition_CI ,  concentration  ))



data_residuals =
  df_for_plot %>%
  select(-intercept, -slope) %>%
  filter(prior=="none") %>%
  nest(data = -dataset) %>%
  mutate(rlm_results = map(
    data,
    ~ .x %>%
      rename(logit_mean = `.median_(Intercept)`) %>%
      MASS::rlm(mean ~ logit_mean , data = .)
  )) %>%
  mutate(residuals = map(
    rlm_results,
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
    list(residuals, weights, dataset),
    ~ enframe(..1) %>%
      mutate(weights = ..2) %>%
      mutate(name = as.numeric(name)) %>%
      ggplot(aes(name, value)) +
      geom_hline(yintercept=0,linetype="dashed") +
      geom_point(size = 0.5) +
      geom_smooth(se=FALSE, span=1, mapping = aes(weight = weights), size=0.5) +
      # stat_smooth(method=function(formula,data,weights=weight) {MASS::rlm(formula,
      #                                                              data,
      #                                                              weights=weight,
      #                                                              method="MM")},
      #             fullrange=TRUE) +
      ggside::geom_ysidedensity() +
      #ggtitle(..3) +
      multipanel_theme +
      theme(axis.title.y = element_blank(), axis.title.x = element_blank())
  ))

plot_residuals =
  data_residuals %>%
  pull(plot) %>%
  wrap_plots(nrow=1)

# Plot no prior
plot_no_prior =
  data_residuals %>%
  unnest(data) %>%
  filter(prior=="none") %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4, size = 0.5 ) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4, size = 0.5) +
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
  facet_wrap(prior ~ dataset, scales = "free_y", nrow=1) +
  scale_color_brewer(palette="Set1") +
  guides(color="none") +
  xlab("Category logit mean") +
  ylab("Category log-concentration")+
  multipanel_theme

plot_prior =
  df_for_plot %>%

  filter(prior!="none") %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4, size = 0.5) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4, size = 0.5) +
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
  facet_wrap(prior ~ dataset, scales = "free_y", nrow=2) +
  scale_color_brewer(palette="Set1") +
  guides(color="none") +
  xlab("Category logit mean") +
  ylab("Category log-concentration")+
  multipanel_theme +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )



 p = ( plot_no_prior / plot_residuals / plot_prior ) +
  plot_layout(guides = "collect", heights = c(1, 1, 2) )  &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'))


 ggsave(
   "dev/article_figures/mean_concentration_association_plot.pdf",
   plot = p,
   units = c("mm"),
   width = 183 ,
   height = 183/7*4 ,
   limitsize = FALSE
 )

 ggsave(
   "dev/article_figures/mean_concentration_association_plot.png",
   plot = p,
   units = c("mm"),
   width = 183 ,
   height = 183/7*4 ,
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
