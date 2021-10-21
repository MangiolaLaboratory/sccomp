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
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/15bcd45b7899b2d72ae46619a3e3cdf8fd6c0e5d/ggplot_theme_multipanel")

  set.seed(42)

# Load multipanel_theme
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

job({
  estimate_UVM =
    readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
    rename(type = `Sample Type`) %>%
    sccomp_glm(
      formula = ~ type,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    ) %>%
    saveRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds")

})

job({
  readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
    mutate(ICB_Response = case_when(
      ICB_Response %in% c("ICB_PR", "ICB_SD") ~ "Responder",
      ICB_Response %in% c("ICB_PD") ~ "Non-responder",
    )) %>%

    tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(ICB_Response) )  |>
    sccomp_glm(
      formula = ~ ICB_Response,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    ) %>%
    saveRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds")
})

job({

  # Analysis about
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/#SD2
  # While each patient showed changes in cluster frequencies between baseline and
  # post-treatment samples (Table S1), there were no consistent changes when aggregating
  # all samples (Table S1).
  readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    separate(res_time, c("res", "time"), sep="_") |>
    sccomp_glm(
      formula = ~ time + res,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    ) %>%
    saveRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds")
})

job({

  # Analysis about
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/#SD2
  # While each patient showed changes in cluster frequencies between baseline and
  # post-treatment samples (Table S1), there were no consistent changes when aggregating
  # all samples (Table S1).
  readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
    separate(res_time, c("res", "time"), sep="_") |>
    sccomp_glm(
      formula = ~ res + time,
      sample, cell_type,
      approximate_posterior_inference = FALSE,
      variance_association = TRUE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    ) %>%
    saveRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds")
})

# https://www.nature.com/articles/s41467-021-21783-3#Abs1
# we find an expansion of Tregs suggesting the early establishment of an immuno-suppressive environment


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
      Sample, cell_type ,
      approximate_posterior_inference = FALSE,
      variance_association = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2))
    ) %>%
    saveRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds")
})



job({
  data_estimates =

    # Import data multi beta
    tribble(
      ~dataset, ~data,


      # Comparing Metastasis vs Primary
      "UVM",
      readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds") %>%
        mutate(heterogeneity_prob_H0 = 1) %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = type) %>%

        # Not tested
        mutate(significant_in_article = FALSE),

      "renal_cell_carcinoma",
      readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds")%>%
        mutate(heterogeneity_prob_H0 = 1) %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = ICB_Response) %>%

        # https://www.sciencedirect.com/science/article/pii/S1535610821001173#fig3
        # The proportions of these cell populations were stable between patients and across ICB exposure status
        mutate(significant_in_article = FALSE),

      # Here factor of interest is subtype== triple-negative versus the rest
      "bc_cells",
      readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
        mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type))))%>%
        unnest(count_data) %>%
        mutate(type = if_else(type =="TRUE", "TNBC", "Other")) %>%
        mutate(factor_of_interest = type) %>%

        # Consistent with the enrichment of TILs and CD8+ T cells in TNBC
        # the T cell clusters IFIT1/c6, LAG3/c8 and MKI67/c11 made up a higher proportion in the TNBC samples
        mutate(significant_in_article = cell_type %in% c("T cells CD8+")),

      "COVID",
      readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds")%>%
        unnest(count_data) %>%
        mutate(is_critical = if_else(is_critical, "Critical", "Non-critical")) %>%
        mutate(factor_of_interest = is_critical) %>%

        # we observed that patients with critical COVID-19 showed a striking depletion of basal cells
        # (P < 1.00 × 10−9) and a strong enrichment for neutrophils (P < 1.00 × 10−9)
        # compared to both patients with moderate COVID-19 and controls
        mutate(significant_in_article = cell_type %in% c("Basal", "Neu")),


      "melanoma_time",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds")%>%
        unnest(count_data) %>%
        mutate(factor_of_interest = time) %>%

        # While each patient showed changes in cluster frequencies between baseline and
        # post-treatment samples (Table S1), there were no consistent changes when aggregating
        # all samples (Table S1).
        mutate(significant_in_article = FALSE),


      "melanoma_responder",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds") %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = res) %>%

        # Here time is pre post checkpoint therapy
        # B cells, CD8+ and CD4+ memory T cells were enriched in responder
        # myeloid cells in non-responder
        mutate(significant_in_article = cell_type %in% c("G1-B cells", "G10-Memory T cells", "G3-Monocytes/Macrophages")),

      "BRCA1_breast",
      readRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds") %>%
        rename(cell_type = CellTypesFinal) %>%
        mutate(heterogeneity_prob_H0 = 1) %>%
        unnest(count_data) %>%
        rename(sample = Sample) %>%
        mutate(factor_of_interest = ptime) %>%

        #
        mutate(significant_in_article = cell_type %in% c("Fb3","Hsp", "Fb9", "Fb2","Fb6", "Fb4", "Hs","cDC1", "CTLs", "CyclingT", "Avd","MdC1", "Mo1", "Tregs", "Mo3","Fb1" ,  "CD83"   ))
    ) %>%

    # Calculate probability
    mutate(data = map(
      data,
      ~ .x %>%
        with_groups(sample, ~ mutate(.x, proportion = (count+1)/sum(count+1)))
    )) %>%


    # Filter significant
    mutate(data = map(
      data,
      ~ .x %>% #filter(.x, heterogeneity_prob_H0 < 0.05 | composition_prob_H0 < 0.05 ) %>%
        mutate(label = case_when(
          composition_prob_H0 < 0.05 ~ "Diff abundant",
          composition_prob_H0 > 0.05 & heterogeneity_prob_H0 < 0.05  ~ "Diff heterogeneous",
          TRUE ~ "Non-significant"
        ) %>% factor(levels = c("Diff abundant", "Diff heterogeneous", "Non-significant"))
        )
    )) %>%

    filter(map_int(data, ~ nrow(.x))>0)
})

job({
  library(tidySingleCellExperiment)
  counts  =
    dir("dev/data_integration/", pattern = "^UMAP_", full.names = TRUE) %>%
    grep("UMAP_oligo", ., invert = TRUE, value = TRUE) %>%
    grep("COVID", ., invert = TRUE, value = TRUE) %>%
    c("dev/data_integration/s41587-020-0602-4_COVID_19.rds") %>%
    map(~ readRDS(.x) %>%
          tidyseurat::as_tibble() %>%
          select(cell, sample, cell_type, UMAP_1, UMAP_2) %>%
          mutate(file = .x)
    ) %>%
    purrr::reduce(bind_rows)
})


# cool_palette = c("#d03161","#ee8080","#bfd8d1","#178a94","#2b374b")
# cool_palette = c("#9cb097", "#e9d6b6", "#8f7b63", "#4c474b", "#415346")
# cool_palette = c("#977b51", "#6a8688", "#714246",  "#233a3c", "#472749", "#d3c1bc", "#9cb097", "#e9d6b6", "#8f7b63", "#4c474b", "#415346")
cool_palette = c("#b58b4c", "#74a6aa", "#a15259",  "#37666a", "#79477c", "#cb9f93", "#9bd18e", "#eece97", "#8f7b63", "#4c474b", "#415346")
color_palette_link = cool_palette %>% setNames(c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))
scales::show_col(cool_palette)

estimate_plots =
  data_estimates %>%

  mutate(data = map(
    data,
    ~ .x %>%
      filter(label %>% is.na %>% `!`) %>%
      unite("color", c(label, significant_in_article)) %>%
      mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))) %>%
      filter(!color %in% c("Non-significant_FALSE", "Diff abundant_TRUE"))

  )) %>%

  filter(map_int(data, ~ nrow(.x))>0) %>%

  # Create plot
  mutate(estimate_plot = purrr::map2(
    data, dataset,
    ~ {

      .x$estimate = .x %>% pull(2)

      if(.y == "BRCA1_breast")
        ggplot() +

        geom_point(aes(factor_of_interest, proportion, color=outlier), size = 0.3, data = .x, height = 0) +
        geom_smooth(
          aes(factor_of_interest, proportion, fill=color), color="black", size = 0.5,
          method="lm",
          data = .x %>% filter(!outlier)
        ) +
        facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=2) +
        scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
        scale_color_manual(values = c("black", "#e11f28")) +
        scale_fill_manual(values = color_palette_link) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Biological condition") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme

      else if(.y == "melanoma_responder")
        ggplot() +
        geom_boxplot(
          aes(factor_of_interest, proportion, fill=color),
          outlier.shape = NA,
          data = .x %>% filter(!outlier), lwd =0.5, fatten = 0.5
        ) +
        geom_point(
          position = position_jitter(width = 0.2, height = 0),
          aes(factor_of_interest, proportion, color=outlier, shape = time),
          size = 0.7, data = .x
        ) +
        facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=1) +
        scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
        scale_color_manual(values = c("black", "#e11f28")) +
        scale_fill_manual(values = color_palette_link) +
        #scale_shape_manual(values = c(1, 19)) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Biological condition") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme

      else
        ggplot() +
        geom_boxplot(
          aes(factor_of_interest, proportion, fill=color),
          outlier.shape = NA,
          data = .x %>% filter(!outlier), lwd =0.5, fatten = 0.5
        ) +
        geom_jitter(aes(factor_of_interest, proportion, color=outlier), size = 0.3, data = .x, height = 0) +
        facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=1) +
        scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
        scale_color_manual(values = c("black", "#e11f28")) +
        scale_fill_manual(values = color_palette_link) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Biological condition") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme
    }

  ))


source("~/PostDoc/oligo_breast//functions.R")

plot_df =
  counts %>%
  with_groups(file, ~ sample_n(.x, 10000)) %>%
  left_join(
    data_estimates %>% mutate(data = map(data, ~ select(.x, sample, cell_type, label, significant_in_article))) %>% unnest(data)
  ) %>% filter(label %>% is.na %>% `!`) %>%
  unite("color", c(label, significant_in_article)) %>%
  mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))) %>%
  mutate(color = case_when(
    color != "Non-significant_FALSE" ~ color
  )) %>%

  # Filter non intersting datasets
  nest(data = -dataset) %>%
  mutate(n = map_int(data, ~ filter(.x, !is.na(color)) %>% nrow())) %>%
  filter(n > 0) %>%

  # Create UMAP
  mutate(UMAP_plot = map2(
    data, dataset,
    ~ ggplot() +
      geom_point(aes(UMAP_1, UMAP_2), data=filter(.x, is.na(color)), fill="grey", size=0.5, alpha=0.5, shape=21, stroke=0) +
      geom_point(aes(UMAP_1, UMAP_2, fill=color), data=filter(.x, !is.na(color)), size=0.5, alpha=1, shape=21, stroke=0) +
      ggrepel::geom_text_repel(
        data= 	get_labels_clusters(
          mutate(.x, cell_type = if_else(!is.na(color), cell_type, "")),
          c(cell_type),
          UMAP_1,
          UMAP_2
        ) ,
        aes(UMAP_1, UMAP_2, label = cell_type), size = 2.5
      ) +
      #facet_wrap(~dataset, scale="free") +
      scale_fill_manual(values = color_palette_link) +
      #scale_fill_brewer(palette="Set1") +
      # scale_fill_manual(values = friendly_cols, na.value = "grey") +
      scale_color_manual(values = friendly_cols, na.value = "grey") +
      guides( fill = "none" ) +
      theme_UMAP +
      theme(title = element_text(size = 7),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()
          ) +
      ggtitle(.y)
  )) %>%

  # add estimate_plots
  left_join( select(estimate_plots, dataset, estimate_plot)   )

# Statistics on novel findings
data_estimates %>%
  mutate(data = map(data, ~ distinct(.x, cell_type, label, significant_in_article))) %>%
  unnest(data) %>%
  unite("color", c(label, significant_in_article)) %>%
  mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))) %>%

  add_count(dataset, name = "dataset_size") %>%
  count(color, dataset, dataset_size) %>%

  nest(data = -dataset) %>%
  mutate(sum_interesting = map_int(data, ~ .x %>% filter(!is.na(color)) %>% pull(n) %>% sum)) %>%
  unnest(data) %>%
  group_by(color) %>%
  summarise(n=sum(n)) %>%
  mutate(prop = n/sum(n))

# How many discrepancies included outliers
data_estimates %>%
  mutate(data = map(data, ~ distinct(.x, cell_type, label, significant_in_article))) %>%
  unnest(data) %>%
  unite("color", c(label, significant_in_article))  %>%

  # Join outlers
  left_join(
    data_estimates %>%
      mutate(number_of_cell_types = map_int(data, ~ distinct(.x, cell_type) %>% nrow())) %>%
      mutate(outliers = map(
        data,
        ~ .x %>%
          distinct(cell_type, outlier)
      )) %>%
      unnest(outliers) %>% select(cell_type, outlier, dataset)
  ) %>%
  filter(color=="Diff abundant_FALSE")

# The comparison between sccomp estimation and the estimation from the selected studies revealed that
# 15% of the calls that disagreed included one or more outliers.
data_estimates %>%
  mutate(data = map(data, ~ distinct(.x, cell_type, label, significant_in_article))) %>%
  unnest(data) %>%
  unite("color", c(label, significant_in_article))  %>%

  # Join outlers
  left_join(
    data_estimates %>%
      mutate(number_of_cell_types = map_int(data, ~ distinct(.x, cell_type) %>% nrow())) %>%
      mutate(outliers = map(
        data,
        ~ .x %>%
          distinct(cell_type, outlier)
      )) %>%
      unnest(outliers) %>% select(cell_type, outlier, dataset)
  ) %>%
  filter(color %in% c("Diff abundant_FALSE", "Non-significant_TRUE")) %>% count(outlier)

plot_novel_results =

  data_estimates %>%
    mutate(data = map(data, ~ distinct(.x, cell_type, label, significant_in_article))) %>%
    unnest(data) %>%
    unite("color", c(label, significant_in_article)) %>%
    mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))) %>%
    mutate(color = case_when(
      color != "Non-significant_FALSE" ~ color
    )) %>%
    add_count(dataset, name = "dataset_size") %>%
    count(color, dataset, dataset_size) %>%

    nest(data = -dataset) %>%
    mutate(sum_interesting = map_int(data, ~ .x %>% filter(!is.na(color)) %>% pull(n) %>% sum)) %>%
    unnest(data) %>%

    ggplot() +
    geom_bar(
      aes(forcats::fct_reorder(dataset, sum_interesting, .desc = TRUE ),  n, fill=fct_rev(factor(color,  exclude = NULL)  )),
      stat = "identity", position = "stack",lwd =0.5, fatten = 0.5
    )+
    scale_fill_manual(values = color_palette_link) +
    xlab("Dataset") +
    ylab("Number of cell-types") +
  multipanel_theme +
    theme( axis.text.x = element_text(angle=30, hjust = 1))

# We identified and quarantined outlier observations in all datasets,
# with an 19% of cell types containing one or more outliers.
data_estimates %>%
  mutate(outliers = map(
    data,
    ~ .x %>%
      distinct(cell_type, outlier)
  )) %>%
  unnest(outliers) %>% filter(outlier)

plot_outliers =
  data_estimates %>%
  mutate(number_of_cell_types = map_int(data, ~ distinct(.x, cell_type) %>% nrow())) %>%
  mutate(outliers = map(
    data,
    ~ .x %>%
      filter(outlier) %>%
      mutate(significant = composition_prob_H0<0.05) %>%
      count(outlier, significant) %>%
      mutate(nn = sum(n))
  )) %>%
  unnest(outliers) %>%

  ggplot() +
  geom_bar(
    aes(forcats::fct_reorder(dataset, nn, .desc = TRUE ),  n, fill=significant  ),
    stat = "identity", position = "stack", lwd =0.5, fatten = 0.5
  )+
  geom_point(aes(dataset, number_of_cell_types) ) +
  scale_fill_manual(values = c("grey", "#e11f28")) +
  xlab("Dataset") +
  ylab("Count of outliers") +
  guides(fill="none") +
  multipanel_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust = 1))


# Mean variance association
job({
  data_mean_variance_association =

    # Import data multi beta
    tribble(
      ~dataset, ~data,


      # Comparing Metastasis vs Primary
      "UVM",
      readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds"),

      "renal_cell_carcinoma",
      readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds"),

      # Here factor of interest is subtype== triple-negative versus the rest
      "bc_cells",
      readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") ,

      "COVID",
      readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds"),


      "melanoma_time",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds"),


      "melanoma_responder",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds") ,

      "BRCA1_breast",
      readRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds")
    )
})

# Plor trends
slopes_df =
  data_mean_variance_association %>%
  mutate(data = map(
    data,
    ~ .x %>%
      attr("mean_concentration_association") %>%
      t() %>%
      as.data.frame %>%
      setNames(c("intercept", "slope", "standard_deviation"))
  )) %>%
  unnest(data)

plot_association_all =
  data_mean_variance_association %>%
  mutate(data = map(data,  ~ .x %>% select(composition_CI,concentration)  )) %>%
  unnest(data) %>%
  unnest(c(composition_CI ,  concentration  )) %>%
  ggplot(aes(`.median_(Intercept)`, -mean)) +
  geom_errorbar(aes(ymin = -`2.5%`, ymax=-`97.5%`, color=dataset),  alpha = 0.4) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
  geom_point(size=0.1) +
  geom_abline(aes(intercept =  -intercept, slope = -slope, color = dataset), linetype = "dashed", data = slopes_df) +
  scale_color_brewer(palette="Set1") +
  xlab("Category logit mean") +
  ylab("Category log-overdispersion") +
  multipanel_theme



# size 1 + 5
p =
  (
  ( plot_novel_results | plot_outliers | plot_association_all ) /
((
  ( plot_df$UMAP_plot[[4]] + plot_df$estimate_plot[[4]] + plot_df$UMAP_plot[[6]] + plot_df$estimate_plot[[6]] +  plot_layout(widths = c(1, 2.5/2, 1, 2.5/2)) ) / # size 1 + 4

    ( plot_df$UMAP_plot[[5]] + plot_df$estimate_plot[[5]] +  plot_layout(widths = c(1, 3.5)) ) / # size 1 + 9

    ( plot_df$UMAP_plot[[3]] + plot_df$estimate_plot[[3]] + plot_df$UMAP_plot[[2]] + plot_df$estimate_plot[[2]] +  plot_layout(widths = c(1, 2, 1, 0.50)) ) / # size 1 + 4

    ( plot_df$UMAP_plot[[1]] + plot_df$estimate_plot[[1]]  +  plot_layout(widths = c(1, 3.5)) )
) +

  # Style
  plot_layout( nrow=4)
  ) )+
  plot_layout(guides = 'collect', height = c(1, 7)) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'))

ggsave(
  "dev/article_figures/novel_results_plot.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 220 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/novel_results_plot.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 220 ,
  limitsize = FALSE
)
