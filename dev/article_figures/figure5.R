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
library(ggsignif)
library(scales)

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

S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)


prior_mean_variable_association = list(intercept = c(5, 5), slope = c(0,  5), standard_deviation = c(5,5))

# job({
#   estimate_UVM =
#     readRDS("dev/data_integration/UVM_single_cell/counts.rds")  |>
#     rename(type = `Sample Type`) %>%
#     sccomp_glm(
#       formula = ~ type,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       variance_association = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#     ) %>%
#     saveRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds")
#
# })
#
# job({
#   readRDS("dev/data_integration/SCP1288_renal_cell_carcinoma.rds")  |>
#     mutate(ICB_Response = case_when(
#       ICB_Response %in% c("ICB_PR", "ICB_SD") ~ "Responder",
#       ICB_Response %in% c("ICB_PD") ~ "Non-responder",
#     )) %>%
#
#     tidyseurat::filter(!is.na(sample) & !is.na(cell_type) & !is.na(ICB_Response) )  |>
#     sccomp_glm(
#       formula = ~ ICB_Response,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       prior_mean_variable_association = prior_mean_variable_association
#     ) %>%
#     saveRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds")
# })
#
# job({
#
#   # Analysis about
#   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/#SD2
#   # While each patient showed changes in cluster frequencies between baseline and
#   # post-treatment samples (Table S1), there were no consistent changes when aggregating
#   # all samples (Table S1).
#   readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
#     separate(res_time, c("res", "time"), sep="_") |>
#     sccomp_glm(
#       formula = ~ time + res,
#       formula_variability = ~ time,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       variance_association = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#     ) %>%
#     saveRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds")
# })
#
# job({
#
#   # Analysis about
#   # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/#SD2
#   # While each patient showed changes in cluster frequencies between baseline and
#   # post-treatment samples (Table S1), there were no consistent changes when aggregating
#   # all samples (Table S1).
#   readRDS("dev/data_integration/GSE120575_melanoma.rds")  |>
#     separate(res_time, c("res", "time"), sep="_") |>
#     sccomp_glm(
#       formula = ~ res + time,
#       formula_variability = ~ res,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       variance_association = TRUE,
#       prior_mean_variable_association = prior_mean_variable_association
#     ) %>%
#     saveRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds")
# })

# https://www.nature.com/articles/s41467-021-21783-3#Abs1
# we find an expansion of Tregs suggesting the early establishment of an immuno-suppressive environment


job({
  data_estimates =

    # Import data multi beta
    tribble(
      ~dataset, ~data,


      # Comparing Metastasis vs Primary
      "UVM",
      readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds") %>%
        mutate(v_pH0= 1, v_FDR=1) %>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%

        unnest(count_data) %>%
        mutate(factor_of_interest = type) %>%

        # Not tested
        mutate(significant_in_article = FALSE)  %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_typePrimary`) %>%
        mutate(v_FDR = v_FDR_typePrimary)  ,

      "renal_cell_carcinoma",
      readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds") %>%
        mutate(v_pH0 = 1, v_FDR = 1) %>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = ICB_Response) %>%

        # https://www.sciencedirect.com/science/article/pii/S1535610821001173#fig3
        # The proportions of these cell populations were stable between patients and across ICB exposure status
        mutate(significant_in_article = FALSE) %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_ICB_ResponseResponder`) %>%
        mutate(v_FDR = `v_FDR_ICB_ResponseResponder`),

      # Here factor of interest is subtype== triple-negative versus the rest
      "bc_cells",
      readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type))))%>%
        unnest(count_data) %>%
        #mutate(type = if_else(type =="TRUE", "TNBC", "Other")) %>%
        mutate(factor_of_interest = type) %>%

        # Consistent with the enrichment of TILs and CD8+ T cells in TNBC
        # the T cell clusters IFIT1/c6, LAG3/c8 and MKI67/c11 made up a higher proportion in the TNBC samples
        mutate(significant_in_article = cell_type %in% c("T cells CD8+", "T_cells_c6_IFIT1", "T_cells_c11_MKI67", "T_cells_c8_CD8+_LAG3", "Myeloid_c1_LAM1_FABP5")) %>%

        # Readapt to ols format
        mutate(c_FDR = pmin(`c_FDR_typeER+`, `c_FDR_typeHER2+`)) %>%
        mutate(v_FDR = pmin(`v_FDR_typeER+`, `v_FDR_typeHER2+`))        ,

      "COVID",
      readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds")%>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        unnest(count_data) %>%
        mutate(is_critical = if_else(is_critical, "Critical", "Non-critical")) %>%
        mutate(factor_of_interest = is_critical) %>%

        # we observed that patients with critical COVID-19 showed a striking depletion of basal cells
        # (P < 1.00 × 10−9) and a strong enrichment for neutrophils (P < 1.00 × 10−9)
        # compared to both patients with moderate COVID-19 and controls
        mutate(significant_in_article = cell_type %in% c("Basal", "Neu"))  %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_is_criticalTRUE`) %>%
        mutate(v_FDR = v_FDR_is_criticalTRUE)   ,


      "melanoma_time",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds")%>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = time) %>%

        # While each patient showed changes in cluster frequencies between baseline and
        # post-treatment samples (Table S1), there were no consistent changes when aggregating
        # all samples (Table S1).
        mutate(significant_in_article = FALSE)  %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_timePre`) %>%
        mutate(v_FDR = v_FDR_timePre),


      "melanoma_responder",
      readRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds") %>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        unnest(count_data) %>%
        mutate(factor_of_interest = res) %>%

        # Here time is pre post checkpoint therapy
        # B cells, CD8+ and CD4+ memory T cells were enriched in responder
        # myeloid cells in non-responder
        mutate(significant_in_article = cell_type %in% c("G1-B cells", "G10-Memory T cells", "G3-Monocytes/Macrophages")) %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_resResponder`) %>%
        mutate(v_FDR = v_FDR_resResponder),

      "BRCA1_breast",

      readRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds") %>%
       # rename(cell_type = CellTypesFinal) %>%
        mutate(v_FDR = 1) %>%
        pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
        unnest(count_data) %>%
        #rename(sample = Sample) %>%
        mutate(factor_of_interest = ptime) %>%

        #
        mutate(significant_in_article = cell_type %in% c("Fb3","Hsp", "Fb9", "Fb2","Fb6", "Fb4", "Hs","cDC1", "CTLs", "CyclingT", "Avd","MdC1", "Mo1", "Tregs", "Mo3","Fb1" ,  "CD83"   )) %>%

        # Readapt to ols format
        mutate(c_FDR = `c_FDR_ptime`) %>%
        mutate(v_FDR = `v_FDR_ptime`)
    ) %>%

    # Calculate probability
    mutate(data = map(
      data,
      ~ .x %>%
        with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)))
    )) %>%


    # Filter significant
    mutate(data = map(
      data,
      ~ .x %>% #filter(.x, variability_prob_H0 < 0.05 | composition_prob_H0 < 0.05 ) %>%
        mutate(label = case_when(
          c_FDR < 0.025 & v_FDR < 0.025 ~ "Diff both",
          c_FDR < 0.025 ~ "Diff abundant",
          v_FDR < 0.025  ~ "Diff heterogeneous",
          TRUE ~ "Non-significant"
        ) %>% factor(levels = c("Diff abundant", "Diff heterogeneous", "Diff both",  "Non-significant"))
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
    grep("GSE115189", ., invert=TRUE, value=TRUE) %>%
    c("dev/data_integration/s41587-020-0602-4_COVID_19.rds") %>%
    imap(~ readRDS(.x) %>% {print(.y); (.)} %>%
          tidyseurat::as_tibble() %>%

          # if brca change cell type column
           when(
             .x=="dev/data_integration//UMAP_SCP1039_bc_cells.rds" ~ mutate(., cell_type = celltype_subset) ,
             ~ (.)
            ) %>%

          select(.cell, sample, cell_type, UMAP_1, UMAP_2) %>%
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


add_multi_shades = function(.data, .palette_col, .color_col){

  .palette_col = enquo(.palette_col)
  .color_col = enquo(.color_col)

  palettes = c("Reds", "Greens", "Purples", "Blues", "RdPu", "Greys", "YlOrBr")

  palettes_df =
    tibble(!!.palette_col := levels(pull(!!.data, !!.palette_col))) %>%
    mutate(.palette = palettes[1:n()])

  .data %>%
    nest(data__ = -!!.palette_col) %>%
    left_join(palettes_df, by=quo_name(.palette_col)) %>%

    mutate(data__ = map2(
      data__, .palette,
      ~ .x %>%
        nest(data__ = -!!.color_col) %>%
        mutate(.color = gradient_n_pal(brewer_pal(palette=.y)(9))(seq(0, 1, length.out=n())) ) %>%
        unnest(data__)
    )) %>%
    unnest(data__)


}

estimate_plots =
  data_estimates %>%

  mutate(data = map(
    data,
    ~ .x %>%
      filter(label %>% is.na %>% `!`) %>%
      unite("color", c(label, significant_in_article)) %>%
      mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Diff both_TRUE", "Diff both_FALSE", "Non-significant_FALSE"))) %>%
      filter(!color %in% c("Non-significant_FALSE"))

  )) %>%

  filter(map_int(data, ~ nrow(.x))>0) %>%

  # Create plot
  mutate(estimate_plot = purrr::map2(
    data, dataset,
    ~ {

      .x$estimate = .x %>% select(contains("c_effect")) %>%  pull(1)

      if(.y == "BRCA1_breast"){
        .x = .x %>% filter(color != "Diff abundant_TRUE")

        ggplot() +

        geom_point(aes(factor_of_interest, proportion, color=outlier), size = 0.3, data = .x, height = 0) +
        geom_smooth(
          aes(factor_of_interest, proportion, fill=color), color="black", size = 0.2,
          method="lm", alpha=1,
          data = .x %>% filter(!outlier)
        ) +
        facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=2) +
        scale_x_continuous(breaks = c(min(.x$factor_of_interest), max(.x$factor_of_interest))) +
          scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +

        scale_color_manual(values = c("black", "#e11f28")) +
          scale_fill_manual(values =  c(
            "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
            "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
            "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
            "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
            "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
          )) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Pseudo-time") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      }
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
        #scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
        scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +

        scale_color_manual(values = c("black", "#e11f28")) +
        scale_fill_manual(values =  c(
          "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
          "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
          "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
          "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
          "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
        )) +
        #scale_shape_manual(values = c(1, 19)) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Biological condition") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme+
        theme(axis.text.x =  element_text(angle=20, hjust = 1)) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )

      else
        ggplot() +
        geom_boxplot(
          aes(factor_of_interest, proportion, fill=color),
          outlier.shape = NA,
          data = .x %>% filter(!outlier), lwd =0.5, fatten = 0.5
        ) +
        geom_jitter(aes(factor_of_interest, proportion, color=outlier), size = 0.3, data = .x, height = 0) +
        facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=1) +
        #scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
        scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +

        scale_color_manual(values = c("black", "#e11f28")) +
      scale_fill_manual(values =  c(
        "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
        "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
        "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
        "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
        "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
      )) +
        #scale_fill_brewer(palette="Set1") +
        xlab("Biological condition") +
        ylab("Cell-group proportion (decimal)") +
        guides(fill="none", color="none") +
        multipanel_theme+
        theme(axis.text.x =  element_text(angle=20, hjust = 1)) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
    }

  ))


source("~/PostDoc/oligo_breast//functions.R")

plot_df =
  counts %>%
  #with_groups(file, ~ sample_n(.x, 10000)) %>%
  inner_join(
    data_estimates %>% mutate(data = map(data, ~ select(.x, sample, cell_type, label, significant_in_article))) %>% unnest(data),
    by = c("sample", "cell_type")
  ) %>%
  #filter(label %>% is.na %>% `!`) %>%
  unite("color", c(label, significant_in_article)) %>%
  mutate(color = case_when(
    color != "Non-significant_FALSE" ~ color
  )) %>%

  # Filter non intersting datasets
  nest(data = -dataset) %>%
  mutate(n = map_int(data, ~ filter(.x, !is.na(color)) %>% nrow())) %>%
  filter(n > 0 | dataset %>% is.na) %>%
  unnest(data) %>%

  # Recover dataset - this is result of bad programming
  nest(data = -file) %>%
  mutate(data = map(
    data,
    ~ {
      ds = filter(.x, dataset %>% is.na %>% `!`) %>% pull(dataset) %>% unique %>% as.character
      mutate(.x, dataset = ds)
    }
  )) %>%
  unnest(data) %>%
  mutate(color = case_when(
    !(is.na(color) | color=="NA_NA") ~ color
  )) %>%
  mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Diff both_TRUE", "Diff both_FALSE", "Non-significant_FALSE"))) %>%


  # Create UMAP
  nest(data = -dataset) %>%

  # Subset
  mutate(data = map2(
    data, dataset,
    ~ when(.x, .y != "bc_cells" ~ sample_n(., 10000), ~ (.))
  )) %>%

  # lot
  mutate(UMAP_plot = map2(
    data, dataset,
    ~ {

      .x = .x %>%
        # Add colouring
        mutate(color = factor(color, levels = c( "Non-significant_TRUE" , "Diff abundant_FALSE" ,  "Diff heterogeneous_FALSE" ,   "Diff both_FALSE",  "Diff abundant_TRUE" ))) %>%
        add_multi_shades(color, cell_type)

      ggplot() +
      geom_point(aes(UMAP_1, UMAP_2), data=filter(.x, is.na(color)), fill="grey", size=0.5, alpha=0.5, shape=21, stroke=0) +
      geom_point(aes(UMAP_1, UMAP_2), data=filter(.x, !is.na(color)), fill=filter(.x, !is.na(color)) %>% pull(.color), size=0.5, alpha=1, shape=21, stroke=0) +
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
      #guides( fill = "none" ) +
      theme_UMAP +
      theme(title = element_text(size = 7),
            axis.line = element_blank(),
            axis.ticks = element_blank()
          )
    }
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
  mutate(dataset = case_when(
    dataset == "BRCA1_breast" ~ "BRCA1",
    dataset == "renal_cell_carcinoma" ~ "RCC",
    dataset == "bc_cells" ~ "BRCA",
    dataset == "melanoma_time" ~ "SKCM time",
    dataset == "melanoma_responder" ~ "SKCM resp.",
    TRUE ~ dataset
  )) %>%
    nest(data = -dataset) %>%
    mutate(sum_interesting = map_int(data, ~ .x %>% filter(!is.na(color)) %>% pull(n) %>% sum)) %>%
    unnest(data) %>%

    ggplot() +
    geom_bar(
      aes(forcats::fct_reorder(dataset, sum_interesting, .desc = TRUE ),  n, fill=fct_rev(factor(color,  exclude = NULL)  )),
      stat = "identity", position = "stack",lwd =0.5, fatten = 0.5
    )+
    scale_fill_manual(values = c(
      "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
      "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
      "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
      "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
      "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
    )) +
    xlab("Dataset") +
    ylab("Number of cell-types") +
  multipanel_theme +
    theme( axis.text.x = element_text(angle=30, hjust = 1))  +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

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
  mutate(dataset = case_when(
    dataset == "BRCA1_breast" ~ "BRCA1",
    dataset == "renal_cell_carcinoma" ~ "RCC",
    dataset == "bc_cells" ~ "BRCA",
    dataset == "melanoma_time" ~ "SKCM time",
    dataset == "melanoma_responder" ~ "SKCM resp.",
    TRUE ~ dataset
  )) %>%
  mutate(number_of_cell_types = map_int(data, ~ distinct(.x, cell_type) %>% nrow())) %>%
  mutate(outliers = map(
    data,
    ~ .x %>%
      filter(outlier) %>%
      mutate(significant = c_FDR<0.025) %>%
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
  theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust = 1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# # Mean variance association
# job({
#   data_mean_variance_association =
#
#     # Import data multi beta
#     tribble(
#       ~dataset, ~data,
#
#
#       # Comparing Metastasis vs Primary
#       "UVM",
#       readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma_no_heterogeneity.rds"),
#
#       "renal_cell_carcinoma",
#       readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma_based_on_ICB_Response.rds"),
#
#       # Here factor of interest is subtype== triple-negative versus the rest
#       "bc_cells",
#       readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") ,
#
#       "COVID",
#       readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds"),
#
#
#       "melanoma_time",
#       readRDS("dev/data_integration/estimate_GSE120575_melanoma_pre_vs_post_treatment.rds"),
#
#
#       "melanoma_responder",
#       readRDS("dev/data_integration/estimate_GSE120575_melanoma_responders_vs_not.rds") ,
#
#       "BRCA1_breast",
#       readRDS("dev/data_integration/estimate_BRCA1_s41467-021-21783-3.rds")
#     )
# })

# # Plor trends
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
#
# plot_association_all =
#   data_mean_variance_association %>%
#   mutate(data = map(data,  ~ .x %>% select(composition_CI,concentration)  )) %>%
#   unnest(data) %>%
#   unnest(c(composition_CI ,  concentration  )) %>%
#   ggplot(aes(`.median_(Intercept)`, mean)) +
#   geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`, color=dataset),  alpha = 0.4) +
#   geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`, color=dataset), alpha = 0.4) +
#   geom_point(size=0.1) +
#   geom_abline(aes(intercept =  intercept, slope = slope, color = dataset), linetype = "dashed", data = slopes_df) +
#   scale_color_brewer(palette="Set1") +
#   xlab("Category logit mean") +
#   ylab("Category log-concentration") +
#   multipanel_theme


# Analysis of breast
plot_BRCA_UMAP =
  counts %>%
  #with_groups(file, ~ sample_n(.x, 10000)) %>%
  inner_join(
    data_estimates %>% mutate(data = map(data, ~ select(.x, sample, cell_type, label, significant_in_article))) %>% unnest(data)
  ) %>%
  #filter(label %>% is.na %>% `!`) %>%
  unite("color", c(label, significant_in_article)) %>%
  mutate(color = case_when(
    color != "Non-significant_FALSE" ~ color
  )) %>%

  # Filter non intersting datasets
  nest(data = -dataset) %>%
  mutate(n = map_int(data, ~ filter(.x, !is.na(color)) %>% nrow())) %>%
  filter(n > 0 | dataset %>% is.na) %>%
  unnest(data) %>%

  # Recover dataset - this is result of bad programming
  nest(data = -file) %>%
  mutate(data = map(
    data,
    ~ {
      ds = filter(.x, dataset %>% is.na %>% `!`) %>% pull(dataset) %>% unique %>% as.character
      mutate(.x, dataset = ds)
    }
  )) %>%
  unnest(data) %>%
  mutate(color = case_when(
    !(is.na(color) | color=="NA_NA") ~ color
  )) %>%
  mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE"))) %>%


  # Create UMAP
  nest(data = -dataset) %>%
  filter(dataset=="bc_cells") %>%
  unnest(data) %>%

  left_join(
    readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>% unnest(count_data) %>% distinct(sample, type )
  ) %>%

  left_join(
    readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
      distinct(cell_type, `c_FDR` , parameter,  v_FDR ) %>%
      pivot_wider(names_from = parameter, values_from = c(c_FDR ,    v_FDR)  )
  ) %>%


  with_groups(sample, ~ sample_n(.x, 3000, replace = TRUE)) %>%

  # Forget about cancer
  mutate(color = case_when(!grepl("Cancer", cell_type) ~ color )) %>%

  # Rename cell types
  left_join(
    read_csv("dev/brca_cell_type_abbreviations.csv") %>%
      setNames(c("cell_type", "cell_type_pretty"))
  ) %>%

  # Add colouring
  mutate(color = factor(color, levels = c( "Non-significant_TRUE" , "Diff abundant_FALSE" ,  "Diff heterogeneous_FALSE" ,   "Diff both_FALSE",  "Diff abundant_TRUE" ))) %>%
  add_multi_shades(color, cell_type) %>%

  # Plot
 { .x = (.)
   ggplot() +
      geom_point(aes(UMAP_1, UMAP_2), data=filter(.x, is.na(color)), fill="grey", size=0.5, alpha=0.5, shape=21, stroke=0, position = "jitter") +
      geom_point(aes(UMAP_1, UMAP_2), fill=filter(.x, !is.na(color)) %>% pull(.color), data=filter(.x, !is.na(color)), size=0.5, alpha=1, shape=21, stroke=0, position = "jitter") +
      ggrepel::geom_text_repel(
        data= 	get_labels_clusters(
          mutate(.x, cell_type_pretty =
                   case_when(
                     is.na(color) ~"",
                     type == "TNBC" ~"",
                     type == "HER2+" & (`c_FDR_typeHER2+`<0.025 | `v_FDR_typeHER2+`<0.025) ~cell_type_pretty,
                     type == "ER+" & (`c_FDR_typeER+`<0.025 | `v_FDR_typeER+`<0.025) ~cell_type_pretty,
                     TRUE ~""
                   )
                 ),
          c(cell_type_pretty, type),
          UMAP_1,
          UMAP_2
        ) ,
        aes(UMAP_1, UMAP_2, label = cell_type_pretty), size = 2, max.overlaps = 100

      ) +
     facet_wrap(~type) +
      #facet_wrap(~dataset, scale="free") +
      scale_fill_manual(values = color_palette_link) +
      #scale_fill_brewer(palette="Set1") +
      # scale_fill_manual(values = friendly_cols, na.value = "grey") +
      scale_color_manual(values = friendly_cols, na.value = "grey") +
      #guides( fill = "none" ) +
      theme_UMAP +
      theme(title = element_text(size = 7),
            axis.line = element_blank(),
            axis.ticks = element_blank()
      )
 }


#
# data_estimates %>%
#
#   mutate(data = map(
#     data,
#     ~ .x %>%
#       filter(label %>% is.na %>% `!`) %>%
#       unite("color", c(label, significant_in_article)) %>%
#       mutate(color = fct_relevel(color, c("Diff abundant_FALSE", "Diff heterogeneous_FALSE", "Non-significant_TRUE", "Diff abundant_TRUE", "Non-significant_FALSE"))) %>%
#       filter(!color %in% c("Non-significant_FALSE", "Diff abundant_TRUE"))
#
#   )) %>%
#
#   filter(map_int(data, ~ nrow(.x))>0) %>%
#   filter(dataset=="bc_cells") %>%
#   unnest(data) %>%
#   # Create plot
#   {
#       .x = (.)
#       .x$estimate = .x %>% pull(3)
#
#       ggplot() +
#       geom_boxplot(
#         aes(factor_of_interest, proportion, fill=color),
#         outlier.shape = NA,
#         data = .x %>% filter(!outlier), lwd =0.5, fatten = 0.5
#       ) +
#       geom_jitter(aes(factor_of_interest, proportion, color=outlier), size = 0.3, data = .x, height = 0) +
#       facet_wrap(~ fct_reorder(cell_type, abs(estimate), .desc = TRUE), scale="free_y", nrow=1) +
#       scale_y_continuous(trans="logit",labels = dropLeadingZero  ) +
#       scale_color_manual(values = c("black", "#e11f28")) +
#       scale_fill_manual(values = color_palette_link) +
#       #scale_fill_brewer(palette="Set1") +
#       xlab("Biological condition") +
#       ylab("Cell-group proportion (decimal)") +
#       guides(fill="none", color="none") +
#       multipanel_theme+
#       theme(axis.text.x =  element_text(angle=20, hjust = 1))
#     }




BRCA_tnbc_vs_all = readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds")
# BRCA_tnbc_vs_all %>%
#   ggplot(aes(`composition_effect_typeER+`, `variability_effect_typeER+`)) +
#   geom_point()  +
#   geom_smooth(method="lm")
#
# BRCA_tnbc_vs_all %>%
#   ggplot(aes(`composition_effect_typeER+`, `variability_effect_typeHER2+`)) +
#   geom_point() +
#   geom_smooth(method="lm")


# Generated

data_proportion =
  BRCA_tnbc_vs_all |>
  pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
  unnest(count_data) |>
  select(cell_type, sample, outlier, count, type, contains("FDR")) |>
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) %>%


  mutate(significance_c = case_when(
    `c_FDR_typeHER2+` < 0.025 & `c_FDR_typeER+` < 0.025 ~ "her2+ er+",
    `c_FDR_typeHER2+` < 0.025 ~ "her2+",
    `c_FDR_typeER+` < 0.025 ~ "er+",
    TRUE ~ "none"
  )) %>%
  mutate(significance_v = case_when(
    `v_FDR_typeHER2+` < 0.025 & `v_FDR_typeER+` < 0.025 ~ "her2+ er+",
    `v_FDR_typeHER2+` < 0.025 ~ "her2+",
    `v_FDR_typeER+` < 0.025 ~ "er+",
    TRUE ~ "none"
  ))  %>%
  left_join(
    plot_df %>% filter(dataset == "bc_cells") %>% unnest(data) %>% distinct(cell_type,  color )
  )  %>%
  # Forget about cancer
  mutate(color = case_when(!grepl("Cancer", cell_type) ~ color )) %>%

  # Rename cell types
  left_join(
    read_csv("dev/brca_cell_type_abbreviations.csv") %>%
      setNames(c("cell_type", "cell_type_pretty"))
  ) %>%

  mutate(c_FDR = pmin(`c_FDR_typeHER2+`, `c_FDR_typeER+`))

simulated_proportion =
  BRCA_tnbc_vs_all |>
  replicate_data(number_of_draws = 100) |>
  left_join(data_proportion |> distinct(type, sample)) %>%

  # Rename cell types
  left_join(
    read_csv("dev/brca_cell_type_abbreviations.csv") %>%
      setNames(c("cell_type", "cell_type_pretty"))
  )

folded_power = function(x, p=1/3){
  x^p - (1-x)^p
}



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

plot_brca_boxplot =
  ggplot() +
  stat_summary(
    aes(type, (generated_proportions)),
    fun.data = calc_boxplot_stat, geom="boxplot",
    fatten = 0.5, lwd=0.2,
    data =
      simulated_proportion %>%

      # Filter uanitles because of limits
      inner_join( filter(data_proportion, !is.na(color)) %>% distinct(cell_type, c_FDR, color)) ,
    color="blue"

  )+
 # geom_jitter(aes(type, (generated_proportions)), color="blue" ,alpha=0.2, size = 0.6, data = simulated_proportion) +

  geom_boxplot(
    aes(type, (proportion), fill=color),
    outlier.shape = NA,
    fatten = 0.5, lwd=0.2,
    data = filter(data_proportion, !outlier) %>% filter(!is.na(color))
  ) +
  geom_jitter(aes(type, (proportion), color=outlier), size = 0.2, data = data_proportion  %>% filter(!is.na(color)), height = 0) +

  #geom_signif(comparisons = list(c("HER2+", "ER+")), map_signif_level=TRUE) +

  facet_wrap(~ fct_reorder(cell_type_pretty, as.numeric(color)), scale="free_y", nrow=2) +
  #scale_y_continuous(trans="logit") +

  scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +

  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values =  c(
    "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
    "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
    "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
    "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
    "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
  )) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  multipanel_theme +
  theme(strip.background =element_rect(fill="white", colour = NA), legend.position = "bottom")  +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=30, hjust = 1)
  )

# ggplot() +
#
#   geom_boxplot(
#     aes(type, folded_power(generated_proportions)),
#     outlier.shape = NA, alpha=0.2,
#     data = simulated_proportion, color="blue"
#   ) +
#   # geom_jitter(aes(type, folded_power(generated_proportions)), color="blue" ,alpha=0.2, size = 0.6, data = simulated_proportion) +
#
#   geom_boxplot(
#     aes(type, folded_power(proportion), fill=significance_v),
#     outlier.shape = NA,
#     data = filter(data_proportion, !outlier)
#   ) +
#   geom_jitter(aes(type, folded_power(proportion), color=outlier), size = 1, data = data_proportion, height = 0) +
#
#   facet_wrap(~ interaction(cell_type), scale="free_y") +
#   #scale_y_continuous(trans="logit") +
#   scale_color_manual(values = c("black", "#e11f28")) +
#   #scale_fill_manual(values = c("white", "#E2D379")) +
#   xlab("Biological condition") +
#   ylab("Cell-group proportion") +
#   multipanel_theme +
#   theme(strip.background =element_rect(fill="white"), legend.position = "bottom")  +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )

# BRCA_tnbc_vs_all |>
#
#   select(cell_type, composition_CI, contains("composition_effect"), contains("composition_pH0")) %>%
#   setNames(gsub("composition_effect", ".median", colnames(.))) %>%
#   setNames(gsub("composition_pH0", "pH0", colnames(.))) %>%
#   unnest(composition_CI) |>
#   select(-contains("Intercept")) %>%
#   pivot_longer(c(3:10), names_sep ="_", names_to = c("quantile", "subtype"), values_to = c("value")) %>%
#   pivot_wider(names_from = quantile, values_from = value) %>%
#   ggplot(aes(x=`.median`, y=fct_reorder(cell_type, .median), label=cell_type)) +
#   geom_vline(xintercept = 0.2, colour="grey") +
#   geom_vline(xintercept = -0.2, colour="grey") +
#   geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=pH0<0.025)) +
#   geom_point() +
#   facet_wrap(~subtype) +
#   scale_color_brewer(palette = "Set1") +
#   xlab("Credible interval of the slope") +
#   ylab("Cell group") +
#   theme_bw() +
#   theme(legend.position = "bottom")
#
# BRCA_tnbc_vs_all |>
#   select(cell_type, variability_CI, contains("variability_effect"), contains("variability_pH0")) %>%
#   setNames(gsub("variability_effect", ".median", colnames(.))) %>%
#   setNames(gsub("variability_pH0", "pH0", colnames(.))) %>%
#   unnest(variability_CI) |>
#   select(-contains("Intercept")) %>%
#   pivot_longer(c(3:10), names_sep ="_", names_to = c("quantile", "subtype"), values_to = c("value")) %>%
#   pivot_wider(names_from = quantile, values_from = value) %>%
#   ggplot(aes(x=`.median`, y=fct_reorder(cell_type, .median))) +
#   geom_vline(xintercept = 0.2, colour="grey") +
#   geom_vline(xintercept = -0.2, colour="grey") +
#   geom_errorbar(aes(xmin=`.lower`, xmax=`.upper`, color=pH0<0.025)) +
#   geom_point() +
#   facet_wrap(~subtype) +
#   scale_color_brewer(palette = "Set1") +
#   xlab("Credible interval of the slope") +
#   ylab("Cell group") +
#   theme_bw() +
#   theme(legend.position = "bottom")

plot_BRCA_summary_er =
  BRCA_tnbc_vs_all %>%
  filter(parameter == "typeER+") %>%
  ggplot(aes(c_effect, -v_effect, label=cell_type)) +
  geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_errorbar(aes(xmin=c_lower, xmax=c_upper, color=c_FDR<0.025, alpha=c_FDR<0.025), size=0.2) +
  geom_errorbar(aes(ymin=-v_lower, ymax=-v_upper, color=v_FDR<0.025, alpha=v_FDR<0.025), size=0.2) +

  geom_point(size=0.2)  +
  annotate("text", x = 0, y = 3.5, label = "Variable", size=2) +
  annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
  scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
  scale_alpha_manual(values = c(0.4, 1)) +
  #ggtitle("ER+") +
  multipanel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_BRCA_summary_her2 =
  BRCA_tnbc_vs_all %>%
  filter(parameter == "typeHER2+") %>%
  ggplot(aes(c_effect, -v_effect, label=cell_type)) +
  geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
  geom_errorbar(aes(xmin=c_lower, xmax=c_upper, color=c_FDR<0.025, alpha=c_FDR<0.025), size=0.2) +
  geom_errorbar(aes(ymin=-v_lower, ymax=-v_upper, color=v_FDR<0.025, alpha=v_FDR<0.025), size=0.2) +

  geom_point(size=0.2)  +
  annotate("text", x = 0, y = 3.5, label = "Variable", size=2) +
  annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
  scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
  scale_alpha_manual(values = c(0.4, 1)) +
  #ggtitle("ER+") +
  multipanel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


 # bc_subtypes =
#   readRDS("dev/data_integration/SCP1039_bc_cells.rds")  |>
#   #mutate(type = subtype=="TNBC") %>%
#   # mutate(type = factor(subtype, levels = c("TNBC", "HER2+", "ER+"))) %>%
#   mutate(cell_type = celltype_subset) %>%
#   nest(data = -subtype) %>%
#   mutate(result = map(
#     data,
#     ~ sccomp_glm(
#       .x,
#       formula = ~ 1,
#       sample, cell_type,
#       approximate_posterior_inference = FALSE,
#       variance_association = FALSE,
#       prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(5.06983, 8.549324))
#     )
#   ))
#
# bc_subtypes %>% saveRDS("dev/bc_subtypes.rds")
bc_subtypes = readRDS("dev/bc_subtypes.rds")

fair_cols <- c("#BF1B0B", "#66ADE5", "#378805")
library(shades)

plot_type_associations_lines =
  bc_subtypes %>%
    mutate(mean_concentration_association = map(result, ~ .x %>% attr("mean_concentration_association") %>% t %>% as.data.frame())) %>%
    unnest(mean_concentration_association) %>%
    unnest(result) %>%
    select(-data) %>%
    unnest(composition_CI) %>%
    unnest(concentration) %>%
    ggplot(aes(`.median_(Intercept)`, mean, color=subtype)) +
    geom_point(size=0.2) +
    geom_abline(
      mapping = aes(intercept = `prec_coeff[1]`, slope = `prec_coeff[2]`, color=subtype),
      linetype = "dashed", size=0.3
    ) +
    scale_color_manual(values = cool_palette) +
    xlab("Multinomial-logit mean") +
    ylab("Log concentration") +
    multipanel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_type_associations_densities =
  bc_subtypes %>%
  mutate(intercept_assoc = map(
    result,
    ~ attr(.x, "fit") %>%
      tidybayes::spread_draws(prec_coeff[A])
  )) %>%
  select(1, 4) %>%
  unnest(intercept_assoc) %>%
  ggplot(aes(prec_coeff,  fill=subtype)) +
  ggdist:: stat_halfeye(alpha = 0.6, size=0.3) +
  #geom_density()+
  facet_wrap(~A, scale="free_x") +
  scale_fill_manual(values = cool_palette) +
  xlab("Mean-concentration association") +
  ylab("density") +
  guides(fill="none") +
  multipanel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


p=plot_BRCA_UMAP /
  (plot_BRCA_summary_er | plot_BRCA_summary_her2 | plot_type_associations_lines | plot_type_associations_densities) /
  plot_brca_boxplot /
  ( plot_df$UMAP_plot[[1]] + plot_df$estimate_plot[[1]]  +  plot_layout(widths = c(1, 3.5)) ) /
 # (  plot_df$UMAP_plot[[2]] + plot_df$estimate_plot[[2]] +  plot_df$UMAP_plot[[3]] + plot_df$estimate_plot[[3]] + plot_layout(widths = c(1, 0.5, 1, 2)) ) /
  # ( plot_df$UMAP_plot[[4]] + plot_df$estimate_plot[[4]] +
      (plot_novel_results + plot_outliers +  plot_spacer() + plot_layout(widths = c(1,1,2))) +
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), plot.tag = element_text(size = 10))



ggsave(
  "dev/article_figures/novel_results_plot.png",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 230 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/novel_results_plot.pdf",
  plot = p,
  units = c("mm"),
  width = 183 ,
  height = 220 ,
  limitsize = FALSE
)

# Only Cancer for supplementary
data_proportion_cancer =
  BRCA_tnbc_vs_all |>
  pivot_wider(names_from = parameter, values_from = c(contains("c_"), contains("v_"))) %>%
  unnest(count_data) |>
  select(cell_type, sample, outlier, count, type, contains("FDR")) |>
  with_groups(sample, ~ mutate(.x, proportion = (count)/sum(count)) ) %>%


  mutate(significance_c = case_when(
    `c_FDR_typeHER2+` < 0.025 & `c_FDR_typeER+` < 0.025 ~ "her2+ er+",
    `c_FDR_typeHER2+` < 0.025 ~ "her2+",
    `c_FDR_typeER+` < 0.025 ~ "er+",
    TRUE ~ "none"
  )) %>%
  mutate(significance_v = case_when(
    `v_FDR_typeHER2+` < 0.025 & `v_FDR_typeER+` < 0.025 ~ "her2+ er+",
    `v_FDR_typeHER2+` < 0.025 ~ "her2+",
    `v_FDR_typeER+` < 0.025 ~ "er+",
    TRUE ~ "none"
  ))  %>%
  left_join(
    plot_df %>% filter(dataset == "bc_cells") %>% unnest(data) %>% distinct(cell_type,  color )
  )  %>%
  # Only Cncer for supplementary
  filter( grepl("Cancer", cell_type) ) %>%

  # Rename cell types
  left_join(
    read_csv("dev/brca_cell_type_abbreviations.csv") %>%
      setNames(c("cell_type", "cell_type_pretty"))
  ) %>%

  mutate(c_FDR = pmin(`c_FDR_typeHER2+`, `c_FDR_typeER+`)) %>%
  mutate(color = "Diff abundant_TRUE")

plot_brca_boxplot_cancer =
  ggplot() +
  stat_summary(
    aes(type, (generated_proportions)),
    fun.data = calc_boxplot_stat, geom="boxplot",
    fatten = 0.5, lwd=0.2,
    data =
      simulated_proportion %>%

      # Filter uanitles because of limits
      inner_join( filter(data_proportion_cancer, !is.na(color)) %>% distinct(cell_type, c_FDR, color)) ,
    color="blue"

  )+

  geom_boxplot(
    aes(type, (proportion), fill=color),
    outlier.shape = NA,
    fatten = 0.5, lwd=0.2,
    data = filter(data_proportion_cancer, !outlier) %>% filter(!is.na(color))
  ) +
  geom_jitter(aes(type, (proportion), color=outlier), size = 0.2, data = data_proportion_cancer  %>% filter(!is.na(color)), height = 0) +


  facet_wrap(~ fct_reorder(cell_type_pretty, as.numeric(color)), scale="free_y", nrow=1) +

  scale_y_continuous(trans="S_sqrt", labels = dropLeadingZero) +

  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values =  c(
    "Non-significant_TRUE" = brewer_pal(palette="Reds")(9)[7],
    "Diff abundant_FALSE" = brewer_pal(palette="Greens")(9)[7],
    "Diff heterogeneous_FALSE"= brewer_pal(palette="Purples")(9)[7] ,
    "Diff both_FALSE"= brewer_pal(palette="Blues")(9)[7],
    "Diff abundant_TRUE" = brewer_pal(palette="Greys")(9)[7]
  )) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  multipanel_theme +
  theme(strip.background =element_rect(fill="white", colour = NA), legend.position = "bottom")  +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=30, hjust = 1)
  )


p_other =
(  plot_brca_boxplot_cancer /
(
  (plot_df$UMAP_plot[[2]] + ggtitle(plot_df$dataset[[2]])) + plot_df$estimate_plot[[2]] +
    plot_layout(widths = c(1, 1))
) /
  (
    (plot_df$UMAP_plot[[3]] + ggtitle(plot_df$dataset[[3]])) + plot_df$estimate_plot[[3]] +
      plot_layout(widths = c(1, 1))
  ) /
(
   (plot_df$UMAP_plot[[5]] + ggtitle(plot_df$dataset[[5]])) + plot_df$estimate_plot[[5]]  +
    plot_layout(widths = c(1,1))
)
)+
  plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), plot.tag = element_text(size = 10))


ggsave(
  "dev/article_figures/novel_results_other_plot.png",
  plot = p_other,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/novel_results_other_plot.pdf",
  plot = p_other,
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)



# # Other mean trend variability trend association
# data_estimates %>%
#   slice(4) %>%
#   pull(data) %>%
#   .[[1]] %>%
#   unnest(c(composition_CI, variability_CI), names_repair = "unique") %>%
#   ggplot(aes(`composition_effect_is_criticalTRUE`, `variability_effect_is_criticalTRUE`, label=cell_type)) +
#   geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   geom_errorbar(aes(xmin=`.lower_is_criticalTRUE...5`, xmax=`.upper_is_criticalTRUE...8`, color=`composition_FDR_is_criticalTRUE`<0.025, alpha=`composition_FDR_is_criticalTRUE`<0.025), size=0.2) +
#   geom_errorbar(aes(ymin=`.lower_is_criticalTRUE...14`, ymax=`.upper_is_criticalTRUE...17`, color=`variability_FDR_is_criticalTRUE`<0.025, alpha=`variability_FDR_is_criticalTRUE`<0.025), size=0.2) +
#
#   geom_point(aes(alpha=`c_FDR_is_criticalTRUE`<0.025), size=0.2)  +
#   annotate("text", x = 0, y = -3.5, label = "Variable", size=2) +
#   annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
#   scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
#   scale_alpha_manual(values = c(0.4, 1)) +
#   #ggtitle("ER+") +
#   multipanel_theme +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )
#
# data_estimates %>%
#   slice(5) %>%
#   pull(data) %>%
#   .[[1]] %>%
#   unnest(c(composition_CI, variability_CI), names_repair = "unique") %>%
#   ggplot(aes(`composition_effect_timePre`, `variability_effect_timePre`, label=cell_type)) +
#   geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   # geom_errorbar(aes(xmin=`.lower_is_criticalTRUE...5`, xmax=`.upper_typeER+...11`, color=`c_FDR_is_criticalTRUE`<0.025, alpha=`composition_FDR_typeER+`<0.025), size=0.2) +
#   # geom_errorbar(aes(ymin=`.lower_typeER+...21`, ymax=`.upper_typeER+...25`, color=`v_FDR_is_criticalTRUE`<0.025, alpha=`variability_FDR_typeER+`<0.025), size=0.2) +
#   #
#   geom_point(aes(alpha=`composition_effect_timePre`<0.025), size=0.2)  +
#   annotate("text", x = 0, y = -3.5, label = "Variable", size=2) +
#   annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
#   scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
#   scale_alpha_manual(values = c(0.4, 1)) +
#   #ggtitle("ER+") +
#   multipanel_theme +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )
#
#
# data_estimates %>%
#   slice(6) %>%
#   pull(data) %>%
#   .[[1]] %>%
#   unnest(c(composition_CI, variability_CI), names_repair = "unique") %>%
#   ggplot(aes(`composition_effect_resResponder`, `variability_effect_resResponder`, label=cell_type)) +
#   geom_vline(xintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   geom_hline(yintercept = c(-0.2, 0.2), colour="grey", linetype="dashed", size=0.3) +
#   # geom_errorbar(aes(xmin=`.lower_is_criticalTRUE...5`, xmax=`.upper_typeER+...11`, color=`c_FDR_is_criticalTRUE`<0.025, alpha=`composition_FDR_typeER+`<0.025), size=0.2) +
#   # geom_errorbar(aes(ymin=`.lower_typeER+...21`, ymax=`.upper_typeER+...25`, color=`v_FDR_is_criticalTRUE`<0.025, alpha=`variability_FDR_typeER+`<0.025), size=0.2) +
#   #
#   geom_point(aes(alpha=`composition_effect_resResponder`<0.025), size=0.2)  +
#   annotate("text", x = 0, y = -3.5, label = "Variable", size=2) +
#   annotate("text", x = 5, y = 0, label = "Abundant", size=2, angle=270) +
#   scale_color_manual(values = c("#D3D3D3", "#E41A1C")) +
#   scale_alpha_manual(values = c(0.4, 1)) +
#   #ggtitle("ER+") +
#   multipanel_theme +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )
