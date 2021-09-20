library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)



estimate_oligo = readRDS("dev/data_integration/estimate_GSE115189_SCP345_SCP424_SCP591_SRR11038995_SRR7244582_10x6K_10x8K.rds")
estimate_oligo %>% attr("mean_concentration_association")
# [1]  5.6260004 -0.6940178
# prec_sd  = 0.3312485

estimate_UVM = readRDS("dev/data_integration/estimate_GSE139829_uveal_melanoma.rds")
estimate_UVM %>% attr("mean_concentration_association")
# [1]  3.0423138 -0.6920534
# prec_sd  = 0.1987102

estimate_renal_cell_carcinoma = readRDS("dev/data_integration/estimate_SCP1288_renal_cell_carcinoma.rds")
estimate_renal_cell_carcinoma %>% attr("mean_concentration_association")
# [1]  4.1541595    -0.7367941
# prec_sd  =  0.5060364

estimate_bc_cells =
  readRDS("dev/data_integration/estimate_SCP1039_bc_cells.rds") %>%
  mutate(count_data  = map(count_data , ~mutate(.x, type = as.character(type))))
estimate_bc_cells %>% attr("mean_concentration_association")
# [1]  3.2800052 -0.7575131
# prec_sd  = 1.0363162

estimate_COVID = readRDS("dev/data_integration/estimate_s41587-020-0602-4_COVID_19.rds")
estimate_COVID %>% attr("mean_concentration_association")
# [1]  3.8746815    -0.8179763
# prec_sd  = 0.6306725

estimate_melanoma = readRDS("dev/data_integration/estimate_GSE120575_melanoma.rds")
estimate_melanoma %>% attr("mean_concentration_association")
# [1]  2.2340179    -0.7466973
# prec_sd  = 0.6198534


# Plor trends
size_df =

    list(
      estimate_oligo = estimate_oligo %>% nrow(),
      estimate_UVM = estimate_UVM  %>% nrow(),
      estimate_renal_cell_carcinoma = estimate_renal_cell_carcinoma  %>% nrow(),
      estimate_bc_cells = estimate_bc_cells  %>% nrow(),
      estimate_COVID = estimate_COVID  %>% nrow(),
      estimate_melanoma = estimate_melanoma  %>% nrow()
    ) %>%
      enframe("dataset", "number_of_cell_types") %>%
      unnest(number_of_cell_types)

plot_outlier =
  list(
  estimate_oligo %>% mutate(dataset = "estimate_oligo"),
  estimate_UVM %>% mutate(dataset = "estimate_UVM"),
  estimate_renal_cell_carcinoma %>% mutate(dataset = "estimate_renal_cell_carcinoma"),
  estimate_bc_cells %>% mutate(dataset = "estimate_bc_cells"),
  estimate_COVID %>% mutate(dataset = "estimate_COVID"),
  estimate_melanoma %>% mutate(dataset = "estimate_melanoma")
) %>%
  reduce(bind_rows) %>%
  unnest(count_data ) %>%
  filter(outlier) %>%
  count(dataset, significant) %>%
  with_groups(dataset, ~ mutate(.x, nn = sum(n))) %>%
  mutate(dataset = gsub("estimate_", "", dataset)) %>%

  ggplot() +
  geom_bar(aes(forcats::fct_reorder(dataset, nn, .desc = TRUE ),  n, fill=significant  ), stat = "identity", position = "stack")+
  geom_point(aes(dataset, number_of_cell_types), data = size_df %>% mutate(dataset = gsub("estimate_", "", dataset))  ) +
  scale_fill_manual(values = c("grey", "#e11f28")) +
  xlab("Dataset") +
  ylab("Count of outliers") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle=30, hjust = 1))


ggsave(
  "dev/outliers_across_datasets_WEHI_seminar_2.pdf",
  plot = plot_outlier,
  units = c("mm"),
  width = 80 ,
  height = 183 ,
  limitsize = FALSE
)
