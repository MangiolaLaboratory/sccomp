library(tidyverse)
library(tidyseurat)
library(sccomp)
library(job)
library(patchwork)
library(forcats)

debugonce(sccomp:::estimate_multi_beta_binomial_glm)

counts_obj =
  readRDS("dev/data_integration/s41587-020-0602-4_COVID_19.rds")  |>
  mutate(is_critical = severity=="critical")

job({
  estimate_COVID =
    counts_obj %>%
    sccomp:::sccomp_glm(
      formula = ~ is_critical,
      sample, cell_type,
      approximate_posterior_inference = F,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      variance_association = TRUE
    )
})


estimate_COVID %>%
  filter(significant)

data_for_plot =
  estimate_COVID %>%
  tidyr::unnest(outliers) %>%
  left_join(counts_obj %>% distinct(sample, is_critical), by = c("sample")) %>%
  group_by(sample) %>%
  mutate(proportion = (count+1)/sum(count+1)) %>%
  ungroup(sample) %>%
  mutate(variance = row_number() %in% c(2,11,24,14))

ggplot() +
  geom_violin(
    aes(is_critical, proportion, fill=significant),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(is_critical, proportion, color=outlier), size = 1, data = data_for_plot, height = 0) +
  facet_wrap(~ fct_reorder(cell_type, abs(.median_is_criticalTRUE), .desc = TRUE), scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

ggsave(
  "dev/variance_analysis_WEHI_talk.pdf",
  units = c("mm"),
  width = 183 ,
  height = 183 ,
  limitsize = FALSE
)


data_for_plot_2 =
  counts_obj %>%
  count(sample, is_critical, cell_type, name = "count") %>%
  group_by(sample) %>%
  mutate(proportion = (count+1)/sum(count+1)) %>%
  ungroup(sample) %>%
  mutate(outlier = cell_type == "pDC" & count <=10)

  ggplot() +
  geom_boxplot(
    aes(is_critical, proportion, fill=cell_type=="pDC"),
    outlier.shape = NA,
    data = data_for_plot_2
  ) +
  geom_jitter(aes(is_critical, proportion, color=outlier), size = 1, data = data_for_plot_2) +
  facet_wrap(~ cell_type, scale="free_y") +
  scale_y_continuous(trans="logit") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))


  data_for_plot_3 =
    counts_obj %>%
    count(sample, is_critical, cell_type, name = "count") %>%
    mutate(count = if_else(cell_type == "pDC", as.double(0), as.double(count))) %>%
    group_by(sample) %>%
    mutate(proportion = (count+1)/sum(count+1)) %>%
    ungroup(sample) %>%
    mutate(outlier = cell_type == "pDC" & count <=10)

  ggplot() +
    geom_boxplot(
      aes(is_critical, proportion, fill=cell_type=="pDC"),
      outlier.shape = NA,
      data = data_for_plot_3
    ) +
    geom_jitter(aes(is_critical, proportion, color=outlier), size = 1, data = data_for_plot_3) +
    facet_wrap(~ cell_type, scale="free_y") +
    scale_y_continuous(trans="logit") +
    scale_color_manual(values = c("black", "#e11f28")) +
    scale_fill_manual(values = c("white", "#E2D379")) +
    xlab("Biological condition") +
    ylab("Cell-group proportion") +
    theme_bw() +
    theme(strip.background =element_rect(fill="white"))
