
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(patchwork)
library(magrittr)
library(scales)

data("seurat_obj")
data("sce_obj")
data("counts_obj")

logit_pretty_trans <- function(){


  if (find.package("functional", quiet = TRUE) %>% length %>% equals(0)) {
    message("Installing functional needed for analyses")
    install.packages("functional", repos = "https://cloud.r-project.org")
  }

  trans <- qlogis
  inv <- plogis

  trans_new("logit",
            transform = trans,
            inverse = inv,
            breaks = functional::Compose(trans, extended_breaks(), inv),
            format = label_scientific(digits = 2)
  )
}


res =
  counts_obj %>%
  sccomp_glm( ~ type, sample, cell_group, count)


out = ggplot() +
  geom_boxplot(
    aes(type, proportion, fill=significant),
    outlier.shape = NA,
    data = data_for_plot %>% filter(!outlier)
  ) +
  geom_jitter(aes(type, proportion, color=outlier), size = 0.6, data = data_for_plot) +
  facet_wrap(~ interaction(cell_group), scale="free_y") +
  scale_y_continuous(trans="logit_pretty") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_fill_manual(values = c("white", "#E2D379")) +
  xlab("Biological condition") +
  ylab("Cell-group proportion") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white")) +
  theme(legend.position = "bottom", strip.background = element_blank()) +
  theme(axis.text.x = element_text(angle=30, vjust = 1, hjust = 1))

ci = res %>%
  ggplot(aes(x=`.median_typecancer`, y=fct_reorder(cell_group, .median_typecancer))) +
  geom_vline(xintercept = 0, colour="grey") +
  geom_errorbar(aes(xmin=`.lower_typecancer`, xmax=`.upper_typecancer`, color=significant)) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  xlab("Credible interval of the slope") +
  ylab("Cell group") +
  theme(legend.position = "bottom")

out + ci + plot_layout( guides = "collect") +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = 'bottom', plot.margin = margin(0, 0, 0, 0, "pt"), text = element_text(size=8))

ggsave(
  "dev/figure_draft_article.png",
  units = c("mm"),
  width = 183 ,
  height = 110 ,
  limitsize = FALSE
)
