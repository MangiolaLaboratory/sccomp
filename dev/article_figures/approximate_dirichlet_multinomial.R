library(tidyverse)
library(extraDistr)
library(sccomp)
library(job)
library(tidybayes)
library(glue)
library(patchwork)
library(rstan)
library(GGally)

# Load multipanel_theme
source("https://gist.githubusercontent.com/stemangiola/fc67b08101df7d550683a5100106561c/raw/305ed9ba2af815fdef3214b9e6171008d5917400/ggplot_theme_multipanel")

set.seed(42)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

dirichlet_multinomial_counts_1 =
  rdirmnom(1000, 1000, c(1,3,10,20)/5) %>%
  as_tibble(rownames = "sample") %>%
  gather(category, count, -sample) %>%
  mutate(count = as.integer(count))



distro_plot_1 =
  dirichlet_multinomial_counts_1 %>%
  ggplot(aes(category, count)) +
  geom_jitter(height = 0,  size=0.01, alpha=0.5) +
  multipanel_theme

# Fit
# job({
estimate_1 =
  dirichlet_multinomial_counts_1 %>%
  sccomp_glm(
    ~1, sample, category, count,
    check_outliers = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(5,5)),
    verbose=TRUE
  )
# })

#cool_palette = c("#b58b4c", "#74a6aa", "#a15259",  "#37666a", "#79477c", "#cb9f93", "#9bd18e", "#eece97", "#8f7b63", "#4c474b", "#415346")
cool_palette = c("#e29a3b", "#3a9997", "#a15259",  "#cc374a", "#94379b", "#ff8d73", "#7be05f", "#ffc571", "#724c24", "#4c474b", "#236633")


estimate_1 %>% saveRDS("dev/article_figures/estimate_1_for_dirichlet_approximation.rds")
library(GGally)
bind_rows(
  readRDS("dev/article_figures/estimate_1_for_dirichlet_approximation.rds")  %>%
    attr("fit") %>%
    tidybayes::gather_draws(beta[C,M]) %>%
    mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>%
    mutate(class="dependent") %>%
    ungroup() %>%
    filter(.draw<=4000) %>%
    filter(M < 4) ,
  readRDS("dev/article_figures/estimate_1_fit_for_dirichlet_approximation_with_just_beta_binomial.rds")  %>%
    tidybayes::gather_draws(beta[C,M]) %>%
    mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>%
    ungroup() %>%
    filter(.draw<=4000) %>%
    mutate(class="independent")
) %>%

  bind_rows(
    fit %>% gather_draws(beta[C, M]) %>% filter(.chain!=1 & M != 4) %>%
      mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>%
      filter(.draw<=4000) %>%
      mutate(class="Dirichlet-mult")
  ) %>%

  spread(M, .value) %>%
  GGally::ggpairs(
    columns = 7:10, aes(color = class, fill=class),
    lower = list(continuous = function(data, mapping, method="lm", ...){
      p <- ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.5,    size=0.1, shape=21, stroke=0) +
        geom_smooth(method=method, ...)
      p
    }),
    diag = list(continuous = function(data, mapping, ...){
      p <- ggplot(data = data, mapping = mapping) +
        geom_density(alpha = 0, size=0.5)
      p
    }),
    upper = list(continuous = wrap("cor", size = 2))
    ) +
  scale_color_manual(values =c("dependent" = cool_palette[3], "independent" = cool_palette[2], "Dirichlet-mult" =cool_palette[1])) +
  scale_fill_manual(values =c("dependent" = cool_palette[3], "independent" = cool_palette[2], "Dirichlet-mult" =cool_palette[1])) +
  multipanel_theme +
  theme(text = element_text(size=8)) +
  theme(  strip.background =element_rect(fill="white", color="white")   ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

plot_estimate =
  estimate_1 %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`), color="#cc6666", alpha = 0.8, width=0) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), color="#cc6666", alpha = 0.8, width=0) +
  geom_point(size=0.2) +
  geom_abline(intercept = 5.7496330, slope = -0.9650953, linetype = "dashed", color="grey") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  multipanel_theme

# Simulate data
fit_object_1 = attr(estimate_1, "fit")
generated_1 =  estimate_1 %>% replicate_data()


both_data_sets_1 =
  estimate_1  %>%
  select(category, count_data ) %>%
  unnest(count_data ) %>%
  mutate(distribution = "Simplex beta binomial") %>%
  bind_rows(
    generated_1 %>%
      rename(count = generated_counts ) %>%
      mutate(distribution = "Dirichlet-multinomial")
  )



# Plot

distro_plot_2 =
  both_data_sets_1 %>%
  ggplot(aes(category, count, color = distribution,  size=distribution)) +
  geom_jitter(height = 0,  alpha=0.5) +
  scale_color_manual(values = c("black", "#cc6666")) +
  scale_size_discrete(range = c(0.05, 0.5)) +
  guides(color="none", fill="none", size="none") +
  multipanel_theme

distro_violin_1 =
  both_data_sets_1 %>%
  ggplot(aes(category, count + 1, fill = distribution)) +
  geom_split_violin(scale = "width") +
  scale_fill_manual(values = c("white", "#cc6666")) +
  #scale_size_discrete(range = c(0.1, 1)) +
  scale_y_log10() +
  guides(color="none", fill="none", size="none") +
  multipanel_theme


# Dirichlet multinomial upper bound
dirichlet_multinomial_counts_2 =
  rdirmnom(1000, 1000, c(1,3,10,50)/20) %>%
  as_tibble(rownames = "sample") %>%
  gather(category, count, -sample) %>%
  mutate(count = as.integer(count))


# distro_plot_3 =
#   dirichlet_multinomial_counts_2 %>%
#   ggplot(aes(category, count)) +
#   geom_jitter(height = 0,  size=0.01, alpha=0.2) +
#   multipanel_theme

# Fit
# job({
estimate_2 =
  dirichlet_multinomial_counts_2 %>%
  sccomp_glm(
    ~1, sample, category, count,
    check_outliers = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(5,5)),
    verbose=TRUE
  )
# })

estimate_2 %>% saveRDS("dev/article_figures/estimate_2_for_dirichlet_approximation.rds")

bind_rows(
  estimate_2 %>% attr("fit") %>% gather_draws(beta[C,M]) %>% mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>% mutate(class="dependent"),
  readRDS("dev/article_figures/estimate_2_fit_for_dirichlet_approximation_with_just_beta_binomial.rds")  %>%
    gather_draws(beta[C,M]) %>% mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>% mutate(class="independent")
) %>%
  ungroup() %>%
  ggplot(aes(.value, color=class)) +
  geom_density() +
  facet_wrap(~M)

readRDS("dev/article_figures/estimate_2_for_dirichlet_approximation.rds") %>% attr("fit") %>% pairs(pars="beta")
readRDS("dev/article_figures/estimate_2_fit_for_dirichlet_approximation_with_just_beta_binomial.rds") %>% pairs(pars="beta")

# Simulate data
fit_object_2 = attr(estimate_2, "fit")
generated_2 =  estimate_2 %>% replicate_data()

bind_rows(
readRDS("dev/article_figures/estimate_2_for_dirichlet_approximation.rds")  %>% attr("fit") %>% tidybayes::gather_draws(beta[C,M]) %>% mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>% mutate(class="dependent"),
readRDS("dev/article_figures/estimate_2_fit_for_dirichlet_approximation_with_just_beta_binomial.rds")  %>%
  tidybayes::gather_draws(beta[C,M]) %>% mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>% mutate(class="independent")
) %>%
  ungroup() %>%
  ggplot(aes(.value, color=class)) +
  geom_density() +
  facet_wrap(~M)

bind_rows(
  readRDS("dev/article_figures/estimate_2_for_dirichlet_approximation.rds")  %>%
    attr("fit") %>%
    tidybayes::gather_draws(beta[C,M]) %>%
    mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>%
    mutate(class="dependent"),
  readRDS("dev/article_figures/estimate_2_fit_for_dirichlet_approximation_with_just_beta_binomial.rds")  %>%
    tidybayes::gather_draws(beta[C,M]) %>% mutate(.value = scale(.value, scale=F) %>% as.numeric()) %>% mutate(class="independent")
) %>%
  ungroup() %>%
  spread(M, .value) %>%
  GGally::ggpairs(columns = 7:10, aes(color = class))

readRDS("dev/article_figures/estimate_2_for_dirichlet_approximation.rds")  %>% attr("fit") %>% pairs(pars="beta")
readRDS("dev/article_figures/estimate_2_fit_for_dirichlet_approximation_with_just_beta_binomial.rds") %>% pairs(pars="beta")



both_data_sets_2 =
  estimate_2  %>%
  select(category, count_data ) %>%
  unnest(count_data ) %>%
  mutate(distribution = "Simplex beta binomial") %>%
  bind_rows(
    generated_2 %>%
      rename(count = generated_counts ) %>%
      mutate(distribution = "Dirichlet-multinomial")
  )

# Plot
distro_plot_4 =
  both_data_sets_2 %>%
  ggplot(aes(category, count, color = distribution,  size=distribution)) +
  geom_jitter(height = 0, alpha=0.5) +
  scale_color_manual(values = c("black", "#cc6666")) +
  scale_size_discrete(range = c(0.05, 0.5)) +
  guides(color="none", fill="none", size="none") +
  multipanel_theme +
  theme(legend.position = "none")

distro_violin_2 =
  both_data_sets_2 %>%
  ggplot(aes(category, count + 1, fill = distribution)) +
  geom_split_violin(scale = "width") +
  scale_fill_manual(values = c("white", "#cc6666")) +
  #scale_size_discrete(range = c(0.1, 1)) +
  scale_y_log10() +
  multipanel_theme




plot_pairs =
  both_data_sets_1 %>%

  # Plot
  select(sample, category, count, distribution) %>%
  mutate(distribution = if_else(distribution == "Simplex beta binomial", "sbb", "dm")) %>%
  spread(category, count ) %>%
  ggpairs(
    3:6,
    mapping=ggplot2::aes(colour = distribution, fill=distribution) ,
    lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1)) ,
    diag = list(discrete="barDiag", continuous = wrap("densityDiag", alpha=0 )),
    upper = list(continuous = wrap("cor", size = 2))
  ) +
  scale_color_manual(values = c("black", "#cc6666")) +
  multipanel_theme +
  theme(text = element_text(size=8)) +
  theme(  strip.background =element_rect(fill="white", color="white")   ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

saveRDS(estimate_1, "dev/article_figures/estimate_1.rds")
saveRDS(estimate_2, "dev/article_figures/estimate_2.rds")

# plot_pairs_dirichlet_multinomial =
#   dirichlet_multinomial_counts %>%
#   spread(category, count) %>%
#   ggpairs(2:5,  lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1))   ) +
#   multipanel_theme +
#   theme(text = element_text(size=8)) +
#   theme(  strip.background =element_rect(fill="white", color="white")   ) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   ggtitle("Dirichlet-multinomial")
#



p1 =
  (
    (distro_plot_1 + plot_estimate + distro_plot_2 + distro_violin_1 ) |
      wrap_elements(ggmatrix_gtable(plot_pairs)) |
      (distro_plot_4 / distro_violin_2)
  ) +
  # Style
  plot_layout(guides = 'collect', width  = c(2, 4, 1)) + plot_annotation(tag_levels = c('A')) &
  theme( plot.margin = margin(0, 0, 0, 0, "pt"), legend.position = "bottom")


ggsave(
  "dev/article_figures/approximation_dirichlet_multinomial.pdf",
  plot = p1,
  units = c("mm"),
  width = 183 ,
  height = 120 ,
  limitsize = FALSE
)

ggsave(
  "dev/article_figures/approximation_dirichlet_multinomial.png",
  plot = p1,
  units = c("mm"),
  width = 183 ,
  height = 120 ,
  limitsize = FALSE
)

