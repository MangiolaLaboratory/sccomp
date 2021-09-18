library(extraDistr)
library(sccomp)
library(job)
library(tidybayes)
library(glue)
library(patchwork)
library(rstan)
library(GGally)

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

dirichlet_multinomial_counts =
  rdirmnom(1000, 1000, c(1,3,10,20)/5) %>%
  as_tibble(rownames = "sample") %>%
  gather(category, count, -sample) %>%
  mutate(count = as.integer(count))



distro_plot_1 =
  dirichlet_multinomial_counts %>%
  ggplot(aes(category, count)) +
  geom_jitter(height = 0, size = 0.1) +
  theme_bw()

# Fit
job({
  estimate =
  dirichlet_multinomial_counts %>%
  sccomp_glm(
    ~1, sample, category, count,
    check_outliers = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
    verbose=TRUE
  )
})

plot_estimate =
  estimate %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`), color="#4DAF4A", alpha = 0.8, width=0) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), color="#4DAF4A", alpha = 0.8, width=0) +
  geom_point() +
  geom_abline(intercept = 5.7496330, slope = -0.9650953, linetype = "dashed", color="grey") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw()

# Simulate data
job({
  fit_object =
  dirichlet_multinomial_counts %>%
  sccomp:::estimate_multi_beta_binomial_glm(
    ~1, sample, category, count,
    check_outliers = FALSE,
    approximate_posterior_inference = FALSE,
    prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
    verbose=TRUE
  )
})


# rng =  rstan::gqs(
#   sccomp:::stanmodels$glm_multi_beta_binomial_generate_date,
#   #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
#   draws =  as.matrix(fit_object$fit),
#   data = fit_object$data_for_model
# )

both_data_sets =estimate  %>%
  select(category, generated_data ) %>%
  unnest(generated_data ) %>%
  distinct() %>%
  rename(count = generated_counts ) %>%
  mutate(distribution = "Simplex beta binomial") %>%
  bind_rows(
    dirichlet_multinomial_counts %>%
      mutate(distribution = "Dirichlet-multinomial")
  )


# both_data_sets =
#   sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
#   filter(.draw == 1) %>%
#   rename(sample = N, count = .value) %>%
#   mutate(category = glue("V{M}"), sample = as.character(sample)) %>%
#   mutate(distribution = "Simplex beta binomial") %>%
#   bind_rows(
#     dirichlet_multinomial_counts %>%
#       mutate(distribution = "Dirichlet-multinomial")
#   )

plot_pairs_beta_binomial =
  sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
  filter(.draw == 1) %>%
  rename(sample = N, count = .value) %>%
  mutate(category = glue("V{M}"), sample = as.character(sample)) %>%
  select(sample, category, count) %>%

  # Afjust counts
  with_groups(sample, ~ mutate(.x, sample_total = sum(count))) %>%
  mutate(multiplier = 1000 / sample_total) %>%
  mutate(count = as.integer(count*multiplier)) %>%
  select(-multiplier, -sample_total) %>%

  # Plot
  spread(category, count) %>%
  ggpairs(2:5,  lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1))   ) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  theme(  strip.background =element_rect(fill="white", color="white")   ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Simplex-bounded Beta-binomial")


plot_pairs_dirichlet_multinomial =
  dirichlet_multinomial_counts %>%
  spread(category, count) %>%
  ggpairs(2:5,  lower = list(continuous = wrap("points", alpha = 0.3,    size=0.1))   ) +
  theme_bw() +
  theme(text = element_text(size=8)) +
  theme(  strip.background =element_rect(fill="white", color="white")   ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("Dirichlet-multinomial")


plot_pairs = wrap_elements(ggmatrix_gtable(plot_pairs_beta_binomial)) + wrap_elements(ggmatrix_gtable(plot_pairs_dirichlet_multinomial))

ggsave(
  "dev/plot_pairs_WEHI_seminar.pdf",
  units = c("mm"),
  width = 183 ,
  height = 183/2 ,
  limitsize = FALSE
)

# Plot

distro_plot_2 =
  both_data_sets %>%
  ggplot(aes(category, count, color = distribution,  size=distribution)) +
  geom_jitter(height = 0) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0.1, 1)) +
  theme_bw() +
  theme(legend.position = "none")

distro_violin_1 =
  both_data_sets %>%
  ggplot(aes(category, count + 1, fill = distribution)) +
  geom_split_violin() +
  scale_fill_manual(values = c("white", "#e11f28")) +
  #scale_size_discrete(range = c(0.1, 1)) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none")


p1 = distro_plot_1 + plot_estimate + distro_plot_2 + distro_violin_1

ggsave(
  "dev/approximation_dirichlet_multinomial_WEHI_seminar.pdf",
  units = c("mm"),
  width = 183 ,
  height = 120 ,
  limitsize = FALSE
)


# Dirichlet multinomial upper bound
dirichlet_multinomial_counts =
  rdirmnom(1000, 1000, c(1,3,10,50)/20) %>%
  as_tibble(rownames = "sample") %>%
  gather(category, count, -sample) %>%
  mutate(count = as.integer(count))


distro_plot_3 =
  dirichlet_multinomial_counts %>%
  ggplot(aes(category, count)) +
  geom_jitter(height = 0, size = 0.1) +
  theme_bw()

# Fit
job({
  estimate =
    dirichlet_multinomial_counts %>%
    sccomp_glm(
      ~1, sample, category, count,
      check_outliers = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      verbose=TRUE
    )
})

plot_estimate_2 =
  estimate %>%
  unnest(concentration) %>%
  ggplot(aes(`.median_(Intercept)`, mean)) +
  geom_errorbar(aes(ymin = `2.5%`, ymax=`97.5%`), color="#4DAF4A", alpha = 0.8, width=0) +
  geom_errorbar(aes(xmin = `.lower_(Intercept)`, xmax=`.upper_(Intercept)`), color="#4DAF4A", alpha = 0.8, width=0) +
  geom_point() +
  geom_abline(intercept = 5.7496330, slope = -0.9650953, linetype = "dashed", color="grey") +
  xlab("Category rate") +
  ylab("Category log-concentration") +
  theme_bw()

# Simulate data
job({
  fit_object =
    dirichlet_multinomial_counts %>%
    sccomp:::estimate_multi_beta_binomial_glm(
      ~1, sample, category, count,
      check_outliers = FALSE,
      approximate_posterior_inference = FALSE,
      prior_mean_variable_association = list(intercept = c(0, 5), slope = c(0,  5), standard_deviation = c(0, 2)),
      verbose=TRUE
    )
})

rng =  rstan::gqs(
  sccomp:::stanmodels$glm_multi_beta_binomial_generate_date,
  #rstan::stan_model("inst/stan/glm_multi_beta_binomial_generate_date.stan"),
  draws =  as.matrix(fit_object$fit),
  data = fit_object$data_for_model
)

both_data_sets =
  sccomp:::draws_to_tibble_x_y(rng, "counts", "N", "M") %>%
  filter(.draw == 1) %>%
  rename(sample = N, count = .value) %>%
  mutate(category = glue("V{M}"), sample = as.character(sample)) %>%
  mutate(distribution = "Simplex beta binomial") %>%
  bind_rows(
    dirichlet_multinomial_counts %>%
      mutate(distribution = "Dirichlet-multinomial")
  )

# Plot
distro_plot_4 =
  both_data_sets %>%
  ggplot(aes(category, count, color = distribution,  size=distribution)) +
  geom_jitter(height = 0) +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0.1, 1)) +
  theme_bw() +
  theme(legend.position = "none")

distro_violin_2 =
  both_data_sets %>%
  ggplot(aes(category, count + 1, fill = distribution)) +
  geom_split_violin() +
  scale_fill_manual(values = c("white", "#e11f28")) +
  #scale_size_discrete(range = c(0.1, 1)) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none")


p2 = distro_plot_3 + plot_estimate_2 + distro_plot_4 + distro_violin_2

ggsave(
  "dev/approximation_dirichlet_multinomial_WEHI_seminar_2.pdf",
  plot = p2,
  units = c("mm"),
  width = 183 ,
  height = 120 ,
  limitsize = FALSE
)

