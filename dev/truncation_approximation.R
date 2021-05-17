
library(tidyverse)
library(furrr)
library(rstan)
library(TailRank)
plan(multisession, workers = 20)
mo <- stan_model(model_code = 'data{ int N; int exposure[N]; int y[N];} parameters{ real<lower=0, upper=1> mu;  real<upper = 10> precision; } model{  target += beta_binomial_lpmf(y | exposure, mu * exp(precision), (1.0 - mu) * exp(precision) );   precision ~ gamma( 4, 2); mu ~ beta(2,2); }' )

shrinkage =
  expand_grid(
    mu = seq(0,1, length.out=102) %>% tail(-1) %>% head(-1),
    precision = seq(1, 6, length.out=100) %>% exp()
  ) %>%

  # Put plateau to shape
  rowwise() %>%
  mutate(correction_factor = 1.0/min(mu, 1-mu) * 0.5) %>%
  ungroup() %>%
  mutate(precision = precision * correction_factor) %>%
  mutate(
    alpha =  mu * precision,
    beta = (1-mu) * precision
  ) %>%
  mutate(fit = future_map2(alpha, beta, ~ {

    exposure = 1000
    y =
      rbb(1000, exposure, .x, .y) %>%
      enframe() %>%
      filter(value <= qbb(0.975, exposure, .x, .y)) %>%
      filter(value >= qbb(0.025, exposure, .x, .y)) %>%
      pull(value)

    summary(sampling(
      mo,
      data=list(y = y, N = length(y), exposure = rep(exposure, length(y))),
      chains=1
      #,
      #iter = 250, warmup = 150
    ))$summary[1:2,"50%"] %>%
      enframe() %>%
      spread(name, value) %>%
      setNames(c("mu_truncated", "precision_truncated")) %>%
      mutate(precision_truncated = exp(precision_truncated))

  })) %>%
  unnest(fit) %>%
  mutate(shrink = precision / precision_truncated)

saveRDS(shrinkage, "dev/shrinkage_positive_control.rds")
saveRDS(shrinkage, "dev/shrinkage.rds")



library(viridis)
my_theme =
  theme_bw() +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_line(size = 0.1),
    text = element_text(size=12),
    legend.position="bottom",
    aspect.ratio=1,
    strip.background = element_blank(),
    axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
    axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
  )
(

  shrinkage %>%
    ggplot(aes(mu, precision)) +
    geom_tile(aes(fill = shrink)) +
    scale_fill_viridis(	option="magma") +
    my_theme

  ) %>%
  ggsave(plot = .,
         "dev/truncatin_approximation.pdf",
         useDingbats=FALSE,
         units = c("mm"),
         width = 183
  )


ggplot(shrinkage %>% head(10000) , aes(log(precision), shrink)) +
  geom_line(aes(color = mu)) + facet_wrap(~mu) + scale_y_log10()
