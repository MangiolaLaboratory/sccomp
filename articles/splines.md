# Smooth (spline) terms in sccomp

Abstract

This vignette introduces `brms`-style smooth terms (`s()`, `t2()`) in
`sccomp`. Smooths let you model non-linear effects of a continuous
covariate — pseudotime, age, dose, time-since-treatment — on cell-type
proportions, while keeping the rest of the `sccomp` modelling pipeline
(compositional likelihood, variability modelling, outlier detection,
random effects) unchanged. We start from the bundled `counts_obj`
dataset, add a synthetic continuous covariate, fit a smooth, predict on
a dense grid, and visualise both the fitted curves and the underlying
basis decomposition.

## Why smooth terms?

`sccomp` fits a (logit-)linear model on the unconstrained scale and then
maps it to the simplex via softmax. If you only have categorical
covariates this is all you need. But with a **continuous covariate** —
pseudotime, age, dose, days since treatment — the relationship between
covariate and proportion is rarely linear. A spline (`s()` from `mgcv`)
replaces a single linear column with a small, flexible basis that lets
the curve bend, while a penalty keeps the bend honest.

`sccomp` adopts the same syntax as `brms` / `mgcv`:

``` r

~ type + s(pseudotime, k = 5)        # smooth of one variable
~ type + t2(pseudotime, dose)        # tensor product of two variables
```

Under the hood every smooth is split into two pieces:

- an unpenalised **null-space** block `Xf` (the linear “trend” part),
  and
- a penalised **wiggly** block `Xr` (the bend-y part).

`Xf` joins the design matrix as ordinary fixed-effect columns; `Xr`
joins one of `sccomp`’s random-effect slots and is governed by a single
smoothing parameter `sds`. This is exactly how `brms` handles smooths,
and it means the Stan model does **not** change — only the R-side data
assembly.

## What is supported

| Feature | Status |
|----|----|
| `s(x, k = ...)` in `formula_composition` | yes (1 slot) |
| `s(x, z, bs = "fs", k = ...)` (factor smooth) | yes (3 slots: 1 wiggly + 2 per-level deviation blocks) |
| Multiple smooth terms in the same formula | yes (each penalty block consumes one RE slot) |
| Mix smooths with parametric terms (`~ type + s(x)`) | yes |
| Mix smooths with explicit `(... | g)` random effects | yes (total ≤ 4 slots) |
| Smooths in `formula_variability` | not yet |
| `t2(x, y, ...)` (tensor-product smooth) | fits, but [`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md) not yet supported |
| `s(x, by = factor)` (independent smooths per level) | works but consumes one slot **per level** |

The 4-slot budget comes from the existing `sccomp` Stan model. If your
formula needs more than 4 random-effect slots in total (explicit
`( | g)` clauses + smooth penalty blocks),
[`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
will stop with an informative error. Note that multi-penalty smooth
types like `bs = "fs"` emit several blocks; in mgcv / brms these
correspond to several `sds_*` parameters.

## Worked example

``` r

library(dplyr)
library(tidyr)
library(ggplot2)
library(sccomp)
data("counts_obj")
```

### 1. Synthetic continuous covariate

`counts_obj` has no continuous covariate, so we attach one — a
*pseudotime* between 0 and 6 — that varies across samples but is
constant within a sample (a covariate must always be sample-level in
`sccomp`).

``` r

set.seed(1)
samples <- levels(counts_obj$sample)
pt_map <- tibble(
  sample     = samples,
  pseudotime = seq(0, 6, length.out = length(samples)) +
                 rnorm(length(samples), 0, 0.05)
)

counts_pt <- counts_obj |> left_join(pt_map, by = "sample")
counts_pt |> distinct(sample, type, pseudotime) |> head()
#> # A tibble: 6 × 3
#>   sample       type   pseudotime
#>   <chr>        <fct>       <dbl>
#> 1 10x_6K       benign    -0.0313
#> 2 10x_8K       benign     0.325 
#> 3 GSE115189    benign     0.590 
#> 4 SCP345_580   benign     1.03  
#> 5 SCP345_860   benign     1.28  
#> 6 SCP424_pbmc1 benign     1.54
```

### 2. Fit `~ type + s(pseudotime, k = 5)`

``` r

fit <- counts_pt |>
  sccomp_estimate(
    formula_composition = ~ type + s(pseudotime, k = 5),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    verbose             = FALSE,
    max_sampling_iterations = 2000
  )
```

Look at the design-matrix column names that
[`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
prints. You should see one new fixed column

    s(pseudotime, k = 5)__lin1

and three penalised basis columns living in an RE slot

    s(pseudotime, k = 5)___basis01
    s(pseudotime, k = 5)___basis02
    s(pseudotime, k = 5)___basis03

The `__lin*` columns are the null-space block (`Xf`); the `___basis*`
columns are the wiggly block (`Xr`). The single smoothing parameter
`sds` for this smooth controls how much the `___basis*` coefficients are
pulled toward zero — exactly the role of `sds_*` in `brms`.

### 3. Tabular output

[`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md)
works unchanged. It returns one row per (cell_group, parameter), and
smooth basis columns appear like any other coefficient:

``` r

fit |>
  sccomp_test() |>
  select(cell_group, parameter, c_lower, c_effect, c_upper, c_FDR) |>
  head(10)
#> Joining with `by = join_by(cell_group, M, parameter)`
#> sccomp model
#> ============
#> 
#> Model specifications:
#>   Family: multi_beta_binomial 
#>   Composition formula: ~type + s(pseudotime, k = 5) 
#>   Variability formula: ~1 
#>   Inference method: pathfinder 
#> 
#> Data: Samples: 20   Cell groups: 36 
#> 
#> Column prefixes: c_ -> composition parameters  v_ -> variability parameters
#> 
#> Convergence diagnostics:
#>   For each parameter, n_eff is the effective sample size and R_k_hat is the potential
#>   scale reduction factor on split chains (at convergence, R_k_hat = 1).
#> 
#> # A tibble: 10 × 6
#>    cell_group parameter                      c_lower c_effect c_upper    c_FDR
#>    <chr>      <chr>                            <dbl>    <dbl>   <dbl>    <dbl>
#>  1 B1         (Intercept)                      0.918   1.19     1.46  0       
#>  2 B1         typecancer                      -0.859  -0.559   -0.229 0.000625
#>  3 B1         s(pseudotime, k = 5)__lin1      -0.445  -0.0367   0.321 0.222   
#>  4 B1         s(pseudotime, k = 5)___basis01  -1.90   -0.653    0.172 0.0585  
#>  5 B1         s(pseudotime, k = 5)___basis02  -1.56   -0.338    0.399 0.151   
#>  6 B1         s(pseudotime, k = 5)___basis03  -0.944  -0.242    0.419 0.143   
#>  7 B2         (Intercept)                      0.534   0.849    1.15  0       
#>  8 B2         typecancer                      -1.16   -0.811   -0.445 0       
#>  9 B2         s(pseudotime, k = 5)__lin1      -0.256   0.211    0.603 0.0684  
#> 10 B2         s(pseudotime, k = 5)___basis01  -2.22   -1.10    -0.165 0.0113
```

`c_effect` on the `___basis*` rows is rarely interpreted directly — the
basis is internal — but the linear `__lin1` row is a clean “is there an
overall monotone trend in pseudotime for this cell group?” test, and the
basis-level `c_FDR` values let you spot non-linear effects.

### 4. Predict on a dense grid

To draw the fitted curve, predict at a regular grid of `pseudotime`
values (remember to also supply any other covariates the formula needs —
here, `type`):

``` r

grid <- expand.grid(
  pseudotime = seq(min(pt_map$pseudotime), max(pt_map$pseudotime), length.out = 60),
  type       = levels(counts_obj$type)
) |>
  as_tibble() |>
  mutate(sample = sprintf("grid_%03d", row_number()))

pred <- fit |>
  sccomp_predict(new_data = grid, number_of_draws = 200)
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 1.238 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.

head(pred)
#> # A tibble: 6 × 10
#>   pseudotime type   sample   cell_group proportion_mean proportion_lower
#>        <dbl> <fct>  <chr>    <chr>                <dbl>            <dbl>
#> 1    -0.0313 benign grid_001 B1                 0.0612           0.0436 
#> 2    -0.0313 benign grid_001 B2                 0.0298           0.0184 
#> 3    -0.0313 benign grid_001 B3                 0.00812          0.00501
#> 4    -0.0313 benign grid_001 BM                 0.00845          0.00549
#> 5    -0.0313 benign grid_001 CD4 1              0.0256           0.0188 
#> 6    -0.0313 benign grid_001 CD4 2              0.0370           0.0257 
#> # ℹ 4 more variables: proportion_upper <dbl>, unconstrained_mean <dbl>,
#> #   unconstrained_lower <dbl>, unconstrained_upper <dbl>
```

[`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md)
already returns the `new_data` columns (`pseudotime`, `type`, `sample`)
alongside `proportion_mean` and the 95% credible interval
`proportion_lower` / `proportion_upper`.

### 5. Plot: fitted curves over the raw proportions

We show the six cell groups with the largest pseudotime range in the
posterior mean.

``` r

top_cgs <- pred |>
  group_by(cell_group) |>
  summarise(rng = max(proportion_mean) - min(proportion_mean)) |>
  arrange(desc(rng)) |>
  slice_head(n = 6) |>
  pull(cell_group)

raw <- counts_pt |>
  group_by(sample) |>
  mutate(prop = count / sum(count)) |>
  ungroup() |>
  filter(cell_group %in% top_cgs)

pred |>
  filter(cell_group %in% top_cgs) |>
  ggplot(aes(pseudotime, proportion_mean, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = proportion_lower, ymax = proportion_upper),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  geom_point(data = raw, aes(y = prop), alpha = 0.5, size = 1.2) +
  facet_wrap(~ cell_group, scales = "free_y", ncol = 3) +
  labs(x = "pseudotime", y = "proportion",
       title = "Fitted compositional curves vs pseudotime") +
  theme_minimal(base_size = 10)
```

![](splines_files/figure-html/plot-curves-1.png)

The blue/red curves are posterior means, ribbons are 95% credible
intervals, points are observed sample proportions. Note the wiggle —
e.g. CD8 1 peaks around pseudotime ≈ 4 — which a purely linear
`~ type + pseudotime` model could not capture.

### 6. Same model with a smaller basis (`k = 3`)

The same parametric + smooth layout works with a lower `k`. Fewer basis
functions mean a smoother curve and less risk of over-fitting when
sample size is modest — here `k = 3` yields one unpenalised linear
column and a single wiggly basis column (instead of three with `k = 5`).

``` r

fit_k3 <- counts_pt |>
  sccomp_estimate(
    formula_composition = ~ type + s(pseudotime, k = 3),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    verbose             = FALSE,
    max_sampling_iterations = 2000
  )
```

The printed design matrix should list `(Intercept)`, `typecancer`, one
`__lin` column, and a single `___basis01` column for the smooth.

``` r

pred_k3 <- fit_k3 |>
  sccomp_predict(new_data = grid, number_of_draws = 200)
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 1.239 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.
```

``` r

pred_k3 |>
  filter(cell_group %in% top_cgs) |>
  ggplot(aes(pseudotime, proportion_mean, colour = type, fill = type)) +
  geom_ribbon(aes(ymin = proportion_lower, ymax = proportion_upper),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  geom_point(data = raw, aes(y = prop, colour = type), alpha = 0.5, size = 1.2) +
  facet_wrap(~ cell_group, scales = "free_y", ncol = 3) +
  labs(x = "pseudotime", y = "proportion",
       title = "Fitted curves: ~ type + s(pseudotime, k = 3)",
       subtitle = "One wiggly basis column — smoother than k = 5") +
  theme_minimal(base_size = 10)
```

![](splines_files/figure-html/plot-curves-k3-1.png)

Compared with the `k = 5` figure above, the credible bands are often a
little tighter and the lines less wiggly: the model has less capacity to
bend, which is appropriate when you want a gentle pseudotime trend on
top of the `type` effect.

## What `s()` actually decomposes into

The smooth metadata used to fit the model is stored on the result, so we
can re-evaluate the basis at any new pseudotime value. This is exactly
what
[`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md)
does internally; we can also do it by hand to *see* the basis.

``` r

smooth_meta    <- sccomp:::get_smooth_results(fit)
smooth_specs   <- smooth_meta$smooth_specs
smooth_re_objs <- smooth_meta$smooth_re_objs

grid_basis <- tibble(
  pseudotime = seq(min(pt_map$pseudotime), max(pt_map$pseudotime),
                   length.out = 200),
  type = levels(counts_obj$type)[1]  # any constant works; smooth ignores it
)

basis <- sccomp:::predict_smooth_at_newdata(
  smooth_specs[[1]], smooth_re_objs[[1]], grid_basis
)
# Plain s() has one penalised block — for multi-penalty smooths
# (e.g. bs = "fs", t2(...)) `Xr_blocks` would contain several matrices.
Xr_single <- basis$Xr_blocks[[1]]

basis_long <- bind_rows(
  as.data.frame(basis$Xf) |>
    setNames(sprintf("lin%d", seq_len(ncol(basis$Xf)))) |>
    mutate(pseudotime = grid_basis$pseudotime) |>
    pivot_longer(-pseudotime, names_to = "column", values_to = "value") |>
    mutate(block = "Xf  (unpenalised, joins fixed effects)"),
  as.data.frame(Xr_single) |>
    setNames(sprintf("basis%d", seq_len(ncol(Xr_single)))) |>
    mutate(pseudotime = grid_basis$pseudotime) |>
    pivot_longer(-pseudotime, names_to = "column", values_to = "value") |>
    mutate(block = "Xr  (penalised, joins one RE slot, shared sds)")
)

ggplot(basis_long, aes(pseudotime, value, colour = column)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ block, scales = "free_y", ncol = 1) +
  labs(x = "pseudotime", y = "basis value",
       title = "Decomposition of s(pseudotime, k = 5)",
       subtitle = "Linear null-space (Xf) is a fixed effect; wiggly basis (Xr) is one RE slot") +
  theme_minimal(base_size = 10)
```

![](splines_files/figure-html/basis-decomp-1.png)

The top panel is the **null-space**: for a default thin-plate `s()` this
is just the linear term in `pseudotime`. The bottom panel is the
**wiggly basis**: each of these columns gets one coefficient per cell
group, and a *single* smoothing standard deviation `sds` per cell group
controls how much wiggle is allowed. Increasing `k` would add more
wiggly basis columns; the penalty would then keep the extra capacity
from over-fitting.

## Factor smooths: per-group curves with shared structure

A common ask is: “fit a smooth per tissue (or per condition, per donor,
…) but partial-pool across groups”. The `mgcv` shorthand for this is
`s(x, factor, bs = "fs")` — a *factor smooth*. It builds one common
wiggly curve and a per-level deviation, governed jointly by a small set
of smoothing parameters (three, in mgcv’s standard reparameterisation:
the common wiggly part and two per-level deviation blocks). In `sccomp`
each of these blocks occupies its own RE slot.

### A worked example

We extend `counts_obj` with a synthetic `tissue` factor (3 levels) on
top of the existing `pseudotime`:

``` r

set.seed(11)
samples <- levels(counts_obj$sample)
covs <- tibble::tibble(
  sample     = samples,
  pseudotime = seq(0, 6, length.out = length(samples)) +
                 rnorm(length(samples), 0, 0.05),
  tissue     = factor(rep(c("blood", "lymph", "tumor"),
                          length.out = length(samples)))
)
counts_tissue <- counts_obj |> dplyr::left_join(covs, by = "sample")
counts_tissue |> dplyr::distinct(sample, type, pseudotime, tissue) |> head()
#> # A tibble: 6 × 4
#>   sample       type   pseudotime tissue
#>   <chr>        <fct>       <dbl> <fct> 
#> 1 10x_6K       benign    -0.0296 blood 
#> 2 10x_8K       benign     0.317  lymph 
#> 3 GSE115189    benign     0.556  tumor 
#> 4 SCP345_580   benign     0.879  blood 
#> 5 SCP345_860   benign     1.32   lymph 
#> 6 SCP424_pbmc1 benign     1.53   tumor
```

Fit `~ s(pseudotime, tissue, bs = "fs", k = 5)`:

``` r

fit_fs <- counts_tissue |>
  sccomp_estimate(
    formula_composition = ~ s(pseudotime, tissue, bs = "fs", k = 5),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    verbose             = FALSE,
    max_sampling_iterations = 2000
  )

# This one smooth occupies 3 RE slots (one per mgcv penalty block);
# you can confirm via:
attr(fit_fs, "model_input")$ncol_X_random_eff
#> [1] 9 3 3 0
```

The three non-zero entries correspond to the wiggly main block (9
columns for `k = 5` × 3 tissues) and two smaller per-level deviation
blocks.

### Predict and plot per-tissue curves

[`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md)
takes a grid that crosses `pseudotime` with `tissue`:

``` r

grid_fs <- expand.grid(
  pseudotime = seq(min(covs$pseudotime), max(covs$pseudotime),
                   length.out = 60),
  tissue     = levels(covs$tissue)
) |>
  as_tibble() |>
  mutate(sample = sprintf("fs_grid_%03d", row_number()))

pred_fs <- fit_fs |>
  sccomp_predict(new_data = grid_fs, number_of_draws = 200)
#> Running standalone generated quantities after 1 MCMC chain, with 1 thread(s) per chain...
#> 
#> Chain 1  Elapsed Time: 1.863 seconds (Generated Quantities) 
#> Chain 1 finished in 0.0 seconds.

head(pred_fs)
#> # A tibble: 6 × 10
#>   pseudotime tissue sample      cell_group proportion_mean proportion_lower
#>        <dbl> <fct>  <chr>       <chr>                <dbl>            <dbl>
#> 1    -0.0296 blood  fs_grid_001 B1                 0.0657           0.0399 
#> 2    -0.0296 blood  fs_grid_001 B2                 0.0240           0.0120 
#> 3    -0.0296 blood  fs_grid_001 B3                 0.0103           0.00674
#> 4    -0.0296 blood  fs_grid_001 BM                 0.00680          0.00420
#> 5    -0.0296 blood  fs_grid_001 CD4 1              0.0280           0.0183 
#> 6    -0.0296 blood  fs_grid_001 CD4 2              0.0420           0.0276 
#> # ℹ 4 more variables: proportion_upper <dbl>, unconstrained_mean <dbl>,
#> #   unconstrained_lower <dbl>, unconstrained_upper <dbl>
```

Pick the most pseudotime-responsive cell groups and plot one curve per
tissue:

``` r

top_cgs_fs <- pred_fs |>
  group_by(cell_group) |>
  summarise(rng = max(proportion_mean) - min(proportion_mean)) |>
  arrange(desc(rng)) |>
  slice_head(n = 6) |>
  pull(cell_group)

raw_fs <- counts_tissue |>
  group_by(sample) |>
  mutate(prop = count / sum(count)) |>
  ungroup() |>
  filter(cell_group %in% top_cgs_fs)

pred_fs |>
  filter(cell_group %in% top_cgs_fs) |>
  ggplot(aes(pseudotime, proportion_mean,
             colour = tissue, fill = tissue)) +
  geom_ribbon(aes(ymin = proportion_lower, ymax = proportion_upper),
              alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  geom_point(data = raw_fs, aes(y = prop, colour = tissue),
             alpha = 0.5, size = 1.2) +
  facet_wrap(~ cell_group, scales = "free_y", ncol = 3) +
  labs(x = "pseudotime", y = "proportion",
       title = "Factor smooth s(pseudotime, tissue, bs = 'fs', k = 5)",
       subtitle = "Per-tissue curves sharing one factor-smooth structure") +
  theme_minimal(base_size = 10)
```

![](splines_files/figure-html/fs-plot-1.png)

The three coloured curves are the same smooth applied to different
tissues, allowed to deviate at each pseudotime point. Where the
deviations are small (e.g. CD8 1) the curves track each other; where
they are large (e.g. CD8 2) one tissue diverges. Whether the deviations
are “real” is expressed by the credible-interval ribbons — overlap means
the data don’t support a difference.

### fs vs `s(x, by = tissue)`

`s(x, by = tissue)` is an alternative spelling that *also* gives a
per-tissue curve, but it treats each tissue’s smooth as independent:
each level gets its own smoothing parameter, and the curves can have
totally unrelated shape. With L tissue levels this consumes L sccomp
slots (one per level). Use `bs = "fs"` when you want a single shared
smoothness scale and partial pooling across levels; use `by =` when the
tissues are clearly heterogeneous and a shared smoothness assumption
would be misleading.

## Multiple smooths and mixing with random effects

[`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
is happy with several smooths, and you can also keep explicit
`(... | g)` random-effect clauses:

``` r

~ type +
  s(pseudotime, k = 5) +
  s(age,        k = 4) +
  (1 | donor)
```

Each `s()` / `t2()` term consumes one of the four available RE slots.
The remaining slots are shared with any explicit `( | g)` clauses, so
the total must stay ≤ 4. If you exceed it,
[`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
will stop with a message naming the smooth and the slot budget.

## Caveats and current limitations

- **No smooths in `formula_variability`** yet.
  [`sccomp_estimate()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_estimate.md)
  errors out if you try.
- **`t2()` tensor smooths** fit fine but
  [`sccomp_predict()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_predict.md)
  is not yet supported for them:
  [`mgcv::PredictMat()`](https://rdrr.io/pkg/mgcv/man/smoothCon.html)
  applies the `t2` constraint differently from `smoothCon()`, producing
  a basis that doesn’t line up with the fitted coefficients. The predict
  call errors out cleanly. Use `s(x, z, bs = "fs")` (shown above) when
  you want a smooth indexed by a factor, or stack univariate
  `s(x) + s(y)` smooths for additive multivariate effects.
- **`s(x, by = factor)`** decomposes into one independent smooth per
  factor level, each costing a separate RE slot. With many levels you
  will run out of the 4-slot budget. For partial-pooling across factor
  levels use `bs = "fs"` (shown above) — it fits in 3 slots regardless
  of the number of levels.
- The smooth columns are part of
  [`sccomp_test()`](https://mangiolalaboratory.github.io/sccomp/reference/sccomp_test.md)’s
  output. The basis-level `c_effect` / `c_FDR` are diagnostics, not
  biological effect sizes; interpret the *curve* (predicted proportions)
  rather than individual basis coefficients.

## Session info

    #> R version 4.6.1 (2026-06-24)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.4 LTS
    #> 
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #> 
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    #> 
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> other attached packages:
    #> [1] sccomp_2.1.34     instantiate_0.2.3 ggplot2_4.0.3     tidyr_1.3.2      
    #> [5] dplyr_1.2.1      
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] tidyselect_1.2.1            farver_2.1.2               
    #>  [3] S7_0.2.2                    fastmap_1.2.0              
    #>  [5] SingleCellExperiment_1.34.0 tensorA_0.36.2.1           
    #>  [7] digest_0.6.39               lifecycle_1.0.5            
    #>  [9] processx_3.9.0              magrittr_2.0.5             
    #> [11] posterior_1.7.0             compiler_4.6.1             
    #> [13] rlang_1.3.0                 sass_0.4.10                
    #> [15] tools_4.6.1                 utf8_1.2.6                 
    #> [17] yaml_2.3.12                 data.table_1.18.4          
    #> [19] knitr_1.51                  labeling_0.4.3             
    #> [21] S4Arrays_1.12.0             htmlwidgets_1.6.4          
    #> [23] DelayedArray_0.38.2         RColorBrewer_1.1-3         
    #> [25] cmdstanr_0.9.0              abind_1.4-8                
    #> [27] withr_3.0.3                 purrr_1.2.2                
    #> [29] BiocGenerics_0.58.1         desc_1.4.3                 
    #> [31] grid_4.6.1                  stats4_4.6.1               
    #> [33] scales_1.4.0                SummarizedExperiment_1.42.0
    #> [35] cli_3.6.6                   rmarkdown_2.31             
    #> [37] crayon_1.5.3                ragg_1.5.2                 
    #> [39] generics_0.1.4              otel_0.2.0                 
    #> [41] tzdb_0.5.0                  cachem_1.1.0               
    #> [43] stringr_1.6.0               splines_4.6.1              
    #> [45] parallel_4.6.1              XVector_0.52.0             
    #> [47] matrixStats_1.5.0           vctrs_0.7.3                
    #> [49] Matrix_1.7-5                jsonlite_2.0.0             
    #> [51] callr_3.8.0                 hms_1.1.4                  
    #> [53] patchwork_1.3.2             IRanges_2.46.0             
    #> [55] S4Vectors_0.50.1            ggrepel_0.9.8              
    #> [57] systemfonts_1.3.2           jquerylib_0.1.4            
    #> [59] glue_1.8.1                  ggside_0.4.1               
    #> [61] pkgdown_2.2.1               ps_1.9.3                   
    #> [63] distributional_0.8.1        stringi_1.8.7              
    #> [65] gtable_0.3.6                GenomicRanges_1.64.0       
    #> [67] tibble_3.3.1                pillar_1.11.1              
    #> [69] htmltools_0.5.9             Seqinfo_1.2.0              
    #> [71] R6_2.6.1                    textshaping_1.0.5          
    #> [73] evaluate_1.0.5              lattice_0.22-9             
    #> [75] Biobase_2.72.0              readr_2.2.0                
    #> [77] backports_1.5.1             bslib_0.11.0               
    #> [79] Rcpp_1.1.2                  SparseArray_1.12.2         
    #> [81] nlme_3.1-169                checkmate_2.3.4            
    #> [83] mgcv_1.9-4                  xfun_0.60                  
    #> [85] fs_2.1.0                    MatrixGenerics_1.24.0      
    #> [87] forcats_1.0.1               prettydoc_0.4.1            
    #> [89] pkgconfig_2.0.3
