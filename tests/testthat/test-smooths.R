## Tests for brms-style smooth terms (s(), t2()) in sccomp formulas.
##
## The unit tests at the top do not need cmdstan: they only exercise the
## R-side machinery in `R/smooths.R` and the patched `parse_formula()`.
## The end-to-end smoke test does need cmdstan and is gated by
## `skip_cmdstan()`.

library(dplyr)

# -----------------------------------------------------------------------------
# Pure-R unit tests
# -----------------------------------------------------------------------------

test_that("has_smooth_terms detects s() / t2() in a formula", {
  expect_false(sccomp:::has_smooth_terms(~ 1))
  expect_false(sccomp:::has_smooth_terms(~ type + age))
  expect_true(sccomp:::has_smooth_terms(~ type + s(age)))
  expect_true(sccomp:::has_smooth_terms(~ s(x, k = 5) + s(y, k = 3)))
  expect_true(sccomp:::has_smooth_terms(~ t2(x, z)))
})


test_that("strip_smooth_terms removes s()/t2() and preserves other terms", {
  strip   <- sccomp:::strip_smooth_terms
  labels  <- function(fm) attr(stats::terms(fm), "term.labels")
  has_int <- function(fm) attr(stats::terms(fm), "intercept") == 1L
  
  # ---- No-op cases --------------------------------------------------------
  
  # NULL passes through unchanged.
  expect_null(strip(NULL))
  
  # Smooth-free formulas are unchanged in term content and intercept setting.
  fm_plain <- ~ type + age + (1 | sample)
  expect_equal(labels(strip(fm_plain)), labels(fm_plain))
  expect_true(has_int(strip(fm_plain)))
  
  fm_noint <- ~ 0 + type + age
  expect_equal(labels(strip(fm_noint)), labels(fm_noint))
  expect_false(has_int(strip(fm_noint)))
  
  # `~ 1` and `~ 0` round-trip.
  expect_true (has_int(strip(~ 1)))
  expect_false(has_int(strip(~ 0)))
  
  # ---- Removes s() / t2(), keeps the rest ---------------------------------
  
  fm <- ~ type + s(pseudotime, k = 5)
  expect_false(sccomp:::has_smooth_terms(strip(fm)))
  expect_equal(labels(strip(fm)), "type")
  expect_true(has_int(strip(fm)))
  
  # `+ 0` is preserved when smooths are stripped.
  fm_noint_s <- ~ 0 + type + s(pseudotime)
  expect_equal(labels(strip(fm_noint_s)), "type")
  expect_false(has_int(strip(fm_noint_s)))
  
  # Smooth-only formula collapses to `~ 1` / `~ 0`.
  expect_equal(labels(strip(~ s(x))), character(0))
  expect_true (has_int(strip(~ s(x))))
  expect_equal(labels(strip(~ 0 + s(x))), character(0))
  expect_false(has_int(strip(~ 0 + s(x))))
  
  # `t2()` is removed just like `s()`.
  expect_equal(labels(strip(~ type + t2(x, z, k = c(3, 4)))), "type")
  
  # Multiple smooths and a `by =` argument: all smooth calls vanish, parametric
  # terms remain. Variables that only appeared inside `s(..., by = ...)` are
  # gone because they were not parametric to begin with.
  fm_many <- ~ type + age + s(x, k = 4) + s(y, by = condition) + t2(u, v)
  expect_equal(labels(strip(fm_many)), c("type", "age"))
  
  # ---- Preserves random-effect clauses verbatim ---------------------------
  
  # Single bar term sits alongside the parametric terms.
  fm_re <- ~ type + s(pseudotime) + (1 | sample)
  expect_equal(labels(strip(fm_re)), c("type", "1 | sample"))
  
  # Multiple bar terms, including multi-LHS bars, are all kept.
  fm_re_multi <- ~ 0 + type + s(x) + (1 | sample) + (age | dataset)
  expect_setequal(
    labels(strip(fm_re_multi)),
    c("type", "1 | sample", "age | dataset")
  )
  expect_false(has_int(strip(fm_re_multi)))
  
  # ---- Preserves transforms and interactions ------------------------------
  
  # Function transforms and interactions survive intact.
  fm_tx <- ~ type * age + I(age^2) + log(exposure) + s(pseudotime)
  out_tx <- strip(fm_tx)
  # `type * age` expands to type + age + type:age in term.labels; check the
  # main effects and the interaction are all retained, and no smooth label
  # remains.
  expect_true (all(c("type", "age", "type:age", "I(age^2)", "log(exposure)")
                   %in% labels(out_tx)))
  expect_false(any(grepl("^s\\(",  labels(out_tx))))
  expect_false(any(grepl("^t2\\(", labels(out_tx))))
  
  # Explicit interaction term with `:` is preserved.
  expect_setequal(
    labels(strip(~ a:b + s(x))),
    "a:b"
  )
})


test_that("parse_formula expands smooth specials to underlying variables", {
  expect_equal(
    sccomp:::parse_formula(~ type + s(pseudotime, k = 5)),
    c("type", "pseudotime")
  )
  # `s(x, by = condition)` should expose both x and condition
  expect_setequal(
    sccomp:::parse_formula(~ s(x, by = condition)),
    c("x", "condition")
  )
  # Random-effect bars are still excluded
  expect_equal(
    sccomp:::parse_formula(~ type + s(age, k = 4) + (1 | dataset)),
    c("type", "age")
  )
})


test_that("parse_formula_smooths produces Xf + Xr with brms shapes (tp, 1 block)", {
  skip_if_not_installed("mgcv")
  set.seed(1)
  dat <- data.frame(x = sort(stats::runif(50, 0, 10)))
  res <- sccomp:::parse_formula_smooths(~ s(x, bs = "tp", k = 8), dat)
  
  expect_length(res$Xf_list, 1L)
  expect_length(res$Xr_list, 1L)              # 1 penalty -> 1 block
  expect_equal(res$Xr_to_smooth, 1L)
  expect_length(res$smooth_specs, 1L)
  expect_length(res$smooth_re_objs, 1L)
  expect_equal(res$smooth_labels, "s(x, bs = \"tp\", k = 8)")
  
  Xf <- res$Xf_list[[1]]
  Xr <- res$Xr_list[[1]]
  expect_equal(nrow(Xf), nrow(dat))
  expect_equal(nrow(Xr), nrow(dat))
  # tp + absorb.cons: 1 unpenalised column + (k-2) wiggly columns
  expect_equal(ncol(Xf), 1L)
  expect_equal(ncol(Xr), 8L - 2L)
  expect_true(all(grepl("___basis", colnames(Xr))))
})


test_that("parse_formula_smooths flattens t2 / fs into multiple Xr blocks", {
  skip_if_not_installed("mgcv")
  set.seed(3)
  
  # t2: 3 penalties -> 3 random blocks
  dat_t <- data.frame(
    x = sort(stats::runif(40, 0, 10)),
    y = sort(stats::runif(40, 0, 5))
  )
  res_t <- sccomp:::parse_formula_smooths(~ t2(x, y, k = c(4, 3)), dat_t)
  expect_length(res_t$Xr_list, 3L)
  expect_equal(res_t$Xr_to_smooth, c(1L, 1L, 1L))
  # All blocks share the same parent smooth label
  expect_length(res_t$smooth_labels, 1L)
  # Multi-block: column names are disambiguated with __b<b>___basis
  expect_true(any(grepl("__b1___basis", colnames(res_t$Xr_list[[1]]))))
  expect_true(any(grepl("__b2___basis", colnames(res_t$Xr_list[[2]]))))
  expect_true(any(grepl("__b3___basis", colnames(res_t$Xr_list[[3]]))))
  
  # fs: factor smooth with 3 levels and k=5 -> 3 penalties -> 3 blocks
  dat_f <- data.frame(
    x = seq(0, 6, length.out = 24),
    z = factor(rep(c("a", "b", "c"), length.out = 24))
  )
  res_f <- sccomp:::parse_formula_smooths(~ s(x, z, bs = "fs", k = 5), dat_f)
  expect_length(res_f$Xr_list, 3L)
  expect_equal(res_f$Xr_to_smooth, c(1L, 1L, 1L))
  # Sum of block widths matches the raw fs basis dimension (k * n_levels = 15)
  expect_equal(sum(vapply(res_f$Xr_list, ncol, integer(1))), 5L * 3L)
  # fs has no null-space block (everything is penalised)
  expect_equal(ncol(res_f$Xf_list[[1]]), 0L)
})


test_that("parse_formula_smooths returns identity for smooth-free formulas", {
  dat <- data.frame(type = letters[1:5], age = 1:5)
  res <- sccomp:::parse_formula_smooths(~ type + age, dat)
  expect_length(res$Xf_list, 0L)
  expect_length(res$Xr_list, 0L)
  expect_length(res$smooth_specs, 0L)
  # Parametric formula identical to input
  expect_equal(deparse(res$parametric_formula), deparse(~ type + age))
})


test_that("predict_smooth_at_newdata reproduces Xr at training data (tp)", {
  skip_if_not_installed("mgcv")
  set.seed(2)
  dat <- data.frame(x = sort(stats::runif(40, 0, 10)))
  res <- sccomp:::parse_formula_smooths(~ s(x, bs = "tp", k = 7), dat)
  
  sm     <- res$smooth_specs[[1]]
  re_obj <- res$smooth_re_objs[[1]]
  Xr_fit <- res$Xr_list[[1]]
  Xf_fit <- res$Xf_list[[1]]
  
  pred <- sccomp:::predict_smooth_at_newdata(sm, re_obj, dat)
  
  expect_length(pred$Xr_blocks, 1L)
  expect_equal(dim(pred$Xr_blocks[[1]]), dim(Xr_fit))
  expect_equal(dim(pred$Xf), dim(Xf_fit))
  
  # Same columns, same rows: round-trip should be numerically identical
  expect_equal(unname(pred$Xr_blocks[[1]]), unname(Xr_fit), tolerance = 1e-10)
  expect_equal(unname(pred$Xf),            unname(Xf_fit), tolerance = 1e-10)
})


test_that("predict_smooth_at_newdata reproduces all Xr blocks at training data (fs)", {
  skip_if_not_installed("mgcv")
  set.seed(4)
  
  dat <- data.frame(
    x = seq(0, 6, length.out = 24),
    z = factor(rep(c("a", "b", "c"), length.out = 24))
  )
  res <- sccomp:::parse_formula_smooths(~ s(x, z, bs = "fs", k = 5), dat)
  sm     <- res$smooth_specs[[1]]
  re_obj <- res$smooth_re_objs[[1]]
  
  pred <- sccomp:::predict_smooth_at_newdata(sm, re_obj, dat)
  expect_length(pred$Xr_blocks, length(res$Xr_list))
  
  for (b in seq_along(pred$Xr_blocks)) {
    expect_equal(
      unname(pred$Xr_blocks[[b]]),
      unname(res$Xr_list[[b]]),
      tolerance = 1e-10
    )
  }
  expect_equal(unname(pred$Xf), unname(res$Xf_list[[1]]), tolerance = 1e-10)
})


test_that("predict_smooth_at_newdata errors cleanly on t2()", {
  skip_if_not_installed("mgcv")
  dat <- data.frame(
    x = sort(stats::runif(30, 0, 10)),
    y = sort(stats::runif(30, 0, 5))
  )
  res <- sccomp:::parse_formula_smooths(~ t2(x, y, k = c(4, 3)), dat)
  expect_error(
    sccomp:::predict_smooth_at_newdata(
      res$smooth_specs[[1]], res$smooth_re_objs[[1]], dat
    ),
    regexp = "t2\\(\\) smooths is not yet supported"
  )
})


test_that("build_smooth_slot produces a sccomp RE-slot structure", {
  Xr <- matrix(stats::rnorm(20), nrow = 5, ncol = 4,
               dimnames = list(NULL, paste0("c", 1:4)))
  slot <- sccomp:::build_smooth_slot(Xr, "s(x)")
  expect_equal(slot$ncol, 4L)
  expect_equal(slot$n_factors, 1L)
  expect_equal(slot$n_groups, 4L)
  expect_equal(dim(slot$gfi), c(1L, 4L))
  expect_equal(rownames(slot$gfi), "s(x)")
  expect_equal(as.vector(slot$gfi), 1:4)
})


# -----------------------------------------------------------------------------
# End-to-end smoke test (requires cmdstan)
# -----------------------------------------------------------------------------

test_that("sccomp_estimate works with a single s() smooth on a real dataset", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("seurat_obj")
  
  fit <- sccomp_estimate(
    seurat_obj,
    formula_composition = ~ s(continuous_covariate, k = 4),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )
  
  mi <- attr(fit, "model_input")
  sr <- sccomp:::get_smooth_results(fit)
  
  # Smooth metadata lives on `model_input` (R-only, not a Stan list element).
  expect_false(is.null(sr))
  expect_equal(length(sr$smooth_specs), 1L)
  expect_equal(length(sr$smooth_re_objs), 1L)
  expect_equal(length(sr$smooth_labels), 1L)
  # Not a Stan list element, and not duplicated on the fitted tibble.
  expect_null(mi$smooth_results)
  expect_null(attr(fit, "smooth_results"))
  
  # Model input slot 1 (no explicit RE clauses) carries the wiggly basis.
  expect_gt(mi$ncol_X_random_eff[1], 0L)
  expect_equal(mi$how_many_factors_in_random_design[1], 1L)
  
  # Fixed-effect design has the unpenalised null-space column too.
  expect_true(any(grepl("__lin", colnames(mi$X))))
})


test_that("fs factor smooth fits and predicts end-to-end (multi-block)", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("counts_obj")
  
  # Attach a per-sample pseudotime and a 3-level tissue factor (both sample-
  # level). counts_obj has 20 samples; we just assign tissues cyclically.
  set.seed(7)
  samples <- levels(counts_obj$sample)
  covs <- tibble::tibble(
    sample     = samples,
    pseudotime = seq(0, 6, length.out = length(samples)) +
                   rnorm(length(samples), 0, 0.05),
    tissue     = factor(rep(c("blood", "lymph", "tumor"),
                            length.out = length(samples)))
  )
  counts_pt <- counts_obj |> dplyr::left_join(covs, by = "sample")
  
  fit <- sccomp_estimate(
    counts_pt,
    formula_composition = ~ s(pseudotime, tissue, bs = "fs", k = 5),
    formula_variability = ~ 1,
    sample              = "sample",
    cell_group          = "cell_group",
    abundance           = "count",
    cores               = 1,
    inference_method    = "pathfinder",
    max_sampling_iterations = 200,
    verbose = FALSE
  )
  
  mi <- attr(fit, "model_input")
  sr <- sccomp:::get_smooth_results(fit)
  
  # One mgcv smooth -> 3 penalty blocks -> 3 RE slots occupied
  expect_equal(length(sr$smooth_specs), 1L)
  expect_equal(sum(mi$ncol_X_random_eff > 0L), 3L)
  
  # All wiggly columns landed in the RE slots; no `__lin` cols (fs has empty
  # null-space)
  expect_false(any(grepl("__lin", colnames(mi$X))))
  
  # Predict on a grid that spans pseudotime x tissue. Must succeed (this is
  # the failure mode the multi-block refactor fixes).
  grid <- expand.grid(
    pseudotime = seq(min(covs$pseudotime), max(covs$pseudotime),
                     length.out = 20),
    tissue     = levels(covs$tissue)
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(sample = sprintf("grid_%03d", dplyr::row_number()))
  
  pred <- sccomp_predict(fit, new_data = grid, number_of_draws = 100)
  expect_true(all(c("pseudotime", "tissue", "proportion_mean") %in% names(pred)))
  expect_equal(nrow(pred), nrow(grid) * length(unique(counts_obj$cell_group)))
  # No NAs in the posterior summary
  expect_true(all(is.finite(pred$proportion_mean)))
})


test_that("smooth slot count beyond N_RE_SLOTS errors cleanly", {
  skip_if_not_installed("mgcv")
  skip_cmdstan()
  data("seurat_obj")
  
  # 4 smooths + 1 RE clause needs 5 slots; budget is 4. Use the same
  # continuous_covariate four times under different `k` (parser doesn't
  # de-duplicate intentionally — each call is a fresh basis).
  expect_error(
    sccomp_estimate(
      seurat_obj,
      formula_composition =
        ~ s(continuous_covariate, k = 3) +
          s(continuous_covariate, k = 4) +
          s(continuous_covariate, k = 5) +
          s(continuous_covariate, k = 6) +
          (1 | group__),
      formula_variability = ~ 1,
      sample = "sample", cell_group = "cell_group",
      cores = 1, inference_method = "pathfinder",
      max_sampling_iterations = 50, verbose = FALSE
    ),
    regexp = "random-effect slot",
    fixed  = FALSE
  )
})


test_that("smooths in formula_variability are rejected", {
  skip_cmdstan()
  data("seurat_obj")
  
  expect_error(
    sccomp_estimate(
      seurat_obj,
      formula_composition = ~ 1,
      formula_variability = ~ s(continuous_covariate),
      sample = "sample", cell_group = "cell_group",
      cores = 1, inference_method = "pathfinder",
      max_sampling_iterations = 50, verbose = FALSE
    ),
    regexp = "smooth terms .* `formula_variability`"
  )
})
