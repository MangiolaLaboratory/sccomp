## Paired pre-fits: `my_estimate_with_variance` (association on) vs
## `my_estimate_exclude_mean_variability_association` (association off in the
## prior) — same composition and variability formulas. Defined in `test-plot.R`
## when CmdStan is available.

test_that("exclude_mean_variability_association is stored and passed to Stan", {
  skip_cmdstan()

  ex <- attr(my_estimate_exclude_mean_variability_association, "model_input")$exclude_mean_variability_association
  expect_true(isTRUE(ex) || identical(as.integer(ex), 1L))

  def <- attr(my_estimate_with_variance, "model_input")$exclude_mean_variability_association
  expect_true(isFALSE(def) || identical(as.integer(def), 0L))
})

test_that("exclude_mean_variability_association: prec_sd and alpha Stan summaries are finite", {
  skip_cmdstan()

  fit_ex <- attr(my_estimate_exclude_mean_variability_association, "fit")
  fit_def <- attr(my_estimate_with_variance, "fit")

  for (nm in c("prec_sd", "prec_slope_1", "prec_intercept_1")) {
    s_ex <- fit_ex$summary(variables = nm)
    s_def <- fit_def$summary(variables = nm)
    expect_true(all(is.finite(s_ex$mean)))
    expect_true(all(is.finite(s_def$mean)))
    rh <- intersect(names(s_ex), c("rhat", "Rhat"))
    if (length(rh) == 1L) {
      expect_true(all(s_ex[[rh]] < 1.05, na.rm = TRUE))
      expect_true(all(s_def[[rh]] < 1.05, na.rm = TRUE))
    }
  }

  sa_ex <- fit_ex$summary(variables = "alpha")
  sa_def <- fit_def$summary(variables = "alpha")
  expect_true(all(is.finite(sa_ex$mean)))
  expect_true(all(is.finite(sa_def$mean)))
})

test_that("exclude_mean_variability_association: sccomp_test and 2D interval plot build cleanly", {
  skip_cmdstan()

  tested_ex <- my_estimate_exclude_mean_variability_association |> sccomp_test()
  tested_def <- my_estimate_with_variance |> sccomp_test()

  expect_true(all(is.finite(tested_ex$v_effect), na.rm = TRUE))
  expect_true(all(is.finite(tested_ex$v_lower), na.rm = TRUE))
  expect_true(all(is.finite(tested_ex$v_upper), na.rm = TRUE))
  expect_true(all(is.finite(tested_def$v_effect), na.rm = TRUE))

  p_ex <- sccomp_plot_intervals_2D(tested_ex, add_marginal_density = FALSE)
  p_def <- sccomp_plot_intervals_2D(tested_def, add_marginal_density = FALSE)

  expect_s3_class(p_ex, "ggplot")
  expect_s3_class(p_def, "ggplot")

  gb_ex <- ggplot2::ggplot_build(p_ex)
  gb_def <- ggplot2::ggplot_build(p_def)

  y_spans_ex <- vapply(
    gb_ex$layout$panel_params,
    function(p) diff(p$y.range),
    numeric(1)
  )
  y_spans_def <- vapply(
    gb_def$layout$panel_params,
    function(p) diff(p$y.range),
    numeric(1)
  )

  expect_true(all(is.finite(y_spans_ex)))
  expect_true(all(is.finite(y_spans_def)))
  expect_true(all(y_spans_ex > 0))
  expect_true(all(y_spans_def > 0))
})

test_that("exclude_mean_variability_association: v_* intervals match raw alpha and don't inflate", {
  # Regression guard: when the mean-variability association is excluded, the
  # Stan model fixes the slope contribution to zero in the variability prior
  # but leaves `prec_slope_*` sampled from its (wide) prior. If
  # `compute_alpha_normalised_draws()` ever resumes applying the
  # `alpha - beta * prec_slope` correction in this mode, every v_* interval
  # picks up prior-sized noise scaled by |beta| (5-10x inflation on real data).
  # We tie the v_* width to the raw alpha 95% width so any regression is
  # caught immediately.
  skip_cmdstan()

  fit_ex <- attr(my_estimate_exclude_mean_variability_association, "fit")
  tested_ex <- my_estimate_exclude_mean_variability_association |> sccomp_test()
  tested_def <- my_estimate_with_variance |> sccomp_test()

  alpha_summary <- fit_ex$summary(variables = "alpha")
  alpha_95_width <- 3.92 * mean(alpha_summary$sd, na.rm = TRUE)

  v_width_ex <- mean(tested_ex$v_upper - tested_ex$v_lower, na.rm = TRUE)
  v_width_def <- mean(tested_def$v_upper - tested_def$v_lower, na.rm = TRUE)

  expect_lt(abs(v_width_ex - alpha_95_width) / alpha_95_width, 0.25)
  expect_lt(v_width_ex / max(v_width_def, .Machine$double.eps), 3)
})
