# Hypothesis testing for smooth (spline) terms.
#
# `sccomp_test_smooth()` tests whether the mean composition differs between two
# positions (or two intervals) of a continuous covariate entered as a smooth
# `s()` term. It is a thin wrapper that:
#   1. predicts the *mean* linear predictor (`mu_unconstrained`, the logit-scale
#      mean — no beta-binomial overdispersion) at the requested covariate
#      positions via the existing replicate-data / generated-quantities pass,
#   2. forms the per-draw contrast `mu(to) - mu(from)` for each cell group, and
#   3. feeds those contrast draws into the same `draws_to_statistics()` engine
#      used by `sccomp_test()`, so the returned `c_pH0` / `c_FDR` mean exactly
#      what they mean elsewhere in sccomp.
#
# The contrast lives on the logit (unconstrained) scale, which is the scale
# `test_composition_above_logit_fold_change` is calibrated for.

#' Test differences along a smooth term
#'
#' @description
#' Tests whether the mean cell-group composition differs between two positions
#' (or two intervals) of a continuous covariate modelled with a smooth `s()`
#' term. The test is a contrast of the model's mean prediction on the logit
#' scale, evaluated at the two positions and summarised with the same
#' probability-of-null (`pH0`) and false-discovery-rate (`FDR`) machinery as
#' [sccomp_test()].
#'
#' @details
#' Because the smooth curve can be non-linear, "difference across an interval"
#' is ambiguous. Two comparison modes are provided:
#' \itemize{
#'   \item \code{"endpoints"} (default): \code{from} and \code{to} are single
#'     values; the contrast is \eqn{\mu(\code{to}) - \mu(\code{from})}.
#'   \item \code{"average"}: \code{from} and \code{to} are length-2 intervals
#'     \code{c(lo, hi)}; the contrast is the average of the curve over the
#'     \code{to} interval minus the average over the \code{from} interval
#'     (each interval sampled at \code{resolution} points).
#' }
#'
#' Any covariate not being varied cancels out of the contrast because the
#' logit-scale predictor is additive; it therefore does not matter what value
#' those covariates take, as long as it is held constant. For factor smooths
#' (\code{s(x, g, bs = "fs")}) or by-factor smooths (\code{s(x, by = g)}) the
#' grouping factor selects \emph{which} curve is tested; set it via \code{at},
#' e.g. \code{at = list(tissue = "tumor")}.
#'
#' The test is about the \strong{mean} composition: it uses the expected
#' linear predictor (\code{mu_unconstrained}), not the overdispersed
#' beta-binomial realisation. It therefore answers "does the expected
#' composition differ between these covariate positions?", not "would an
#' individual sample differ?".
#'
#' @param fit The result of [sccomp_estimate()] with a smooth term in
#'   \code{formula_composition}.
#' @param smooth Character. The continuous covariate inside the smooth to test
#'   (e.g. \code{"pseudotime"}). Can be \code{NULL} when the model has exactly
#'   one smooth term over one continuous covariate, in which case it is
#'   auto-detected.
#' @param from Numeric. The reference position. A single value for
#'   \code{comparison = "endpoints"}, or a length-2 interval \code{c(lo, hi)}
#'   for \code{comparison = "average"}.
#' @param to Numeric. The comparison position, same shape as \code{from}.
#' @param comparison One of \code{"endpoints"} or \code{"average"}. See details.
#' @param at Optional named list fixing the value of other covariates (e.g.
#'   the grouping factor of a factor smooth). Covariates not supplied are held
#'   at their first observed value (which cancels for non-grouping covariates).
#' @param resolution Integer. Number of grid points used to approximate each
#'   interval average when \code{comparison = "average"}.
#' @param test_composition_above_logit_fold_change Positive numeric. Effect
#'   threshold for the hypothesis test, on the logit scale — identical meaning
#'   to the argument of [sccomp_test()].
#' @param percent_false_positive Numeric in (0, 100). Used for the credible
#'   interval width, as in [sccomp_test()].
#' @param number_of_draws Integer. Number of posterior draws used for the
#'   prediction pass.
#' @param mcmc_seed Integer. Seed for the generated-quantities pass.
#' @param robust Logical. Currently unused placeholder for API parity with
#'   [sccomp_predict()]; the effect is always the posterior mean of the
#'   contrast.
#'
#' @return A tibble with one row per cell group:
#' \itemize{
#'   \item \code{cell_group} — the cell group tested.
#'   \item \code{smooth} — the continuous covariate tested.
#'   \item \code{from}, \code{to} — the compared positions (as labels).
#'   \item \code{c_lower}, \code{c_effect}, \code{c_upper} — 95% CI and posterior
#'     mean of the logit-scale contrast \eqn{\mu(\code{to}) - \mu(\code{from})}.
#'   \item \code{c_pH0} — probability the effect is within the null region.
#'   \item \code{c_FDR} — false-discovery rate across cell groups.
#' }
#'
#' @examples
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'     # add a continuous covariate
#'     counts_obj$pseudotime <- as.numeric(factor(counts_obj$sample))
#'
#'     fit <- sccomp_estimate(
#'       counts_obj,
#'       ~ s(pseudotime, k = 4), ~ 1, "sample", "cell_group", "count",
#'       cores = 1
#'     )
#'
#'     fit |> sccomp_test_smooth(from = 2, to = 8)
#'   }
#' }
#'
#' @export
sccomp_test_smooth <- function(fit,
                               smooth = NULL,
                               from,
                               to,
                               comparison = c("endpoints", "average"),
                               at = NULL,
                               resolution = 20L,
                               test_composition_above_logit_fold_change = 0.1,
                               percent_false_positive = 5,
                               number_of_draws = 500,
                               mcmc_seed = sample_seed(),
                               robust = FALSE) {
  check_and_install_cmdstanr()
  UseMethod("sccomp_test_smooth", fit)
}

#' @importFrom dplyr tibble bind_rows left_join group_by summarise ungroup arrange mutate select distinct pull
#' @importFrom tidyr pivot_wider nest unnest
#' @importFrom rlang enquo quo_name sym
#' @export
sccomp_test_smooth.sccomp_tbl <- function(fit,
                                          smooth = NULL,
                                          from,
                                          to,
                                          comparison = c("endpoints", "average"),
                                          at = NULL,
                                          resolution = 20L,
                                          test_composition_above_logit_fold_change = 0.1,
                                          percent_false_positive = 5,
                                          number_of_draws = 500,
                                          mcmc_seed = sample_seed(),
                                          robust = FALSE) {

  # Define the variables as NULL to avoid CRAN NOTES
  M <- NULL
  .chain <- NULL
  .iteration <- NULL
  .draw <- NULL
  .value <- NULL
  .smooth_group <- NULL
  parameter <- NULL

  comparison   <- match.arg(comparison)
  .sample      <- attr(fit, ".sample")
  .cell_group  <- attr(fit, ".cell_group")
  model_input  <- attr(fit, "model_input")
  smooth_specs <- get_smooth_results(fit)$smooth_specs

  if (is.null(smooth_specs) || length(smooth_specs) == 0)
    stop("sccomp says: this fit has no smooth (s()/t2()) terms; use sccomp_test() for parametric contrasts.")

  # ---- Resolve which continuous covariate to test -------------------------
  # For each smooth, the continuous variable(s) are `term` minus the grouping
  # factor `fterm` (present only for factor smooths).
  smooth_continuous <- lapply(smooth_specs, function(sm) setdiff(sm$term, sm$fterm))
  all_continuous    <- unique(unlist(smooth_continuous))

  if (is.null(smooth)) {
    if (length(all_continuous) != 1)
      stop(sprintf(
        "sccomp says: the model has several smooth covariates (%s); specify one via `smooth = `.",
        paste(all_continuous, collapse = ", ")
      ))
    smooth <- all_continuous
  } else if (!smooth %in% all_continuous) {
    stop(sprintf(
      "sccomp says: `%s` is not a continuous smooth covariate in this model (available: %s).",
      smooth, paste(all_continuous, collapse = ", ")
    ))
  }

  # ---- Validate from/to shapes against comparison mode --------------------
  if (comparison == "endpoints") {
    if (length(from) != 1L || length(to) != 1L)
      stop("sccomp says: for comparison = 'endpoints', `from` and `to` must be single values.")
    from_grid <- from
    to_grid   <- to
    from_label <- as.character(from)
    to_label   <- as.character(to)
  } else {
    if (length(from) != 2L || length(to) != 2L)
      stop("sccomp says: for comparison = 'average', `from` and `to` must be length-2 intervals c(lo, hi).")
    from_grid <- seq(from[1], from[2], length.out = resolution)
    to_grid   <- seq(to[1], to[2], length.out = resolution)
    from_label <- sprintf("[%s, %s]", from[1], from[2])
    to_label   <- sprintf("[%s, %s]", to[1], to[2])
  }

  # ---- Build new_data at the requested positions --------------------------
  new_data <- build_smooth_test_newdata(
    fit, smooth, from_grid, to_grid, at, .sample
  )

  # ---- Predict the mean linear predictor (logit scale) --------------------
  rng <- replicate_data(
    fit,
    formula_composition = NULL,       # use the fitted model
    formula_variability = ~ 1,
    new_data            = new_data,
    number_of_draws     = number_of_draws,
    mcmc_seed           = mcmc_seed
  )

  sample_col   <- quo_name(.sample)
  sample_names <- new_data |> pull(!!.sample)

  # Per-draw logit-scale mean (mu_unconstrained), keeping chain/iteration so
  # the diagnostics in draws_to_statistics remain meaningful.
  draws <-
    rng |>
    draws_to_tibble_x_y("mu_unconstrained", "M", "N") |>
    # Map N -> sample name -> group ("from"/"to")
    left_join(
      tibble(
        N     = seq_along(sample_names),
        !!sample_col := sample_names
      ),
      by = "N"
    ) |>
    left_join(
      new_data |> select(!!.sample, .smooth_group),
      by = sample_col
    ) |>
    # Map M -> cell group name
    left_join(
      tibble(
        M = seq_len(ncol(model_input$y)),
        !!quo_name(.cell_group) := colnames(model_input$y)
      ),
      by = "M"
    )

  # ---- Form the contrast per (cell_group, draw) ---------------------------
  # Average over grid points within each group (a no-op for endpoints, where
  # each group has one point), then difference to - from.
  contrast_draws <-
    draws |>
    group_by(!!.cell_group, M, .chain, .iteration, .draw, .smooth_group) |>
    summarise(.value = mean(.value), .groups = "drop") |>
    pivot_wider(names_from = .smooth_group, values_from = .value) |>
    mutate(.value = to - from) |>
    mutate(parameter = sprintf("%s: %s -> %s", smooth, from_label, to_label)) |>
    select(!!.cell_group, M, parameter, .chain, .iteration, .draw, .value)

  # ---- Reuse the standard pH0 / FDR engine --------------------------------
  result <-
    contrast_draws |>
    draws_to_statistics(
      percent_false_positive / 100,
      test_composition_above_logit_fold_change,
      !!.cell_group,
      "c_"
    ) |>
    select(-M, -parameter) |>
    mutate(smooth = smooth, from = from_label, to = to_label) |>
    select(!!.cell_group, smooth, from, to, everything())

  result
}


#' Build the sample-wise new_data for a smooth test
#'
#' Two groups of rows ("from" and "to"), each covering the requested grid of
#' the continuous covariate. All other covariates in the composition formula
#' are held at a constant value (from `at`, else the first observed value), so
#' they cancel in the contrast; a grouping factor supplied via `at` selects
#' which curve of a factor / by-factor smooth is tested.
#'
#' @keywords internal
#' @noRd
#' @importFrom dplyr tibble bind_rows distinct mutate row_number all_of
build_smooth_test_newdata <- function(fit, smooth, from_grid, to_grid, at, .sample) {

  sample_col  <- quo_name(.sample)
  count_data  <- attr(fit, "count_data")
  formula_composition <- attr(fit, "formula_composition")

  covariates <- parse_formula(formula_composition)
  other_covs <- setdiff(covariates, smooth)

  # Defaults for non-varied covariates: `at` override, else first observed.
  ref_row <- count_data |> distinct(!!.sample, .keep_all = TRUE) |> head(1)
  cov_value <- function(v) {
    if (!is.null(at) && v %in% names(at)) return(at[[v]])
    if (v %in% colnames(count_data)) return(ref_row[[v]])
    stop(sprintf("sccomp says: covariate `%s` needed by the model is missing from the data and not supplied via `at`.", v))
  }
  other_values <- lapply(other_covs, cov_value)
  names(other_values) <- other_covs

  make_group <- function(x_grid, group_label) {
    df <- tibble(!!smooth := x_grid, .smooth_group = group_label)
    for (v in other_covs) df[[v]] <- other_values[[v]]
    df
  }

  new_data <-
    bind_rows(
      make_group(from_grid, "from"),
      make_group(to_grid,   "to")
    ) |>
    mutate(!!sample_col := sprintf("smoothtest_%03d", row_number()))

  new_data
}
