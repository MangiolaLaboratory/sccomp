
#' Plot 1D Intervals for Cell-group Effects
#'
#' This function creates a series of 1D interval plots for cell-group effects, highlighting significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test. A value of 0.2 correspond to a change in cell proportion of 10% for a cell type with baseline proportion of 50%. That is, a cell type goes from 45% to 50%. When the baseline proportion is closer to 0 or 1 this effect thrshold has consistent value in the logit uncontrained scale.
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "pH0".
#' @param factor Optional character string selecting one model factor to plot. If provided, plots are restricted to that factor plus `(Intercept)`.
#' @param sort_by Character vector indicating how to sort taxa. Options are "none" (default), "effect" (by effect size), "significance" (by FDR/pH0), or "alphabetical".
#' @importFrom patchwork wrap_plots
#' @importFrom forcats fct_reorder fct_inorder
#' @importFrom tidyr drop_na
#'
#' @export
#'
#' @return A combined plot of 1D interval plots.
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~1,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1
#'     ) |>
#'     sccomp_test()
#'
#'   # Example usage:
#'   my_plot = sccomp_plot_intervals_1D(estimate, sort_by = "effect")
#'
#'   }
#' }
#'
#'
sccomp_plot_intervals_1D = function(
    .data,
    factor = NULL,
    significance_threshold = 0.05,
    test_composition_above_logit_fold_change = .data |> attr("test_composition_above_logit_fold_change"),
    show_fdr_message = TRUE,
  significance_statistic = c("pH0", "FDR"),
  sort_by = c("none", "effect", "significance", "alphabetical")
) {
  significance_statistic <- match.arg(significance_statistic)
  sort_by <- match.arg(sort_by)

  # Define the variables as NULL to avoid CRAN NOTES
  parameter <- NULL
  estimate <- NULL
  value <- NULL
  pH0 <- NULL
  FDR <- NULL
  effect <- NULL

  .cell_group = attr(.data, ".cell_group")

  # Check if test have been done
  if(.data |> select(ends_with("FDR")) |> ncol() |> equals(0))
    stop("sccomp says: to produce plots, you need to run the function sccomp_test() on your estimates.")

  .data <- subset_results_by_factor(.data, factor, keep_intercept = TRUE)

  plot_list =
    .data |>

    # Reshape data
    select(-contains("n_eff"), -contains("R_k_hat"), -contains("_ess"), -contains("_rhat")) |>
    pivot_longer(c(contains("c_"), contains("v_")), names_sep = "_", names_to = c("which", "estimate")) |>
    pivot_wider(names_from = estimate, values_from = value) |>

    # Nest data by parameter and which
    nest(data = -c(parameter, which)) |>
    mutate(plot = pmap(
      list(data, which, parameter),
      function(plot_data, plot_which, plot_param) {
        # Check if there are any statistics to plot
        if(plot_data |> filter(!is.na(effect)) |> nrow() |> equals(0))
          return(NA)

        # Choose color variable and legend
        if (significance_statistic == "FDR") {
          color_aes <- aes(xmin = lower, xmax = upper, color = FDR < significance_threshold)
          color_scale <- scale_color_manual(values = c("grey40", "red"))
          legend_title <- "FDR < significance_threshold"
        } else {
          color_aes <- aes(xmin = lower, xmax = upper, color = pH0 < significance_threshold)
          color_scale <- scale_color_manual(values = c("grey40", "red"))
          legend_title <- "pH0 < significance_threshold"
        }

        # Determine y-axis variable based on sort_by
        # Use string-based column selection instead of NSE
        if (sort_by == "none") {
          # No sorting
          plot_data$y_var <- plot_data[[.cell_group]]
        } else if (sort_by == "effect") {
          # Sort by absolute effect size
          plot_data$y_var <- fct_reorder(plot_data[[.cell_group]], abs(plot_data$effect))
        } else if (sort_by == "significance") {
          # Sort by significance
          if (significance_statistic == "FDR") {
            plot_data$y_var <- fct_reorder(plot_data[[.cell_group]], plot_data$FDR, .desc = TRUE)
          } else {
            plot_data$y_var <- fct_reorder(plot_data[[.cell_group]], plot_data$pH0, .desc = TRUE)
          }
        } else if (sort_by == "alphabetical") {
          # Alphabetical sorting
          plot_data <- plot_data %>% arrange(.data[[.cell_group]])
          plot_data$y_var <- fct_inorder(plot_data[[.cell_group]])
        }

        ggplot(plot_data, aes(x = effect, y = y_var)) +
          geom_vline(xintercept = test_composition_above_logit_fold_change, colour = "grey") +
          geom_vline(xintercept = -test_composition_above_logit_fold_change, colour = "grey") +
          geom_errorbar(color_aes) +
          geom_point() +
          color_scale +
          xlab("Credible interval of the effect") +
          ylab("Cell group") +
          ggtitle(sprintf("%s %s", plot_which, plot_param)) +
          sccomp_theme() +
          theme(legend.position = "bottom") +
          guides(color = guide_legend(title = legend_title))
      }
    )) %>%

    # Filter out NA plots
    filter(!is.na(plot)) |>
    pull(plot)

  # Combine all individual plots into one plot
  combined_plot <- plot_list |>
    wrap_plots(ncol = plot_list |> length() |> sqrt() |> ceiling())

  # Only show the FDR message if significance_statistic == "FDR" and show_fdr_message is TRUE
  if (significance_statistic == "FDR" && show_fdr_message) {
    combined_plot <- combined_plot + theme(plot.caption = ggplot2::element_text(hjust = 0))
    combined_plot <- combined_plot + patchwork::plot_annotation(
      caption = paste(
        "Bayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)",
        "\nFDR-significant populations may cross fold change thresholds because Bayesian FDR considers posterior probabilities rather than p-values.",
        "\nThe method sorts null hypothesis probabilities in ascending order and calculates cumulative averages for robust false discovery control.",
        sep = ""
      )
    )
  }
  combined_plot
}

#' Plot 2D Intervals for Mean-Variance Association
#'
#' This function creates a 2D interval plot for mean-variance association, handling both single and bimodal models.
#' It highlights significant differences based on a given significance threshold.
#'
#' @param .data Data frame containing the main data.
#' @param significance_threshold Numeric value specifying the significance threshold for highlighting differences. Default is 0.05.
#' @param test_composition_above_logit_fold_change A positive integer. It is the effect threshold used for the hypothesis test.
#' @param show_fdr_message Logical. Whether to show the Bayesian FDR interpretation message on the plot. Default is TRUE.
#' @param significance_statistic Character vector indicating which statistic to highlight. Default is "pH0".
#' @param factor Optional character string selecting one model factor to plot. If provided, plots are restricted to that factor plus `(Intercept)`.
#' @param add_marginal_density Logical. Whether to add marginal density plots on adjusted panels. Default is TRUE.
#' @param omit_ci Logical. Whether to omit credible interval error bars. Default is FALSE.
#'
#' @importFrom dplyr filter arrange mutate if_else row_number bind_rows distinct slice pull with_groups
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_errorbar geom_point geom_line geom_blank aes xlab ylab facet_wrap theme_bw theme labs guides guide_legend scale_color_manual scale_alpha_manual scale_fill_manual element_rect element_blank element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom stringr str_detect
#' @importFrom ggside geom_ysidedensity theme_ggside_void scale_ysidex_continuous
#'
#' @export
#'
#' @return A ggplot object representing the 2D interval plot.
#'
#' @examples
#'
#' print("cmdstanr is needed to run this example.")
#'
#' \donttest{
#'   if (instantiate::stan_cmdstan_exists()) {
#'     data("counts_obj")
#'
#'     estimate <- sccomp_estimate(
#'       counts_obj,
#'       ~ type,
#'       ~type,
#'       "sample",
#'       "cell_group",
#'       "count",
#'       cores = 1,
#'       bimodal_mean_variability_association = TRUE
#'     ) |>
#'     sccomp_test()
#'
#'     # Example usage:
#'     my_plot = sccomp_plot_intervals_2D(estimate)
#'
#'   }
#' }
#'
sccomp_plot_intervals_2D <- function(
    .data,
    factor = NULL,
    significance_threshold = 0.05,
    test_composition_above_logit_fold_change =
      .data |> attr("test_composition_above_logit_fold_change"),
    show_fdr_message = TRUE,
    significance_statistic = c("pH0", "FDR"),
    add_marginal_density = TRUE,
    omit_ci = FALSE
) {

  significance_statistic <- match.arg(significance_statistic)

  # Locals declared as NULL to silence R CMD check NOTEs about "no visible binding"
  # for tidyverse NSE references (dplyr/ggplot evaluate these as column names).
  v_effect <- NULL
  parameter <- NULL
  c_effect <- NULL
  c_lower <- NULL
  c_upper <- NULL
  c_FDR <- NULL
  v_lower <- NULL
  v_upper <- NULL
  v_FDR <- NULL
  cell_type_label <- NULL
  pH0 <- NULL
  FDR <- NULL
  c_pH0 <- NULL
  v_pH0 <- NULL
  component <- NULL
  assigned_component <- NULL
  v_value <- NULL

  .cell_group <- attr(.data, ".cell_group")

  # The plot relies on FDR / pH0 columns that only exist after sccomp_test();
  # fail fast with a user-meaningful message rather than later with a cryptic NSE error.
  if(.data |> select(ends_with("FDR")) |> ncol() == 0)
    stop("sccomp says: you need to run sccomp_test() first.")

  # Intercept is kept even when subsetting to a single factor because the Intercept
  # panel acts as the visual baseline against which factor-specific effects are read.
  .data <- subset_results_by_factor(.data, factor, keep_intercept = TRUE)

  fit <- attr(.data, "fit")

  # Determine model topology *before* pulling Stan summaries: `prec_*_2` variables
  # are only declared in the bimodal Stan program, so blindly requesting them
  # would error on single-component fits.
  bimodal_flag <- attr(.data, "model_input")$bimodal_mean_variability_association
  # no check needed; bimodal_flag should always be present
  bimodal_flag <- isTRUE(as.logical(bimodal_flag))

  # Posterior summaries of the mean-variability regression coefficients
  #   v_pred = -(prec_intercept + prec_slope * c_effect)
  # These are pulled once and reused per-parameter via `param_idx`.
  prec_intercept_1_summary <- fit$summary("prec_intercept_1")
  prec_slope_1_summary <- fit$summary("prec_slope_1")

  # Only parameters that actually have a v_effect (i.e. are estimated under the
  # variability design Xa) can be plotted on the 2D space.
  param_names <- .data |>
    filter(!is.na(v_effect)) |>
    distinct(parameter) |>
    pull("parameter")

  # Map each parameter name to its column index in the variability design matrix
  # `Xa` because Stan stores `prec_*` indexed by Xa column position, not by name.
  param_idx <- match(param_names, colnames(attr(.data, "model_input")$Xa))
  if (any(is.na(param_idx))) {
    stop("sccomp says: could not map selected parameters to model coefficients.")
  }

  # `params_list` flattens the Stan posterior summaries into one entry per parameter,
  # carrying only the scalars (`intercept`, `slope`, or per-component variants)
  # required by downstream geometry, so the heavy `fit$summary` tibbles are not
  # re-scanned inside per-row mutate() / lapply() calls.
  if (!bimodal_flag) {
    params_list <- lapply(seq_along(param_names), function(a) {
      param_name <- param_names[a]
      idx <- param_idx[a]
      list(
        parameter = param_name,
        intercept = prec_intercept_1_summary$mean[idx],
        slope = prec_slope_1_summary$mean[idx]
      )
    })

    # Print the fitted lines so users can sanity-check the regression underlying
    # the visual adjustment (especially useful when comparing across fits).
    message("=== Single Model Parameters ===")
    for(i in 1:length(params_list)) {
      p <- params_list[[i]]
      message(sprintf("\n%s:", p$parameter))
      message(sprintf("  v = -(%.3f + %.3f \u00d7 c)", p$intercept, p$slope))
    }
    message("")

  } else {
    prec_intercept_2_summary <- fit$summary("prec_intercept_2")
    prec_slope_2_summary <- fit$summary("prec_slope_2")
    # `mix_p` is the posterior mixing weight of component 1; reported in the
    # caption so users can judge which component dominates the population.
    mix_p <- fit$summary("mix_p") |> pull(mean)

    params_list <- lapply(seq_along(param_names), function(a) {
      param_name <- param_names[a]
      idx <- param_idx[a]

      list(
        parameter = param_name,
        intercept_1 = prec_intercept_1_summary$mean[idx],
        slope_1 = prec_slope_1_summary$mean[idx],
        slope_2 = prec_slope_2_summary$mean[idx],
        intercept_2 = prec_intercept_2_summary$mean[idx]
      )
    })

    message("=== Bimodal Model Parameters ===")
    for(i in 1:length(params_list)) {
      p <- params_list[[i]]
      message(sprintf("\n%s:", p$parameter))
      message(sprintf("  Component 1: v = -(%.3f + %.3f \u00d7 c)", p$intercept_1, p$slope_1))
      message(sprintf("  Component 2: v = -(%.3f + %.3f \u00d7 c)", p$intercept_2, p$slope_2))
    }
    message("")
  }

  # The 2D plot juxtaposes two views of v_effect for each parameter:
  #   - "raw":      v_effect with the mean-variability association reintroduced,
  #                 i.e. the data as it would look *without* the model's correction;
  #                 the fitted regression line should pass through this cloud.
  #   - "adjusted": v_effect with the association removed (what sccomp uses for
  #                 inference). Cells off-axis here are the genuine outliers.
  # Since `.data$v_effect` is already the adjusted form (it comes from the
  # `alpha_normalised` draws), we *reverse the adjustment* to construct "raw"
  # and we use v_effect verbatim for "adjusted".
  if (!bimodal_flag) {
    .data_raw_list <- lapply(params_list, function(params) {
      .data %>%
        filter(parameter == params$parameter) %>%
        mutate(
          v_effect = v_effect - params$slope * c_effect,
          v_lower = v_lower - params$slope * c_effect,
          v_upper = v_upper - params$slope * c_effect,
          parameter = paste0(params$parameter, ", raw")
        )
    })
  } else {
    # Bimodal case: each cell type's "raw" position is generated by *one* of the
    # two mixture components, but the latent assignment is not directly observed.
    # We hard-assign by nearest predicted v_effect (smallest residual against
    # each component's regression line), then reverse the adjustment using only
    # that component's slope. This produces a visually coherent "raw" scatter
    # in which each point lies near its parent regression line rather than at
    # an arbitrary average of the two.
    .data_raw_list <- lapply(params_list, function(params) {
      .data %>%
        filter(parameter == params$parameter) %>%
        rowwise() %>%
        mutate(
          raw_v_comp1 = v_effect - params$slope_1 * c_effect,
          raw_v_comp2 = v_effect - params$slope_2 * c_effect,
          pred_comp1 = -(params$intercept_1 + params$slope_1 * c_effect),
          pred_comp2 = -(params$intercept_2 + params$slope_2 * c_effect),
          assigned_component = if_else(
            abs(raw_v_comp1 - pred_comp1) < abs(raw_v_comp2 - pred_comp2), 1, 2
          ),
          slope_to_use = if_else(assigned_component == 1, params$slope_1, params$slope_2),
          v_effect = v_effect - slope_to_use * c_effect,
          v_lower = v_lower - slope_to_use * c_effect,
          v_upper = v_upper - slope_to_use * c_effect,
          parameter = paste0(params$parameter, ", raw")
        ) %>%
        ungroup() %>%
        select(-raw_v_comp1, -raw_v_comp2, -pred_comp1, -pred_comp2, -slope_to_use)
    })
  }

  .data_raw <- bind_rows(.data_raw_list)

  # "adjusted" rows are the original v_effect, only the parameter label changes
  # so facet_wrap can place them in a separate panel.
  .data_adjusted_list <- lapply(params_list, function(params) {
    .data %>%
      filter(parameter == params$parameter) %>%
      mutate(parameter = paste0(params$parameter, ", adjusted"))
  })
  .data_adjusted <- bind_rows(.data_adjusted_list)

  .data_plot <- bind_rows(.data_raw, .data_adjusted)

  # Interleave "raw" before "adjusted" for each parameter so the facets read
  # left-to-right as "before -> after the model's adjustment". param_order is
  # also reused by the helper to ensure stable panel order across patchwork and
  # facet_wrap callers.
  param_order <- c()
  for(p in params_list) {
    param_order <- c(param_order, paste0(p$parameter, ", raw"), paste0(p$parameter, ", adjusted"))
  }

  .data_plot$parameter <- factor(.data_plot$parameter, levels = param_order)

  # Label only the top-3 cell groups by significance per panel, and only in
  # ", adjusted" panels because the "raw" panel is meant to show the underlying
  # association (labels there would clutter without adding analytical value).
  # Two sequential passes — first by c_FDR (abundance), then by v_FDR
  # (variability) — guarantee a cell is labeled if it is significant on
  # *either* axis, while preserving the abundance label when both apply.
  .data_plot <- .data_plot %>%
    filter(!is.na(v_effect)) %>%
    with_groups(
      parameter,
      ~ .x %>%
        arrange(c_FDR) %>%
        mutate(
          cell_type_label = if_else(
            row_number() <= 3 &
              c_FDR < significance_threshold &
              str_detect(parameter, ", adjusted$"),
            !!.cell_group,
            ""
          )
        )
    ) %>%
    with_groups(
      parameter,
      ~ .x %>%
        arrange(v_FDR) %>%
        mutate(
          cell_type_label = if_else(
            row_number() <= 3 &
              v_FDR < significance_threshold &
              str_detect(parameter, ", adjusted$") &
              cell_type_label == "",
            !!.cell_group,
            cell_type_label
          )
        )
    )

  # Encode significance in the errorbar aesthetics so the plot is readable
  # without consulting the underlying table. Significance is only meaningful on
  # ", adjusted" panels — raw panels are unconditionally desaturated to convey
  # "descriptive, not inferential". FDR (Stephens) uses red, pH0 (posterior
  # tail probability) uses blue, mirroring their distinct meanings.
  if (significance_statistic == "FDR") {
    color_c_aes <- aes(
      xmin = c_lower, xmax = c_upper,
      color = c_FDR < significance_threshold & str_detect(parameter, ", adjusted$"),
      alpha = c_FDR < significance_threshold & str_detect(parameter, ", adjusted$")
    )
    color_v_aes <- aes(
      ymin = v_lower, ymax = v_upper,
      color = v_FDR < significance_threshold & str_detect(parameter, ", adjusted$"),
      alpha = v_FDR < significance_threshold & str_detect(parameter, ", adjusted$")
    )
    color_scale <- scale_color_manual(values = c("#D3D3D3", "#E41A1C"))
    alpha_scale <- scale_alpha_manual(values = c(0.4, 1))
    legend_title <- "FDR < significance_threshold"
  } else {
    color_c_aes <- aes(
      xmin = c_lower, xmax = c_upper,
      color = c_pH0 < significance_threshold & str_detect(parameter, ", adjusted$"),
      alpha = c_pH0 < significance_threshold & str_detect(parameter, ", adjusted$")
    )
    color_v_aes <- aes(
      ymin = v_lower, ymax = v_upper,
      color = v_pH0 < significance_threshold & str_detect(parameter, ", adjusted$"),
      alpha = v_pH0 < significance_threshold & str_detect(parameter, ", adjusted$")
    )
    color_scale <- scale_color_manual(values = c("#D3D3D3", "#377EB8"))
    alpha_scale <- scale_alpha_manual(values = c(0.4, 1))
    legend_title <- "pH0 < significance_threshold"
  }

  # Prepare regression line data based on model type
  if (!bimodal_flag) {
    regression_data_all <- lapply(params_list, function(params) {
      raw_param <- paste0(params$parameter, ", raw")
      param_data <- .data_plot %>% filter(parameter == raw_param)
      if(nrow(param_data) == 0) return(NULL)

      c_range <- range(param_data$c_effect, na.rm = TRUE)
      c_seq <- seq(c_range[1], c_range[2], length.out = 100)
      v_pred <- -(params$intercept + params$slope * c_seq)

      data.frame(
        c_effect = c_seq,
        v_effect = v_pred,
        parameter = raw_param,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()

    adjusted_lines_all <- lapply(params_list, function(params) {
      adj_param <- paste0(params$parameter, ", adjusted")
      adj_data <- .data_plot %>% filter(parameter == adj_param)
      if(nrow(adj_data) == 0) return(NULL)

      c_range_adj <- range(adj_data$c_effect, na.rm = TRUE)
      mean_adj <- mean(adj_data$v_effect, na.rm = TRUE)

      data.frame(
        c_effect = c_range_adj,
        v_effect = rep(mean_adj, 2),
        parameter = adj_param,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()

  } else {
    # Bimodal: two regression lines and two horizontal references per panel.
    # The `component` column lets the plot helper colour them distinctly.
    regression_data_all <- lapply(params_list, function(params) {
      raw_param <- paste0(params$parameter, ", raw")
      param_data <- .data_plot %>% filter(parameter == raw_param)
      if(nrow(param_data) == 0) return(NULL)

      c_range <- range(param_data$c_effect, na.rm = TRUE)
      c_seq <- seq(c_range[1], c_range[2], length.out = 100)

      v_pred_1 <- -(params$intercept_1 + params$slope_1 * c_seq)
      v_pred_2 <- -(params$intercept_2 + params$slope_2 * c_seq)

      bind_rows(
        data.frame(
          c_effect = c_seq, v_effect = v_pred_1,
          parameter = raw_param, component = "Component 1"
        ),
        data.frame(
          c_effect = c_seq, v_effect = v_pred_2,
          parameter = raw_param, component = "Component 2"
        )
      )
    }) %>% bind_rows()

    adjusted_lines_all <- lapply(params_list, function(params) {
      adj_param <- paste0(params$parameter, ", adjusted")
      c_range_adj <- range(.data_plot$c_effect[.data_plot$parameter == adj_param], na.rm = TRUE)

      bind_rows(
        data.frame(
          c_effect = c_range_adj, v_effect = rep(-params$intercept_1, 2),
          parameter = adj_param, component = "Component 1"
        ),
        data.frame(
          c_effect = c_range_adj, v_effect = rep(-params$intercept_2, 2),
          parameter = adj_param, component = "Component 2"
        )
      )
    }) %>% bind_rows()
  }

  # Caption is reserved for FDR mode: users that opt into FDR-based inference
  # benefit from the explicit citation/interpretation reminder; pH0 users are
  # assumed to already understand the posterior-tail interpretation.
  if (significance_statistic == "FDR" && show_fdr_message) {
    if (!bimodal_flag) {
      caption_text <- paste(
        "Single model with parameter-specific slopes",
        "\nBlue line shows mean-variability relationship for each parameter",
        "\n'Adjusted' panels show variability after removing parameter-specific mean-variability association",
        "\nMarginal density shows posterior distribution of variability intercept parameters",
        "\nBayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)",
        sep = ""
      )
    } else {
      caption_text <- paste(
        sprintf("Bimodal model: mix_p = %.3f (Component 1 weight)", mix_p),
        "\nSolid blue line = Component 1 | Dashed orange line = Component 2",
        "\n'Adjusted' panels show variability after removing mean-variability association",
        "\nMarginal density shows posterior distribution of variability intercept parameters",
        "\nBayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)",
        sep = ""
      )
    }
  } else {
    caption_text <- NULL
  }

  # Both modes (with and without marginal density) produce a single faceted plot; the only difference is
  # whether ggside's y-side area layer is added on top. We rely on ggside
  # rather than patchwork-of-per-panel-plots so axis titles are not duplicated
  # and panel widths stay grid-aligned.
  build_2d_interval_plot <- function(plot_data, regression_data, adjusted_lines, density_data = NULL) {

    # Force the canonical factor order on every input: the data preparation
    # above does it on `.data_plot`, but `regression_data` / `adjusted_lines` /
    # `density_data` come from separate pipelines and must agree on level order
    # for facet_wrap to keep panels in sync.
    plot_data <- plot_data %>%
      mutate(parameter = factor(as.character(parameter), levels = param_order))

    if (!is.null(regression_data) && nrow(regression_data) > 0) {
      regression_data <- regression_data %>%
        mutate(parameter = factor(as.character(parameter), levels = param_order))
    }

    if (!is.null(adjusted_lines) && nrow(adjusted_lines) > 0) {
      adjusted_lines <- adjusted_lines %>%
        mutate(parameter = factor(as.character(parameter), levels = param_order))
    }

    if (!is.null(density_data) && nrow(density_data) > 0) {
      density_data <- density_data %>%
        mutate(parameter = factor(as.character(parameter), levels = param_order))
    }

    # Decision boundary at ± test_composition_above_logit_fold_change: cells
    # outside this box are the only ones a hypothesis test could call significant.
    p <- ggplot(plot_data, aes(c_effect, v_effect)) +
      geom_vline(
        xintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
        colour = "grey", linetype = "dashed", linewidth = 0.3
      ) +
      geom_hline(
        yintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
        colour = "grey", linetype = "dashed", linewidth = 0.3
      )

    # Regression / reference lines.
    # `inherit.aes = FALSE` shields these geoms from the top-level
    # aes(c_effect, v_effect) mapping so the lines use their own data verbatim
    # — important because their tibbles have a different row schema (no CI
    # columns, no cell_group, etc.) than `plot_data`.
    # In bimodal mode the two components get distinct visual treatment so the
    # eye can separate them without consulting the caption.
    if (!bimodal_flag) {
      if(!is.null(regression_data) && nrow(regression_data) > 0) {
        p <- p + geom_line(data = regression_data, mapping = aes(c_effect, v_effect),
                           color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
      }
      if(!is.null(adjusted_lines) && nrow(adjusted_lines) > 0) {
        p <- p + geom_line(data = adjusted_lines, mapping = aes(c_effect, v_effect),
                           color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
      }
    } else {
      if(!is.null(regression_data) && nrow(regression_data) > 0) {
        p <- p +
          geom_line(data = regression_data %>% filter(component == "Component 1"),
                    mapping = aes(c_effect, v_effect), color = "#0072B2",
                    linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
          geom_line(data = regression_data %>% filter(component == "Component 2"),
                    mapping = aes(c_effect, v_effect), color = "#D55E00",
                    linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
      }
      if(!is.null(adjusted_lines) && nrow(adjusted_lines) > 0) {
        p <- p +
          geom_line(data = adjusted_lines %>% filter(component == "Component 1"),
                    mapping = aes(c_effect, v_effect), color = "#0072B2",
                    linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
          geom_line(data = adjusted_lines %>% filter(component == "Component 2"),
                    mapping = aes(c_effect, v_effect), color = "#D55E00",
                    linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
      }
    }

    # Credible intervals. When omitted, we substitute invisible `geom_blank`
    # layers at the CI bounds so the plot's coordinate system is still trained
    # by the same data range — keeping axes identical whether errorbars are
    # drawn or not (otherwise the no-CI plot would zoom in onto just the points).
    if (!omit_ci) {
      p <- p +
        geom_errorbar(color_c_aes, linewidth = 0.2) +
        geom_errorbar(color_v_aes, linewidth = 0.2) +
        color_scale +
        alpha_scale +
        guides(color = guide_legend(title = legend_title), alpha = "none")
    } else {
      p <- p +
        geom_blank(aes(x = c_lower, y = v_lower)) +
        geom_blank(aes(x = c_upper, y = v_upper))
    }

    # Points are drawn *after* the errorbars/lines so they sit on top and
    # remain readable when CIs overlap; labels go last for highest z-order.
    # Note `geom_text_repel` uses `-v_effect` (y-axis flip) because labels are
    # placed relative to the inverted visual reading where higher v means
    # higher dispersion — this matches users' mental model of "outlier upward".
    p <- p +
      geom_point(size = 0.2) +
      geom_text_repel(
        aes(c_effect, -v_effect, label = cell_type_label),
        size = 2.5,
        data = plot_data %>% filter(cell_type_label != ""),
        max.overlaps = 20
      ) +
      xlab("c_effect (Abundance effect)") +
      ylab("v_effect (Variability effect)") +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank()
      )

    # Marginal posterior densities, drawn via ggside so they participate in the
    # same facet system as the main scatter (no patchwork → no duplicated axis
    # titles, no panel-width misalignment). We pass raw posterior draws (one
    # row per draw, with `parameter` set to the panel they belong to) and let
    # `geom_ysidedensity` compute the kernel density per facet. ggside lays
    # the result on the y-side, mapping the data's `y` aesthetic to the main
    # panel's y-axis. `scales = "free"` propagates so each panel's density
    # rescales to fit its own y-range.
    if (!is.null(density_data) && nrow(density_data) > 0) {
      density_aes <- if (bimodal_flag) {
        aes(y = v_value, fill = component)
      } else {
        aes(y = v_value)
      }
      p <- p +
        geom_ysidedensity(
          data = density_data,
          mapping = density_aes,
          alpha = 0.5,
          position = "identity",
          inherit.aes = FALSE
        ) +
        # Drop ggside's own axis decorations: the density value axis is
        # uninformative on its own (relative scale per panel) and would
        # otherwise clutter the strip.
        theme_ggside_void() +
        scale_ysidex_continuous(expand = c(0, 0))

      if (bimodal_flag) {
        p <- p + scale_fill_manual(values = c("Component 1" = "#0072B2", "Component 2" = "#D55E00"))
      }
    }

    # `scales = "free"` is intentional: raw and adjusted panels live on
    # different natural scales (raw v_effect can span much wider than adjusted
    # residuals), and the same applies across parameters.
    p + facet_wrap(~ parameter, scales = "free", ncol = 2)
  }

  # Posterior draws of `-prec_intercept_*[param_idx]` for every parameter,
  # replicated across both that parameter's facets (", raw" and ", adjusted").
  # The intercept density is a property of the parameter, not of the visual
  # mode, so showing it on every facet is both honest and what ggside requires
  # to attach a side layer per panel (omitting a facet's density would leave
  # that side panel empty, breaking visual rhythm).
  # The sign flip mirrors the parameterisation v = -(prec_intercept + slope·c)
  # so the density lies on the same axis orientation as v_effect. We pass raw
  # draws (not pre-computed densities) because ggside's density geom computes
  # the kernel per facet itself.
  density_data_all <- if (add_marginal_density) {
    bind_rows(lapply(seq_along(param_names), function(i) {
      param_name <- param_names[i]
      idx <- param_idx[i]

      component_vars <- if (bimodal_flag) {
        c("Component 1" = paste0("prec_intercept_1[", idx, "]"),
          "Component 2" = paste0("prec_intercept_2[", idx, "]"))
      } else {
        c("Component 1" = paste0("prec_intercept_1[", idx, "]"))
      }

      bind_rows(lapply(seq_along(component_vars), function(j) {
        var_name <- component_vars[[j]]
        values <- as.vector(fit$draws(variables = var_name, format = "draws_df")[[var_name]])
        draws_per_param <- tibble(
          component = names(component_vars)[j],
          v_value = -values
        )
        bind_rows(
          draws_per_param %>% mutate(parameter = paste0(param_name, ", raw")),
          draws_per_param %>% mutate(parameter = paste0(param_name, ", adjusted"))
        )
      }))
    }))
  } else {
    NULL
  }

  p <- build_2d_interval_plot(
    .data_plot,
    regression_data_all,
    adjusted_lines_all,
    density_data = density_data_all
  )

  if (!is.null(caption_text)) {
    p <- p +
      theme(plot.caption = element_text(hjust = 0, size = 9)) +
      labs(caption = caption_text)
  }

  p
}

#' Soft-deprecated aliases (call [sccomp_plot_intervals_1D()] / [sccomp_plot_intervals_2D()] instead).
#'
#' @importFrom lifecycle deprecate_soft
#' @export
#' @noRd
plot_1D_intervals <- function(...) {
  deprecate_soft("2.1.29", "plot_1D_intervals()", "sccomp_plot_intervals_1D()")
  sccomp_plot_intervals_1D(...)
}

#' @export
#' @noRd
plot_2D_intervals <- function(...) {
  deprecate_soft("2.1.29", "plot_2D_intervals()", "sccomp_plot_intervals_2D()")
  sccomp_plot_intervals_2D(...)
}