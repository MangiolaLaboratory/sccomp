
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
#' @importFrom forcats fct_reorder
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
#'
#' @importFrom dplyr filter arrange mutate if_else row_number bind_rows distinct slice pull with_groups
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_errorbar geom_point geom_line geom_area aes facet_wrap theme_bw theme labs guides guide_legend scale_color_manual scale_alpha_manual scale_fill_manual scale_y_continuous coord_flip theme_void element_rect element_text margin
#' @importFrom ggrepel geom_text_repel
#' @importFrom stringr str_detect
#' @importFrom patchwork plot_annotation wrap_plots plot_layout
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
    add_marginal_density = TRUE
) {

  significance_statistic <- match.arg(significance_statistic)

  # Define variables to avoid CRAN NOTES
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

  .cell_group <- attr(.data, ".cell_group")

  # Check if test has been done
  if(.data |> select(ends_with("FDR")) |> ncol() == 0)
    stop("sccomp says: you need to run sccomp_test() first.")

  .data <- subset_results_by_factor(.data, factor, keep_intercept = TRUE)

  # Extract fitted model and mean-variability regression coefficients
  fit <- attr(.data, "fit")
  prec_intercept_1_summary <- fit$summary("prec_intercept_1")
  prec_slope_1_summary <- fit$summary("prec_slope_1")
  prec_intercept_2_summary <- tryCatch(
    fit$summary("prec_intercept_2"),
    error = function(e) tibble()
  )
  prec_slope_2_summary <- tryCatch(
    fit$summary("prec_slope_2"),
    error = function(e) tibble()
  )

  param_names <- .data |>
    filter(!is.na(v_effect)) |>
    distinct(parameter) |>
    pull(parameter)

  param_idx <- match(param_names, colnames(attr(.data, "model_input")$Xa))
  if (any(is.na(param_idx))) {
    stop("sccomp says: could not map selected parameters to model coefficients.")
  }

  # Derive model type from stored model metadata
  bimodal_flag <- attr(.data, "model_input")$bimodal_mean_variability_association
  if (is.null(bimodal_flag)) {
    stop("sccomp says: cannot infer model type because `bimodal_mean_variability_association` is missing from model metadata.")
  }
  bimodal_flag <- isTRUE(as.logical(bimodal_flag))

  # Extract parameters based on model type
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

    cat("=== Single Model Parameters ===\n")
    for(i in 1:length(params_list)) {
      p <- params_list[[i]]
      cat(sprintf("\n%s:\n", p$parameter))
      cat(sprintf("  v = -(%.3f + %.3f × c)\n", p$intercept, p$slope))
    }
    cat("\n")

  } else {
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

    cat("=== Bimodal Model Parameters ===\n")
    for(i in 1:length(params_list)) {
      p <- params_list[[i]]
      cat(sprintf("\n%s:\n", p$parameter))
      cat(sprintf("  Component 1: v = -(%.3f + %.3f × c)\n", p$intercept_1, p$slope_1))
      cat(sprintf("  Component 2: v = -(%.3f + %.3f × c)\n", p$intercept_2, p$slope_2))
    }
    cat("\n")
  }

  # v_effect already comes from alpha_normalised (computed in R from draws)
  # "raw" panel: ADD BACK entanglement to show raw alpha
  # "adjusted" panel: USE v_effect AS-IS

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

  # Adjusted panel: v_effect as-is (already from alpha_normalised)
  .data_adjusted_list <- lapply(params_list, function(params) {
    .data %>%
      filter(parameter == params$parameter) %>%
      mutate(parameter = paste0(params$parameter, ", adjusted"))
  })
  .data_adjusted <- bind_rows(.data_adjusted_list)

  .data_plot <- bind_rows(.data_raw, .data_adjusted)

  # Set parameter factor levels
  param_order <- c()
  for(p in params_list) {
    param_order <- c(param_order, paste0(p$parameter, ", raw"), paste0(p$parameter, ", adjusted"))
  }

  .data_plot$parameter <- factor(.data_plot$parameter, levels = param_order)

  # Add labels for significant cell groups
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
            !!sym(.cell_group),
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
            !!sym(.cell_group),
            cell_type_label
          )
        )
    )

  # Choose color aesthetics based on significance statistic
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

  # Add caption based on model type
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

  # Add marginal density plots if requested
  if (add_marginal_density) {

    plot_list <- lapply(param_order, function(param) {

      param_data <- .data_plot %>% filter(parameter == param)
      if(nrow(param_data) == 0) return(NULL)

      # Create main plot
      p_param <- ggplot(param_data, aes(c_effect, v_effect)) +
        geom_vline(
          xintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
          colour = "grey", linetype = "dashed", linewidth = 0.3
        ) +
        geom_hline(
          yintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
          colour = "grey", linetype = "dashed", linewidth = 0.3
        )

      # Add regression lines
      if (!bimodal_flag) {
        reg_data <- regression_data_all %>% filter(parameter == param)
        if(!is.null(reg_data) && nrow(reg_data) > 0) {
          p_param <- p_param +
            geom_line(data = reg_data, mapping = aes(c_effect, v_effect),
                      color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
        }

        adj_line <- adjusted_lines_all %>% filter(parameter == param)
        if(!is.null(adj_line) && nrow(adj_line) > 0) {
          p_param <- p_param +
            geom_line(data = adj_line, mapping = aes(c_effect, v_effect),
                      color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
        }

      } else {
        reg_data <- regression_data_all %>% filter(parameter == param)
        if(!is.null(reg_data) && nrow(reg_data) > 0) {
          p_param <- p_param +
            geom_line(data = reg_data %>% filter(component == "Component 1"),
                      mapping = aes(c_effect, v_effect), color = "#0072B2",
                      linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
            geom_line(data = reg_data %>% filter(component == "Component 2"),
                      mapping = aes(c_effect, v_effect), color = "#D55E00",
                      linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
        }

        adj_line <- adjusted_lines_all %>% filter(parameter == param)
        if(!is.null(adj_line) && nrow(adj_line) > 0) {
          p_param <- p_param +
            geom_line(data = adj_line %>% filter(component == "Component 1"),
                      mapping = aes(c_effect, v_effect), color = "#0072B2",
                      linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
            geom_line(data = adj_line %>% filter(component == "Component 2"),
                      mapping = aes(c_effect, v_effect), color = "#D55E00",
                      linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
        }
      }

      # Add error bars, points, and labels
      p_param <- p_param +
        geom_errorbar(color_c_aes, linewidth = 0.2) +
        geom_errorbar(color_v_aes, linewidth = 0.2) +
        geom_point(size = 0.2) +
        geom_text_repel(
          aes(c_effect, -v_effect, label = cell_type_label),
          size = 2.5,
          data = param_data %>% filter(cell_type_label != ""),
          max.overlaps = 20
        ) +
        color_scale +
        alpha_scale +
        xlab("c_effect (Abundance effect)") +
        ylab("v_effect (Variability effect)") +
        ggtitle(param) +
        theme_bw() +
        theme(
          legend.position = "bottom",
          strip.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank()
        ) +
        guides(color = guide_legend(title = legend_title), alpha = "none")

      # Add marginal density for adjusted panels (not Intercept)
      if (str_detect(param, ", adjusted$") && !str_detect(param, "Intercept")) {

        if (!bimodal_flag) {
          param_idx <- which(sapply(params_list, function(p) paste0(p$parameter, ", adjusted") == param))

          if (length(param_idx) > 0) {
            intercept_var_name <- paste0("prec_intercept_1[", param_idx, "]")

            tryCatch({
              intercept_draws <- fit$draws(variables = intercept_var_name, format = "draws_df")
              intercept_values <- as.vector(intercept_draws[[intercept_var_name]])

              dens <- density(-intercept_values, na.rm = TRUE)
              dens_df <- data.frame(x = dens$x, y = dens$y)


              y_range <- ggplot_build(p_param)$layout$panel_params[[1]]$y.range

              p_density <- ggplot(dens_df, aes(x = x, y = y)) +
                geom_area(alpha = 0.5, position = "identity") +
                geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
                coord_flip(xlim = y_range) +
                scale_y_continuous(expand = c(0, 0)) +
                xlab("Posterior Probability") +
                theme_void() +
                theme(
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 6),
                  axis.title.y = element_text(angle = 90, size = 7, vjust = 0.5)
                )

              p_combined <- p_param + p_density +
                plot_layout(ncol = 2, widths = c(5, 0.6), guides = "collect") &
                theme(legend.position = "bottom")

              return(p_combined)
            }, error = function(e) {
              warning(sprintf("Could not extract intercept draws for %s: %s", intercept_var_name, e$message))
              return(p_param)
            })
          }

        } else {
          param_idx <- which(sapply(params_list, function(p) paste0(p$parameter, ", adjusted") == param))

          if (length(param_idx) > 0) {
            intercept1_var_name <- paste0("prec_intercept_1[", param_idx, "]")
            intercept2_var_name <- paste0("prec_intercept_2[", param_idx, "]")

            tryCatch({
              intercept1_draws <- fit$draws(variables = intercept1_var_name, format = "draws_df")
              intercept2_draws <- fit$draws(variables = intercept2_var_name, format = "draws_df")

              intercept1_values <- as.vector(intercept1_draws[[intercept1_var_name]])
              intercept2_values <- as.vector(intercept2_draws[[intercept2_var_name]])

              dens1 <- density(-intercept1_values, na.rm = TRUE)
              dens2 <- density(-intercept2_values, na.rm = TRUE)

              dens_df <- bind_rows(
                data.frame(x = dens1$x, y = dens1$y, component = "Component 1"),
                data.frame(x = dens2$x, y = dens2$y, component = "Component 2")
              )


              y_range <- ggplot_build(p_param)$layout$panel_params[[1]]$y.range

              p_density <- ggplot(dens_df, aes(x = x, y = y, fill = component)) +
                geom_area(alpha = 0.5, position = "identity") +
                geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
                scale_fill_manual(values = c("Component 1" = "#0072B2", "Component 2" = "#D55E00")) +
                coord_flip(xlim = y_range) +
                scale_y_continuous(expand = c(0, 0)) +
                xlab("Posterior Probability") +
                theme_void() +
                theme(
                  plot.margin = margin(t = 0, r = 0, b = 0, l = 5),
                  legend.position = "none",
                  axis.title.y = element_text(angle = 90, size = 7, vjust = 0.5)
                )

              p_combined <- p_param + p_density +
                plot_layout(ncol = 2, widths = c(5, 0.6), guides = "collect") &
                theme(legend.position = "bottom")

              return(p_combined)
            }, error = function(e) {
              warning(sprintf("Could not extract intercept draws: %s", e$message))
              return(p_param)
            })
          }
        }
      }

      return(p_param)
    })


    plot_list <- plot_list[!sapply(plot_list, is.null)]
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)

    if (!is.null(caption_text)) {
      combined_plot <- combined_plot +
        plot_annotation(
          caption = caption_text,
          theme = theme(plot.caption = element_text(hjust = 0, size = 9))
        )
    }

    return(combined_plot)

  } else {
    # Return faceted plot without marginal densities
    p <- ggplot(.data_plot, aes(c_effect, v_effect)) +
      geom_vline(
        xintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
        colour = "grey", linetype = "dashed", linewidth = 0.3
      ) +
      geom_hline(
        yintercept = c(-test_composition_above_logit_fold_change, test_composition_above_logit_fold_change),
        colour = "grey", linetype = "dashed", linewidth = 0.3
      )

    # Add regression lines
    if (!bimodal_flag) {
      if(!is.null(regression_data_all) && nrow(regression_data_all) > 0) {
        p <- p + geom_line(data = regression_data_all, mapping = aes(c_effect, v_effect),
                           color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
      }
      if(!is.null(adjusted_lines_all) && nrow(adjusted_lines_all) > 0) {
        p <- p + geom_line(data = adjusted_lines_all, mapping = aes(c_effect, v_effect),
                           color = "#0072B2", linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE)
      }
    } else {
      if(!is.null(regression_data_all) && nrow(regression_data_all) > 0) {
        p <- p +
          geom_line(data = regression_data_all %>% filter(component == "Component 1"),
                    mapping = aes(c_effect, v_effect), color = "#0072B2",
                    linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
          geom_line(data = regression_data_all %>% filter(component == "Component 2"),
                    mapping = aes(c_effect, v_effect), color = "#D55E00",
                    linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
      }
      if(!is.null(adjusted_lines_all) && nrow(adjusted_lines_all) > 0) {
        p <- p +
          geom_line(data = adjusted_lines_all %>% filter(component == "Component 1"),
                    mapping = aes(c_effect, v_effect), color = "#0072B2",
                    linewidth = 0.5, alpha = 0.8, inherit.aes = FALSE) +
          geom_line(data = adjusted_lines_all %>% filter(component == "Component 2"),
                    mapping = aes(c_effect, v_effect), color = "#D55E00",
                    linewidth = 0.5, alpha = 0.8, linetype = "dashed", inherit.aes = FALSE)
      }
    }

    p <- p +
      geom_errorbar(color_c_aes, linewidth = 0.2) +
      geom_errorbar(color_v_aes, linewidth = 0.2) +
      geom_point(size = 0.2) +
      geom_text_repel(
        aes(c_effect, -v_effect, label = cell_type_label),
        size = 2.5,
        data = .data_plot %>% filter(cell_type_label != ""),
        max.overlaps = 20
      ) +
      color_scale +
      alpha_scale +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      xlab("c_effect (Abundance effect)") +
      ylab("v_effect (Variability effect)") +
      theme_bw() +
      theme(
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(title = legend_title), alpha = "none")

    if (!is.null(caption_text)) {
      p <- p +
        theme(plot.caption = element_text(hjust = 0, size = 9)) +
        labs(caption = caption_text)
    }

    return(p)
  }
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