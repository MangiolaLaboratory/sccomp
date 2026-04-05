#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(utils)
})

parse_args <- function(args) {
  out <- list(
    repo = NULL,
    label = NULL,
    out = NULL,
    base_lib = NULL,
    n_datasets = 8L,
    n_samples_per_group = 20L,
    inference_method = "pathfinder",
    max_sampling_iterations = 4000L
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    value <- if (i < length(args)) args[[i + 1L]] else NULL

    if (key == "--repo") out$repo <- value
    if (key == "--label") out$label <- value
    if (key == "--out") out$out <- value
    if (key == "--base-lib") out$base_lib <- value
    if (key == "--n-datasets") out$n_datasets <- as.integer(value)
    if (key == "--n-samples-per-group") out$n_samples_per_group <- as.integer(value)
    if (key == "--inference-method") out$inference_method <- value
    if (key == "--max-sampling-iterations") out$max_sampling_iterations <- as.integer(value)

    i <- i + 2L
  }

  if (is.null(out$repo) || is.null(out$label) || is.null(out$out)) {
    stop(
      "Required arguments: --repo <path> --label <version_label> --out <output_csv>\n",
      "Optional: --n-datasets <int> --n-samples-per-group <int> ",
      "--inference-method <pathfinder|hmc|variational> --max-sampling-iterations <int>"
    )
  }

  out
}

softmax <- function(x) {
  z <- x - max(x)
  exp(z) / sum(exp(z))
}

simulate_composition_data <- function(seed, n_samples_per_group = 20L) {
  set.seed(seed)

  cell_groups <- paste0("cg", seq_len(10))
  p_control <- c(0.26, 0.18, 0.15, 0.11, 0.08, 0.07, 0.05, 0.04, 0.035, 0.025)
  log_fc <- c(0.8, -0.6, 0.0, 0.0, 0.4, -0.3, 0.0, 0.0, 0.0, 0.0)
  p_treated <- softmax(log(p_control) + log_fc)

  make_one_sample <- function(sample_id, condition) {
    total_count <- rpois(1, lambda = 2000) + 400L
    probs <- if (condition == "control") p_control else p_treated
    counts <- as.integer(rmultinom(1, size = total_count, prob = probs))

    data.frame(
      sample = sample_id,
      condition = condition,
      cell_group = cell_groups,
      count = counts,
      stringsAsFactors = FALSE
    )
  }

  control_ids <- paste0("ctrl_", seq_len(n_samples_per_group))
  treated_ids <- paste0("trt_", seq_len(n_samples_per_group))

  control_df <- do.call(rbind, lapply(control_ids, make_one_sample, condition = "control"))
  treated_df <- do.call(rbind, lapply(treated_ids, make_one_sample, condition = "treated"))
  out <- rbind(control_df, treated_df)
  out$condition <- factor(out$condition, levels = c("control", "treated"))
  out
}

run_single_dataset <- function(data, dataset_seed, fit_seed, inference_method, max_sampling_iterations, use_new_api) {
  if (use_new_api) {
    fit <- sccomp::sccomp_estimate(
      data,
      formula_composition = ~ condition,
      formula_variability = ~ 1,
      sample = "sample",
      cell_group = "cell_group",
      abundance = "count",
      cores = 1,
      verbose = FALSE,
      inference_method = inference_method,
      mcmc_seed = fit_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = FALSE
    )
  } else {
    # Keep symbols visible for linters while passing tidy-eval columns to old API.
    sample <- NULL
    cell_group <- NULL
    count <- NULL

    fit <- sccomp::sccomp_estimate(
      data,
      formula_composition = ~ condition,
      formula_variability = ~ 1,
      .sample = sample,
      .cell_group = cell_group,
      .count = count,
      cores = 1,
      verbose = FALSE,
      inference_method = inference_method,
      mcmc_seed = fit_seed,
      max_sampling_iterations = max_sampling_iterations,
      pass_fit = FALSE
    )
  }

  tested <- sccomp::sccomp_test(fit, pass_fit = FALSE)
  composition_rows <- tested[grepl("condition", tested$parameter, fixed = TRUE), , drop = FALSE]

  data.frame(
    dataset_seed = dataset_seed,
    fit_seed = fit_seed,
    n_composition_tests = nrow(composition_rows),
    n_significant_fdr_0_05 = sum(composition_rows$c_FDR < 0.05, na.rm = TRUE),
    n_significant_fdr_0_10 = sum(composition_rows$c_FDR < 0.10, na.rm = TRUE),
    min_c_FDR = suppressWarnings(min(composition_rows$c_FDR, na.rm = TRUE)),
    median_c_FDR = suppressWarnings(stats::median(composition_rows$c_FDR, na.rm = TRUE)),
    mean_abs_effect = suppressWarnings(mean(abs(composition_rows$c_effect), na.rm = TRUE))
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
repo <- normalizePath(args$repo, mustWork = TRUE)
label <- args$label
out_csv <- args$out

message("Installing sccomp from: ", repo)
lib_path <- tempfile("sccomp-lib-")
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)

install_status <- system2(
  command = file.path(R.home("bin"), "R"),
  args = c("CMD", "INSTALL", "-l", shQuote(lib_path), shQuote(repo)),
  stdout = TRUE,
  stderr = TRUE
)
if (!is.null(attr(install_status, "status")) && attr(install_status, "status") != 0) {
  stop("Package installation failed for ", repo, "\n", paste(install_status, collapse = "\n"))
}

.libPaths(c(lib_path, .libPaths()))
if (!is.null(args$base_lib)) {
  .libPaths(c(lib_path, normalizePath(args$base_lib, mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages(library(sccomp))
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  message("cmdstanr version in use: ", as.character(utils::packageVersion("cmdstanr")))
}

# Avoid cross-version cached model collisions (important for old/new comparisons).
if ("clear_stan_model_cache" %in% getNamespaceExports("sccomp")) {
  try(sccomp::clear_stan_model_cache(), silent = TRUE)
}
cache_dir <- path.expand("~/.sccomp_models")
if (dir.exists(cache_dir)) {
  unlink(cache_dir, recursive = TRUE, force = TRUE)
}

use_new_api <- "sample" %in% names(formals(sccomp::sccomp_estimate))

dataset_seeds <- seq_len(args$n_datasets) + 1000L
fit_seeds <- seq_len(args$n_datasets) + 5000L

rows <- vector("list", length(dataset_seeds))
for (i in seq_along(dataset_seeds)) {
  d_seed <- dataset_seeds[[i]]
  f_seed <- fit_seeds[[i]]
  message(sprintf("[%s] dataset %d/%d (data seed=%d, fit seed=%d)",
                  label, i, length(dataset_seeds), d_seed, f_seed))

  dat <- simulate_composition_data(d_seed, n_samples_per_group = args$n_samples_per_group)
  rows[[i]] <- run_single_dataset(
    dat,
    dataset_seed = d_seed,
    fit_seed = f_seed,
    inference_method = args$inference_method,
    max_sampling_iterations = args$max_sampling_iterations,
    use_new_api = use_new_api
  )
}

res <- do.call(rbind, rows)
res$version_label <- label
res$inference_method <- args$inference_method
res$n_samples_per_group <- args$n_samples_per_group
res$max_sampling_iterations <- args$max_sampling_iterations

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(res, out_csv, row.names = FALSE)

summary_row <- data.frame(
  version_label = label,
  inference_method = args$inference_method,
  mean_n_significant_fdr_0_05 = mean(res$n_significant_fdr_0_05),
  median_n_significant_fdr_0_05 = stats::median(res$n_significant_fdr_0_05),
  mean_min_c_FDR = mean(res$min_c_FDR),
  mean_median_c_FDR = mean(res$median_c_FDR),
  mean_abs_effect = mean(res$mean_abs_effect)
)

message("\nBenchmark summary:")
print(summary_row, row.names = FALSE)
message("\nSaved results to: ", out_csv)
