% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{sccomp_estimate}
\alias{sccomp_estimate}
\title{Main Function for SCCOMP Estimate}
\usage{
sccomp_estimate(
  .data,
  formula_composition = ~1,
  formula_variability = ~1,
  .sample,
  .cell_group,
  .abundance = NULL,
  cores = detectCores(),
  bimodal_mean_variability_association = FALSE,
  percent_false_positive = 5,
  inference_method = "pathfinder",
  prior_mean = list(intercept = c(0, 1), coefficients = c(0, 1)),
  prior_overdispersion_mean_association = list(intercept = c(5, 2), slope = c(0, 0.6),
    standard_deviation = c(10, 20)),
  .sample_cell_group_pairs_to_exclude = NULL,
  output_directory = "sccomp_draws_files",
  verbose = TRUE,
  enable_loo = FALSE,
  noise_model = "multi_beta_binomial",
  exclude_priors = FALSE,
  use_data = TRUE,
  mcmc_seed = sample(1e+05, 1),
  max_sampling_iterations = 20000,
  pass_fit = TRUE,
  sig_figs = 9,
  ...,
  .count = NULL,
  approximate_posterior_inference = NULL,
  variational_inference = NULL
)
}
\arguments{
\item{.data}{A tibble including cell_group name column, sample name column,
abundance column (counts or proportions), and factor columns.}

\item{formula_composition}{A formula describing the model for differential abundance.}

\item{formula_variability}{A formula describing the model for differential variability.}

\item{.sample}{A column name as a symbol for the sample identifier.}

\item{.cell_group}{A column name as a symbol for the cell-group identifier.}

\item{.abundance}{A column name as a symbol for the cell-group abundance, which can be counts (> 0) or proportions (between 0 and 1, summing to 1 across \code{.cell_group}).}

\item{cores}{Number of cores to use for parallel calculations.}

\item{bimodal_mean_variability_association}{Logical, whether to model mean-variability as bimodal.}

\item{percent_false_positive}{A real number between 0 and 100 for outlier identification.}

\item{inference_method}{Character string specifying the inference method to use ('pathfinder', 'hmc', or 'variational').}

\item{prior_mean}{A list specifying prior knowledge about the mean distribution, including intercept and coefficients.}

\item{prior_overdispersion_mean_association}{A list specifying prior knowledge about mean/variability association.}

\item{.sample_cell_group_pairs_to_exclude}{A column name indicating sample/cell-group pairs to exclude.}

\item{output_directory}{A character string specifying the output directory for Stan draws.}

\item{verbose}{Logical, whether to print progression details.}

\item{enable_loo}{Logical, whether to enable model comparison using the LOO package.}

\item{noise_model}{A character string specifying the noise model (e.g., 'multi_beta_binomial').}

\item{exclude_priors}{Logical, whether to run a prior-free model.}

\item{use_data}{Logical, whether to run the model data-free.}

\item{mcmc_seed}{An integer seed for MCMC reproducibility.}

\item{max_sampling_iterations}{Integer to limit the maximum number of iterations for large datasets.}

\item{pass_fit}{Logical, whether to include the Stan fit as an attribute in the output.}

\item{sig_figs}{Number of significant figures to use for Stan model output. Default is 9.}

\item{...}{Additional arguments passed to the \code{cmdstanr::sample} function.}

\item{.count}{DEPRECATED. Use .abundance instead.}

\item{approximate_posterior_inference}{DEPRECATED. Use inference_method instead.}

\item{variational_inference}{DEPRECATED. Use inference_method instead.}
}
\value{
A tibble (\code{tbl}) with the following columns:
\itemize{
\item cell_group - The cell groups being tested.
\item parameter - The parameter being estimated from the design matrix described by the input \code{formula_composition} and \code{formula_variability}.
\item factor - The covariate factor in the formula, if applicable (e.g., not present for Intercept or contrasts).
\item c_lower - Lower (2.5\%) quantile of the posterior distribution for a composition (c) parameter.
\item c_effect - Mean of the posterior distribution for a composition (c) parameter.
\item c_upper - Upper (97.5\%) quantile of the posterior distribution for a composition (c) parameter.
\item c_pH0 - Probability of the null hypothesis (no difference) for a composition (c). This is not a p-value.
\item c_FDR - False-discovery rate of the null hypothesis for a composition (c).
\item c_n_eff - Effective sample size for a composition (c) parameter.
\item c_R_k_hat - R statistic for a composition (c) parameter, should be within 0.05 of 1.0.
\item v_lower - Lower (2.5\%) quantile of the posterior distribution for a variability (v) parameter.
\item v_effect - Mean of the posterior distribution for a variability (v) parameter.
\item v_upper - Upper (97.5\%) quantile of the posterior distribution for a variability (v) parameter.
\item v_pH0 - Probability of the null hypothesis for a variability (v).
\item v_FDR - False-discovery rate of the null hypothesis for a variability (v).
\item v_n_eff - Effective sample size for a variability (v) parameter.
\item v_R_k_hat - R statistic for a variability (v) parameter.
\item count_data - Nested input count data.
}
}
\description{
The \code{sccomp_estimate} function performs linear modeling on a table of cell counts or proportions,
which includes a cell-group identifier, sample identifier, abundance (counts or proportions), and factors
(continuous or discrete). The user can define a linear model using an R formula,
where the first factor is the factor of interest. Alternatively, \code{sccomp} accepts
single-cell data containers (e.g., Seurat, SingleCellExperiment, cell metadata, or
group-size) and derives the count data from cell metadata.
}
\examples{

print("cmdstanr is needed to run this example.")
# Note: Before running the example, ensure that the 'cmdstanr' package is installed:
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev/", getOption("repos")))

\donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimate <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      sample,
      cell_group,
      count,
      cores = 1
    )
    
   # Note! 
   # If counts are available, do not use proportion.
   # Using proportion ignores the high uncertainty of low counts
   
   estimate_proportion <- sccomp_estimate(
      counts_obj,
      ~ type,
      ~1,
      sample,
      cell_group,
      proportion,
      cores = 1
    )
    
  }
}

}
