% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{sccomp_remove_unwanted_variation}
\alias{sccomp_remove_unwanted_variation}
\title{sccomp_remove_unwanted_variation}
\usage{
sccomp_remove_unwanted_variation(
  .data,
  formula_composition = ~1,
  formula_variability = NULL,
  cores = detectCores()
)
}
\arguments{
\item{.data}{A tibble. The result of sccomp_estimate.}

\item{formula_composition}{A formula. The formula describing the model for differential abundance, for example ~treatment. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.}

\item{formula_variability}{A formula. The formula describing the model for differential variability, for example ~treatment. In most cases, if differentially variability is of interest, the formula should only include the factor of interest as a large anount of data is needed to define variability depending to each factors. This formula can be a sub-formula of your estimated model; in this case all other factor will be factored out.#' @param cores Integer, the number of cores to be used for parallel calculations.}

\item{cores}{Integer, the number of cores to be used for parallel calculations.}
}
\value{
A tibble (\code{tbl}) with the following columns:
\itemize{
\item \strong{sample} - A character column representing the sample name for which data was adjusted.
\item \strong{cell_group} - A character column representing the cell group being tested.
\item \strong{adjusted_proportion} - A numeric column representing the adjusted proportion after removing unwanted variation.
\item \strong{adjusted_counts} - A numeric column representing the adjusted counts after removing unwanted variation.
\item \strong{logit_residuals} - A numeric column representing the logit residuals calculated after adjustment.
}
}
\description{
This function uses the model to remove unwanted variation from a dataset using the estimated of the model. For example if you fit your data with this formula \code{~ factor_1 + factor_2} and use this formula to remove unwanted variation \code{~ factor_1}, the \code{factor_2} will be factored out.
}
\examples{

message("Use the following example after having installed install.packages(\"cmdstanr\", repos = c(\"https://stan-dev.r-universe.dev/\", getOption(\"repos\")))")

\donttest{
  if (instantiate::stan_cmdstan_exists()) {
    data("counts_obj")

    estimates = sccomp_estimate(
      counts_obj,
      ~ type, ~1, sample, cell_group, count,
      cores = 1
    ) |>
    sccomp_remove_unwanted_variation()
  }
}

}
