% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{plot.sccomp_tbl}
\alias{plot.sccomp_tbl}
\title{plot}
\usage{
\method{plot}{sccomp_tbl}(.data, significance_threshold = 0.025)
}
\arguments{
\item{.data}{A tibble including a cell_group name column | sample name column | read counts column | factor columns | Pvalue column | a significance column}

\item{significance_threshold}{A real. FDR threshold for labelling significant cell-groups.}
}
\value{
A \code{ggplot}
}
\description{
This function plots a summary of the results of the model.
}
\examples{

data("counts_obj")

estimate =
  sccomp_estimate(
  counts_obj ,
   ~ type, ~1, sample, cell_group, count,
    approximate_posterior_inference = "all",
    cores = 1
  )

# estimate |> plot()

}