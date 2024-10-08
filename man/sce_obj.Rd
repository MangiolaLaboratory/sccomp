% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{sce_obj}
\alias{sce_obj}
\title{sce_obj}
\format{
A \code{SingleCellExperiment} object with the following structure:
\itemize{
\item \strong{assays}: Two assays: counts (raw RNA counts) and logcounts (log-transformed counts).
\item \strong{rowData}: No additional row-level metadata is present.
\item \strong{colData}: Metadata for each cell, including six fields: sample, type, nFeature_RNA, ident, and others.
\item \strong{dim}: 1 feature and 106,297 cells.
\item \strong{colnames}: Cell identifiers for all 106,297 cells.
}
}
\usage{
data(sce_obj)
}
\value{
A \code{SingleCellExperiment} object containing single-cell RNA expression data.
}
\description{
Example \code{SingleCellExperiment} object containing gene expression data for 106,297 cells across two assays: counts and logcounts. The object includes metadata and assay data for RNA expression, which can be used directly in differential analysis functions like \code{sccomp_glm}.
}
\keyword{datasets}
