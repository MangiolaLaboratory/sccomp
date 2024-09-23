#' counts_obj
#'
#' @description
#' A tidy example dataset containing cell counts per cell group (cluster), sample, and phenotype for differential analysis. This dataset represents the counts of cells in various phenotypes and cell groups across multiple samples.
#'
#' @return A `tibble` representing cell counts per cluster, with columns for sample, type, phenotype, cell group, and counts.
#' 
#' @format A tidy data frame with the following columns:
#' \itemize{
#'   \item \strong{sample}: Factor, representing the sample identifier.
#'   \item \strong{type}: Factor, indicating the sample type (e.g., benign, cancerous).
#'   \item \strong{phenotype}: Factor, representing the cell phenotype (e.g., B_cell, HSC, etc.).
#'   \item \strong{count}: Integer, representing the number of cells for each cell group within each sample.
#'   \item \strong{cell_group}: Factor, representing the cell group (e.g., BM, B1, Dm, etc.).
#' }
#' @usage data(counts_obj)
#' 
#' @importFrom utils data
"counts_obj"


#' seurat_obj
#'
#' @description
#' Example `Seurat` object containing gene expression data for 106,297 cells across a single assay. The object includes RNA counts and data layers, but no variable features are defined. This dataset can be directly used with functions like `sccomp_glm` for differential abundance analysis.
#'
#' @return A `Seurat` object containing single-cell RNA expression data.
#'
#' @format A `Seurat` object with the following structure:
#' \itemize{
#'   \item \strong{assays}: Contains gene expression data. The active assay is `RNA`, with 1 feature and no variable features.
#'   \item \strong{layers}: Two layers: counts and data, representing raw and processed RNA expression values, respectively.
#'   \item \strong{samples}: 106,297 samples (cells) within the RNA assay.
#' }
#' @usage data(seurat_obj)
#' 
#' @importFrom utils data
"seurat_obj"


#' sce_obj
#'
#' @description
#' Example `SingleCellExperiment` object containing gene expression data for 106,297 cells across two assays: counts and logcounts. The object includes metadata and assay data for RNA expression, which can be used directly in differential analysis functions like `sccomp_glm`.
#'
#' @return A `SingleCellExperiment` object containing single-cell RNA expression data.
#'
#' @format A `SingleCellExperiment` object with the following structure:
#' \itemize{
#'   \item \strong{assays}: Two assays: counts (raw RNA counts) and logcounts (log-transformed counts).
#'   \item \strong{rowData}: No additional row-level metadata is present.
#'   \item \strong{colData}: Metadata for each cell, including six fields: sample, type, nFeature_RNA, ident, and others.
#'   \item \strong{dim}: 1 feature and 106,297 cells.
#'   \item \strong{colnames}: Cell identifiers for all 106,297 cells.
#' }
#' @usage data(sce_obj)
#' 
#' @importFrom utils data
"sce_obj"


#' multipanel_theme
#'
#' @description
#' A custom `ggplot2` theme used for creating publication-quality multi-panel plots. This theme modifies the appearance of plots by adjusting text sizes, spacing between panels, and axis formatting, ensuring better readability for complex figures.
#'
#' @return A `ggplot2` theme object.
#'
#' @format A `ggplot2` theme with the following adjustments:
#' \itemize{
#'   \item \strong{text}: Font size adjustments for plot titles, axis labels, and legend text.
#'   \item \strong{panel.spacing}: Adjusts the spacing between panels in multi-panel plots.
#'   \item \strong{axis.text}: Customises axis text appearance for better readability.
#' }
#' @usage data(multipanel_theme)
"multipanel_theme"