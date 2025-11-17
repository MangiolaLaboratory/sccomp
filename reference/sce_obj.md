# sce_obj

Example `SingleCellExperiment` object containing gene expression data
for 106,297 cells across two assays: counts and logcounts. The object
includes metadata and assay data for RNA expression, which can be used
directly in differential analysis functions like `sccomp_glm`.

## Usage

``` r
data(sce_obj)
```

## Format

A `SingleCellExperiment` object with the following structure:

- **assays**: Two assays: counts (raw RNA counts) and logcounts
  (log-transformed counts).

- **rowData**: No additional row-level metadata is present.

- **colData**: Metadata for each cell, including six fields: sample,
  type, nFeature_RNA, ident, and others.

- **dim**: 1 feature and 106,297 cells.

- **colnames**: Cell identifiers for all 106,297 cells.

## Value

A `SingleCellExperiment` object containing single-cell RNA expression
data.
