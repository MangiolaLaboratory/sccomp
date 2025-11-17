# seurat_obj

Example `Seurat` object containing gene expression data for 106,297
cells across a single assay. The object includes RNA counts and data
layers, but no variable features are defined. This dataset can be
directly used with functions like `sccomp_glm` for differential
abundance analysis.

## Usage

``` r
data(seurat_obj)
```

## Format

A `Seurat` object with the following structure:

- **assays**: Contains gene expression data. The active assay is `RNA`,
  with 1 feature and no variable features.

- **layers**: Two layers: counts and data, representing raw and
  processed RNA expression values, respectively.

- **samples**: 106,297 samples (cells) within the RNA assay.

## Value

A `Seurat` object containing single-cell RNA expression data.
