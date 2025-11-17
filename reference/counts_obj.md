# counts_obj

A tidy example dataset containing cell counts per cell group (cluster),
sample, and phenotype for differential analysis. This dataset represents
the counts of cells in various phenotypes and cell groups across
multiple samples.

## Usage

``` r
data(counts_obj)
```

## Format

A tidy data frame with the following columns:

- **sample**: Factor, representing the sample identifier.

- **type**: Factor, indicating the sample type (e.g., benign,
  cancerous).

- **phenotype**: Factor, representing the cell phenotype (e.g., B_cell,
  HSC, etc.).

- **count**: Integer, representing the number of cells for each cell
  group within each sample.

- **cell_group**: Factor, representing the cell group (e.g., BM, B1, Dm,
  etc.).

## Value

A `tibble` representing cell counts per cluster, with columns for
sample, type, phenotype, cell group, and counts.
