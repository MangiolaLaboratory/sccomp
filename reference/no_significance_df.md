# no_significance_df

A small example dataset containing cell counts across samples,
conditions, and cell groups. This dataset is used to demonstrate the use
of `sccomp` functions in scenarios where there is no significant
difference in cell composition between conditions.

## Usage

``` r
data(no_significance_df)
```

## Format

A tibble with the following columns:

- **sample**: Character. Identifier for each sample.

- **condition**: Character. Experimental condition or group (e.g., "X"
  or "Y").

- **cell_group**: Character. Cell group or cell type (e.g., "A", "B").

- **count**: Numeric. Count of cells in the given sample, condition, and
  cell group.

## Value

A tibble with 34 rows and 4 columns: `sample`, `condition`,
`cell_group`, and `count`.

## Examples

``` r
data(no_significance_df)
head(no_significance_df)
#> # A tibble: 6 × 4
#>   sample condition cell_group count
#>   <chr>  <chr>     <chr>      <dbl>
#> 1 1      X         A            713
#> 2 12     X         A            556
#> 3 12     X         B             43
#> 4 16     X         A            164
#> 5 17     X         A            579
#> 6 17     X         B             21
```
