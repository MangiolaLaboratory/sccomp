# Plot Scatterplot of Cell-group Proportion

This function creates a scatterplot of cell-group proportions,
optionally highlighting significant differences based on a given
significance threshold.

## Usage

``` r
plot_scatterplot(
  .data,
  data_proportion,
  factor_of_interest,
  .cell_group,
  .sample,
  significance_threshold = 0.05,
  my_theme
)
```

## Arguments

- .data:

  Data frame containing the main data.

- data_proportion:

  Data frame containing proportions of cell groups.

- factor_of_interest:

  A factor indicating the biological condition of interest.

- .cell_group:

  The cell group to be analysed.

- .sample:

  The sample identifier.

- significance_threshold:

  Numeric value specifying the significance threshold for highlighting
  differences. Default is 0.025.

- my_theme:

  A ggplot2 theme object to be applied to the plot.

## Value

A ggplot object representing the scatterplot.

## Examples

``` r
# Example usage:
# plot_scatterplot(.data, data_proportion, "condition", "cell_group", "sample", 0.025, theme_minimal())
```
