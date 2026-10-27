# Handle missing values in a data frame

An NA means the covariate is unknown for that row, so it resolves to the
average of the fitted covariate: evenly over the levels of a factor, and
to the mean of a continuous one.

## Usage

``` r
handle_missing_values(data)
```

## Arguments

- data:

  A data frame holding the formula covariates, with categorical columns
  already declared as factors over the fitted levels and continuous
  columns already z-scored against the fitted rows

## Value

A data frame with missing values handled:

- categorical: NA becomes a real "NA" level, last in the ordering, which
  the caller spreads evenly over the declared levels

- continuous: NA becomes 0
