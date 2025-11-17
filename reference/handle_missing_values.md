# Handle missing values in a data frame

Handle missing values in a data frame

## Usage

``` r
handle_missing_values(data)
```

## Arguments

- data:

  A data frame to process

## Value

A data frame with missing values handled:

- For categorical variables (factor/character): NA values are converted
  to a factor level "NA"

- For numeric variables: NA values are replaced with the mean of the
  column
