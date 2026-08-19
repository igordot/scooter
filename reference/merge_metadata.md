# Function to merge two metadata tables together.

`x` and `y` must not share column names other than `cell` -
[`stop()`](https://rdrr.io/r/base/stop.html)s if they do, rather than
letting `left_join()` silently rename the collision to
`<col>.x`/`<col>.y`.

## Usage

``` r
merge_metadata(x, y, log_file = NULL)
```

## Arguments

- x:

  A Seurat object or tibble with cell IDs (a `cell` column, or as
  rownames).

- y:

  A tibble with cell IDs (a `cell` column, or as rownames).

- log_file:

  A log filename.

## Value

A metadata file merged on cell identifiers.
