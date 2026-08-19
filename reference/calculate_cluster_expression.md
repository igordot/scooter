# Per-cluster average expression (non-log space).

Writes `expression-mean-<label>.csv`, and additionally
`expression-mean-<label>-per-sample.csv` when `x` has more than one
sample.

## Usage

``` r
calculate_cluster_expression(x, label)
```

## Arguments

- x:

  Seurat object with an identity set, such as via
  [`set_identity()`](https://igordot.github.io/scooter/reference/set_identity.md).

- label:

  Suffix for the output filenames.

## Value

`x`, invisibly. Called for its file output.
