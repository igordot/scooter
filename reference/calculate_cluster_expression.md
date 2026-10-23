# Per-cluster average expression (non-log space) and summed raw counts.

Writes `expression-mean-<label>.csv` from
[`Seurat::AverageExpression()`](https://satijalab.org/seurat/reference/AverageExpression.html).
Writes `expression-raw-<label>.csv` from
[`Seurat::AggregateExpression()`](https://satijalab.org/seurat/reference/AggregateExpression.html).
When `x` has more than one sample, it also writes a `-per-sample`
version of each table.

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
