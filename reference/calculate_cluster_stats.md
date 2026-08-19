# Per-cluster cell counts and a joined metadata+embeddings table.

Writes `metadata-<label>.csv.gz` (cell metadata joined with the
tSNE/UMAP embeddings, skipped if neither reduction is present) and
`summary-<label>.csv` (cell counts and percentages per cluster, split
further by sample into `summary-<label>-per-sample.csv` when that
differs).

## Usage

``` r
calculate_cluster_stats(x, label)
```

## Arguments

- x:

  Seurat object with an identity set, such as via
  [`set_identity()`](https://igordot.github.io/scooter/reference/set_identity.md).

- label:

  Suffix for the output filenames.

## Value

`x`, invisibly. Called for its file output.
