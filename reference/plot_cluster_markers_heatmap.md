# Heatmap of the top cluster markers, at several top-N cutoffs.

Heatmap of the top cluster markers, at several top-N cutoffs.

## Usage

``` r
plot_cluster_markers_heatmap(x, markers_tbl, num_genes, filename_base)
```

## Arguments

- x:

  Seurat object with an identity set.

- markers_tbl:

  Marker table, such as from
  [`calculate_cluster_markers()`](https://igordot.github.io/scooter/reference/calculate_cluster_markers.md) -
  either the standard `log2FC`/`cluster` columns, or the pairwise
  `log2FC_min`/`cluster` columns.

- num_genes:

  Top-N genes per cluster to include, one heatmap per value.

- filename_base:

  Path prefix for the output files.

## Value

Nothing; called for its file output.
