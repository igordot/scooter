# Plot a gene list in multiple ways: UMAP, dot plot, violin, and per-cluster bar plot.

Plot a gene list in multiple ways: UMAP, dot plot, violin, and
per-cluster bar plot.

## Usage

``` r
plot_cluster_markers_top(
  x,
  genes,
  filename_base,
  color_scheme = get_color_scheme("clusters")
)
```

## Arguments

- x:

  Seurat object with an identity set.

- genes:

  Genes to plot.

- filename_base:

  Path prefix for the output files.

- color_scheme:

  (optional) Named vector of cluster colors, used by the violin and bar
  plots. Defaults to `get_color_scheme("clusters")`.

## Value

Nothing; called for its file output.
