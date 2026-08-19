# Plot a UMAP colored by one or more cluster/grouping variables.

Plot a UMAP colored by one or more cluster/grouping variables.

## Usage

``` r
plot_dr_umap_clusters(
  x,
  cluster_var,
  color_scheme = get_color_scheme("clusters")
)
```

## Arguments

- x:

  Seurat object with a "umap" reduction.

- cluster_var:

  Metadata column(s) to color by, as a single string (comma-separated
  for more than one).

- color_scheme:

  (optional) Named vector of colors. Defaults to
  `get_color_scheme("clusters")`.

## Value

Nothing; called for its file output, one `dr-umap-<cluster_var>.png` per
variable.
