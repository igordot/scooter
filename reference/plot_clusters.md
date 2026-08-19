# Plot a cluster resolution on tSNE and UMAP.

A single cluster says nothing, and past the color scheme's length the
colors start repeating, so nothing is plotted outside
`2:length(color_scheme)` clusters.

## Usage

``` r
plot_clusters(
  x,
  resolution,
  filename_base,
  color_scheme = get_color_scheme("clusters")
)
```

## Arguments

- x:

  Seurat object.

- resolution:

  Metadata column (or bare resolution value, resolved via
  [`check_identity_column()`](https://igordot.github.io/scooter/reference/check_identity_column.md))
  to set as the identity before plotting.

- filename_base:

  Path prefix for the two plots (`-tsne`/`-umap` are appended by
  [`plot_dr_group()`](https://igordot.github.io/scooter/reference/plot_dr_group.md)).

- color_scheme:

  (optional) Named vector of cluster colors. Defaults to
  `get_color_scheme("clusters")`.

## Value

`x` with `resolution` set as its identity.
