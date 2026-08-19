# Scatter plot of a reduction, colored by a grouping variable.

The categorical counterpart of
[`plot_dr_feature()`](https://igordot.github.io/scooter/reference/plot_dr_feature.md):
cells are colored by a metadata column (or by the current identities)
rather than by a continuous feature. Cells are shuffled so no one group
is drawn on top of the others.

## Usage

``` r
plot_dr_group(
  x,
  reduction = "umap",
  group_by = NULL,
  color_scheme = NULL,
  pt_size = NULL,
  na_value = "grey50",
  file_prefix = NULL,
  file_format = c("png", "pdf"),
  width = NULL,
  height = 6
)
```

## Arguments

- x:

  Seurat object.

- reduction:

  Dimensionality reduction to plot on (e.g. "umap", "tsne", "pca").

- group_by:

  (optional) Metadata column to color by. `NULL` uses the current
  identities.

- color_scheme:

  (optional) Vector of colors. Defaults to the "clusters" scheme.

- pt_size:

  (optional) Point size. Defaults to `get_dr_point_size(x)`.

- na_value:

  Color for cells with no value in `group_by`.

- file_prefix:

  (optional) Path to save to.

- file_format:

  File formats to write when `file_prefix` is given.

- width:

  (optional) Saved width in inches. `NULL` measures the plot's own
  legend and widens to fit it - see
  [`save_dr_plot()`](https://igordot.github.io/scooter/reference/save_dr_plot.md).

- height:

  Saved height in inches.

## Value

A ggplot object.
