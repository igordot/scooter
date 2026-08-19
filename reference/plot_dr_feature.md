# Plot a single feature overlaid on a dimensionality reduction.

One point per cell, colored by expression/abundance of `feature` (a
gene, antibody, or any other value `FeaturePlot()` can fetch). Unlike
passing a multi-color ramp straight to `FeaturePlot(cols = )`, the ramp
is applied afterward with `scale_color_gradientn()` so the color scale
stays continuous - see the "cols" note below.

## Usage

``` r
plot_dr_feature(
  x,
  feature,
  assay = NULL,
  reduction = "umap",
  cells = NULL,
  color_scheme = NULL,
  pt_size = NULL
)
```

## Arguments

- x:

  Seurat object.

- feature:

  Feature name (gene, antibody, or metadata column) to plot.

- assay:

  (optional) Assay to fetch `feature` from. `NULL` uses
  `DefaultAssay(x)`.

- reduction:

  Dimensionality reduction to plot on (e.g. "umap", "tsne", "pca").

- cells:

  (optional) Cells to plot, in the given order (controls z-order of
  overlapping points). Defaults to a random shuffle of all cells so no
  single cell/group dominates due to plot order.

- color_scheme:

  (optional) Vector of colors for the low-to-high expression gradient.

- pt_size:

  (optional) Point size. Defaults to `get_dr_point_size(x)`.

## Value

A ggplot object.

## Seurat's `FeaturePlot(cols = )` quirk

As of Seurat 5, passing a `cols` vector of anything other than exactly 2
colors makes `FeaturePlot()` bin the continuous values into 2 groups
before coloring (`cut(x, breaks = 2)`), which turns a smooth expression
gradient into a flat low/high split. Seurat 4 only did this when `cols`
had exactly 2 colors, so a long custom ramp used to pass through
untouched. This function sidesteps the quirk (and the version
difference) entirely by never passing a multi-color `cols` to
`FeaturePlot()` - the ramp is layered on afterward instead.
