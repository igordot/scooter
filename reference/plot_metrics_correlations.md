# Scatter plots of the per-cell QC metrics against each other.

One panel per pair in `pairs`, colored by sample. No legend, since this
only ever runs on a single sample - see
[`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md).
Pairs are a curated list rather than every combination of a metric
vector, since that would make the panel count quadratic in the number of
metrics (5 metrics -\> 10 panels via
[`combn()`](https://rdrr.io/r/utils/combn.html)) and most pairs are not
biologically meaningful together.

## Usage

``` r
plot_metrics_correlations(
  x,
  pairs = list(c("total_counts", "detected_genes"), c("pct_mito", "detected_genes"),
    c("pct_mito", "pct_ribo")),
  group_by = "orig.ident",
  color_scheme = NULL,
  file = NULL,
  width = 16,
  height = 5
)
```

## Arguments

- x:

  Seurat object.

- pairs:

  A list of length-2 character vectors, `c(feature1, feature2)`, one per
  panel.

- group_by:

  Metadata column to color the points by.

- color_scheme:

  (optional) Named vector of colors. Defaults to the "samples" scheme,
  named by the sorted levels of `group_by`.

- file:

  Path to save the plot to. `NULL` returns it without writing. Parent
  directories are created as needed.

- width, height:

  Saved size in inches. The default is wide and short: one square panel
  per pair of metrics, side by side.

## Value

A combined plot, one panel per pair.
