# Violin plots of the per-cell QC metrics.

One violin panel per metric, grouped by sample. The metrics are metadata
columns, so `layer` is set explicitly: `VlnPlot()` otherwise searches
for the "data" layer, which does not exist before `NormalizeData()` has
been run.

## Usage

``` r
plot_metrics_distribution(
  x,
  metrics,
  group_by = "orig.ident",
  color_scheme = NULL,
  file = NULL,
  width = 12,
  height = 6
)
```

## Arguments

- x:

  Seurat object.

- metrics:

  Metadata columns to plot, one panel each. No default - the caller
  decides, since the right set depends on context (e.g. a single sample
  can afford `pct_ribo`/`pct_hb` alongside
  `detected_genes`/`total_counts`/`pct_mito`; a multi-sample plot
  usually can't). A `pct_*` metric that is under 1% in every cell is
  skipped (a flat near-zero violin says nothing) - its range is reported
  with a message instead.

- group_by:

  Metadata column to group the violins by.

- color_scheme:

  (optional) Named vector of colors. Defaults to the "samples" scheme,
  named by the sorted levels of `group_by`.

- file:

  Path to save the plot to. `NULL` returns it without writing. Parent
  directories are created as needed.

- width, height:

  Saved size in inches. Widen it when there are many samples.

## Value

A combined plot, one panel per plotted metric.
