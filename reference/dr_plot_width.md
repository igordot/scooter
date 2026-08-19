# Width for a reduction plot, widened so many legend labels stay legible.

A label-count estimate made without building the plot.
[`plot_dr_group()`](https://igordot.github.io/scooter/reference/plot_dr_group.md)
now sizes itself from the actual rendered legend instead (see
[`save_dr_plot()`](https://igordot.github.io/scooter/reference/save_dr_plot.md)),
so this has no callers left in the package - kept as a cheap standalone
estimate for callers that don't have a plot object yet.

## Usage

``` r
dr_plot_width(x, group_by = "orig.ident")
```

## Arguments

- x:

  Seurat object, or the number of legend labels.

- group_by:

  Metadata column whose distinct values become the plot's legend.
  Ignored if `x` is already a number.

## Value

Numeric width in inches.
