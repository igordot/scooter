# Save a reduction plot, widened to fit its own legend.

Unlike
[`dr_plot_width()`](https://igordot.github.io/scooter/reference/dr_plot_width.md)
(a label-count guess made before the plot exists), this measures the
legend actually rendered on `plot` - so long label text widens the saved
image even when the label count is small. Requires an open graphics
device to get accurate text metrics, so it opens and closes a throwaway
one around the measurement.

## Usage

``` r
save_dr_plot(
  plot,
  file_prefix,
  file_format = c("png", "pdf"),
  width = NULL,
  height = 6,
  panel_width = height
)
```

## Arguments

- plot:

  A ggplot object, such as from
  [`plot_dr_group()`](https://igordot.github.io/scooter/reference/plot_dr_group.md).
  Forwarded straight to `ggsave(plot = )`.

- file_prefix:

  Path to save to without an extension.

- file_format:

  File formats to write.

- width:

  (optional) Saved width in inches. `NULL` measures `plot`'s legend and
  adds it to `panel_width`.

- height:

  Saved height in inches.

- panel_width:

  Width of the plot panel itself, excluding the legend, in inches. Only
  used when `width` is `NULL`.

## Value

`plot`, invisibly.
