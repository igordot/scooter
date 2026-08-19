# Report the size of an object's slots, largest first.

Works on any S4 object or list, not just Seurat objects. It descends
into S4 slots and list elements, bounded by `max_depth`. So a Seurat
object, an assay pulled from one (`x[["RNA"]]`), a `DimReduc`, or a
plain nested list all get the same recursive breakdown.

## Usage

``` r
profile_object_size(x, max_depth = 3)
```

## Arguments

- x:

  Any object - a Seurat object, a piece pulled from one (an assay via
  `x[["RNA"]]`, a reduction via `x[["pca"]]`, ...), a plain list, or
  anything else
  [`object.size()`](https://rdrr.io/r/utils/object.size.html) accepts.

- max_depth:

  How many levels to descend before reporting a subtree's total size
  instead of continuing to break it down.

## Value

A tibble with one row per slot/element, largest first: `path` (its
position, `" > "`-separated) and `size` (human-readable, e.g. "1.2 Mb").

## Details

An `Assay5`'s `layers` slot is itself a list. So
`layers > counts`/`layers > data`/`layers > scale.data` are sized
separately, not lumped into one `layers` total.

A matrix, data frame, or atomic vector is neither S4 nor a list. Passing
one of these just reports its own size, the same as calling
[`object.size()`](https://rdrr.io/r/utils/object.size.html) directly.
There is nothing to descend into.
