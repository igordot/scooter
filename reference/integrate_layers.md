# Integrate the layers of a Seurat object.

Runs the Seurat 5 `IntegrateLayers()` workflow. It adds what calling
`IntegrateLayers()` directly leaves for the caller to handle:

## Usage

``` r
integrate_layers(
  x,
  assay = "RNA",
  int_reduction = "cca",
  num_dim,
  k_anchor = 10,
  k_weight = 100,
  log_file = NULL,
  ...
)
```

## Arguments

- x:

  Seurat object with split layers and a "pca" reduction.

- assay:

  Assay to integrate.

- int_reduction:

  Integration method: "cca", "rpca", or "harmony".

- num_dim:

  Number of dimensions to use (5-50). For "cca"/"rpca", this value is
  passed as `dims = 1:num_dim` to
  `CCAIntegration()`/`RPCAIntegration()`, which both take a real `dims`
  argument. `HarmonyIntegration()` has no `dims` parameter. It runs on
  the full `"pca"` embedding as-is, through `Embeddings(orig)`,
  unsliced. So for `int_reduction = "harmony"`, `num_dim` only takes
  effect indirectly, through however many components `orig.reduction`
  ("pca") already has. The normal
  [`integrate_seurat_object()`](https://igordot.github.io/scooter/reference/integrate_seurat_object.md)
  pipeline is unaffected:
  [`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)
  already computes that "pca" reduction with exactly `num_dim`
  components. A direct call is different. If `x`'s "pca" reduction
  already has more components than `num_dim`, then `num_dim` has no
  effect on the harmony path.

- k_anchor:

  Neighbors to use when picking anchors. Ignored by "harmony". A value
  larger than Seurat's default of 5 strengthens the alignment.

- k_weight:

  Neighbors to use when weighting the corrections. Must be smaller than
  the number of cells in the smallest batch. Reduced to 25 when the
  smallest batch has fewer than 100 cells. Batches smaller than 25 cells
  are an error.

- log_file:

  Filename for the log file.

- ...:

  Extra arguments for whichever method `IntegrateLayers()` dispatches
  to - e.g. `theta`/`lambda`/`sigma`/`tau`/`nclust` for "harmony", or
  `k.filter`/`sample.tree`/ `normalization.method` for "cca"/"rpca". An
  unrecognized name errors (via
  [`rlang::check_dots_used()`](https://rlang.r-lib.org/reference/check_dots_used.html))
  rather than being silently dropped by the method function's own `...`.

## Value

A Seurat object with rejoined layers and an `int_reduction`-named
reduction (`"cca"`, `"rpca"`, or `"harmony"`).

## Details

- Clamps `k_weight` to fit the smallest batch. Seurat's own default of
  100 errors outright below that size. See `k_weight` below.

- Excludes `dims`/`k.anchor` for `"harmony"`. `HarmonyIntegration()`
  silently ignores both instead of erroring, so passing them would rely
  on that silent drop. See `num_dim`/`k_anchor` below.

- Forces `future::plan(future::sequential)`. Benchmarks show
  multisession gives no gain here, and is a documented liability at this
  scale elsewhere in this package.

- Rejoins the split layers with `JoinLayers()` afterward.
  `IntegrateLayers()` itself leaves this undone.

Unlike the older anchor workflow, this function creates no "integrated"
assay. Clustering and visualization use the corrected embedding
directly, instead of re-running PCA. This function expects
[`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)
to have run first. It is the one step of
[`integrate_seurat_object()`](https://igordot.github.io/scooter/reference/integrate_seurat_object.md)
that calls `IntegrateLayers()` directly.
