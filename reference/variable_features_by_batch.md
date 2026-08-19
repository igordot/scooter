# Variable features of each batch of a layer-split object.

After
[`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)
the per-batch variable features are stored as
`vf_<method>_<layer>_variable` columns in the assay meta data rather
than being reachable through `VariableFeatures()`, which returns the
combined set.

## Usage

``` r
variable_features_by_batch(x, assay = "RNA")
```

## Arguments

- x:

  Seurat object with split layers.

- assay:

  Assay to read the variable features from.

## Value

A named list of variable gene vectors, one element per batch.

## Details

[`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md)
now clears its own assay's `vf_*` columns before recomputing them. But
an object built before that guard existed can still carry a whole
history of them: each sample's own pre-merge run from
[`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md),
and the merge step's own post-merge preview run. `DietSeurat()` clears
layers, reductions, and graphs, but never assay meta data.
[`merge()`](https://rdrr.io/r/base/merge.html) disambiguates colliding
names with a numeric suffix, for example `vf_vst_counts.1_variable`.

Matching every `_variable$` column would pick up that whole history
along with the batches actually being asked about. This function instead
anchors on `Layers(x[[assay]], search = "counts")`, the object's current
split state. That keeps stale columns out, regardless of which
generation of the object this is.
