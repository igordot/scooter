# Integrate a merged Seurat object across batches.

The integrate step of the pipeline in one call. It plays the same role
for a merged object that
[`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)
plays for a single sample.

## Usage

``` r
integrate_seurat_object(
  x,
  assay = "RNA",
  int_reduction = "cca",
  num_dim,
  batch_var = "orig.ident",
  num_variable_genes = 3000,
  num_neighbors = 30,
  k_anchor = 10,
  k_weight = 100,
  log_file = NULL,
  ...
)
```

## Arguments

- x:

  Merged Seurat object, such as the output of
  [`merge_seurat_objects()`](https://igordot.github.io/scooter/reference/merge_seurat_objects.md).

- assay:

  Assay to integrate. Passed straight through to
  [`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md).

- int_reduction:

  Integration method: "cca", "rpca", or "harmony".

- num_dim:

  Number of dimensions to use (5-50).

- batch_var:

  Metadata column identifying the batch.

- num_variable_genes:

  Number of variable features per batch.

- num_neighbors:

  Neighbors for UMAP.

- k_anchor:

  Neighbors to use when picking anchors. Passed straight through to
  [`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md) -
  see there for defaults/caveats (ignored by "harmony").

- k_weight:

  Neighbors to use when weighting the corrections. Passed straight
  through to
  [`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md) -
  see there for defaults/caveats (also ignored by "harmony").

- log_file:

  Filename for the log file.

- ...:

  Passed straight through to
  [`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md) -
  extra method-specific tuning arguments (harmony's
  `theta`/`lambda`/`sigma`/..., or cca/rpca's
  `k.filter`/`sample.tree`/...).

## Value

A Seurat object with rejoined layers, an `int_reduction`-named reduction
(`"cca"`, `"rpca"`, or `"harmony"`), and "tsne"/"umap" reductions
computed from it. Also writes the variable-gene overlap, QC
distribution, integrated-embedding, and tSNE/UMAP plots to the working
directory.

## Details

The function runs these steps:

1.  Split the assay layers by batch and prepare them
    ([`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)).

2.  Report how much the per-batch variable genes overlap
    ([`plot_var_genes_euler()`](https://igordot.github.io/scooter/reference/plot_var_genes_euler.md),
    [`plot_var_genes_upset()`](https://igordot.github.io/scooter/reference/plot_var_genes_upset.md)).

3.  Correct the batch effect
    ([`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md)).

4.  Compute a first-pass tSNE/UMAP on the corrected embedding.

The tSNE/UMAP from step 4 is a preview, capped at `num_dim`. It differs
from the real one
[`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md)
computes later over its own `num_dim`.
[`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md)
has no default `reduction`. Pass it this function's `int_reduction` name
explicitly. Otherwise it may use the uncorrected "pca" reduction that
[`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)
also leaves on the object.
