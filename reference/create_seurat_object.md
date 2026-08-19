# Create a usable Seurat object from a sample.

Everything needed to take a sample from a counts file to an object ready
for clustering, in one call.

## Usage

``` r
create_seurat_object(
  path,
  sample_name,
  num_mads = 3,
  min_genes = 500,
  max_genes = NULL,
  min_counts = NULL,
  max_counts = NULL,
  max_mt = 10,
  normalization_method = "log",
  num_variable_genes = 3000,
  num_pcs = 50,
  num_dim = 30,
  log_file = NULL
)
```

## Arguments

- path:

  Path to a 10x directory or a flat counts file.

- sample_name:

  Sample/library name, used as the cell ID prefix.

- num_mads, min_genes, max_genes, min_counts, max_counts, max_mt:

  Quality cutoffs, passed to
  [`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md)
  — see there for how the outlier stage and the fixed cutoffs combine.

- normalization_method:

  Normalization method, passed to
  [`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md)
  (as its own `method` argument): "log" or "sct".

- num_variable_genes:

  Number of variable features.

- num_pcs:

  Principal components to compute. Reduced automatically when the object
  is too small to support that many.

- num_dim:

  Number of dimensions for the first-pass tSNE/UMAP, capped at
  `num_pcs`.

- log_file:

  Filename for the log file.

## Value

A Seurat object with normalized data, variable features, scaled data,
and "pca", "tsne", and "umap" reductions — ready for
[`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md).
Also writes the variable-feature, PCA, tSNE, and UMAP diagnostic plots
(via
[`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md),
[`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md),
[`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md),
and
[`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md))
and the metadata+embeddings table (via
[`save_metadata()`](https://igordot.github.io/scooter/reference/save_metadata.md),
as `metadata.csv.gz`) to the working directory.

## Details

The function runs these steps:

1.  Read the counts
    ([`read_counts_file()`](https://igordot.github.io/scooter/reference/read_counts_file.md)).

2.  Build the object
    ([`initialize_seurat_object()`](https://igordot.github.io/scooter/reference/initialize_seurat_object.md)).

3.  Attach any antibody capture data as an "ADT" assay
    ([`add_seurat_assay()`](https://igordot.github.io/scooter/reference/add_seurat_assay.md)).

4.  Apply the quality cutoffs
    ([`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md)).

5.  Normalize, select variable features, and scale
    ([`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md),
    which does all three for either method).

6.  Run PCA
    ([`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md)).

7.  Compute a first-pass tSNE/UMAP
    ([`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md),
    [`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md)).

This function attaches the ADT assay *before* filtering, on purpose.
[`add_seurat_assay()`](https://igordot.github.io/scooter/reference/add_seurat_assay.md)
restricts the object to the barcodes both matrices share. Attaching ADT
after filtering would apply the cutoffs to a different set of cells.

PCA, and the tSNE/UMAP computed from it, runs on whichever assay
normalization left active. `SCTransform()` creates and activates "SCT".
The log path leaves "RNA" active. So `dr_assay` is decided directly from
`normalization_method`, not read back off the object. Reading it back
off `DefaultAssay()` would reintroduce the implicit dependency an
explicit `assay` argument exists to avoid.

The tSNE/UMAP computed here are a first-pass look, capped at `num_dim`
dimensions.
[`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md)
recomputes both once a clustering resolution's own `num_dim` is chosen.
So this pair is not the final one used for the actual clustering plots.
