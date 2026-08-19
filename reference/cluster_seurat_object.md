# Reduce dimensions, cluster, and plot the result.

The clustering step of the pipeline in one call: run tSNE and UMAP on an
existing "pca" reduction
([`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md),
[`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md)),
build the SNN graph and cluster over a range of resolutions
([`calculate_clusters()`](https://igordot.github.io/scooter/reference/calculate_clusters.md)),
and plot each reduction by cluster.

## Usage

``` r
cluster_seurat_object(
  x,
  assay = "RNA",
  reduction,
  num_dim,
  num_neighbors = 30,
  metadata_file = NULL,
  log_file = NULL
)
```

## Arguments

- x:

  Seurat object with `reduction` computed, such as the output of
  [`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)
  (which only ever produces "pca") or
  [`integrate_seurat_object()`](https://igordot.github.io/scooter/reference/integrate_seurat_object.md)
  (which also produces one named after its integration method: "cca",
  "rpca", or "harmony").

- assay:

  Assay to attribute the tSNE/UMAP reductions to.
  [`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)
  creates its own assay, so it needs no such parameter. This function
  creates none: it only calls
  [`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md)/[`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md),
  which read from an existing assay. So `assay` is explicit here, rather
  than read back off `DefaultAssay(x)`.

- reduction:

  Reduction to compute tSNE/UMAP from and to cluster on - "pca" after
  [`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md),
  or the integration method's own name ("cca"/"rpca"/"harmony") after
  [`integrate_seurat_object()`](https://igordot.github.io/scooter/reference/integrate_seurat_object.md)
  to cluster on the batch-corrected embedding instead of the
  pre-correction one. No default: which reduction is the intended one
  for clustering differs by which pipeline branch got you here, so every
  caller states it explicitly rather than risk silently picking the
  wrong one.

- num_dim:

  Principal components to use, bounded to 5-50.

- num_neighbors:

  Neighbors for UMAP and for the SNN graph.

- metadata_file:

  Path to write the cell metadata and embeddings to, via
  [`save_metadata()`](https://igordot.github.io/scooter/reference/save_metadata.md).
  `NULL` writes nothing.

- log_file:

  Filename for the log file.

## Value

A Seurat object with "tsne" and "umap" reductions and one `res.<x>`
column per retained clustering resolution.

## Details

This function computes both reductions before clustering, since the
cluster plots draw on them.
[`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md)
and
[`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md)
already write their own sample-colored scatter to the working directory,
so this function skips plotting "by sample" again.

The per-resolution cluster plots also write unconditionally, into a
`clusters-resolutions/` subdirectory. This matches every other verb's
standard-output convention.
