# Merge multiple Seurat objects into one.

Merges a list of Seurat objects, rejoins the split layers that Seurat 5
leaves behind, and removes genes that are detected in too few cells.
Features of the non-RNA assays (such as ADT and HTO) are always kept.

## Usage

``` r
merge_seurat_objects(seurat_object_list, min_cells = NULL, log_file = NULL)
```

## Arguments

- seurat_object_list:

  List of Seurat objects and/or analysis directory paths (mixing both is
  fine).

- min_cells:

  Keep genes detected in at least this many cells. If `NULL`, 10 is
  used, or 0.1% of the cells when there are more than 50,000 of them.

- log_file:

  Filename for the log file.

## Value

A merged Seurat object.

## Details

Elements of `seurat_object_list` may be Seurat objects, or paths to an
analysis directory containing a `seurat_obj.qs2` (as written by the
CLI). Which one it is does not matter to the caller: this function
resolves each element to an object first.

Once every input is an object, this function strips each one's own
reductions, scaled data, and clustering columns before merging. A
per-sample PCA, tSNE/UMAP, or clustering — written by
[`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)/[`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md)
— is meaningless once merged.
[`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md)/[`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md)
recompute all of it downstream anyway.

Writes `metrics-distribution.png`, a per-sample QC violin plot, to the
working directory.
