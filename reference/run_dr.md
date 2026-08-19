# Run dimensionality reduction: pca, tsne, or umap

Run dimensionality reduction: pca, tsne, or umap

## Usage

``` r
run_dr(x, dr_method, suffix = "", assay = "RNA", ...)
```

## Arguments

- x:

  A matrix or a Seurat object.

- dr_method:

  Dimensionality reduction method ("pca", "tsne", or "umap").

- suffix:

  Tag inserted into the reduction name and key.

- assay:

  Assay to use (Seurat methods). Always passed explicitly here rather
  than left `NULL`, which would otherwise fall back on whatever
  `DefaultAssay(x)` happens to be.

- ...:

  Passed to the chosen reduction function.

## Value

The reduced data (matrix input) or a Seurat object with the reduction
added. For a Seurat object, also writes the diagnostic plot(s) for the
chosen method to the working directory: PCA gets a sample-colored
scatter and an elbow plot (plus a heatmap past 15 components); tSNE/UMAP
get a sample-colored scatter.

## Details

`run_dr()` is a thin dispatcher: it routes to
[`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md),
[`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md),
or
[`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md)
based on `dr_method`, forwarding `suffix` and anything in `...`. It only
carries the parameters common to every method; method-specific ones
(`num_pcs`, `num_neighbors`, `reduction`, ...) are passed through `...`
to the chosen function. Each of those functions is itself generic: a
plain matrix takes the `.default` method (`irlba`/`Rtsne`/`uwot`), a
Seurat object takes the `.Seurat` method (native
`RunPCA()`/`RunTSNE()`/`RunUMAP()`, storing a `"<method><suffix>"`
reduction).

## See also

run_pca(), run_tsne(), run_umap()
