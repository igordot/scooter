# Transfer labels from a reference Seurat object

Reference-based annotation wrapping Seurat's reference-mapping workflow
([`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html) +
[`Seurat::TransferData()`](https://satijalab.org/seurat/reference/TransferData.html))
to annotate `query` with labels (cell type, or any other categorical
metadata) from an already integrated and annotated reference object.
Optionally also projects `query` onto the reference's UMAP via
[`Seurat::MapQuery()`](https://satijalab.org/seurat/reference/MapQuery.html).

## Usage

``` r
transfer_labels(
  query,
  ref,
  query_label_col,
  ref_label_col,
  query_assay,
  ref_assay,
  normalization,
  ref_reduction = "pca",
  num_dim,
  anchor_reduction = "pcaproject",
  ref_umap = FALSE,
  verbose = TRUE
)
```

## Arguments

- query:

  Query Seurat object.

- ref:

  Reference Seurat object - already integrated and annotated.

- query_label_col:

  Metadata column name added to `query` for the transferred labels.

- ref_label_col:

  Metadata column in `ref` holding the labels to transfer (cell type, or
  any other category).

- query_assay:

  Assay in `query` to use for anchor finding. `"RNA"` (log-normalized)
  or `"SCT"`, depending on how `query` was built. No default, like the
  params below: this and `ref_assay`/`normalization`/ `num_dim` all
  describe how a specific object was built.

- ref_assay:

  Assay in `ref` to use for anchor finding.

- normalization:

  Normalization method shared by `ref` and `query`
  (`"log"`/`"LogNormalize"` or `"sct"`/`"SCT"`).
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html)
  takes one value for both objects, so `ref` and `query` must actually
  share a method.

- ref_reduction:

  Reduction in `ref` to project `query` onto.

- num_dim:

  Number of dimensions to use for anchor finding and label transfer.

- anchor_reduction:

  Dimensional reduction workflow
  [`Seurat::FindTransferAnchors()`](https://satijalab.org/seurat/reference/FindTransferAnchors.html)
  uses to find anchors.

- ref_umap:

  If `TRUE`, also runs
  [`Seurat::MapQuery()`](https://satijalab.org/seurat/reference/MapQuery.html)
  and adds a `"ref.umap"` reduction to `query`, projecting it onto
  `ref`'s UMAP. `ref` needs a UMAP built with `return.model = TRUE`,
  which this function does not check for - `MapQuery()` errors on its
  own if it is missing.

- verbose:

  Print progress messages.

## Value

`query` with `query_label_col` (a factor) added to its metadata, and a
`"ref.umap"` reduction added if `ref_umap = TRUE`. Also writes the
predictions table (`cell`/`predicted.id`/`prediction.score.max`) to
`annotation/annotation-transfer-label-<query_label_col>.csv.gz`, and a
`dr-<reduction>-<query_label_col>.png` UMAP colored by the transferred
labels for each of `"umap"`/`"ref.umap"` already present on `query`.
