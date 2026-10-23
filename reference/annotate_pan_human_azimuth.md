# Annotate with Pan-human Azimuth

Annotates `query` with Pan-human Azimuth, a neural network trained on 23
human tissues. Install the AzimuthAPI package from GitHub
(`satijalab/AzimuthAPI`).

## Usage

``` r
annotate_pan_human_azimuth(query, query_assay = "RNA", verbose = TRUE)
```

## Arguments

- query:

  Query Seurat object.

- query_assay:

  Assay in `query` to annotate.

- verbose:

  Print progress messages.

## Value

`query` with these columns added to its metadata:
`pan_human_azimuth_label` (the most specific label),
`pan_human_azimuth_final_level_confidence` (its confidence),
`pan_human_azimuth_broad`, `pan_human_azimuth_medium`,
`pan_human_azimuth_fine`, `pan_human_azimuth_final_level_labels`,
`pan_human_azimuth_full_consistent_hierarchy` and
`pan_human_azimuth_full_hierarchical_labels`. The function writes the
columns of AzimuthAPI, under their own names and after a `cell` column,
to `annotation/annotation-pan-human-azimuth.csv.gz`. It also writes
`dr-umap-<column>.png` for each metadata column, except those with too
many groups to color.

## Details

Unlike
[`transfer_labels()`](https://igordot.github.io/scooter/reference/transfer_labels.md),
the function needs no reference object.

## References

Pan-human Azimuth website: <https://satijalab.org/pan_human_azimuth>

Sarkar S, Li Z, Molla G, Shenoy A, Zhang B, Collins D, et al.
Organism-scale annotation with Pan-human Azimuth. bioRxiv. 2026 Jul 21.
[doi:10.64898/2026.07.16.738997](https://doi.org/10.64898/2026.07.16.738997)
