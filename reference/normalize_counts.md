# Normalize the counts, select variable features, and scale.

A thin wrapper that picks the right Seurat function for `method` and
fixes the parameters this pipeline uses, so either path hands back a
fully processed object — no separate variance step needed.
`method = "sct"` already worked this way natively: `SCTransform()`
produces the counts, variable features and scaled data in one pass.
`method = "log"` now mirrors that explicitly: `NormalizeData()` followed
by `FindVariableFeatures()` and `ScaleData()`.

## Usage

``` r
normalize_counts(
  x,
  method = "log",
  num_variable_genes = 3000,
  assay = "RNA",
  log_file = NULL
)
```

## Arguments

- x:

  A Seurat object.

- method:

  "log" for `NormalizeData()` + `FindVariableFeatures()` + `ScaleData()`
  (log-normalized to a scale factor of 10,000), or "sct" for
  `SCTransform()`.

- num_variable_genes:

  Number of variable features.

- assay:

  Assay to normalize.

- log_file:

  Log file.

## Value

The processed Seurat object. Also writes the variable-feature plot
(`variance-features.png`) to the working directory.

## Details

No covariates are regressed out (`vars.to.regress = NULL`) for either
path. Germain et al. (2020) found that the "common practice of
regressing out cell covariates, such as the detection rate or proportion
of mitochondrial reads\[,\] nearly always had a negative impact."

On the log path, this function drops any `vf_*` columns already on
`assay`'s feature meta data before recomputing them. An earlier call can
leave these columns behind: each sample's own pre-merge run, or a merge
step's preview run. Left alone, they would just accumulate call after
call. See
[`variable_features_by_batch()`](https://igordot.github.io/scooter/reference/variable_features_by_batch.md)
for what reading an accumulated column list, instead of a current one,
leads to.

`SCTransform()` never writes this kind of per-gene bookkeeping column.
So the sct path has nothing to clean up. The check runs the same way for
both paths; it is a no-op for sct, not something conditioned on
`method`.
