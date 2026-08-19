# Add mitochondrial, ribosomal, and hemoglobin percentages to a Seurat object.

Detected by gene symbol prefix, case-insensitively so the same patterns
match both human (`MT-`, `RPS`/`RPL`, `HBA`/`HBB`) and mouse (`mt-`,
`Rps`/`Rpl`, `Hba-a1`/`Hbb-bs`) naming. `RPS6K*` (ribosomal protein S6
kinase) is excluded from `pct_ribo` - it matches the `RPS` prefix but is
a signaling kinase, not a ribosomal protein.

## Usage

``` r
add_gene_class_percent(x)
```

## Arguments

- x:

  A Seurat object.

## Value

Seurat object with `pct_mito`, `pct_ribo`, and `pct_hb` metadata columns
added.

## Details

High expression levels of mitochondrial genes could be an indicator of
poor sample quality, leading to a high fraction of apoptotic or lysing
cells (see reference below).

`pct_ribo` is normally 15-45%; the exact fraction varies by cell type
and overall cell health. "Although both ribosomal proteins and ribosomal
RNA (rRNA) make up the ribosome complex, ribosomal protein transcripts
are not equivalent to ribosomal RNA (rRNA). Ribosomal protein transcript
detection will not necessarily correlate with either rRNA or
mitochondrial transcripts." (see reference below)

High `pct_hb` is a contamination signal in PBMC and most solid tissue:
per-sample `pct_hb` flags incomplete red blood cell lysis in PBMC preps,
or inadequate perfusion in solid tissue.

## References

`pct_mito`:
<https://kb.10xgenomics.com/s/article/360001086611-Why-do-I-see-a-high-level-of-mitochondrial-gene-expression>

`pct_ribo`:
<https://kb.10xgenomics.com/s/article/218169723-What-fraction-of-reads-map-to-ribosomal-proteins>
