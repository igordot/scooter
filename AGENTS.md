# AGENTS.md

`scooter` is an R package (roxygen2 + testthat, not on CRAN). It loads,
preprocesses, and plots single-cell RNA-seq data on top of Seurat 5 and
6. `inst/extdata/scrna-scooter.R` is a docopt CLI that runs the package
as a pipeline.

Roxygen documents each function. Read the function or `man/<fn>.Rd`
first. This file covers only rules and traps that the code does not
show.

## Rules

- Do not make git commits or write to git history. Leave changes in the
  working tree.
- `NAMESPACE` and `man/` are generated. Run
  `Rscript -e 'roxygen2::roxygenise()'`. Do not edit them.
- Run `air format .` on the whole package after edits.
- Use `usethis` helpers when one fits (`use_package()`, `use_test()`).
- Loose `.png`, `.pdf`, and `.log` files in the repo root are generated
  output. Delete them.
- End each response with a list of the files that changed. If you
  installed the package, give the install path.

``` bash
Rscript -e 'testthat::test_local()'
Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-import.R")'
R CMD build . && R CMD check scooter_*.tar.gz
```

## Writing style

These rules apply to all prose: roxygen, `#` comments, this file, and
your responses.

- Write about 80% of the way toward ASD-STE100 Simplified Technical
  English.
- Use direct language and short sentences. Put one main idea in each
  sentence.
- Prefer familiar words when they preserve the meaning.
- Remove repetition and unnecessary introductory phrases.
- Preserve important technical details, assumptions, conditions,
  exceptions, and uncertainty.
- Do not remove information only to make the answer shorter.
- Prioritize clarity and accuracy over strict compliance with the
  standard.
- Use the domain words that the code already uses. Do not invent
  abstractions or labels for them.
- Function, package, and Seurat names stay as they are.
- Keep `#` comments to one line unless the code is not clear without
  more.
- A `@param` entry says what the parameter is, in one sentence. Do not
  add history or caveats that a caller does not need.

## Seurat 5 and 6

- Do not access assay slots. Use `GetAssayData()`/`SetAssayData()` with
  `layer=`.
- On split layers, `GetAssayData()` stops with an error. `LayerData()`
  only warns and returns the first layer, so it silently drops the other
  samples. Use `LayerData()` only after `JoinLayers()`.
- [`merge()`](https://rdrr.io/r/base/merge.html) leaves layers split
  (`counts.1`, `counts.2`, …). Call
  [`SeuratObject::JoinLayers()`](https://satijalab.github.io/seurat-object/reference/SplitLayers.html)
  after each merge.
- Make reductions with `CreateDimReducObject()`, not
  `new("DimReduc", ...)`.
- `VariableFeatures<-` with `NULL` does nothing on `Assay5`. To clear
  the variable features, use `DietSeurat()`.

## Package vs. CLI

- A package function wraps a Seurat function with the lab’s settings,
  plus the extra steps and plots of that step. Do not add a wrapper that
  only forwards arguments.
- The CLI only sequences package functions, parses arguments, and writes
  the logs. It is also a staging area for unfinished functions. These
  move into the package later.
- If a CLI function has the same name as a package export, move its
  extra behavior into the package function. Then delete the CLI copy. Do
  not rename the CLI copy.
- Functions still to move: `add_adt_assay_qc()` (the only ADT CLR
  normalization), `split_seurat_object()` (target `R/split.R`),
  `plot_adt_qc()`, `split_adt_hto_assay()`, `plot_hto_qc()`.

## Code conventions

- Verbs that touch expression data use S3: `.default` works on a plain
  matrix, and `.Seurat` extracts the data, calls `.default`, and puts
  the result back. Do not add this dispatch if the Seurat function
  already has it (check
  [`methods()`](https://rdrr.io/r/utils/methods.html)).
- In `dr.R`, the `.Seurat` methods call the Seurat functions directly.
  Always pass `assay` explicitly. Do not rely on `DefaultAssay()`.
- Each `.Seurat` method writes its plots to the working directory every
  time. There is no parameter to turn this off.
  `tests/testthat/teardown.R` deletes the files after the tests.
- Names: `read_*` reads from disk, `initialize_*` builds a raw object,
  `create_*` returns a usable object. Write `object` in function names
  and `obj` in variable names.
- Do not rename the file `seurat_obj.qs2`. Existing analysis dirs use
  it.
- Name a single-object first parameter `x`.
- Do not split one operation into two functions only because the input
  type differs.
- Add a parameter only if the function body uses it.
- Cell IDs are `sample_name:barcode`, with the Cell Ranger `-1` removed.
- Keep `log_file = NULL` and log with
  [`write_message()`](https://igordot.github.io/scooter/reference/write_message.md),
  one fact per call. Use
  [`message()`](https://rdrr.io/r/base/message.html) for separator
  banners.
- Get colors from
  [`get_color_scheme()`](https://igordot.github.io/scooter/reference/get_color_scheme.md).
- In new code, prefer intermediate variables to `|>` chains.
- Add a shim to `R/deprecated.R` only if the old name was in a commit
  (`git status` shows `M`, not `??`).

## Traps

- The CLI `create` step calls
  [`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)
  directly, not through
  [`do.call()`](https://rdrr.io/r/base/do.call.html). A top-level
  [`do.call()`](https://rdrr.io/r/base/do.call.html) breaks
  `LogSeuratCommand()`, and
  [`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md)
  then fails.
- The `Sys.sleep(1)` calls between plot writes keep the file times in
  order. Do not remove them.
- `mad_outlier_bounds()` uses all cells on purpose, to match
  `scater::isOutlier()`. Do not add a debris pre-filter.
- Per-sample splits use up to 30 dimensions. Merged and pooled analyses
  use up to 50. Keep both.
