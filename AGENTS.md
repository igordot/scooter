# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`scooter` is an R package (roxygen2 + testthat, not on CRAN) providing a data-agnostic
framework for loading, preprocessing, and plotting single-cell RNA-seq datasets on top of Seurat.

**Function-level behavior lives in roxygen, not here.** Every exported function's own `@details`/
`@param`/`@return` already explains its rationale, defaults, and what it writes — read the function
(or `man/<fn>.Rd`) first. This file covers what roxygen can't: cross-cutting conventions, the CLI
script (`inst/extdata/scrna-scooter.R`, a docopt script with no roxygen of its own), and
historical/deferred notes that don't belong in documentation of current behavior.

## Writing style

Write all prose for this project in Simplified Technical English (STE, per ASD-STE100). This
applies to CLAUDE.md itself, roxygen comments (`@param`/`@details`/`@return`), inline `#` code
comments, and chat responses while working in this repo.

- Write one short sentence per instruction or fact. Keep sentences under 20 words.
- Use active voice. Write "the function reads the file," not "the file is read by the function."
- Use simple verb tenses. Prefer simple present. Avoid compound or complex tenses.
- Use one term for one thing. Do not switch between synonyms for the same concept.
- Avoid long noun strings. Write "the pattern to access the assay data layer," not "the assay
  data layer access pattern."
- Prefer plain, common words over technical synonyms when a plain word means the same thing.
- Technical names are exempt from the vocabulary rule above. Function names, package names, and
  R/Seurat/genomics terms (`GetAssayData()`, layer, assay, `DimReduc`, ...) stay as they are — the
  sentence-structure rules above still apply around them.
- Inline `#` code comments stay to one line. Use more than one line only when the code's behavior
  is not obvious without it. This rule covers inline comments only, not roxygen (`#'`) or file
  headers.
- When you edit this file, skip rename history. Skip facts that the referenced file shows at once.

## Environment (read this first)

**The package targets Seurat v5.** The code uses the `layer=` argument of
`GetAssayData()`/`SetAssayData()`, which SeuratObject 4 does not have. The `Seurat >= 5.3.1` floor
comes from `leiden_method = "igraph"` in `calculate_clusters()`; everything else only needs 5.0.0.

Whatever `R`/`Rscript` is first on `PATH` needs to satisfy the Seurat/SeuratObject versions in
`DESCRIPTION` and have the dev toolchain installed (`devtools`, `roxygen2`, `testthat`, `pkgload`).
Check with `Rscript -e 'packageVersion("Seurat")'` if unsure. Then everything below runs with a
plain `Rscript` call; no env activation step is needed.

`scooter-package.R`'s `@import dplyr`/`@import Seurat` tags mean those packages' functions resolve
everywhere without per-function `@importFrom` tags — don't add them elsewhere. `DESCRIPTION` sets
`Roxygen: list(markdown = TRUE)` — without it, `[fn()]` cross-references render as dead literal
bracket text, and a bare `%`/`{`/`}` in doc prose corrupts Rd parsing. Keep markdown mode on.

### Seurat 5 conventions used here

- **Never touch assay slots directly.** `Assay5` has no `@counts` / `@data` / `@scale.data`; it holds
  a named *layer* list. Read with `GetAssayData(obj, assay=, layer=)`, write with
  `SetAssayData(obj, assay=, layer=, new.data=)`. Both dispatch on `Assay` and `Assay5`.
- **The two accessors disagree on split layers.** `GetAssayData()` is formally superseded by
  `LayerData()`, but that matters less than the disagreement. `GetAssayData()` **aborts**
  (*"doesn't work for multiple layers in v5 assay"*); `LayerData()` only **warns** and returns the
  *first* matching layer, silently dropping every other sample. `save_counts()` is the one function
  using `LayerData()`, and only because it calls `JoinLayers()` first. Without that guard,
  `GetAssayData()`'s hard error is the safer default.
- **`merge()` leaves layers split** (`counts.1`, `counts.2`, …). Any merge must be followed by
  `SeuratObject::JoinLayers()`, which `merge_seurat_objects()` already does.
- **`DimReduc` keys are validated** against `^[a-zA-Z][a-zA-Z0-9]*_$`. Build reductions with
  `CreateDimReducObject()` (which repairs the key) rather than `new("DimReduc", ...)`.
- **Clearing variable features is not possible via `VariableFeatures<-`.** `Assay5` stores them as
  `vf_*` columns in the assay `meta.data`, and assigning `character()`/`NULL` is a silent no-op. Use
  `DietSeurat(obj, layers = c("counts", "data"), dimreducs = NULL, graphs = NULL)` instead — variable
  features are then left to be overwritten by `normalize_counts()`.

## Commands

```bash
# load for interactive work
Rscript -e 'pkgload::load_all(".")'

# full test suite
Rscript -e 'testthat::test_local()'

# a single test file
Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-import.R")'

# regenerate NAMESPACE + man/ after changing roxygen comments
Rscript -e 'roxygen2::roxygenise()'

# format (config in air.toml) — run on the whole package, not just the files you touched; it's fast
air format .

# full check
R CMD build . && R CMD check scooter_*.tar.gz
```

`NAMESPACE` and everything in `man/` are generated — never edit them by hand.

Use `usethis` helpers for package scaffolding when one fits. For example, use `use_package()` to add a
dependency and `use_test()` to add a test file.

**Never make git commits.** Leave all changes in the working tree for review; do not `git commit`,
`git rm`, or otherwise write to git history. Running the tests, `R CMD build`/`check`, and
`R CMD INSTALL` (into whichever R library the environment's `R`/`Rscript` resolves to) is fine and
does not need to be asked about.

**Always end a response with a summary of what changed**: every file edited, created, or deleted,
listed by path. If the package was installed, say so and give the install path. If nothing was
edited, say that too.
Report outcomes (files, install path, test counts). Do not list routine steps such as `air format .`.

Loose `.png`, `.pdf`, and `.log` files in the repo root are generated output. Delete them without
asking. Other file types and `tests/testthat/` files do not follow this rule.

## Architecture

### What the package is for vs. what the CLI is for

Read this before deciding where a function belongs — it settles most "should this move?" questions.

- **The package enhances Seurat.** Each verb wraps a Seurat function with the parameters this lab
  has settled on, plus the extra steps and the plots that belong to that step. The added value *is*
  the extra steps and plots — a thin binding that only forwards arguments is not worth having. So
  plotting and writing a step's outputs are legitimately package work, not something to push out to
  the caller.
- **The CLI groups related package functions into a pipeline.** `inst/extdata/scrna-scooter.R`
  exists to sequence package verbs into `create` / `cluster` / `integrate` / … and to handle
  argument parsing, the analysis directory, and logging. It should call, not implement.
- **The CLI is also a staging area.** New or unfinished functions live there first and move into the
  package once they have settled. Treat a CLI-local function as provisional by default: the question
  is *when* it moves, not *whether*.

**Shadowing is a symptom of a pending move, not a pattern to keep.** When a CLI function has the
same name as a package export, that means the CLI holds extra behavior the package function should
absorb. The fix is to move that behavior into the package function and delete the CLI copy — **not**
to rename the CLI copy around the collision.

### The S3 dual-dispatch pattern

This is the core idiom. Each public verb is a generic with two methods: `.default` runs the
actual algorithm on a plain matrix/data frame and knows nothing about Seurat; `.Seurat` pulls the
relevant matrix/embedding/metadata out of the object, calls the generic again (landing in
`.default`), then wraps the result back into the object. `merge_metadata()` and `filter_cells()` are
clear examples of the split. When adding a function that touches expression data, follow it.

**Do not add dispatch Seurat already provides.** Check the underlying Seurat function's `methods()`
before reaching for this pattern. `NormalizeData()` and `SCTransform()` already dispatch on a matrix
and a Seurat object, so `normalize_counts()` is not dual-dispatch.

**`dr.R` is the exception worth knowing.** The `.Seurat` methods of `run_pca`/`run_tsne`/`run_umap`
call the *native* Seurat functions rather than recursing into `.default`. For a Seurat object,
Seurat computes the reduction. **`assay` is always passed explicitly, never left `NULL`**, so the
reduction never silently depends on `DefaultAssay(x)` at call time. `cluster_seurat_object()` and
`create_seurat_object()` take the same stance.

**Plotting is not optional and is not a parameter.** Every `.Seurat` method in this package writes
its standard output unconditionally, every time, flat into the *working directory* (no
subdirectory). See each function's own roxygen `@return` for what it writes. The test suite runs
from `tests/testthat/`, so this output lands there. `tests/testthat/teardown.R` sweeps up the known
filenames after the whole suite.

### Conventions

[create_seurat_object()] runs the whole load → filter → normalize → reduce → cluster chain in one
call. Prefer it when building an object from scratch; use the individual verbs when you need to
intervene between steps.

**Style: this package follows the [tidyverse style guide](https://style.tidyverse.org/)**, formatted
by `air format .`.

**Roxygen prose stays short and non-repetitive.** A `@param` entry states what the parameter is.
Explain more only when that is not obvious from the name and type alone. Most `@param` entries
should be describable in one clear sentence — resist adding a second sentence, caveat, or backstory
unless a caller actually needs it to call the function correctly (a parameter's history, or why an
*old* design would have been wrong, is not that). Cut any word or clause that doesn't add information the reader doesn't
already have. When several params share one reason (e.g. several have no default because each
describes how a specific object was built), state it once and let the rest give only their own
fact — don't repeat "same reasoning as X" on every one, and don't restate what's already obvious
from what was just said (e.g. no "either form works" right after listing the accepted forms).

**Naming conventions**, applied deliberately — match them when adding functions:

- `read_*` reads from disk; `initialize_*` is the bare constructor (a raw-counts object, not useful
  on its own); `create_*` returns the usable object.
- **`object` in function names, `obj` in variable names.** Functions spell it out
  (`merge_seurat_objects()`); locals stay short (`seurat_obj`, `so`). The serialized file is
  `seurat_obj.qs2` — a *filename*, so never rename it with a blanket `_obj` -> `_object`
  substitution or every existing analysis dir is orphaned.
- **Don't split one operation into two functions just because the input type differs.**
  `merge_seurat_objects()` takes a list that may mix in-memory Seurat objects and analysis-directory
  paths, resolving whichever it got itself.
- `split_layers_by_batch()` (package) splits the assay **layers inside one object** for
  `IntegrateLayers()`. The CLI's `split_seurat_object()` splits a cohort into **separate analysis
  dirs**, each with its own `seurat_obj.qs2`. Different outputs, similar names.

Cross-cutting conventions:

- Adding a metadata column in `initialize_seurat_object()` shifts the hardcoded `ncol()`
  expectations in `test-utils.R`.
- **Cell IDs are `sample_name:barcode`.** The prefix is applied at load time and Seurat objects are
  created with `names.delim = ":"`, so `orig.ident` recovers the sample. Trailing `-1` from Cell
  Ranger is stripped.
- **`log_file = NULL` is threaded through nearly every function** and passed to `write_message()`,
  which messages to console and optionally appends to the log. Keep the parameter. Don't
  reintroduce a bare `write()` to the log. **Each CLI step writes its own log file, named for that
  step** (`create.log`, `merge.log`, `integrate.log`, `cluster.log`, `hto.log`). `identify`,
  `plot umap` and `de` write no step log file yet.
- Colors come from `get_color_scheme("clusters"|"samples")`, not ad-hoc palettes.
- **Name a single-object first parameter `x`** (Seurat object, matrix, or data frame), as base R
  S3 methods do. A list of objects, or a function with two objects, keeps descriptive names.
- **Add a parameter only if the function body uses it.** Do not add one to match sibling functions.
- **Separator banners use `message()`, not `write_message()`.** Write each fact to the log with its
  own `write_message()` call. Do not join facts into one multi-line `glue()`.
- **Prefer intermediate variables to `|>` chains in new code.** This helps debugging. Do not remove
  existing pipes without a request.

### Sample integration

**The integrated reduction is named bare `int_reduction` (`"cca"`/`"rpca"`/`"harmony"`), not
`"integrated.<method>"`.** `CCAIntegration()`/`RPCAIntegration()`/`JointPCAIntegration()` all default
`new.reduction` to the identical `"integrated.dr"`, while `HarmonyIntegration()` defaults to bare
`"harmony"`. So there is no consistent Seurat convention to match, and the bare names are unambiguous
in this codebase.

### Pipeline entry point

The package is normally driven through `inst/extdata/scrna-scooter.R`, a docopt CLI, not from R.
A representative single-sample run:

```bash
./scrna-scooter.R create out-$sample_name --sample_name=$sample_name --path=$sample_path \
  --min_genes=500 --mt=15 && \
./scrna-scooter.R cluster out-$sample_name --num_dim=30 --reduction=pca && \
./scrna-scooter.R plot umap out-$sample_name --genes_csv=genes.core.csv &
```

Each step is a separate process reading and rewriting `out-$sample_name/seurat_obj.qs2` (serialized
with `qs2::qs_save()`/`qs_read()`, not `saveRDS`), which is why they are `&&`-chained. `identify
<dir> --cluster_var=<x>` is the usual next step; `merge`, `integrate`, `de`, and `hto` are optional
branches. Everything runs relative to the analysis dir (the script `setwd()`s into it), but
`--path` and `--genes_csv` are `normalizePath()`d beforehand, so relative paths to those are
resolved against your shell's cwd.

**`create`** — The step `stop()`s if the analysis dir already exists; a rerun needs the old dir
moved or removed first. `cluster` / `plot` / `identify` instead require it to exist. The tSNE/UMAP
computed here are a first-pass preview, separate from the real one `cluster --num_dim=<n>` computes.

- **`add_adt_assay_qc()` is the one CLI function `create` still calls, and it is not just QC.** It
  runs the ADT assay's CLR normalization — nothing else in the pipeline normalizes ADT at all.
  Its cell-matching logic is now inert, because `add_seurat_assay()` already intersected the cells.
- **`filter_cells()` has no targeted-panel exemption.** Set `min_genes`/`min_counts` explicitly for
  your data's scale, or `NULL`. The `inst/extdata` 10x fixture tops out at ~250 genes/cell, so
  fixture-based filter tests need an explicit low `min_genes` or `min_genes = NULL`.

**`cluster`** — The CLI's DR steps call the `scooter::run_*` verbs, not `Seurat::Run*` directly. So
`cluster` and `integrate` UMAP go through `uwot::umap2()` + euclidean, matching the package. **Leiden
returns 1-based labels while SLM returns 0-based**, so the `C01`/`C02`/… relabelling shifts from the
observed minimum rather than assuming either; changing `algorithm` must not break this. The `res.X`
column names are what you later pass as `--cluster_var` to `identify`.

**`identify --cluster_var=<x>`** — works in a `clusters-<var>-clust<n>/` subdirectory (created and
`setwd()`'d into first). `calculate_cluster_markers()` runs standard `wilcox`, then pairwise `wilcox`
too if under 20 clusters. It calls `plot_cluster_markers_heatmap()` and `plot_cluster_markers_top()`
itself, so `identify` never calls those two directly.

**`hto`** — both steps are CLI-local: `split_adt_hto_assay()` demultiplexes and calls `plot_hto_qc()`
itself; the CLI then keeps only `hash.class == "Singlet"` cells and hands off to
`split_seurat_object()`, which writes one `split-<value>/` subdirectory per `library_hash` value
(skipping any with 10 or fewer cells). Each split gets its own fresh normalize/PCA/tSNE/UMAP pass,
because a `DietSeurat()` earlier stripped `scale.data`/reductions. **`split_seurat_object()`'s
`dr_num_dim`/`dr_num_neighbors` params are never reachable from the CLI**: the `hto` call site
doesn't pass them, so `dr_num_dim = 0` skips clustering.

**`plot umap (--genes_csv=<x> | --cluster_var=<x>)`** — docopt makes the two mutually exclusive.
`--genes_csv` calls the CLI-local `plot_dr_umap_genes()`. **It silently switches `DefaultAssay` to
ADT** when the gene list overlaps the ADT assay better than RNA — surprising if a sample has ADT and
you meant RNA. `--cluster_var` calls the package function `plot_dr_umap_clusters()`.

**`de --cluster_var=<x> --group_var=<x>`** — shares the `identify`/`de` block's setup and
`clusters-<var>-clust<n>/` subdirectory. A lone `de` run without `identify` still creates and moves
into it.

**No shadow functions remain in the CLI.** `merge`'s and `split_seurat_object()`'s sequences call
`normalize_counts()` + `run_pca()` from scratch before the tSNE/UMAP preview, since the object
reaching them has no `scale.data`/PCA of its own yet. The two blocks are similar glue, left inline
rather than wrapped in a function. Revisit only if a third call site needs the same sequence.

### Files that are not live code

- **`R/functions-unprocessed.R`** — entirely commented out. Legacy CLI functions kept as a porting
  backlog. Don't treat it as an API; mine it only when reimplementing something with no live
  counterpart yet.
- **`R/deprecated.R`** holds thin wrappers, one per package function actually *renamed* (not
  deleted/absorbed). Each wrapper only forwards arguments — no logic lives here, ever; fix the
  replacement function instead. A function *absorbed* into another has no single replacement to
  point a `deprecate_warn()` at, so it gets no shim.
- **Before adding a shim for a rename, check `git status` on the file the old name lives in.** A
  shim is only warranted if the old name was ever in a *committed* version — otherwise there is no
  possible external caller to protect, and the shim is pure cruft. `M` (modified, previously
  committed) means the name was real and needs one; `??` (untracked, never committed) means it
  didn't.

### Move backlog (CLI staging -> package)

Moving a function means its extra steps and plots become part of the package verb it wraps, and the
CLI copy is **deleted, not renamed** — never rename the CLI copy just to dodge a collision with the
package export it should eventually absorb into.

Recurring pattern worth knowing before tackling what's left: a CLI-only global (like
`colors_clusters`) usually has to become a real parameter first (`color_scheme`, defaulting to
`get_color_scheme("clusters")`); a moved function's own dependencies (like `pheatmap`) may need
adding to `Imports`; and moving a function into the package can activate an existing same-name
shadow that was only ever harmless because nothing called it unqualified from inside the package.

Remaining:

1. `add_adt_assay_qc()` — ADT-specific CLR normalization with no package home.
2. `split_seurat_object()` — the last composite in staging. `R/split.R` is its natural home.
3. The rest of the ADT/HTO block — `plot_adt_qc()`, `split_adt_hto_assay()`, `plot_hto_qc()`, the
   least settled.

### Deferred / planned

- **`suffix` (`run_pca()`/`run_tsne()`/`run_umap()`) never propagates between calls.** It always
  defaults to `""`, and nothing reads it back off a computed reduction. So the output reduction name
  is always literally `"tsne"`/`"umap"`. Calling either twice on one object (e.g. a first-pass
  preview, then the real one) silently overwrites the previous reduction's *data*. The plot filenames
  avoid this by baking assay/reduction/num_dim/group_by into `file =`; the object storage has no
  equivalent. Flagged, not fixed; leave `suffix`'s blank default alone unless asked.

### Non-obvious gotchas

- **The CLI's `create` step calls `create_seurat_object()` directly, not through `do.call()`, on
  purpose.** A `do.call()` at the top level of the call stack makes `SeuratObject`'s
  `LogSeuratCommand()` misattribute commands logged inside `create_seurat_object()` (e.g.
  `FindVariableFeatures()`) to the `do.call()`. The correct `FindVariableFeatures.<assay>` entry
  never gets recorded, so `normalize_counts()`'s later `HVFInfo()` call fails with "Please run
  either 'FindVariableFeatures' or 'SCTransform'". `create_seurat_object()` accepts a length-zero
  `num_dim` and falls back to its own default, so the CLI needs no `do.call()`.
- **`plot_dr_group()` sleeps one second between writes on purpose.** It keeps the modification times
  of the output files distinct so they sort in the order they were produced. Removing it as an
  "optimization" silently scrambles that ordering.
- **Don't build `filter_cells()`'s second `subset()` from only the default assay's features.**
  `subset(features =)` deletes any assay with no feature in that list, so a bare
  `subset(features = names(detected))` built from `rowSums()` on the default assay alone silently
  drops the ADT/HTO assays (their `nCount_*` metadata columns survive, which is what makes it hard to
  spot). Carry the non-default assays' features through, the same guard `merge_seurat_objects()`
  already has.
- **The CLI no longer writes `counts.*.csv.gz` at all.** The `save_counts()` calls in `create`,
  `merge` and `integrate` are commented out on purpose — the tables are rarely used and are
  expensive to produce. Run `scooter::save_counts()` by hand against the saved `.qs2` object when
  one is actually needed. Do not "fix" this by uncommenting them.
- **`cluster.R` does not shadow `Seurat::FindAllMarkers()`.** The package-wide `@import Seurat` makes
  the real one resolve unqualified — don't reintroduce a local copy.
- `FindClusters.default()` names its columns `res.<resolution>` directly — the `seurat_snn_res.`
  prefix only comes from the `.Seurat` method.
- `run_umap.default()` (the matrix path) uses `uwot::umap2()` and **must pin `nn_method = "annoy"`**:
  a direct `umap2()`'s automatic nearest-neighbor choice is not reproducible under `set.seed()`, so
  dropping that argument silently makes every UMAP run produce a different embedding. The `.Seurat`
  path calls `Seurat::RunUMAP(umap.method = "uwot2")`, which stays reproducible via `seed.use`
  **without** the pin (verified on large data).
- **`normalize_counts()`, `cluster_seurat_object()`, and `integrate_layers()` each force
  `future::plan(future::sequential)` unconditionally before doing anything else** (each function's
  own roxygen states the conclusion; the benchmark data behind it lives only here). Benchmarked
  against real objects (1.7k-224k cells), and the answer differs by function:
  - `normalize_counts()`: sequential ties or wins at every size tested; `multisession` is up to 6.5x
    *slower*, not faster, and can fail outright somewhere around 100k cells (`future.globals.maxSize`,
    or a memory-pressure kill — not a fixed cell-count threshold; it depends on the memory available
    to the job).
  - `calculate_clusters()`: `multisession` can be meaningfully faster on mid-size objects, but it
    crashed (worker OOM-killed) at 50k cells in this session's memory allocation. Merged/integrated
    objects are exactly where the pipeline tends to reach this scale, so unconditional sequential was
    kept anyway, trading away the speedup for predictability.
  - `integrate_layers()`: no measured difference either way, up to ~50k cells across 8 batches.
  The guards live inside each function, not in a shared wrapper, because their call graphs don't
  funnel through one. `integrate_layers()` is a standalone entry point for interactive use. The
  guards are deliberate even when nothing upstream sets a plan — defense against a caller having its
  own unrelated multisession plan set. `FindNeighbors()`'s own nearest-neighbor step was never at
  risk: it force-sets its own sequential plan internally and restores the caller's plan on exit.
- **`mad_outlier_bounds()` runs over all cells on purpose.** Debris cells widen the MAD and so
  loosen the upper bound. The bounds match `scater::isOutlier()` exactly. Do not add a debris
  pre-filter. Use `--max_genes` to remove extreme cells.
- **Per-sample splits use up to 30 dimensions. Merged and pooled analyses use up to 50.** Single
  samples are more uniform. Do not unify the two caps.
