#' Integrate a merged Seurat object across batches.
#'
#' The integrate step of the pipeline in one call. It plays the same role for
#' a merged object that [create_seurat_object()] plays for a single sample.
#'
#' The function runs these steps:
#'
#' 1. Split the assay layers by batch and prepare them
#'    ([split_layers_by_batch()]).
#' 2. Report how much the per-batch variable genes overlap
#'    ([plot_var_genes_euler()], [plot_var_genes_upset()]).
#' 3. Correct the batch effect ([integrate_layers()]).
#' 4. Compute a first-pass tSNE/UMAP on the corrected embedding.
#'
#' The tSNE/UMAP from step 4 is a preview, capped at `num_dim`. It differs
#' from the real one [cluster_seurat_object()] computes later over its own
#' `num_dim`. [cluster_seurat_object()] has no default `reduction`. Pass it
#' this function's `int_reduction` name explicitly. Otherwise it may use the
#' uncorrected "pca" reduction that [split_layers_by_batch()] also leaves on
#' the object.
#'
#' @param x Merged Seurat object, such as the output of
#'   [merge_seurat_objects()].
#' @param assay Assay to integrate. Passed straight through to
#'   [integrate_layers()].
#' @param int_reduction Integration method: "cca", "rpca", or "harmony".
#' @param num_dim Number of dimensions to use (5-50).
#' @param batch_var Metadata column identifying the batch.
#' @param num_variable_genes Number of variable features per batch.
#' @param num_neighbors Neighbors for UMAP.
#' @param k_anchor Neighbors to use when picking anchors. Passed straight
#'   through to [integrate_layers()] - see there for defaults/caveats (ignored
#'   by "harmony").
#' @param k_weight Neighbors to use when weighting the corrections. Passed
#'   straight through to [integrate_layers()] - see there for defaults/caveats
#'   (also ignored by "harmony").
#' @param log_file Filename for the log file.
#' @param ... Passed straight through to [integrate_layers()] - extra
#'   method-specific tuning arguments (harmony's `theta`/`lambda`/`sigma`/...,
#'   or cca/rpca's `k.filter`/`sample.tree`/...).
#'
#' @return A Seurat object with rejoined layers, an `int_reduction`-named
#'   reduction (`"cca"`, `"rpca"`, or `"harmony"`), and "tsne"/"umap"
#'   reductions computed from it. Also writes the variable-gene overlap, QC
#'   distribution, integrated-embedding, and tSNE/UMAP plots to the working
#'   directory.
#'
#' @importFrom grDevices dev.off png
#' @importFrom Matrix colSums rowSums
#' @export
integrate_seurat_object <- function(
  x,
  assay = "RNA",
  int_reduction = "cca",
  num_dim,
  batch_var = "orig.ident",
  num_variable_genes = 3000,
  num_neighbors = 30,
  k_anchor = 10,
  k_weight = 100,
  log_file = NULL,
  ...
) {
  message("\n\n ===== integrate ===== \n\n")
  write_message(glue("object cells: {ncol(x)}"), log_file)
  write_message(glue("object genes: {nrow(x)}"), log_file)

  # clean up the object, split its layers by batch, and normalize/scale/PCA
  # each batch
  x <- split_layers_by_batch(
    x,
    batch_var = batch_var,
    num_variable_genes = num_variable_genes,
    num_dim = num_dim,
    log_file = log_file
  )

  var_genes_list <- variable_features_by_batch(x)

  # eulerr and UpSetR return base plots, so they go through a device rather
  # than ggsave()
  euler_plot <- plot_var_genes_euler(var_genes_list)
  if (!is.null(euler_plot)) {
    png(
      "variance-vargenes-euler.png",
      res = 200,
      width = 5,
      height = 5,
      units = "in"
    )
    print(euler_plot)
    dev.off()
  }

  png(
    "variance-vargenes-upset.png",
    res = 200,
    width = 8,
    height = 5,
    units = "in"
  )
  print(plot_var_genes_upset(var_genes_list))
  dev.off()

  # correct the batch effect in the reduced dimensional space (adds an
  # `int_reduction`-named reduction and rejoins the layers)
  message("\n\n ===== integrate_layers() ===== \n\n")
  x <- integrate_layers(
    x,
    num_dim = num_dim,
    int_reduction = int_reduction,
    assay = assay,
    k_anchor = k_anchor,
    k_weight = k_weight,
    log_file = log_file,
    ...
  )

  # encode batch name as factor (also sets alphabetical sample order)
  x@meta.data[[batch_var]] <- factor(x@meta.data[[batch_var]])

  counts <- GetAssayData(x, assay = assay, layer = "counts")
  write_message(
    glue("min cells per gene: {min(Matrix::rowSums(counts > 0))}"),
    log_file
  )
  write_message(
    glue("min genes per cell: {min(Matrix::colSums(counts > 0))}"),
    log_file
  )
  write_message(
    glue("min counts per cell: {round(min(Matrix::colSums(counts)), 3)}"),
    log_file
  )

  # integration always combines 2+ batches, so pct_ribo/pct_hb are left out -
  # not worth the width
  plot_metrics_distribution(
    x,
    metrics = c("detected_genes", "total_counts", "pct_mito"),
    file = "metrics-distribution.png",
    width = metrics_plot_width(x)
  )

  # plot the integrated embedding - IntegrateLayers() already produced the
  # batch-corrected embedding, so scaling and re-running PCA here would discard
  # the correction. Plot both batch_var and orig.ident when they differ (e.g.
  # batching on something coarser than sample).
  for (group_by in unique(c(batch_var, "orig.ident"))) {
    plot_dr_group(
      x,
      reduction = int_reduction,
      group_by = group_by,
      color_scheme = get_color_scheme("samples"),
      file_prefix = glue(
        "dr-{int_reduction}-{DefaultAssay(x)}-{num_dim}-{group_by}"
      ),
      file_format = "png"
    )
  }

  message("\n\n ===== run_tsne() ===== \n\n")
  # use tSNE as a tool to visualize, not for clustering directly on tSNE
  # components
  x <- run_tsne(
    x,
    reduction = int_reduction,
    num_dim = num_dim,
    assay = DefaultAssay(x),
    file_format = "png"
  )
  for (group_by in unique(c(batch_var, "orig.ident"))) {
    plot_dr_group(
      x,
      reduction = "tsne",
      group_by = group_by,
      color_scheme = get_color_scheme("samples"),
      file_prefix = glue(
        "dr-tsne-{DefaultAssay(x)}-{int_reduction}-{num_dim}-{group_by}"
      ),
      file_format = "png"
    )
  }

  message("\n\n ===== run_umap() ===== \n\n")
  x <- run_umap(
    x,
    reduction = int_reduction,
    num_dim = num_dim,
    num_neighbors = num_neighbors,
    assay = DefaultAssay(x),
    file_format = "png"
  )
  for (group_by in unique(c(batch_var, "orig.ident"))) {
    plot_dr_group(
      x,
      reduction = "umap",
      group_by = group_by,
      color_scheme = get_color_scheme("samples"),
      file_prefix = glue(
        "dr-umap-{DefaultAssay(x)}-{int_reduction}-{num_dim}-{group_by}"
      ),
      file_format = "png"
    )
  }

  save_metadata(x)

  x
}

#' Integrate the layers of a Seurat object.
#'
#' Runs the Seurat 5 `IntegrateLayers()` workflow. It adds what calling
#' `IntegrateLayers()` directly leaves for the caller to handle:
#'
#' - Clamps `k_weight` to fit the smallest batch. Seurat's own default of 100
#'   errors outright below that size. See `k_weight` below.
#' - Excludes `dims`/`k.anchor` for `"harmony"`. `HarmonyIntegration()`
#'   silently ignores both instead of erroring, so passing them would rely
#'   on that silent drop. See `num_dim`/`k_anchor` below.
#' - Forces `future::plan(future::sequential)`. Benchmarks show multisession
#'   gives no gain here, and is a documented liability at this scale
#'   elsewhere in this package.
#' - Rejoins the split layers with `JoinLayers()` afterward.
#'   `IntegrateLayers()` itself leaves this undone.
#'
#' Unlike the older anchor workflow, this function creates no "integrated"
#' assay. Clustering and visualization use the corrected embedding directly,
#' instead of re-running PCA. This function expects
#' [split_layers_by_batch()] to have run first. It is the one step of
#' [integrate_seurat_object()] that calls `IntegrateLayers()` directly.
#'
#' @param x Seurat object with split layers and a "pca" reduction.
#' @param assay Assay to integrate.
#' @param int_reduction Integration method: "cca", "rpca", or "harmony".
#' @param num_dim Number of dimensions to use (5-50). For "cca"/"rpca", this
#'   value is passed as `dims = 1:num_dim` to
#'   `CCAIntegration()`/`RPCAIntegration()`, which both take a real `dims`
#'   argument. `HarmonyIntegration()` has no `dims` parameter. It runs on the
#'   full `"pca"` embedding as-is, through `Embeddings(orig)`, unsliced. So
#'   for `int_reduction = "harmony"`, `num_dim` only takes effect indirectly,
#'   through however many components `orig.reduction` ("pca") already has.
#'   The normal [integrate_seurat_object()] pipeline is unaffected:
#'   [split_layers_by_batch()] already computes that "pca" reduction with
#'   exactly `num_dim` components. A direct call is different. If `x`'s
#'   "pca" reduction already has more components than `num_dim`, then
#'   `num_dim` has no effect on the harmony path.
#' @param k_anchor Neighbors to use when picking anchors. Ignored by "harmony".
#'   A value larger than Seurat's default of 5 strengthens the alignment.
#' @param k_weight Neighbors to use when weighting the corrections. Must be
#'   smaller than the number of cells in the smallest batch. Reduced to 25 when
#'   the smallest batch has fewer than 100 cells. Batches smaller than 25 cells
#'   are an error.
#' @param log_file Filename for the log file.
#' @param ... Extra arguments for whichever method `IntegrateLayers()`
#'   dispatches to - e.g. `theta`/`lambda`/`sigma`/`tau`/`nclust` for
#'   "harmony", or `k.filter`/`sample.tree`/ `normalization.method` for
#'   "cca"/"rpca". An unrecognized name errors (via [rlang::check_dots_used()])
#'   rather than being silently dropped by the method function's own `...`.
#'
#' @return A Seurat object with rejoined layers and an `int_reduction`-named
#'   reduction (`"cca"`, `"rpca"`, or `"harmony"`).
#'
#' @importFrom rlang check_installed check_dots_used
#' @export
integrate_layers <- function(
  x,
  assay = "RNA",
  int_reduction = "cca",
  num_dim,
  k_anchor = 10,
  k_weight = 100,
  log_file = NULL,
  ...
) {
  num_dim <- as.integer(num_dim)
  if (num_dim < 5) {
    stop("too few dims: ", num_dim)
  }
  if (num_dim > 50) {
    stop("too many dims: ", num_dim)
  }

  int_methods <- list(
    cca = CCAIntegration,
    rpca = RPCAIntegration,
    harmony = HarmonyIntegration
  )
  if (!int_reduction %in% names(int_methods)) {
    stop(glue("integration reduction type {int_reduction} is not valid"))
  }

  # fail before the expensive steps rather than inside IntegrateLayers()
  if (int_reduction == "harmony") {
    check_installed("harmony", reason = "to integrate with Harmony.")
  }

  # k.weight has to fit inside the smallest batch or IntegrateLayers() errors
  # out read the batch sizes off the layers themselves: the object may have
  # been split by any metadata column, so orig.ident is not necessarily the
  # batch
  counts_layers <- SeuratObject::Layers(x[[assay]], search = "counts")
  smallest_batch <- min(vapply(
    counts_layers,
    function(lyr) length(SeuratObject::Cells(x[[assay]], layer = lyr)),
    numeric(1)
  ))
  if (smallest_batch < 25) {
    stop(glue(
      "smallest batch has only {smallest_batch} cells: too few to integrate"
    ))
  }

  # k.weight has to fit inside the smallest batch, and the default of 100 is
  # tuned for large ones
  if (smallest_batch < 100) {
    k_weight <- 25
    write_message(
      glue(
        "smallest batch is {smallest_batch} cells: using k_weight {k_weight}"
      ),
      log_file
    )
  }

  # IntegrateLayers()'s cca/rpca methods call
  # FindIntegrationAnchors()/FastRPCAIntegration() internally. Both switch to
  # a future_lapply()-parallel path once nbrOfWorkers() > 1. This stays
  # sequential in case a caller has an unrelated multisession plan set.
  # Benchmarks show no benefit up to ~50k cells/8 batches (731 vs 746 sec).
  # See normalize_counts() (R/normalize.R) for where multisession is an
  # active liability instead.
  future::plan(future::sequential)

  int_args <- list(
    object = x,
    method = int_methods[[int_reduction]],
    orig.reduction = "pca",
    new.reduction = int_reduction,
    assay = assay,
    k.weight = k_weight,
    verbose = FALSE
  )
  # harmony takes neither anchors nor dims. HarmonyIntegration() has no
  # `dims` formal at all, so passing one falls into its own unused `...` and
  # is silently dropped instead of erroring. This is the same failure mode
  # as an unmatched `...` argument. See the roxygen note on `num_dim` above
  # for what controls its effective dimensionality instead.
  if (int_reduction != "harmony") {
    int_args$dims <- 1:num_dim
    int_args$k.anchor <- k_anchor
  }

  # extra method-specific tuning (harmony's theta/lambda/sigma/..., or
  # cca/rpca's k.filter/ sample.tree/...) arrives via `...` rather than named
  # params here - there is too much of it, disjoint per method, to hand-copy
  # and keep in sync with Seurat/harmony's own defaults
  int_args <- c(int_args, list(...))

  # IntegrateLayers() always computes a `features` default, even though
  # int_args never sets one. It forwards that default to whichever method
  # function is chosen. HarmonyIntegration() ignores the argument and
  # reports so once per session via rlang::inform(). The message is
  # harmless. int_args cannot prevent it by leaving `features` unset, since
  # the default comes from IntegrateLayers() itself.
  integrated_obj <- suppressMessages(do.call(IntegrateLayers, int_args))

  # catches a `...` argument no method function ever touched: a typo, or an
  # argument that only applies to a different int_reduction. Without this
  # check, that argument would be silently dropped instead of erroring.
  check_dots_used()

  # the per-batch layers are no longer needed once the correction is computed,
  # and leaving them split makes GetAssayData() error out downstream
  integrated_obj <- SeuratObject::JoinLayers(integrated_obj, assay = assay)

  write_message(glue("integrated cells: {ncol(integrated_obj)}"), log_file)
  write_message(glue("integrated genes: {nrow(integrated_obj)}"), log_file)
  write_message(glue("integrated reduction: {int_reduction}"), log_file)

  return(integrated_obj)
}

#' Variable features of each batch of a layer-split object.
#'
#' After [split_layers_by_batch()] the per-batch variable features are stored
#' as `vf_<method>_<layer>_variable` columns in the assay meta data rather than
#' being reachable through `VariableFeatures()`, which returns the combined
#' set.
#'
#' [normalize_counts()] now clears its own assay's `vf_*` columns before
#' recomputing them. But an object built before that guard existed can still
#' carry a whole history of them: each sample's own pre-merge run from
#' [create_seurat_object()], and the merge step's own post-merge preview
#' run. `DietSeurat()` clears layers, reductions, and graphs, but never
#' assay meta data. `merge()` disambiguates colliding names with a numeric
#' suffix, for example `vf_vst_counts.1_variable`.
#'
#' Matching every `_variable$` column would pick up that whole history along
#' with the batches actually being asked about. This function instead
#' anchors on `Layers(x[[assay]], search = "counts")`, the object's current
#' split state. That keeps stale columns out, regardless of which
#' generation of the object this is.
#'
#' @param x Seurat object with split layers.
#' @param assay Assay to read the variable features from.
#'
#' @return A named list of variable gene vectors, one element per batch.
#'
#' @export
variable_features_by_batch <- function(x, assay = "RNA") {
  feature_meta <- x[[assay]][[]]

  # counts.SAMPLE layer names are the current split state's own source of truth
  # for which batches exist right now, unlike the accumulated set of
  # vf_..._variable columns
  count_layers <- SeuratObject::Layers(x[[assay]], search = "counts")
  batch_names <- sub("^counts\\.", "", count_layers)
  batch_names <- setdiff(batch_names, "counts")

  if (!length(batch_names)) {
    stop(
      "no per-batch variable features found: run split_layers_by_batch() first"
    )
  }

  var_genes_list <- lapply(batch_names, function(batch_name) {
    # match on the layer suffix only (method name/prefix varies) so a stale
    # column from an earlier normalize_counts() run - which never shares this
    # exact "counts.<batch>" suffix - can't be mistaken for the current one
    suffix <- paste0("counts.", batch_name, "_variable")
    column_name <- colnames(feature_meta)[endsWith(
      colnames(feature_meta),
      suffix
    )]
    if (length(column_name) != 1) {
      stop(glue(
        "expected exactly one variable-feature column for batch ",
        "'{batch_name}', found {length(column_name)}"
      ))
    }
    rownames(feature_meta)[which(feature_meta[[column_name]])]
  })
  names(var_genes_list) <- batch_names

  return(var_genes_list)
}
