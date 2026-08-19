#' Run dimensionality reduction: pca, tsne, or umap
#'
#' @details
#' `run_dr()` is a thin dispatcher: it routes to `run_pca()`, `run_tsne()`, or
#' `run_umap()` based on `dr_method`, forwarding `suffix` and anything in
#' `...`. It only carries the parameters common to every method;
#' method-specific ones (`num_pcs`, `num_neighbors`, `reduction`, ...) are
#' passed through `...` to the chosen function. Each of those functions is
#' itself generic: a plain matrix takes the `.default` method
#' (`irlba`/`Rtsne`/`uwot`), a Seurat object takes the `.Seurat` method (native
#' `RunPCA()`/`RunTSNE()`/`RunUMAP()`, storing a `"<method><suffix>"`
#' reduction).
#'
#' @param x A matrix or a Seurat object.
#' @param dr_method Dimensionality reduction method ("pca", "tsne", or "umap").
#' @param suffix Tag inserted into the reduction name and key.
#' @param assay Assay to use (Seurat methods). Always passed explicitly here
#'   rather than left `NULL`, which would otherwise fall back on whatever
#'   `DefaultAssay(x)` happens to be.
#' @param ... Passed to the chosen reduction function.
#'
#' @return The reduced data (matrix input) or a Seurat object with the
#'   reduction added. For a Seurat object, also writes the diagnostic plot(s)
#'   for the chosen method to the working directory: PCA gets a sample-colored
#'   scatter and an elbow plot (plus a heatmap past 15 components); tSNE/UMAP
#'   get a sample-colored scatter.
#'
#' @export
#' @seealso run_pca(), run_tsne(), run_umap()
run_dr <- function(x, dr_method, suffix = "", assay = "RNA", ...) {
  switch(
    dr_method,
    pca = run_pca(x, suffix = suffix, assay = assay, ...),
    tsne = run_tsne(x, suffix = suffix, assay = assay, ...),
    umap = run_umap(x, suffix = suffix, assay = assay, ...),
    stop("dr_method must be one of pca, tsne, umap")
  )
}

#' Run PCA
#'
#' The default method runs `irlba::prcomp_irlba()` on a matrix; the Seurat
#' method runs `Seurat::RunPCA()` on the object and stores a `"pca<suffix>"`
#' reduction.
#'
#' @param x A matrix (default) or a Seurat object.
#' @param num_pcs Number of principal components.
#' @param suffix Tag inserted into the reduction/column names.
#' @param assay Assay to use (Seurat method). Always set explicitly.
#' @param features Explicit features to use (Seurat method).
#' @param var_features Use the assay's variable features (Seurat method).
#' @param ... Passed to the underlying function.
#'
#' @return A named list (default) or a Seurat object with the reduction added
#'   (Seurat method). The Seurat method also writes the diagnostic plots to the
#'   working directory: a sample-colored scatter and an elbow plot, plus a
#'   heatmap past 15 components.
#'
#' @export
run_pca <- function(x, ...) {
  UseMethod("run_pca")
}

#' @rdname run_pca
#' @importFrom irlba prcomp_irlba
#' @importFrom rlang check_dots_used
#' @export
run_pca.default <- function(x, num_pcs = 50, suffix = "", ...) {
  # make sure the number of pcs isn't too large
  npcs <- min(num_pcs, nrow(x))

  # calculate PCA
  pca_out <- irlba::prcomp_irlba(x, n = npcs, verbose = FALSE)

  # Extract information
  feature.loadings <- pca_out$rotation
  sdev <- pca_out$sdev
  cell.embeddings <- pca_out$x

  rownames(feature.loadings) <- colnames(x)
  colnames(feature.loadings) <- paste0("PC", suffix, "_", seq_len(npcs))
  rownames(cell.embeddings) <- rownames(x)
  colnames(cell.embeddings) <- paste0("PC", suffix, "_", seq_len(npcs))

  # this method takes no extra params at all - `...` exists only so
  # run_dr()/run_pca() can forward blindly, so anything landing here is always
  # a mistake, never a legitimate pass-through
  check_dots_used()

  return(list(
    feature.loadings = feature.loadings,
    cell.embeddings = cell.embeddings,
    sdev = sdev,
    pca_out = pca_out
  ))
}

#' @rdname run_pca
#' @importFrom rlang check_dots_used
#' @export
run_pca.Seurat <- function(
  x,
  num_pcs = 50,
  suffix = "",
  assay = "RNA",
  features = NULL,
  var_features = FALSE,
  ...
) {
  if (var_features) {
    features <- VariableFeatures(x[[assay]])
  }

  reduction <- glue("pca{suffix}")
  x <- RunPCA(
    x,
    assay = assay,
    features = features,
    npcs = num_pcs,
    reduction.name = reduction,
    reduction.key = glue("PC{suffix}_"),
    verbose = FALSE,
    ...
  )
  # Catches a `...` argument RunPCA() never touched: a typo, or a parameter
  # Seurat renamed. Without this check, it would be silently dropped
  # instead of erroring.
  check_dots_used()

  # read back off the computed reduction rather than `features` - that stays
  # NULL (and RunPCA() falls back to the assay's own VariableFeatures())
  # whenever neither `features` nor `var_features` is set, so it does not
  # always reflect what was actually used
  num_genes_used <- nrow(Loadings(x, reduction = reduction))

  plot_dr_group(
    x,
    reduction = reduction,
    group_by = "orig.ident",
    color_scheme = get_color_scheme("samples"),
    file_prefix = glue("dr-{reduction}-{assay}-{num_genes_used}"),
    file_format = "png",
    width = 10,
    height = 6
  )

  # exploring the primary sources of heterogeneity only makes sense with enough
  # PCs to look at
  if (num_pcs > 15) {
    png(
      glue("variance-{reduction}-heatmap.png"),
      res = 300,
      width = 10,
      height = 16,
      units = "in"
    )
    DimHeatmap(
      x,
      reduction = reduction,
      dims = 1:15,
      nfeatures = 20,
      cells = min(250, ncol(x)),
      fast = TRUE
    )
    dev.off()
  }

  elbow_plot <- ElbowPlot(x, reduction = reduction, ndims = num_pcs)
  ggsave(
    glue("variance-{reduction}-elbow.png"),
    plot = elbow_plot,
    width = 8,
    height = 5,
    units = "in"
  )

  x
}

#' Validate a run_tsne()/run_umap() Seurat-method input specification
#'
#' `reduction`+`num_dim`, `features`, `var_features`, and (UMAP only)
#' `graph` are one-of-several ways to tell `RunTSNE()`/`RunUMAP()` what to
#' reduce. Leaving all of them unset does not reliably error.
#'
#' A stale caller using a renamed parameter passes an unmatched or misnamed
#' argument. `...` silently absorbs it at every layer down to
#' `Rtsne()`/`uwot::umap2()`. So `reduction`, `num_dim`, and `features` stay
#' at their `NULL` defaults with no complaint. The call only errors if that
#' leftover combination happens to be invalid — Seurat's fallback
#' `stop("Unknown way of running tSNE")`. Otherwise it silently runs on the
#' wrong input.
#'
#' This function checks two things: that exactly one way is set, and that
#' whichever one is set actually exists on `x`. So a stale or misnamed
#' argument fails loudly here instead.
#'
#' @param x A Seurat object.
#' @param assay Assay `features`/`var_features` are checked against.
#' @param reduction Reduction name, or `NULL`.
#' @param num_dim Number of dimensions to take from `reduction`, or `NULL`.
#' @param features Explicit feature names, or `NULL`.
#' @param var_features Whether the assay's variable features are being used.
#' @param graph Graph name, or `NULL` (UMAP only - unused for tSNE).
#'
#' @noRd
validate_dr_input <- function(
  x,
  assay,
  reduction = NULL,
  num_dim = NULL,
  features = NULL,
  var_features = FALSE,
  graph = NULL
) {
  if (!assay %in% Assays(x)) {
    stop(glue("assay '{assay}' not found; available: {toString(Assays(x))}"))
  }

  if (!is.null(reduction) || !is.null(num_dim)) {
    if (is.null(reduction) || is.null(num_dim)) {
      stop(
        "reduction and num_dim must be set together, not one without the other"
      )
    }
  }

  ways <- c(
    reduction = !is.null(reduction),
    features = !is.null(features),
    var_features = isTRUE(var_features),
    graph = !is.null(graph)
  )
  if (sum(ways) != 1) {
    stop(glue(
      "specify exactly one of reduction+num_dim, features, var_features, ",
      "or graph ({sum(ways)} given)"
    ))
  }

  if (!is.null(reduction)) {
    if (!reduction %in% Reductions(x)) {
      stop(glue(
        "reduction '{reduction}' not found; ",
        "available: {toString(Reductions(x))}"
      ))
    }
    num_available <- ncol(Embeddings(x[[reduction]]))
    num_dim_num <- suppressWarnings(as.numeric(num_dim))
    if (
      length(num_dim_num) != 1 ||
        is.na(num_dim_num) ||
        num_dim_num < 1 ||
        num_dim_num > num_available
    ) {
      stop(glue(
        "num_dim ({num_dim}) must be a number between 1 and {num_available} ",
        "(the dimensions available in reduction '{reduction}')"
      ))
    }
  }

  if (!is.null(features)) {
    missing_features <- setdiff(features, rownames(x[[assay]]))
    if (length(missing_features)) {
      stop(glue(
        "features not found in assay '{assay}': {toString(missing_features)}"
      ))
    }
  }

  if (isTRUE(var_features) && !length(VariableFeatures(x[[assay]]))) {
    stop(glue(
      "assay '{assay}' has no variable features; run normalize_counts() first"
    ))
  }

  if (!is.null(graph) && !graph %in% SeuratObject::Graphs(x)) {
    stop(glue(
      "graph '{graph}' not found; ",
      "available: {toString(SeuratObject::Graphs(x))}"
    ))
  }

  invisible(NULL)
}

#' Run TSNE
#'
#' The default method runs `Rtsne::Rtsne()` on a matrix; the Seurat method runs
#' `Seurat::RunTSNE()` on the object and stores a `"tsne<suffix>"` reduction.
#'
#' @param x A matrix (default) or a Seurat object.
#' @param seed.use Seed to use.
#' @param dim.embed Number of tSNE embeddings to return.
#' @param suffix Tag inserted into the reduction/column names.
#' @param assay Assay to use (Seurat method). Always set explicitly.
#' @param reduction Existing reduction to take dimensions from (Seurat method).
#' @param num_dim Number of dimensions to take from `reduction` (Seurat method).
#' @param features Explicit features to use (Seurat method).
#' @param var_features Use the assay's variable features (Seurat method).
#' @param file_format File formats to write the scatter plot in (Seurat method).
#' @param ... Passed to the underlying function.
#'
#' @return A tSNE embedding matrix (default) or a Seurat object with the
#'   reduction added. The Seurat method also writes a sample-colored scatter of
#'   the result to the working directory, as `dr.<reduction>.png`/`.pdf` (or
#'   just `.png` if `file_format = "png"`).
#'
#' @export
run_tsne <- function(x, ...) {
  UseMethod("run_tsne")
}

#' @rdname run_tsne
#' @importFrom Rtsne Rtsne
#' @importFrom rlang check_dots_used
#' @export
run_tsne.default <- function(x, seed.use = 1, dim.embed = 2, suffix = "", ...) {
  set.seed(seed = seed.use)
  # run tsne
  tsne.data <- Rtsne(x, dims = dim.embed)$Y

  colnames(x = tsne.data) <- paste0(
    "tSNE",
    suffix,
    "_",
    seq_len(ncol(x = tsne.data))
  )
  rownames(x = tsne.data) <- rownames(x = x)

  # this method takes no extra params at all - see the matching note in
  # run_pca.default()
  check_dots_used()

  return(tsne.data)
}

#' @rdname run_tsne
#' @importFrom rlang check_dots_used
#' @export
run_tsne.Seurat <- function(
  x,
  seed.use = 1,
  dim.embed = 2,
  suffix = "",
  assay = "RNA",
  reduction = NULL,
  num_dim = NULL,
  features = NULL,
  var_features = FALSE,
  file_format = c("png", "pdf"),
  ...
) {
  validate_dr_input(
    x,
    assay = assay,
    reduction = reduction,
    num_dim = num_dim,
    features = features,
    var_features = var_features
  )
  if (var_features) {
    features <- VariableFeatures(x[[assay]])
  }
  # take num_dim dimensions from the reduction; unused inputs stay NULL and
  # RunTSNE picks the single one that is set (it errors if none or several are)
  dims <- if (is.null(num_dim)) NULL else seq_len(as.integer(num_dim))

  dr_name <- glue("tsne{suffix}")
  x <- RunTSNE(
    x,
    reduction = reduction,
    dims = dims,
    features = features,
    dim.embed = dim.embed,
    seed.use = seed.use,
    reduction.name = dr_name,
    reduction.key = glue("tSNE{suffix}_"),
    ...
  )
  # Catches a `...` argument RunTSNE() never touched: a typo, or a
  # parameter Seurat renamed. Without this check, it would be silently
  # dropped instead of erroring.
  check_dots_used()

  group_by <- "orig.ident"
  # a features-driven run (var_features/features set instead of
  # reduction+num_dim) has no assay/reduction/dims to report - that trio is
  # only meaningful on the reduction path
  dr_info <- if (is.null(reduction)) {
    ""
  } else {
    glue("-{assay}-{reduction}-{num_dim}")
  }
  plot_dr_group(
    x,
    reduction = dr_name,
    group_by = group_by,
    color_scheme = get_color_scheme("samples"),
    file_prefix = glue("dr-{dr_name}{dr_info}-{group_by}"),
    file_format = file_format
  )

  x
}

#' Run UMAP
#'
#' The default method runs `uwot::umap2()` on a matrix; the Seurat method runs
#' `Seurat::RunUMAP()` with `umap.method = "uwot2"` (also `uwot::umap2()`) and
#' stores a `"umap<suffix>"` reduction. Defaults match `Seurat::RunUMAP()` (30
#' neighbors, minimum distance 0.3) rather than the lower `uwot` defaults.
#'
#' @param x A matrix (default) or a Seurat object.
#' @param num_neighbors Number of neighbors.
#' @param min_dist Minimum distance between points in the embedding.
#' @param suffix Tag inserted into the reduction/column names.
#' @param assay Assay to use (Seurat method). Always set explicitly.
#' @param reduction Existing reduction to take dimensions from (Seurat method).
#' @param num_dim Number of dimensions to take from `reduction` (Seurat method).
#' @param features Explicit features to use (Seurat method).
#' @param var_features Use the assay's variable features (Seurat method).
#' @param graph Graph to use as input (Seurat method).
#' @param file_format File formats to write the scatter plot in (Seurat method).
#' @param ... Passed to the underlying function.
#'
#' @return A UMAP embedding matrix (default) or a Seurat object with the
#'   reduction added. The Seurat method also writes a sample-colored scatter of
#'   the result to the working directory, as `dr.<reduction>.png`/`.pdf` (or
#'   just `.png` if `file_format = "png"`).
#'
#' @importFrom uwot umap2
#' @export
run_umap <- function(x, ...) {
  UseMethod("run_umap")
}

#' @rdname run_umap
#' @importFrom rlang check_dots_used
#' @export
run_umap.default <- function(
  x,
  num_neighbors = 30,
  min_dist = 0.3,
  suffix = "",
  ...
) {
  # umap2() is uwot's current entry point, with the same arguments as
  # umap(). It defaults to batch = TRUE (embedding reproducible across
  # thread counts) and init_sdev = "range".
  #
  # nn_method must be pinned. umap2()'s automatic choice is not
  # reproducible under set.seed(). "annoy" is reproducible, and is the
  # approximate method that scales to large datasets.
  umap_out <- umap2(
    x,
    n_neighbors = num_neighbors,
    min_dist = min_dist,
    nn_method = "annoy",
    verbose = FALSE
  )

  colnames(umap_out) <- paste0("UMAP", suffix, "_", seq_len(ncol(umap_out)))
  rownames(umap_out) <- rownames(x)

  # this method takes no extra params at all - see the matching note in
  # run_pca.default()
  check_dots_used()

  return(umap_out)
}

#' @rdname run_umap
#' @importFrom rlang check_dots_used
#' @export
run_umap.Seurat <- function(
  x,
  num_neighbors = 30,
  min_dist = 0.3,
  suffix = "",
  assay = "RNA",
  reduction = NULL,
  num_dim = NULL,
  features = NULL,
  var_features = FALSE,
  graph = NULL,
  file_format = c("png", "pdf"),
  ...
) {
  validate_dr_input(
    x,
    assay = assay,
    reduction = reduction,
    num_dim = num_dim,
    features = features,
    var_features = var_features,
    graph = graph
  )
  if (var_features) {
    features <- VariableFeatures(x[[assay]])
  }
  # take num_dim dimensions from the reduction; unused inputs stay NULL and
  # RunUMAP picks the single one that is set (it errors if none or several are)
  dims <- if (is.null(num_dim)) NULL else seq_len(as.integer(num_dim))

  dr_name <- glue("umap{suffix}")

  # umap.method = "uwot2" calls uwot::umap2(). Unlike run_umap.default, the
  # nearest-neighbor method cannot be pinned here: RunUMAP does not expose
  # nn_method for the reduction/features path, so umap2 auto-selects it (HNSW
  # when RcppHNSW is installed, otherwise annoy). For reference, the clustering
  # neighbor step FindNeighbors() defaults to nn.method = "annoy" and
  # annoy.metric = "euclidean".
  x <- RunUMAP(
    x,
    reduction = reduction,
    dims = dims,
    features = features,
    graph = graph,
    assay = assay,
    umap.method = "uwot2",
    n.neighbors = num_neighbors,
    min.dist = min_dist,
    metric = "euclidean",
    seed.use = 1,
    verbose = FALSE,
    reduction.name = dr_name,
    reduction.key = glue("UMAP{suffix}_"),
    ...
  )
  # Catches a `...` argument RunUMAP() never touched: a typo, or a
  # parameter Seurat renamed. Without this check, it would be silently
  # dropped instead of erroring.
  check_dots_used()

  group_by <- "orig.ident"
  # a features-driven run (var_features/features set instead of
  # reduction+num_dim) has no assay/reduction/dims to report - that trio is
  # only meaningful on the reduction path
  dr_info <- if (is.null(reduction)) {
    ""
  } else {
    glue("-{assay}-{reduction}-{num_dim}")
  }
  plot_dr_group(
    x,
    reduction = dr_name,
    group_by = group_by,
    color_scheme = get_color_scheme("samples"),
    file_prefix = glue("dr-{dr_name}{dr_info}-{group_by}"),
    file_format = file_format
  )

  x
}
