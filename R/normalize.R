#' Normalize the counts, select variable features, and scale.
#'
#' A thin wrapper that picks the right Seurat function for `method` and fixes
#' the parameters this pipeline uses, so either path hands back a fully
#' processed object — no separate variance step needed. `method = "sct"`
#' already worked this way natively: `SCTransform()` produces the counts,
#' variable features and scaled data in one pass. `method = "log"` now mirrors
#' that explicitly: `NormalizeData()` followed by `FindVariableFeatures()` and
#' `ScaleData()`.
#'
#' No covariates are regressed out (`vars.to.regress = NULL`) for either path.
#' Germain et al. (2020) found that the "common practice of regressing out cell
#' covariates, such as the detection rate or proportion of mitochondrial
#' reads\[,\] nearly always had a negative impact."
#'
#' On the log path, this function drops any `vf_*` columns already on
#' `assay`'s feature meta data before recomputing them. An earlier call can
#' leave these columns behind: each sample's own pre-merge run, or a merge
#' step's preview run. Left alone, they would just accumulate call after
#' call. See [variable_features_by_batch()] for what reading an accumulated
#' column list, instead of a current one, leads to.
#'
#' `SCTransform()` never writes this kind of per-gene bookkeeping column. So
#' the sct path has nothing to clean up. The check runs the same way for
#' both paths; it is a no-op for sct, not something conditioned on `method`.
#'
#' @param x A Seurat object.
#' @param method "log" for `NormalizeData()` + `FindVariableFeatures()` +
#'   `ScaleData()` (log-normalized to a scale factor of 10,000), or "sct" for
#'   `SCTransform()`.
#' @param num_variable_genes Number of variable features.
#' @param assay Assay to normalize.
#' @param log_file Log file.
#'
#' @return The processed Seurat object. Also writes the variable-feature plot
#'   (`variance-features.png`) to the working directory.
#'
#' @export
normalize_counts <- function(
  x,
  method = "log",
  num_variable_genes = 3000,
  assay = "RNA",
  log_file = NULL
) {
  method <- match.arg(method, c("log", "sct"))

  # Benchmarks show multisession buys nothing here (~tied with sequential,
  # 1.7k-35k cells). It also crashes outright on a large merged object
  # (future.globals.maxSize). Always run sequential, before the switch()
  # below: NormalizeData() can hit this crash too, not just ScaleData().
  future::plan(future::sequential)

  # drop variable-feature columns left by an earlier normalize_counts() call on
  # this assay - see the roxygen note above
  old_vf_cols <- grep("^vf_", colnames(x[[assay]][[]]), value = TRUE)
  if (length(old_vf_cols)) {
    x[[assay]][[old_vf_cols]] <- NULL
  }

  x <- switch(
    method,
    log = NormalizeData(
      x,
      assay = assay,
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    ),
    sct = SCTransform(
      x,
      assay = assay,
      variable.features.n = num_variable_genes,
      verbose = FALSE
    )
  )

  # SCTransform() already selected variable features and scaled; "log" needs
  # the extra steps to match
  if (method == "log") {
    message("\n\n ===== scale variable features ===== \n\n")

    x <- FindVariableFeatures(
      x,
      assay = assay,
      selection.method = "vst",
      nfeatures = num_variable_genes,
      verbose = FALSE
    )

    # features = VariableFeatures(): scaling every gene (the default) is much
    # slower, and nothing downstream reads scale.data for non-variable genes.
    # RunPCA() only uses VariableFeatures(). vars.to.regress = NULL: see the
    # roxygen note above (Germain et al., 2020) for why.
    x <- ScaleData(
      x,
      features = VariableFeatures(x, assay = assay),
      assay = assay,
      vars.to.regress = NULL,
      verbose = FALSE
    )
  }

  # both methods have selected variable features by this point
  hvf_tbl <- HVFInfo(x) |>
    round(3) |>
    rownames_to_column("gene")

  # HVFInfo() names its ranking column differently per method:
  # "variance.standardized" for FindVariableFeatures(selection.method = "vst")
  # (the log path), "residual_variance" for SCTransform() (the sct path) - sort
  # by whichever is actually present
  rank_col <- intersect(
    c("variance.standardized", "residual_variance"),
    colnames(hvf_tbl)
  )[1]
  hvf_tbl <- hvf_tbl |> arrange(-.data[[rank_col]])

  var_plot <- VariableFeaturePlot(x, pt.size = 0.5)
  var_plot <- LabelPoints(
    var_plot,
    points = head(hvf_tbl$gene, 30),
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
  )
  # VariableFeaturePlot()'s x-axis is log10(mean expression). ggplot2
  # objects are lazy: a gene with a legitimate zero mean expression (-Inf on
  # that scale) only triggers the scale_x_log10() warning once the plot
  # renders here, not when the plot object above is built. Suppressing the
  # warning at construction, inside VariableFeaturePlot() itself, would be a
  # no-op.
  suppressWarnings(ggsave(
    "variance-features.png",
    plot = var_plot,
    width = 12,
    height = 5,
    units = "in"
  ))

  x
}
