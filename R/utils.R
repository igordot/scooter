#' Set identity of the Seurat object.
#'
#' Wrapper for SeuratObject::Idents() with extra safety checks.
#'
#' @param x A Seurat object.
#' @param identity_column The name of the identity column to pull from object
#'   metadata.
#
#' @return A Seurat object with an updated identity set.
#'
#' @importFrom SeuratObject Idents
#'
#' @export
set_identity <- function(x, identity_column) {
  so <- x

  # check if the specified identity column is one of meta data columns
  identity_column <- check_identity_column(
    x = x,
    identity_column = identity_column
  )

  # set identities based on specified column
  message("setting identity column: ", identity_column)
  SeuratObject::Idents(so) <- identity_column

  return(so)
}

#' Save the counts matrix as a single table.
#'
#' Layers left split by a `merge()` are rejoined first, since `LayerData()`
#' would otherwise return only the first of them, and warn rather than error
#' while doing it.
#'
#' @param x A Seurat object.
#' @param assay Assay to retrieve the counts from.
#' @param layer Layer to retrieve. The default "data" layer holds the normalized
#'   and log-transformed expression used for visualization and most
#'   differential expression tests.
#' @param digits Number of digits to round the values to. Ignored for a layer
#'   that already holds whole numbers, such as raw counts.
#' @param file Path to write the table to as a csv, or `NULL` to skip writing
#'   and just return the tibble. Has no default: every caller states where the
#'   table goes.
#'
#' @return A tibble with one row per gene.
#'
#' @importFrom data.table fwrite
#' @importFrom methods is
#' @export
save_counts <- function(
  x,
  assay = "RNA",
  layer = "data",
  digits = 3,
  file
) {
  # `file` is only read at the end, so force it here: a missing argument should
  # be reported before densifying the matrix, not after
  force(file)

  # The counts are stored sparse, but as.matrix() densifies them. Long
  # vectors are not supported everywhere: 32-bit builds, and the many .Call
  # interfaces that still index with int. So past 2^31 elements, the dense
  # matrix is unusable, no matter how much memory the machine has.
  counts_dim <- dim(x[[assay]])
  # prod() returns a double, so this does not overflow the very limit it is
  # testing
  if (prod(counts_dim) > .Machine$integer.max) {
    stop(glue(
      "{assay} {layer} matrix is too large to write: {counts_dim[1]} genes x ",
      "{counts_dim[2]} cells is over the 2^31 element limit that many R and ",
      "Matrix operations still assume"
    ))
  }

  # merge() leaves the layers split ("data.1", "data.2", ...). LayerData()
  # then returns just the first one, with a warning, and quietly drops every
  # other sample. Rejoin first. Joining an already-joined assay is a no-op,
  # but it still copies the matrix: only the split case pays that cost.
  if (length(SeuratObject::Layers(x[[assay]], search = layer)) > 1) {
    x <- SeuratObject::JoinLayers(x, assay = assay)
  }

  counts_mat <- SeuratObject::LayerData(x, assay = assay, layer = layer)

  # A counts layer holds whole numbers. Rounding it is a pass over the
  # matrix that changes nothing. Test the values, not the storage type:
  # integer-valued counts are still held as doubles. On a sparse matrix that
  # test only touches the nonzeros; the implicit zeros are whole numbers
  # either way. Rounding here, rather than after as.matrix(), also skips the
  # zeros.
  counts_values <- if (is(counts_mat, "sparseMatrix")) {
    counts_mat@x
  } else {
    counts_mat
  }
  if (any(counts_values != trunc(counts_values), na.rm = TRUE)) {
    counts_mat <- round(counts_mat, digits)
  }

  counts_tbl <- counts_mat |>
    as.matrix() |>
    as_tibble(rownames = "gene") |>
    arrange(.data$gene)

  if (!is.null(file)) {
    # fwrite() rather than write_csv() because these are the largest tables the
    # pipeline emits; it gzips on the ".gz" extension by itself. Thread count
    # comes from data.table::setDTthreads().
    fwrite(counts_tbl, file)
  }

  counts_tbl
}

#' Function to merge two metadata tables together.
#'
#' `x` and `y` must not share column names other than `cell` - `stop()`s if
#' they do, rather than letting `left_join()` silently rename the collision to
#' `<col>.x`/`<col>.y`.
#'
#' @param x A Seurat object or tibble with cell IDs (a `cell` column, or as
#'   rownames).
#' @param y A tibble with cell IDs (a `cell` column, or as rownames).
#' @param log_file A log filename.
#'
#' @return A metadata file merged on cell identifiers.
#'
#' @export
merge_metadata <- function(x, y, log_file = NULL) {
  UseMethod("merge_metadata")
}

#' @export
merge_metadata.default <- function(x, y, log_file = NULL) {
  # add a "cell" column from rownames, unless one already exists
  ensure_cell_column <- function(df) {
    df <- as.data.frame(df)
    if ("cell" %in% colnames(df)) df else rownames_to_column(df, "cell")
  }
  x <- ensure_cell_column(x)
  y <- ensure_cell_column(y)

  # Besides "cell", left_join() silently renames a column name collision to
  # "<col>.x"/"<col>.y" instead of erroring. That collision then surfaces as
  # a confusing missing-column error downstream. stop() here instead, while
  # the cause is still obvious.
  shared_cols <- intersect(
    setdiff(colnames(x), "cell"),
    setdiff(colnames(y), "cell")
  )
  if (length(shared_cols)) {
    stop(glue(
      "x and y share column names besides 'cell': ",
      "{toString(shared_cols)}"
    ))
  }

  # warn rather than error, since a partial overlap (such as filtered vs.
  # unfiltered metadata for the same sample) is a normal caller mistake worth
  # flagging, not a reason to stop
  only_x <- setdiff(x$cell, y$cell)
  only_y <- setdiff(y$cell, x$cell)
  if (length(only_x) || length(only_y)) {
    warning(glue(
      "cells do not fully overlap: {length(only_x)} only in x, ",
      "{length(only_y)} only in y, {length(intersect(x$cell, y$cell))} in both"
    ))
  }

  # compile all cell metadata into a single table
  left_join(x, y, by = "cell")
}

#' @export
merge_metadata.Seurat <- function(x, y, log_file = NULL) {
  x <- x@meta.data

  s_metadata <- merge_metadata(x, y, log_file = log_file) |>
    column_to_rownames("cell")

  return(s_metadata)
}

#' Function to extract data from Seurat object.
#'
#' @param x A Seurat object.
#' @param assay Assay such as RNA.
#' @param slot Slot such as counts. Default is scale.data.
#' @param features Features from assay.
#' @param reduction Character vector of reduction types.
#' @param metadata Boolean. To grab metadata or not
#'
#' @return A metadata file merged on cell identifiers.
#'
#' @importFrom purrr reduce
#' @export
as_data_frame_seurat <- function(
  x,
  assay = NULL,
  slot = NULL,
  features = NULL,
  reduction = NULL,
  metadata = TRUE
) {
  # if metadata, extract metadata
  if (metadata == TRUE) {
    metadata_out <- x@meta.data |>
      rownames_to_column("cell")
  } else {
    metadata_out <- NULL
  }

  if (!is.null(assay)) {
    assay_out <- as_data_frame_seurat_assay(
      x = x,
      assay = assay,
      slot = slot,
      features = features
    )
  } else {
    assay_out <- NULL
  }

  if (!is.null(reduction)) {
    reduction_out <- as_data_frame_seurat_reduction(
      x = x,
      reduction = reduction
    )
  } else {
    reduction_out <- NULL
  }

  data <- list(metadata_out, assay_out, reduction_out)
  idx <- which(sapply(data, function(x) !is.null(x)))

  data <- data[idx]

  # merge the extracted data
  if (length(data) == 1) {
    data_out <- as.data.frame(data)
  } else {
    data_out <- reduce(data, .f = full_join, by = "cell")
  }

  return(data_out)
}

#' Save the cell metadata and the reduction embeddings as a single table.
#'
#' Joins the requested reductions onto the cell metadata on a `cell` column.
#' Reductions that are not present in the object are skipped, so this works
#' before and after the dimensionality reduction steps have run.
#'
#' @param x A Seurat object.
#' @param reduction Reductions to join onto the metadata.
#' @param digits Number of digits to round the embeddings to.
#' @param file Path to write the table to as a csv. Set to `NULL` to skip
#'   writing and just return the tibble.
#'
#' @return A tibble with one row per cell.
#'
#' @export
save_metadata <- function(
  x,
  reduction = c("tsne", "umap"),
  digits = 3,
  file = "metadata.csv.gz"
) {
  reduction <- intersect(reduction, SeuratObject::Reductions(x))

  metadata_tbl <- x@meta.data |>
    as_tibble(rownames = "cell")

  # rounded here rather than through as_data_frame_seurat(), which returns the
  # embeddings as-is
  for (reduction_name in reduction) {
    reduction_tbl <- SeuratObject::Embeddings(x, reduction = reduction_name) |>
      round(digits) |>
      as.data.frame() |>
      rownames_to_column("cell")
    metadata_tbl <- full_join(metadata_tbl, reduction_tbl, by = "cell")
  }

  metadata_tbl <- arrange(metadata_tbl, .data$cell)

  if (!is.null(file)) {
    write_csv(metadata_tbl, file)
  }

  metadata_tbl
}

#' Report the size of an object's slots, largest first.
#'
#' Works on any S4 object or list, not just Seurat objects. It descends into
#' S4 slots and list elements, bounded by `max_depth`. So a Seurat object,
#' an assay pulled from one (`x[["RNA"]]`), a `DimReduc`, or a plain nested
#' list all get the same recursive breakdown.
#'
#' An `Assay5`'s `layers` slot is itself a list. So
#' `layers > counts`/`layers > data`/`layers > scale.data` are sized
#' separately, not lumped into one `layers` total.
#'
#' A matrix, data frame, or atomic vector is neither S4 nor a list. Passing
#' one of these just reports its own size, the same as calling
#' `object.size()` directly. There is nothing to descend into.
#'
#' @param x Any object - a Seurat object, a piece pulled from one (an assay via
#'   `x[["RNA"]]`, a reduction via `x[["pca"]]`, ...), a plain list, or
#'   anything else `object.size()` accepts.
#' @param max_depth How many levels to descend before reporting a subtree's
#'   total size instead of continuing to break it down.
#'
#' @return A tibble with one row per slot/element, largest first: `path` (its
#'   position, `" > "`-separated) and `size` (human-readable, e.g. "1.2 Mb").
#'
#' @export
profile_object_size <- function(x, max_depth = 3) {
  walk <- function(value, path, depth) {
    size <- as.numeric(object.size(value))
    if (
      depth >= max_depth ||
        !(isS4(value) || (is.list(value) && !is.data.frame(value)))
    ) {
      return(setNames(size, path))
    }
    children <- if (isS4(value)) slotNames(value) else names(value)
    if (!length(children)) {
      return(setNames(size, path))
    }
    child_get <- if (isS4(value)) {
      function(nm) slot(value, nm)
    } else {
      function(nm) value[[nm]]
    }
    unlist(lapply(children, function(nm) {
      walk(child_get(nm), paste(c(path, nm), collapse = " > "), depth + 1)
    }))
  }

  sizes <- sort(walk(x, class(x)[1], 0), decreasing = TRUE)
  tibble(
    path = names(sizes),
    size = vapply(
      sizes,
      function(b) format(structure(b, class = "object_size"), units = "auto"),
      character(1)
    )
  )
}

#' Check identity of the Seurat object.
#'
#' @param x A Seurat object.
#' @param identity_column The name of the identity column to pull from object
#'   metadata.
#
#' @return The name of the identity column, potentially corrected if
#'   resolution.
#'
#' @export
check_identity_column <- function(x, identity_column) {
  # check if the grouping variable is one of meta data columns
  if (!(identity_column %in% colnames(x@meta.data))) {
    # check if grouping variable is the resolution value (X.X instead of
    # res.X.X)
    res_column <- stringr::str_c("res.", identity_column)
    if (res_column %in% colnames(x@meta.data)) {
      identity_column <- res_column
    } else {
      stop("unknown grouping variable: ", identity_column)
    }
  }

  return(identity_column)
}

#' Find the single .rds/.qs2 file a path refers to.
#'
#' A directory resolves to the one .rds/.qs2 file inside it. A file path
#' passes through unchanged. This logic is split out of
#' [resolve_seurat_object()] so callers can find the actual file being
#' read, not just the directory a caller happened to pass, without
#' duplicating the directory-search logic.
#'
#' @param x A path to a `.rds`/`.qs2` file, or a directory containing exactly
#'   one such file.
#'
#' @return The path to the `.rds`/`.qs2` file.
#'
#' @noRd
find_seurat_file <- function(x) {
  if (!file.exists(x)) {
    stop(glue("path not found: {x}"))
  }

  if (!dir.exists(x)) {
    return(x)
  }

  candidates <- list.files(
    x,
    pattern = "\\.(qs2|rds)$",
    ignore.case = TRUE,
    full.names = TRUE
  )
  if (length(candidates) == 0) {
    stop(glue("no .rds or .qs2 file found in directory: {x}"))
  }
  if (length(candidates) > 1) {
    stop(glue("multiple .rds/.qs2 files found in directory, expected one: {x}"))
  }
  candidates
}

#' Resolve a Seurat object, given either the object itself or a path to one on
#' disk.
#'
#' Meant for parameters that accept either an in-memory Seurat object or a
#' serialized one - resolving here means callers do not need to check which
#' they got themselves.
#'
#' @param x A Seurat object, a path to one saved as `.rds` or `.qs2`, or a
#'   directory containing exactly one such file.
#'
#' @return A Seurat object.
#'
#' @importFrom methods is
#' @export
resolve_seurat_object <- function(x) {
  if (is(x, "Seurat")) {
    return(x)
  }

  if (!is.character(x) || length(x) != 1) {
    stop(glue(
      "x must be a Seurat object or a single file path, not {class(x)[1]}"
    ))
  }

  x <- find_seurat_file(x)

  obj <- if (grepl("\\.qs2$", x, ignore.case = TRUE)) {
    qs2::qs_read(x, validate_checksum = TRUE, nthreads = 4)
  } else if (grepl("\\.rds$", x, ignore.case = TRUE)) {
    readRDS(x)
  } else {
    stop(glue("unsupported file extension (expected .rds or .qs2): {x}"))
  }

  if (!is(obj, "Seurat")) {
    stop(glue("{x} did not contain a Seurat object"))
  }

  obj
}

#' Small function to write to message and to log file.
#'
#' @param message_str A string to write as a message.
#' @param log_file A log filename.
#
#' @return A message and writes the message to the specified log file.
#'
#' @examples
#' write_message(message_str = "Finished Step 1", log_file = "log.file.txt")
#' @export
write_message <- function(message_str, log_file = NULL) {
  # Small function to write to message and to log file if log file is not null
  message(message_str)
  if (!is.null(log_file)) {
    write(message_str, file = log_file, append = TRUE)
  }
}

as_data_frame_seurat_assay <- function(
  x,
  assay = NULL,
  slot = NULL,
  features = NULL
) {
  # If assay is specified
  if (!is.null(assay)) {
    # and slot is specified
    if (!is.null(slot)) {
      # get data from that assay and layer
      s_obj_assay <- GetAssayData(x, assay = assay, layer = slot)
      if (!is.null(features)) {
        s_obj_assay <- as.data.frame(s_obj_assay[features, , drop = FALSE])
        rownames(s_obj_assay) <- features
      }
      s_obj_assay_out <- as.data.frame(t(s_obj_assay)) |>
        rownames_to_column("cell")
    } else {
      s_obj_assay <- GetAssayData(x, assay = assay, layer = "data")
      if (!is.null(features)) {
        s_obj_assay <- as.data.frame(s_obj_assay[features, , drop = FALSE])
        colnames(s_obj_assay) <- features
      }
      s_obj_assay_out <- as.data.frame(s_obj_assay) |>
        rownames_to_column("cell")
    }
  } else {
    stop("No assay specified")
  }

  return(s_obj_assay_out)
}

as_data_frame_seurat_reduction <- function(x, reduction) {
  # Extract the cell embeddings for each requested reduction.
  reductions_to_save <- lapply(reduction, function(reduction_name) {
    as.data.frame(x@reductions[[reduction_name]]@cell.embeddings)
  })

  # set rownames to columns for easier joining
  reductions_to_save <- lapply(
    reductions_to_save,
    function(x) {
      rownames_to_column(
        .data = x,
        "cell"
      )
    }
  )
  # join lists by the new column
  s_obj_reduction <- reduce(reductions_to_save, .f = left_join, by = "cell")

  return(s_obj_reduction)
}

#' Function to create a color vector.
#'
#' @param seurat_obj A Seurat object.
#' @param group Assay such as RNA.
#'
#' @return A vector of colors.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggsci pal_igv
#' @export
create_color_vect <- function(seurat_obj, group = "orig.ident") {
  # create a vector of colors for the Idents of seurat_obj
  sample_names <- switch(
    class(seurat_obj),
    Seurat = seurat_obj[[group]] |> unique() |> arrange(get(group)),
    data.frame = unique(seurat_obj)
  )

  colors_samples <- c(
    brewer.pal(5, "Set1"),
    brewer.pal(8, "Dark2"),
    pal_igv("default")(51)
  )
  # create a named color scheme to ensure names and colors are in the proper
  # order
  sample_names[] <- lapply(sample_names, as.character)
  colors_samples_named <- colors_samples[1:nrow(sample_names)]
  names(colors_samples_named) <- sample_names[, 1]
  return(colors_samples_named)
}
