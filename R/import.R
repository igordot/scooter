#' Read in Gene Expression and Antibody Capture data from a 10x Genomics Cell
#' Ranger sparse matrix or from a text file.
#'
#' @param path Path to directory containing 10x matrix, or path to a text file.
#' @param sample_name A character that will be used as a prefix for all cell
#'   names.
#' @param log_file Filename for the log file.
#'
#' @return Named list of matrices. One matrix for each data type. Currently the
#'   only two data types are 'Gene Expression' and 'Antibody Capture'
#'
#' @importFrom stringr str_subset str_c
#' @importFrom data.table fread
#' @export
read_counts_file <- function(path, sample_name, log_file = NULL) {
  if (file.exists(path) && !dir.exists(path)) {
    counts_df <- fread(path, stringsAsFactors = FALSE, data.table = FALSE)
    # data.table doesn't read in rownames assuming first column is the
    # gene/feature names
    counts_df <- counts_df |> column_to_rownames(colnames(counts_df)[1])

    # sometimes there is an 'unmapped' row
    unmapped <- str_detect(toupper(rownames(counts_df)), "UNMAPPED")
    if (any(unmapped)) {
      counts_df <- counts_df[!unmapped, ]
    }

    # if the cell names all have -1 at the end, remove the -1
    cell_names <- colnames(counts_df)
    if (all(str_detect(string = cell_names, pattern = "\\-1$"))) {
      cell_names <- as.vector(as.character(sapply(
        X = cell_names,
        function(x) stringr::str_sub(x, end = (nchar(x) - 2))
      )))
    }
    colnames(counts_df) <- cell_names

    # determine data type based on number of features
    if (nrow(counts_df) < 1000) {
      counts_matrix <- list("Antibody Capture" = counts_df)
    } else {
      counts_matrix <- list("Gene Expression" = counts_df)
    }
  } else {
    if (!dir.exists(path)) {
      stop(glue("{path} is not a file and is not a directory"))
    }

    # directories should contain matrix.mtx (or matrix.mtx.gz) files
    mtx_paths <- list.files(
      path = path,
      pattern = "matrix.mtx",
      full.names = TRUE,
      recursive = TRUE
    )

    # matrix.mtx file should be in filtered_feature (cellranger count) or
    # sample_feature (cellranger multi) directory
    mtx_paths <- str_subset(
      mtx_paths,
      "filtered_gene_bc_mat|filtered_feature_bc_mat|sample_feature_bc_mat"
    )

    if (length(mtx_paths) > 0) {
      data_dir <- dirname(mtx_paths[1])

      message_str <- glue("counts dir: {data_dir}")
      write_message(message_str, log_file)

      counts_matrix <- import_mtx(data_dir)
    } else {
      # no filtered_feature_bc_matrix-style mtx dir - fall back to an h5 counts
      # matrix
      h5_paths <- list.files(
        path = path,
        pattern = "^filtered_feature_bc_matrix\\.h5$",
        full.names = TRUE,
        recursive = TRUE
      )
      if (!length(h5_paths)) {
        h5_paths <- list.files(
          path = path,
          pattern = "\\.h5$",
          full.names = TRUE,
          recursive = TRUE
        )
      }
      if (length(h5_paths) > 1) {
        stop(glue(
          "multiple h5 files found in {path}, none named ",
          "filtered_feature_bc_matrix.h5:\n",
          "{str_c(h5_paths, collapse = '\n')}"
        ))
      }
      if (!length(h5_paths)) {
        stop(glue("{path} does not contain matrix.mtx or an h5 file"))
      }

      message_str <- glue("counts matrix h5: {h5_paths}")
      write_message(message_str, log_file)

      counts_matrix <- Read10X_h5(h5_paths)
      if (!is.list(counts_matrix)) {
        counts_matrix <- list("Gene Expression" = counts_matrix)
      }
    }
  }

  # rename the column(cell) names by adding a sample name prefix to every
  # barcode
  counts_out <- list()
  for (i in 1:length(counts_matrix)) {
    current_mat <- counts_matrix[[i]]
    colnames(current_mat) <- str_c(sample_name, ":", colnames(current_mat))
    counts_out[[i]] <- current_mat
  }
  names(counts_out) <- names(counts_matrix)

  return(counts_out)
}

#' Create a new Seurat object from a matrix.
#'
#' For an RNA assay the per-cell QC metrics are added at the same time:
#' `pct_mito`/`pct_ribo`/ `pct_hb` (via [add_gene_class_percent()]), plus
#' `detected_genes` and `total_counts` as clearer aliases of
#' `nFeature_RNA`/`nCount_RNA`, and `sample_name` as a clearer alias of
#' `orig.ident`. All six are metadata columns as soon as the object exists,
#' so [plot_metrics_distribution()] (which takes its `metrics` explicitly, no
#' default) and [plot_metrics_correlations()] can use any of them without
#' further preparation.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param assay Seurat assay to add the data to.
#' @param min_cells Include genes/features detected in at least this many cells.
#' @param min_genes Include cells where at least this many genes/features are
#'   detected.
#' @param project Project name for Seurat object.
#' @param log_file Filename for the logfile.
#'
#' @return Seurat object.
#'
#' @importFrom Matrix rowSums
#' @export
initialize_seurat_object <- function(
  counts_matrix,
  assay = "RNA",
  min_cells = 1,
  min_genes = 1,
  project = "proj",
  log_file = NULL
) {
  # check that the size of the input matrix is reasonable
  if (ncol(counts_matrix) < 10) {
    stop(glue("matrix contains too few cells: {ncol(counts_matrix)}"))
  }

  # CreateSeuratObject() only accepts dgCMatrix and warns while coercing
  # anything else, and the import functions supply other classes (readMM()
  # returns dgTMatrix, flat files a base matrix)
  counts_matrix <- as.sparse(counts_matrix)

  # remove genes with very few counts
  counts_matrix <- counts_matrix[Matrix::rowSums(counts_matrix) > 0, ]

  write_message(glue("input cells: {ncol(counts_matrix)}"), log_file)
  write_message(glue("input genes: {nrow(counts_matrix)}"), log_file)

  # Create seurat object
  so <- CreateSeuratObject(
    counts = counts_matrix,
    project = "proj",
    assay = assay,
    names.field = 1,
    names.delim = ":",
    min.cells = min_cells,
    min.features = min_genes
  )

  # add the per-cell QC metrics. `detected_genes`/`total_counts`/`sample_name`
  # are clearer names for Seurat's `nFeature_RNA`/`nCount_RNA`/`orig.ident`,
  # and are what plot_metrics_correlations() looks for by default (and what
  # callers of plot_metrics_distribution() pass as `metrics`), so the object
  # comes out ready to plot rather than needing the rename at every call
  # site. The originals are left in place: filter_cells() still reads those.
  if (assay == "RNA") {
    so$sample_name <- so$orig.ident
    so$detected_genes <- so$nFeature_RNA
    so$total_counts <- so$nCount_RNA
    so <- add_gene_class_percent(so)
  }

  return(so)
}

#' Add assay to Seurat object.
#'
#' @param x Seurat object.
#' @param assay Seurat assay to add the matrix to.
#' @param counts_matrix Raw counts matrix.
#' @param log_file Filename for the log file.
#'
#' @return Seurat object of cells found in both the existing object and new
#'   data.
#'
#' @importFrom methods is
#' @export
add_seurat_assay <- function(x, assay, counts_matrix, log_file = NULL) {
  if (!is(x, "Seurat")) {
    stop(glue("{x} is not a Seurat object. Cannot add Assay"))
  }

  if (assay %in% names(x)) {
    stop(glue("{assay} already exists in the Seurat object"))
  }

  # use cells that are found in both antibody capture and RNA
  cells_to_use <- intersect(colnames(x), colnames(counts_matrix))

  if (length(x) != length(cells_to_use)) {
    message_str <- glue(
      "{ncol(x) - length(cells_to_use)} cells in seurat object are not in ",
      "counts matrix"
    )
    write_message(message_str, log_file)
  }

  if (ncol(counts_matrix) != length(cells_to_use)) {
    message_str <- glue(
      "{ncol(x) - ncol(counts_matrix)} cells in counts matrix not in scrna ",
      "matrix"
    )
    write_message(message_str, log_file)
  }

  # subset counts by joint cell barcodes keep it sparse: CreateAssay5Object()
  # warns while coercing anything that is not a dgCMatrix
  counts_matrix <- as.sparse(counts_matrix[, cells_to_use, drop = FALSE])
  x <- subset(x, cells = cells_to_use)

  # add assay (v5 assay, to match the object created by
  # initialize_seurat_object())
  x[[assay]] <- SeuratObject::CreateAssay5Object(counts = counts_matrix)

  return(x)
}

#' Read in 10x Genomics Cell Ranger Matrix Market format data.
#'
#' @param data_path Path to directory that holds the files output from 10x.
#' @param gene_column The column with the gene names.
#' @param log_file Filename for the log file.
#'
#' @return Named list of matrices. One matrix for each data type as specified
#'   in the third column of the features.tsv file. As of Oct 3rd 2019, the two
#'   options are `Gene Expression` and `Antibody Capture`
#'
#' @import utils
#' @importFrom Matrix readMM t
#' @importFrom stringr str_detect str_sub
#' @export
import_mtx <- function(data_path, gene_column = 2, log_file = NULL) {
  # Heavily sourced from Seurat caveat the name of the files have to be
  # features.tsv.gz, matrix.mtx.gz, barcodes.tsv.gz

  # check if the directory exists
  if (!dir.exists(paths = data_path)) {
    stop(glue("dir {data_path} does not exist"))
  }

  feature_names <- read.delim(
    file = file.path(data_path, "features.tsv.gz"),
    header = FALSE,
    stringsAsFactors = FALSE
  )

  data <- readMM(file = file.path(data_path, "matrix.mtx.gz"))

  # set the rownames to unique gene names
  rownames(data) <- make.unique(feature_names[, gene_column])

  # read in cell barcodes
  cell_names <- readLines(file.path(data_path, "barcodes.tsv.gz"))

  # if the cell names all have -1 at the end, remove the -1
  if (all(str_detect(string = cell_names, pattern = "\\-1$"))) {
    cell_names <- as.vector(as.character(sapply(
      X = cell_names,
      function(x) stringr::str_sub(x, end = (nchar(x) - 2))
    )))
  }

  # set the barcodes as colnames of the gene expression data
  colnames(data) <- cell_names

  # a third column indicates that there may be more than one datatype
  if (ncol(feature_names) > 2) {
    # get the different datatypes
    data_types <- factor(feature_names$V3)
    data_levels <- levels(data_types)

    if (length(data_levels) > 1) {
      for (data_level in data_levels) {
        write_message(
          glue("{data_level}: {sum(data_types == data_level)} features"),
          log_file
        )
      }
    }

    # split the data into a list with a matrix for each datatype
    data <- lapply(
      X = data_levels,
      FUN = function(l) {
        return(data[data_types == l, ])
      }
    )

    # name the elements in the list by datatype
    names(data) <- data_levels
  } else {
    data <- list("Gene Expression" = data)
  }

  return(data)
}

#' Add mitochondrial, ribosomal, and hemoglobin percentages to a Seurat object.
#'
#' Detected by gene symbol prefix, case-insensitively so the same patterns
#' match both human (`MT-`, `RPS`/`RPL`, `HBA`/`HBB`) and mouse (`mt-`,
#' `Rps`/`Rpl`, `Hba-a1`/`Hbb-bs`) naming. `RPS6K*` (ribosomal protein S6
#' kinase) is excluded from `pct_ribo` - it matches the `RPS` prefix but is a
#' signaling kinase, not a ribosomal protein.
#'
#' High expression levels of mitochondrial genes could be an indicator of poor
#' sample quality, leading to a high fraction of apoptotic or lysing cells (see
#' reference below).
#'
#' `pct_ribo` is normally 15-45%; the exact fraction varies by cell type and
#' overall cell health. "Although both ribosomal proteins and ribosomal RNA
#' (rRNA) make up the ribosome complex, ribosomal protein transcripts are not
#' equivalent to ribosomal RNA (rRNA). Ribosomal protein transcript detection
#' will not necessarily correlate with either rRNA or mitochondrial
#' transcripts." (see reference below)
#'
#' High `pct_hb` is a contamination signal in PBMC and most solid tissue:
#' per-sample `pct_hb` flags incomplete red blood cell lysis in PBMC preps, or
#' inadequate perfusion in solid tissue.
#'
#' @param x A Seurat object.
#'
#' @return Seurat object with `pct_mito`, `pct_ribo`, and `pct_hb` metadata
#'   columns added.
#'
#' @references
#' `pct_mito`:
#' <https://kb.10xgenomics.com/s/article/360001086611-Why-do-I-see-a-high-level-of-mitochondrial-gene-expression>
#'
#' `pct_ribo`:
#' <https://kb.10xgenomics.com/s/article/218169723-What-fraction-of-reads-map-to-ribosomal-proteins>
#'
#' @importFrom Matrix colSums
#' @export
add_gene_class_percent <- function(x) {
  counts <- GetAssayData(x, assay = "RNA", layer = "counts")
  gene_names <- rownames(counts)
  total_counts <- Matrix::colSums(counts)

  is_mito <- grepl("^mt-", gene_names, ignore.case = TRUE)
  is_ribo <- grepl("^rp[sl][0-9]", gene_names, ignore.case = TRUE) &
    !grepl("^rps6k", gene_names, ignore.case = TRUE)
  is_hb <- grepl(
    "^hb[abdegmqz][0-9]*$|^hb[ab]-",
    gene_names,
    ignore.case = TRUE
  )

  pct <- function(is_gene_set) {
    round(
      100 * Matrix::colSums(counts[is_gene_set, , drop = FALSE]) / total_counts,
      3
    )
  }
  x$pct_mito <- pct(is_mito)
  x$pct_ribo <- pct(is_ribo)
  x$pct_hb <- pct(is_hb)

  x
}
