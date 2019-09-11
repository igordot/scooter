#' Create a new Seurat object from a matrix.
#'
#' @param counts_matrix A matrix of raw counts.
#' @param min_cells X.
#' @param min_genes X.
#' @param out_dir Directory where all the output files will be saved.
#' @param color_scheme A character vector of colors for plots (per sample/library).
#'
#' @return Seurat object.
#'
#' @examples
#' create_seurat_obj(counts_matrix = counts_matrix)
#'
#' @import Matrix dplyr readr tibble
#' @importFrom glue glue
#' @importFrom Seurat CreateSeuratObject AddMetaData GenePlot
#' @export
create_seurat_obj <- function(counts_matrix, min_cells = 10, min_genes = 100, out_dir = ".", color_scheme = NULL) {

  # create output directories
  if (!dir.exists(out_dir)) dir.create(out_dir)
  qc_dir = glue("{out_dir}/qc")
  if (!dir.exists(qc_dir)) dir.create(qc_dir)

  # color scheme
  if (is.null(color_scheme)) color_scheme <- get_color_scheme(type = "samples")

  # check that the size of the input matrix is reasonable
  if (ncol(counts_matrix) < 10) stop(glue("matrix contains too few cells: {ncol(counts_matrix)}"))
  if (nrow(counts_matrix) < 100) stop(glue("matrix contains too few genes: {nrow(counts_matrix)}"))

  # remove genes with very few counts
  counts_matrix <- counts_matrix[Matrix::rowSums(counts_matrix) > 0, ]

  message("input cells: ", ncol(counts_matrix))
  message("input genes: ", nrow(counts_matrix))

  # log to file
  write(glue("input cells: {ncol(counts_matrix)}"), file = "create.log", append = TRUE)
  write(glue("input genes: {nrow(counts_matrix)}"), file = "create.log", append = TRUE)

  # save counts matrix as a csv file (to be consistent with the rest of the tables)
  raw_data <- counts_matrix %>% as.matrix() %>% as.data.frame() %>% rownames_to_column("gene") %>% arrange(gene)
  write_csv(raw_data, path = "counts.raw.csv.gz")

  s_obj <- CreateSeuratObject(
    raw.data = counts_matrix,
    min.cells = min_cells, min.genes = min_genes, project = "proj",
    names.field = 1, names.delim = ":"
  )

  message(glue("usable cells: {ncol(s_obj@data)}"))
  message(glue("usable genes: {nrow(s_obj@data)}"))

  # log to file
  write(glue("usable cells: {ncol(s_obj@data)}"), file = "create.log", append = TRUE)
  write(glue("usable genes: {nrow(s_obj@data)}"), file = "create.log", append = TRUE)

  # nGene and nUMI are automatically calculated for every object by Seurat
  # calculate the percentage of mitochondrial genes here and store it in percent.mito using the AddMetaData
  mt_genes <- grep("^MT-", rownames(s_obj@data), ignore.case = TRUE, value = TRUE)
  percent_mt <- Matrix::colSums(s_obj@raw.data[mt_genes, ]) / Matrix::colSums(s_obj@raw.data)
  percent_mt <- round(percent_mt * 100, digits = 3)

  # add columns to object@meta.data, and is a great place to stash QC stats
  s_obj <- AddMetaData(s_obj, metadata = percent_mt, col.name = "percent.mito")

  # dist_unfilt_plot <- VlnPlot(s_obj, c("nGene", "nUMI", "percent.mito"),
  #                             nCol = 3, group.by = "orig.ident",
  #                             point.size.use = 0.1, x.lab.rot = TRUE, size.title.use = 12, cols.use = color_scheme
  # )
  # ggsave("qc.distribution.unfiltered.png", plot = dist_unfilt_plot, width = 10, height = 5, units = "in")
  plot_distribution(
    so = s_obj, features = c("nGene", "nUMI", "percent.mito"), grouping = "orig.ident",
    color_scheme = color_scheme,
    filename = glue("{qc_dir}/qc.distribution.unfiltered.png"), width = 10, height = 5
  )

  # check for high mitochondrial percentage or low UMI content
  png("qc.correlations.unfiltered.png", res = 200, width = 10, height = 5, units = "in")
  par(mfrow = c(1, 2))
  GenePlot(s_obj, gene1 = "nUMI", gene2 = "percent.mito", cex.use = 0.5, col.use = color_scheme)
  GenePlot(s_obj, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.5, col.use = color_scheme)
  dev.off()

  # check distribution of gene counts and mitochondrial percentage
  # low_quantiles <- c(0.05, 0.02, 0.01, 0.001)
  # high_quantiles <- c(0.95, 0.98, 0.99, 0.999)
  # message("nGene low percentiles:")
  # s_obj@meta.data$nGene %>% quantile(low_quantiles) %>% round(1) %>% print()
  # message(" ")
  # message("nGene high percentiles:")
  # s_obj@meta.data$nGene %>% quantile(high_quantiles) %>% round(1) %>% print()
  # message(" ")
  # message("percent.mito high percentiles:")
  # s_obj@meta.data$percent.mito %>% quantile(high_quantiles) %>% round(1) %>% print()
  # message(" ")

  return(s_obj)
}
