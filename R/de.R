#' Calculate differentially expressed genes within each subpopulation/cluster
#'
#' @param seurat_obj Gene expression data.
#' @param cluster_column Metadata column specifying the groups to split by.
#' @param group_column Metadata column specifying the groups for differential expressin within each split.
#' @param test Statistical method to use.
#' @param out_path Output path.
#' @param write Boolean to save results to disk.
#' @param log_file log file.
#'
#' @return .
#'
#' @import dplyr readr Seurat tibble
#' @export
differential_expression_per_cluster <- function(seurat_obj, cluster_column, group_column, test = "wilcox", out_path = ".", write = TRUE, log_file = NULL) {
  if (!dir.exists(out_path)) {
    stop("out_path does not exist")
  }

  # create a separate sub-directory for differential expression results
  cluster_column <- check_identity_column(seurat_obj, cluster_column)
  de_dir <- glue("{out_path}/differential-expression")
  if (!dir.exists(de_dir)) {
    dir.create(de_dir)
  }

  # results table
  de_all_genes_tbl <- tibble()

  # get DE genes for each cluster
  seurat_obj <- set_identity(seurat_obj = seurat_obj, identity_column = cluster_column)
  clusters <- levels(seurat_obj)
  for (clust_name in clusters) {
    message(glue("calculating DE genes for cluster {clust_name}"))

    # subset to the specific cluster
    clust_obj <- subset(seurat_obj, idents = clust_name)

    # revert back to the differential expression grouping variable labels
    clust_obj <- set_identity(seurat_obj = clust_obj, identity_column = group_column)

    message("cluster cells: ", ncol(clust_obj))
    message("cluster groups: ", paste(levels(clust_obj), collapse = ", "))

    # continue if cluster has multiple groups and more than 10 cells in each group
    if (n_distinct(Idents(clust_obj)) > 1 && min(table(Idents(clust_obj))) > 10) {

      # iterate through sample/library combinations (relevant if more than two)
      group_combinations <- combn(levels(clust_obj), m = 2, simplify = TRUE)
      for (combination_num in 1:ncol(group_combinations)) {

        # determine combination
        g1 <- group_combinations[1, combination_num]
        g2 <- group_combinations[2, combination_num]
        comparison_label <- glue("{g1}-vs-{g2}")
        message(glue("comparison: {clust_name} {g1} vs {g2}"))

        filename_label <- glue("{de_dir}/de.{cluster_column}-{clust_name}.{comparison_label}.{test}")

        # find differentially expressed genes (default Wilcoxon rank sum test)
        de_genes <- FindMarkers(clust_obj,
          ident.1 = g1, ident.2 = g2, assay = "RNA", test.use = test,
          min.pct = 0.1, logfc.threshold = 0, base = 2, fc.name = "log2FC", only.pos = FALSE,
          print.bar = FALSE
        )

        # perform some light filtering and clean up
        de_genes <-
          de_genes %>%
          tibble::rownames_to_column("gene") %>%
          dplyr::mutate(cluster = clust_name, group1 = g1, group2 = g2, de_test = test) %>%
          dplyr::select(cluster, group1, group2, de_test, gene, log2FC, p_val, p_val_adj) %>%
          dplyr::mutate(
            log2FC = round(log2FC, 3),
            p_val = if_else(p_val < 0.00001, p_val, round(p_val, 5)),
            p_val_adj = if_else(p_val_adj < 0.00001, p_val_adj, round(p_val_adj, 5))
          ) %>%
          dplyr::arrange(p_val_adj, p_val)

        message(glue("{comparison_label} num genes: {nrow(de_genes)}"))

        # save stats table
        if (write) {
          write_csv(de_genes, glue("{filename_label}.csv"))
        }
        # add cluster genes to all genes
        de_all_genes_tbl <- bind_rows(de_all_genes_tbl, de_genes)
      }
    } else {
      message("skip cluster: ", clust_name)
    }

    message(" ")
  }

  # save stats table
  if (write) {
    write_csv(de_all_genes_tbl, glue("{de_dir}/de.{cluster_column}.{group_column}.{test}.all.csv"))
    de_sig_genes_tbl <- de_all_genes_tbl %>% dplyr::filter(p_val_adj < 0.01)
    write_csv(de_sig_genes_tbl, glue("{de_dir}/de.{cluster_column}.{group_column}.{test}.sig.csv"))
    de_top_genes_tbl <- de_all_genes_tbl %>%
      dplyr::group_by(cluster, group1, group2) %>%
      dplyr::slice_min(p_val_adj, n = 100) %>%
      dplyr::slice_min(p_val, n = 100)
    write_csv(de_top_genes_tbl, glue("{de_dir}/de.{cluster_column}.{group_column}.{test}.top100.csv"))
  }

  return(de_all_genes_tbl)
}

#' Calculate differential expression between specific groups
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param metadata_column Column in metadata.
#' @param list_groups dataframe of groups to compare in the metadata column.
#' @param log_fc_threshold Log fc threshold.
#' @param min.pct Minimum percentage of cells a gene must be expressed in to be tested.
#' @param test.use Test to use.
#' @param out_path output path.
#' @param write Boolean to write or not.
#' @param log_file log file.
#'
#' @return .
#'
#' @import dplyr
#' @import Seurat
differential_expression_paired <- function(data, metadata, metadata_column, list_groups = NULL, log_fc_threshold = 0.5, min.pct = 0.1, test.use = "wilcox", out_path = ".", write = FALSE, log_file = NULL) {

  # get unique combinations of the factors in the metadata column of interest.
  # Compare all possible combinations
  if (is.null(list_groups)) {
    # get unique combinations
    unique_ids <- unique(unique(metadata[, metadata_column]))[[1]]
    list_groups <- expand.grid(unique_ids,
      unique_ids,
      stringsAsFactors = FALSE
    )
    list_groups <- list_groups %>%
      filter(Var1 != Var2)
    indx <- !duplicated(t(apply(list_groups, 1, sort)))
    list_groups <- list_groups[indx, ]
  } else {
    list_groups <- list_groups
  }

  # If there is a column called cell leave it, if not, make rownames cell
  if ("cell" %in% colnames(metadata)) {
    metadata <- metadata %>% as.tibble()
  } else {
    metadata <- metadata %>% as_tibble(rownames = "cell")
  }

  diff_exp <- list()

  for (j in 1:nrow(list_groups)) {

    # Current comparison
    current_group <- list_groups[j, ]
    group1 <- current_group[1][[1]]
    group2 <- current_group[2][[1]]

    # get the cells for one side
    cell_group1 <- metadata %>%
      filter(get(metadata_column) == group1) %>%
      select("cell")

    # get the cells for the other side
    cell_group2 <- metadata %>%
      filter(get(metadata_column) == group2) %>%
      select("cell")

    message_str <- glue("{group1} versus {group2}")
    write_message(message_str, log_file)

    # Run FindMarkers. If there aren't any genes that pass the logfc
    # cutoff or something, then go to the next comparison
    current_comparison <- tryCatch(
      {
        current_comparison <- Seurat:::FindMarkers.default(
          object = as.matrix(data),
          reduction = NULL,
          slot = "data",
          cells.1 = cell_group1$cell,
          cells.2 = cell_group2$cell,
          logfc.threshold = log_fc_threshold,
          test.use = test.use,
          min.pct = min.pct
        )

        current_comparison_filt <- current_comparison %>%
          select(p_val, avg_logFC, p_val_adj) %>%
          rownames_to_column("gene") %>%
          mutate(group_1 = rep(group1)) %>%
          mutate(group_2 = rep(group2))
      },
      error = function(e) {
        e
      }
    )

    if (inherits(current_comparison, "error")) next

    # write out each comparison
    if (write) {
      write_excel_csv(current_comparison_filt,
        path = glue("{out_path}/diffexp-{metadata_column}-{group1}.{group2}.csv")
      )
    }

    # save data in a list
    diff_exp[[glue("{group1}.{group2}")]] <- current_comparison_filt
  }

  return(diff_exp)
}

#' Calculate differential expression for one group versus all
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param metadata_column Column in metadata.
#' @param log_fc_threshold Log fc threshold.
#' @param min.pct Minimum percentage of cells a gene must be expressed in to be tested.
#' @param test.use Test to use.
#' @param out_path output path.
#' @param write Boolean to write or not.
#' @param log_file log file.
#'
#' @return .
#'
#' @import dplyr
#' @import Seurat
differential_expression_global <- function(data, metadata, metadata_column,
                                           log_fc_threshold = 0.5, min.pct = 0.1,
                                           test.use = "wilcox", out_path = ".",
                                           write = FALSE, log_file = NULL) {


  # If there is a column called cell leave it, if not, make rownames cell
  if ("cell" %in% colnames(metadata)) {
    metadata <- metadata %>%
      as_tibble()
  } else {
    metadata <- metadata %>%
      as_tibble(rownames = "cell")
  }

  diff_exp <- list()

  list_groups <- unique(select(.data = metadata, !!sym(metadata_column)))[[1]]

  for (j in 1:length(list_groups)) {
    current_group <- list_groups[j]

    cell_group1 <- metadata %>%
      filter(get(metadata_column) == current_group) %>%
      select("cell")

    cell_group2 <- metadata %>%
      filter(get(metadata_column) != current_group) %>%
      select("cell")

    message_str <- glue("{current_group} versus all")
    write_message(message_str, log_file)

    current_comparison <- tryCatch(
      {
        current_comparison <- Seurat:::FindMarkers.default(
          object = as.matrix(data),
          reduction = NULL,
          slot = "data",
          cells.1 = cell_group1$cell,
          cells.2 = cell_group2$cell,
          logfc.threshold = log_fc_threshold,
          test.use = test.use,
          min.pct = min.pct
        )

        current_comparison_filt <- current_comparison %>%
          select(p_val, avg_logFC, p_val_adj) %>%
          rownames_to_column("gene") %>%
          mutate(group_1 = rep(current_group)) %>%
          mutate(group_2 = rep("All"))
      },
      error = function(e) {
        e
      }
    )

    if (inherits(current_comparison, "error")) next

    if (write) {
      write_excel_csv(current_comparison_filt,
        path = glue("{out_path}/diffexp-{metadata_column}-{current_group}.All.csv")
      )
    }

    diff_exp[[glue("{current_group}.All")]] <- current_comparison_filt
  }
  return(diff_exp)
}
