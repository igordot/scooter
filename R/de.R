#' Calculate differentially expressed genes within each subpopulation/cluster
#'
#' @param x A Seurat object.
#' @param cluster_column Metadata column specifying the groups to split by.
#' @param group_column Metadata column specifying the groups for differential
#'   expressin within each split.
#' @param test Statistical method to use.
#' @param out_path Output path.
#' @param write Boolean to save results to disk.
#' @param log_file log file.
#'
#' @return .
#'
#' @importFrom readr write_csv
#' @export
differential_expression_per_cluster <- function(
  x,
  cluster_column,
  group_column,
  test = "wilcox",
  out_path = ".",
  write = TRUE,
  log_file = NULL
) {
  if (!dir.exists(out_path)) {
    stop("out_path does not exist")
  }

  # create a separate sub-directory for differential expression results
  cluster_column <- check_identity_column(x, cluster_column)
  de_dir <- glue("{out_path}/differential-expression")
  if (!dir.exists(de_dir)) {
    dir.create(de_dir)
  }

  # results table
  de_all_genes_tbl <- tibble()

  # get DE genes for each cluster
  x <- set_identity(x = x, identity_column = cluster_column)
  clusters <- levels(x)
  for (clust_name in clusters) {
    message(glue("calculating DE genes for cluster {clust_name}"))

    # subset to the specific cluster
    clust_obj <- subset(x, idents = clust_name)

    # revert back to the differential expression grouping variable labels
    clust_obj <- set_identity(x = clust_obj, identity_column = group_column)

    message("cluster cells: ", ncol(clust_obj))
    message("cluster groups: ", paste(levels(clust_obj), collapse = ", "))

    # continue if cluster has multiple groups and more than 10 cells in each
    # group
    if (
      n_distinct(Idents(clust_obj)) > 1 && min(table(Idents(clust_obj))) > 10
    ) {
      # iterate through sample/library combinations (relevant if more than two)
      group_combinations <- combn(levels(clust_obj), m = 2, simplify = TRUE)
      for (combination_num in 1:ncol(group_combinations)) {
        # determine combination
        g1 <- group_combinations[1, combination_num]
        g2 <- group_combinations[2, combination_num]
        comparison_label <- glue("{g1}-vs-{g2}")
        message(glue("comparison: {clust_name} {g1} vs {g2}"))

        filename_label <- glue(
          "{de_dir}/de-{cluster_column}-{clust_name}-{comparison_label}-{test}"
        )

        # find differentially expressed genes (default Wilcoxon rank sum test)
        de_genes <- FindMarkers(
          clust_obj,
          ident.1 = g1,
          ident.2 = g2,
          assay = "RNA",
          test.use = test,
          min.pct = 0.1,
          logfc.threshold = 0,
          base = 2,
          fc.name = "log2FC",
          only.pos = FALSE,
          print.bar = FALSE
        )

        # perform some light filtering and clean up
        de_genes <-
          de_genes |>
          tibble::rownames_to_column("gene") |>
          dplyr::mutate(
            cluster = clust_name,
            group1 = g1,
            group2 = g2,
            de_test = test
          ) |>
          dplyr::select(
            cluster,
            group1,
            group2,
            de_test,
            gene,
            log2FC,
            p_val,
            p_val_adj
          ) |>
          dplyr::mutate(
            log2FC = round(log2FC, 5),
            p_val = if_else(p_val < 0.00001, p_val, round(p_val, 5)),
            p_val_adj = if_else(
              p_val_adj < 0.00001,
              p_val_adj,
              round(p_val_adj, 5)
            )
          ) |>
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
    write_csv(
      de_all_genes_tbl,
      glue("{de_dir}/de-{cluster_column}-{group_column}-{test}-all.csv")
    )
    de_sig_genes_tbl <- de_all_genes_tbl |> dplyr::filter(p_val_adj < 0.05)
    write_csv(
      de_sig_genes_tbl,
      glue("{de_dir}/de-{cluster_column}-{group_column}-{test}-sig.csv")
    )
    de_top_genes_tbl <- de_all_genes_tbl |>
      dplyr::group_by(cluster, group1, group2) |>
      dplyr::slice_min(p_val_adj, n = 100) |>
      dplyr::slice_min(p_val, n = 100)
    write_csv(
      de_top_genes_tbl,
      glue("{de_dir}/de-{cluster_column}-{group_column}-{test}-top100.csv")
    )
  }

  return(de_all_genes_tbl)
}
