
max_scores <- function(scores, method, threshold = 0) {
  # gets max pos or zero
  # idk how to do Igors yet
  scores$unknown <- rep(threshold, nrow(scores))
  scores <- scores %>%
    column_to_rownames("cell")
  scores$module <- colnames(scores)[apply(scores,1,which.max)]
  return(scores)
}

#' Get geneset scores.
#'
#' @param module_tbl geneset table.
#' @param counts_raw Raw counts
#' @param min_cpm .
#' @param limit_pct .
#'
#' @return .
#'
#' @import tidyr
#' @export
geneset_score = function(module_tbl, counts_raw, min_cpm = 0, limit_pct = 1) {
  # perform the cell type enrichment calculation based on rescaled values

  module_list <- module_tbl %>%
    filter(.$gene %in% rownames(counts_raw)) %>%
    with(split(.$gene, celltype))

  if (!is(counts_raw, "matrix")) { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }

  # filter matrix for expressed genes only
  counts_raw = filter_mat_by_cpm(counts_raw = counts_raw, min_cpm = min_cpm)

  # rescale matrix for expressed genes only
  counts_raw_subs = normalize_mat_by_gene(counts_raw = counts_raw[unlist(module_list), ], limit_pct = limit_pct)


  # check if enough genes pass filter
  if (min(lengths(module_list)) < 3) { stop("too few genes per celltype") }

  # calculate average z-score per celltype
  celltype_scores_tbl = tibble()
  for (ct in names(module_list)) {
    celltype_scores_tbl =
      bind_rows(
        celltype_scores_tbl,
        tibble(
          cell = colnames(counts_raw_subs),
          celltype = ct,
          score = colMeans(counts_raw_subs[module_list[[ct]], ])
        )
      )
    ct_scores = colnames(counts_raw_subs)
  }

  celltype_scores_tbl <- celltype_scores_tbl %>%
    spread(celltype, score) %>%
    rename_at(vars(-contains("cell")), list(~paste0(., ".score")))

  return(celltype_scores_tbl)
}

filter_mat_by_cpm = function(counts_raw, min_cpm = 0) {
  # filter matrix by a specified CPM value (higher in at least one sample/column for each gene/row)

  if (!is(counts_raw, "matrix")) { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }
  #if (nrow(counts_raw) < 10000) { stop("expression matrix has too few genes") }

  # expression level equivalent to 1 CPM (1 for 1m, 0.01 for 10k)
  exp_cpm1 = (1 / 1000000) * median(colSums(counts_raw))
  # expression level equivalent to the requested CPM
  min_exp = exp_cpm1 * min_cpm
  # filtered expression matrix
  counts_raw = counts_raw[matrixStats::rowMaxs(counts_raw) > min_exp, ]

  return(counts_raw)

}

rescale_vector = function(x, limit_pct = 1) {
  x / quantile(x, limit_pct)
}

normalize_mat_by_gene = function(counts_raw, limit_pct = 1) {

  if (limit_pct > 1) { stop("percentile should be expressed as a fraction") }
  if (!is(counts_raw, "matrix")) { stop("expression matrix is not a matrix") }
  if (max(counts_raw) < 100) { stop("expression values appear to be log-scaled") }

  counts_raw = apply(counts_raw, MARGIN = 1, FUN = rescale_vector, limit_pct = limit_pct)
  counts_raw = t(counts_raw)
  counts_raw[counts_raw > 1] = 1

  return(counts_raw)

}
