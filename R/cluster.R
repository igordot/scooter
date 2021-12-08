#' Run dimensionality reduction, pca, tse, and umap
#'
#' @param pcs Data.
#' @param num_dim Number of PCs to use for tsne and umap.
#' @param num_neighbors Number of neighbors to use for umap.
#' @param log_file log file.
#' @param res Resolution
#' @param algorithm See Seurat::FindClusters()
#'
#' @return .
#'
#' @import dplyr
#' @import Seurat
#' @export
calculate_clusters <- function(pcs, num_dim, log_file, num_neighbors = 30, res = NULL, algorithm = 3) {
  # TODO: allow UMAP graphs to be used

  message_str <- "========== Seurat::FindNeighbors() =========="
  write_message(message_str, log_file)

  if (num_dim < 5) stop("too few dims: ", num_dim)
  if (num_dim > 50) stop("too many dims: ", num_dim)

  #snn_graph <- Seurat:::FindNeighbors.default(
   # pcs[, 1:num_dim],
   # distance.matrix = FALSE,
   # k.param = num_neighbors,
   # compute.SNN = TRUE,
   # prune.SNN = 1 / 15,
   # nn.eps = 0,
   # force.recalc = TRUE
  #)
  
  t = Seurat:::FindNeighbors.Seurat(seurat_obj, reduction = "pcalognorm")

  snn_graph <- as.Graph(t@graphs$RNA_snn)

  message_str <- "\n\n ========== Seurat::FindClusters() ========== \n\n"
  write_message(message_str, log_file)

  if (is.null(res)) {
    res_range <- seq(0.1, 2.5, 0.1)
    if (nrow(pcs) > 1000) res_range <- c(res_range, 3, 4, 5, 6, 7, 8, 9)
  } else {
    res_range <- res
  }

  clusters <- Seurat:::FindClusters.default(snn_graph,
    algorithm = algorithm,
    resolution = res_range,
    verbose = FALSE
  )
  return(clusters)
}

#' Calculate differential expression for
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param metadata_column Column in metadata.
#' @param list_groups dataframe of groups to compare in the metadata column.
#' @param proj_name Project name as a prefix for the output file.
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
#' @export
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

#' Calculate differential expression for
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
#' @export
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

#' Get cluster averages.
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param group Column in metadata.
#'
#' @return .
#'
#' @import dplyr
#' @import data.table
#' @export
calc_clust_averages <- function(metadata, data, group) {

  # get relevant metadata
  metadata <- metadata %>%
    select("cell", group)

  # manitpulate data to merge with metadata
  # genes are already columns, transpose, so genes are rows
  data <- data %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("cell")

  current_data <- full_join(metadata, data, by = "cell") %>%
    column_to_rownames("cell")

  current_data <- setDT(current_data)
  current_mean <- current_data[, lapply(.SD, mean), by = .(get(group)), .SDcols = 2:ncol(current_data)]

  colnames(current_mean)[which(colnames(current_mean) == "get")] <- group

  current_mean <- as.data.frame(current_mean)

  return(current_mean)
}

FindAllMarkers <- function(
                           object,
                           assay = NULL,
                           features = NULL,
                           logfc.threshold = 0.25,
                           test.use = "wilcox",
                           slot = "data",
                           min.pct = 0.1,
                           min.diff.pct = -Inf,
                           node = NULL,
                           verbose = TRUE,
                           only.pos = FALSE,
                           max.cells.per.ident = Inf,
                           random.seed = 1,
                           latent.vars = NULL,
                           min.cells.feature = 3,
                           min.cells.group = 3,
                           pseudocount.use = 1,
                           return.thresh = 1e-2,
                           ...) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    tree <- Tool(object = object, slot = "BuildClusterTree")
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers(
          object = object,
          assay = assay,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
          features = features,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          slot = slot,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          pseudocount.use = pseudocount.use,
          ...
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (class(x = genes.de[[i]]) == "character") {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c("bimod", "t")) {
        gde <- gde[order(gde$p_val, -gde[, 2]), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = gde.all[, 2] > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE, immediate. = TRUE)
  }
  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " versus all:\n\t", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  return(gde.all)
}
