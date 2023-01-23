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

  # snn_graph <- Seurat:::FindNeighbors.default(
  # pcs[, 1:num_dim],
  # distance.matrix = FALSE,
  # k.param = num_neighbors,
  # compute.SNN = TRUE,
  # prune.SNN = 1 / 15,
  # nn.eps = 0,
  # force.recalc = TRUE
  # )

  t <- Seurat:::FindNeighbors.Seurat(seurat_obj, reduction = "pcalognorm")

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

#' Get cluster averages.
#'
#' @param data Gene expression data.
#' @param metadata Metadata.
#' @param group Column in metadata.
#'
#' @return .
#'
#' @import dplyr
#' @importFrom data.table .SD setDT
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

FindAllMarkers <- function(object,
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
