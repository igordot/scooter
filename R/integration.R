
gather_seurat_objects <- function(seurat_obj_paths, assay){
  # seurat obj paths is a named list

  seurat_obj_list = list()
  var_genes_list = list()

  for (i in 1:length(seurat_obj_paths)) {

    current_seurat_obj = seurat_obj_paths[[i]]
    current_seurat_obj_name = names(seurat_obj_paths[i])
    if(is.null(current_seurat_obj_name)) stop("The list is not named")


    if (!file.exists(current_seurat_obj)) stop(glue("seurat object rds {current_seurat_obj} does not exist"))

    # load seurat object
    seurat_obj_list[[i]] = readRDS(current_seurat_obj)

    # clean up object
    seurat_obj_list[[i]]@assays$RNA@scale.data = matrix()
    seurat_obj_list[[i]]@reductions = list()
    seurat_obj_list[[i]]@meta.data = seurat_obj_list[[i]]@meta.data %>% select(-starts_with("snn_res"))


    var_genes_list[[current_seurat_obj_name]] = VariableFeatures(seurat_obj_list[[i]], assay = assay)

  }
  return(list(seurat_objects = seurat_obj_list,
              variable_genes = var_genes_list))
}

integrate_seurat_log <- function(seurat_obj_list, num_dim){
  message("\n\n ========== Seurat::FindIntegrationAnchors() ========== \n\n")

  # find the integration anchors
  anchors = FindIntegrationAnchors(object.list = seurat_obj_list, anchor.features = 2000, dims = 1:num_dim)
  rm(seurat_obj_list)

  message("\n\n ========== Seurat::IntegrateData() ========== \n\n")

  # integrating all genes may cause issues and may not add any relevant information
  # integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim, features.to.integrate = exp_genes)
  integrated_obj = IntegrateData(anchorset = anchors, dims = 1:num_dim)
  rm(anchors)
  return(integrated_obj)
}

integrate_seurat_sct <- function(seurat_obj_list, num_dim){

  s_obj_features <- SelectIntegrationFeatures(object.list = seurat_obj_list, nfeatures = 3000)

  seurat_obj_list <- lapply(X = seurat_obj_list, FUN = function(x) {
    x <- RunPCA(x, features = s_obj_features, verbose = FALSE) })

  seurat_obj_list <- PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = s_obj_features)

  anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT",
                                    anchor.features = s_obj_features, dims = 1:num_dim, reduction = 'rpca')

  integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

  return(integrated_obj)

}