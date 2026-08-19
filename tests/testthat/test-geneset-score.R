test_that("geneset score can be calculated", {
  counts <- read_counts_file(
    sample_name = "test",
    path = system.file("extdata", "", package = "scooter")
  )
  s_obj <- initialize_seurat_object(
    counts_matrix = counts$`Gene Expression`,
    assay = "RNA",
    log_file = NULL
  )
  s_obj <- add_seurat_assay(
    x = s_obj,
    assay = "ADT",
    counts_matrix = counts$`Antibody Capture`,
    log_file = NULL
  )
  s_obj_filt <- filter_cells(
    x = s_obj,
    log_file = NULL,
    min_counts = 1,
    min_genes = 1,
    max_genes = NULL,
    max_mt = 10
  )

  module_tbl <- data.frame(
    gene = c(
      "CST3",
      "TYROBP",
      "LST1",
      "AIF1",
      "FTL",
      "MALAT1",
      "LTB",
      "IL32",
      "IL7R",
      "CD2",
      "NAPSA",
      "GMFG"
    ),
    celltype = c(
      rep("cell.A", 6),
      rep("cell.B", 6)
    ),
    stringsAsFactors = FALSE
  )

  geneset <- geneset_score(
    counts_raw = as.matrix(GetAssayData(
      s_obj_filt,
      assay = "RNA",
      layer = "counts"
    )),
    module_tbl = module_tbl
  )

  geneset <- max_scores(scores = geneset, method = "test", threshold = 0)
  expect_equal(length(unique(geneset$module)), 2)
})
