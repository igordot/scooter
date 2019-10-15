test_that("geneset score can be calculated", {
  counts <- load_sample_counts_matrix(sample_name = "test",
                                      path = system.file("extdata", "",
                                                         package = "scooter"))
  s_obj <- create_seurat_obj(counts_matrix = counts$`Gene Expression`,
                             assay = "RNA", log_file = NULL)
  s_obj <- add_seurat_assay(seurat_obj = s_obj,
                            assay = "ADT",
                            counts_matrix = counts$`Antibody Capture`,
                            log_file = NULL)
  s_obj_filt <- filter_data(data = s_obj,
                            log_file = NULL,
                            min_genes = NULL,
                            max_genes = NULL,
                            max_mt = 10)

  module_tbl <- data.frame(gene = c("CST3", "TYROBP", "LST1",
                                    "AIF1", "FTL", "MALAT1",
                                    "LTB", "IL32", "IL7R",
                                    "CD2", "NAPSA", "GMFG"),
                           celltype = c(rep("cell.A", 6),
                                        rep("cell.B", 6)),
                           stringsAsFactors =  FALSE)

  geneset <- geneset_score(counts_raw = as.matrix(s_obj_filt@assays$RNA@counts),
                module_tbl = module_tbl)

  geneset <- max_scores(scores = geneset, method = "test", threshold = 0)
  expect_equal(length(unique(geneset$module)), 2)
})