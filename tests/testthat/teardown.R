# the Seurat-object verbs (filter_cells(), normalize_counts(), run_pca(), run_tsne(), run_umap(),
# create_seurat_object(), integrate_seurat_object(), ...) write their plots and metadata tables
# automatically, flat into the working directory - that is the point, not a side effect to
# suppress. During a test run the working directory is tests/testthat/ itself (testthat resets it
# before every test file, which rules out redirecting it once for the whole run via a setup.R), so
# sweep up generated output here rather than letting it accumulate in the source tree. No fixture
# in tests/testthat/ uses these extensions, so matching by extension alone is safe.
generated_files <- list.files(pattern = "\\.(png|pdf|csv|csv\\.gz)$")
unlink(generated_files)
