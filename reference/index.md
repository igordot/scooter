# Package index

## All functions

- [`add_gene_class_percent()`](https://igordot.github.io/scooter/reference/add_gene_class_percent.md)
  : Add mitochondrial, ribosomal, and hemoglobin percentages to a Seurat
  object.
- [`add_seurat_assay()`](https://igordot.github.io/scooter/reference/add_seurat_assay.md)
  : Add assay to Seurat object.
- [`as_data_frame_seurat()`](https://igordot.github.io/scooter/reference/as_data_frame_seurat.md)
  : Function to extract data from Seurat object.
- [`calculate_cluster_expression()`](https://igordot.github.io/scooter/reference/calculate_cluster_expression.md)
  : Per-cluster average expression (non-log space).
- [`calculate_cluster_markers()`](https://igordot.github.io/scooter/reference/calculate_cluster_markers.md)
  : Calculate cluster markers (versus all other clusters, or pairwise)
  and plot them.
- [`calculate_cluster_stats()`](https://igordot.github.io/scooter/reference/calculate_cluster_stats.md)
  : Per-cluster cell counts and a joined metadata+embeddings table.
- [`calculate_clusters()`](https://igordot.github.io/scooter/reference/calculate_clusters.md)
  : Identify clusters of cells by graph-based clustering.
- [`check_identity_column()`](https://igordot.github.io/scooter/reference/check_identity_column.md)
  : Check identity of the Seurat object.
- [`cluster_seurat_object()`](https://igordot.github.io/scooter/reference/cluster_seurat_object.md)
  : Reduce dimensions, cluster, and plot the result.
- [`create_color_vect()`](https://igordot.github.io/scooter/reference/create_color_vect.md)
  : Function to create a color vector.
- [`create_seurat_object()`](https://igordot.github.io/scooter/reference/create_seurat_object.md)
  : Create a usable Seurat object from a sample.
- [`differential_expression_per_cluster()`](https://igordot.github.io/scooter/reference/differential_expression_per_cluster.md)
  : Calculate differentially expressed genes within each
  subpopulation/cluster
- [`dr_plot_width()`](https://igordot.github.io/scooter/reference/dr_plot_width.md)
  : Width for a reduction plot, widened so many legend labels stay
  legible.
- [`filter_cells()`](https://igordot.github.io/scooter/reference/filter_cells.md)
  : Filter cells based on the number of genes, counts, and mitochondrial
  reads.
- [`geneset_score()`](https://igordot.github.io/scooter/reference/geneset_score.md)
  : Get geneset scores.
- [`get_color_scheme()`](https://igordot.github.io/scooter/reference/get_color_scheme.md)
  : Determine the color scheme.
- [`get_dr_point_size()`](https://igordot.github.io/scooter/reference/get_dr_point_size.md)
  : Determine the point size for reduced dimensions scatter plots
  (smaller for larger datasets).
- [`get_test_counts_matrix()`](https://igordot.github.io/scooter/reference/get_test_counts_matrix.md)
  : Get an example counts matrix.
- [`import_mtx()`](https://igordot.github.io/scooter/reference/import_mtx.md)
  : Read in 10x Genomics Cell Ranger Matrix Market format data.
- [`initialize_seurat_object()`](https://igordot.github.io/scooter/reference/initialize_seurat_object.md)
  : Create a new Seurat object from a matrix.
- [`integrate_layers()`](https://igordot.github.io/scooter/reference/integrate_layers.md)
  : Integrate the layers of a Seurat object.
- [`integrate_seurat_object()`](https://igordot.github.io/scooter/reference/integrate_seurat_object.md)
  : Integrate a merged Seurat object across batches.
- [`merge_metadata()`](https://igordot.github.io/scooter/reference/merge_metadata.md)
  : Function to merge two metadata tables together.
- [`merge_seurat_objects()`](https://igordot.github.io/scooter/reference/merge_seurat_objects.md)
  : Merge multiple Seurat objects into one.
- [`normalize_counts()`](https://igordot.github.io/scooter/reference/normalize_counts.md)
  : Normalize the counts, select variable features, and scale.
- [`plot_cluster_markers_heatmap()`](https://igordot.github.io/scooter/reference/plot_cluster_markers_heatmap.md)
  : Heatmap of the top cluster markers, at several top-N cutoffs.
- [`plot_cluster_markers_top()`](https://igordot.github.io/scooter/reference/plot_cluster_markers_top.md)
  : Plot a gene list in multiple ways: UMAP, dot plot, violin, and
  per-cluster bar plot.
- [`plot_clusters()`](https://igordot.github.io/scooter/reference/plot_clusters.md)
  : Plot a cluster resolution on tSNE and UMAP.
- [`plot_distribution()`](https://igordot.github.io/scooter/reference/plot_distribution.md)
  : Plot the distribution of specified features/variables.
- [`plot_dr_feature()`](https://igordot.github.io/scooter/reference/plot_dr_feature.md)
  : Plot a single feature overlaid on a dimensionality reduction.
- [`plot_dr_group()`](https://igordot.github.io/scooter/reference/plot_dr_group.md)
  : Scatter plot of a reduction, colored by a grouping variable.
- [`plot_dr_umap_clusters()`](https://igordot.github.io/scooter/reference/plot_dr_umap_clusters.md)
  : Plot a UMAP colored by one or more cluster/grouping variables.
- [`plot_metrics_correlations()`](https://igordot.github.io/scooter/reference/plot_metrics_correlations.md)
  : Scatter plots of the per-cell QC metrics against each other.
- [`plot_metrics_distribution()`](https://igordot.github.io/scooter/reference/plot_metrics_distribution.md)
  : Violin plots of the per-cell QC metrics.
- [`plot_var_genes_euler()`](https://igordot.github.io/scooter/reference/plot_var_genes_euler.md)
  : Euler diagram of the variable genes shared between batches.
- [`plot_var_genes_upset()`](https://igordot.github.io/scooter/reference/plot_var_genes_upset.md)
  : UpSet plot of the variable genes shared between batches.
- [`profile_object_size()`](https://igordot.github.io/scooter/reference/profile_object_size.md)
  : Report the size of an object's slots, largest first.
- [`read_counts_file()`](https://igordot.github.io/scooter/reference/read_counts_file.md)
  : Read in Gene Expression and Antibody Capture data from a 10x
  Genomics Cell Ranger sparse matrix or from a text file.
- [`resolve_seurat_object()`](https://igordot.github.io/scooter/reference/resolve_seurat_object.md)
  : Resolve a Seurat object, given either the object itself or a path to
  one on disk.
- [`run_dr()`](https://igordot.github.io/scooter/reference/run_dr.md) :
  Run dimensionality reduction: pca, tsne, or umap
- [`run_pca()`](https://igordot.github.io/scooter/reference/run_pca.md)
  : Run PCA
- [`run_tsne()`](https://igordot.github.io/scooter/reference/run_tsne.md)
  : Run TSNE
- [`run_umap()`](https://igordot.github.io/scooter/reference/run_umap.md)
  : Run UMAP
- [`save_counts()`](https://igordot.github.io/scooter/reference/save_counts.md)
  : Save the counts matrix as a single table.
- [`save_dr_plot()`](https://igordot.github.io/scooter/reference/save_dr_plot.md)
  : Save a reduction plot, widened to fit its own legend.
- [`save_metadata()`](https://igordot.github.io/scooter/reference/save_metadata.md)
  : Save the cell metadata and the reduction embeddings as a single
  table.
- [`set_identity()`](https://igordot.github.io/scooter/reference/set_identity.md)
  : Set identity of the Seurat object.
- [`split_layers_by_batch()`](https://igordot.github.io/scooter/reference/split_layers_by_batch.md)
  : Split a Seurat object into per-batch layers and prepare it for
  integration.
- [`transfer_labels()`](https://igordot.github.io/scooter/reference/transfer_labels.md)
  : Transfer labels from a reference Seurat object
- [`variable_features_by_batch()`](https://igordot.github.io/scooter/reference/variable_features_by_batch.md)
  : Variable features of each batch of a layer-split object.
- [`write_message()`](https://igordot.github.io/scooter/reference/write_message.md)
  : Small function to write to message and to log file.
