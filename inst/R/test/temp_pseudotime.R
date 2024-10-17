data <- qs::qread("C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_4/NSM-NS-NC-M/NSM-NS-NC-M.qs")
cds <- SeuratWrappers::as.cell_data_set(data, assay = "SCT")

# load match_clustering function
cds <- match_clustering(data, cds, match_seurat_clustering = FALSE)

cds <- monocle3::learn_graph(cds)


for (partition in seq_along(table(monocle3::partitions(cds)))) {
  message("partition: ", partition)

  # get partition celltype
  partition_cells <- colnames(cds[, monocle3::partitions(cds) == partition])
  if ("mapmycells_supercluster" %in% colnames(cds@colData)) {
    partition_celltype <- names(which.max(table(cds@colData$mapmycells_supercluster[colnames(cds) %in% partition_cells])))
  } else if ("kriegstein.seurat.custom.clusters.mean" %in% colnames(cds@colData)) {
    partition_celltype <- names(which.max(table(cds@colData$kriegstein.seurat.custom.clusters.mean[colnames(cds) %in% partition_cells])))
  } else {
    partition_celltype <- "other"
  }

  if (partition_celltype %in% names(pseudotime_root_markers)) {
    genes_of_interest <- pseudotime_root_markers[[partition_celltype]]
  } else if (TRUE) {
    for (r_label in names(pseudotime_root_markers)) {
      if (grepl(tolower(r_label), partition_celltype)) {
        genes_of_interest <- pseudotime_root_markers[[r_label]]
        break
      }
    }
  } else {
    ## technical default, non-biologically relevant: pick gene with most expressino for this partition,
    genes_of_interest <- names(which.max(Matrix::rowSums(cds@assays@data$counts[colnames(cds) %in% partition_cells])))
  }


  # create df of cell names and gene(s) of interest values
  df <- SeuratObject::FetchData(data, genes_of_interest)
  # get cell name with highest summed expression for gene(s) of interest
  if (length(genes_of_interest %in% rownames(data)) > 1) {
    cell_name <- names(which.max(apply(df[monocle3::partitions(cds) == partition, ], 1, sum)))
  } else if (length(genes_of_interest) == 1) {
    cell_name <- rownames(df)[which.max(df[monocle3::partitions(cds) == partition, ])]
  }
  # calculate pseudo-time from the cell with highest expression
  cds <- monocle3::order_cells(cds, root_cells = cell_name)
  # get cell name at furthest point to invert pseudo-time calculations
  cell_name <- names(which.max(monocle3::pseudotime(cds)[monocle3::partitions(cds) == partition]))
  # recalculate pseudo-time from root cell furthest from cell with highest expression - as trajectory starting point
  cds <- monocle3::order_cells(cds, root_cells = cell_name)

  ### PLOTTING
}




pseudotime_root_markers <- list(
  "Microglia" = c("AIF1"),
  "Astrocyte" = c("VIM", "S100B", "SOX9"),
  "Neuron"    = c("MAP2", "DCX", "NEUROG2"),
  "Dividing"  = c("MKI67"),
  "other"     = c("FOXJ1")
)



# TODO update plot_cells.adjusted in package (doc install git)
# TODO  plot features with trajectory
plot_cells.adjusted(cds,
                    genes = c("S100B", "MAP2", "VIM", "APOE"),
                    label_cell_groups=FALSE,
                    show_trajectory_graph=T)





## TODO could do this per partition, in above for loop
# TODO try graph_test and finding modules and identifying genes on pseudotime trajectory
partition_cds <- cds[, monocle3::partitions(cds) == partition]
plot_cells.adjusted(partition_cds, color_cells_by="partition")
pr_graph_test_res <- monocle3::graph_test(partition_cds, neighbor_graph="knn", cores=8)
# Moran_I value: +1 is genes only locally expressed, 0 = no effect, -1 is anti-spatial-correlation
pr_graph_test_res <- pr_graph_test_res[order(-pr_graph_test_res$morans_I),]
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))





# Gene modules ------------------------------------------------------------
# TODO find gene modules
partition_cds <- monocle3::preprocess_cds(partition_cds, num_dim = 100)
partition_cds <- monocle3::reduce_dimension(partition_cds)
gene_module_df <- monocle3::find_gene_modules(partition_cds[pr_deg_ids,], resolution=1e-2)
# TODO plot modules
cell_group_df <- tibble::tibble(cell=row.names(SummarizedExperiment::colData(partition_cds)),
                                cell_group=monocle3::partitions(cds)[colnames(partition_cds)])
agg_mat <- monocle3::aggregate_gene_expression(partition_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
# TODO visualize modules
plot_cells.adjusted(partition_cds,
           genes=gene_module_df %>% dplyr::filter(module %in% c(8, 28, 33, 37)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)






# Graph test DE -----------------------------------------------------------
partition_cds <- cds[, monocle3::partitions(cds) == 3]

# cds_graph_test_res <- monocle3::graph_test(cds, neighbor_graph="principal_graph", cores=4)
partition_graph_test_res <- monocle3::graph_test(partition_cds, neighbor_graph="principal_graph", cores=4)
# TODO find genes DE based on the trajectory
## TODO try without partition selection
# partition_cds <- monocle3::reduce_dimension(cds = partition_cds, reduction_method = "UMAP")
# partition_cds <- monocle3::cluster_cells(cds = partition_cds, reduction_method = "UMAP")
# partition_cds <- monocle3::learn_graph(partition_cds, use_partition = TRUE)
# pr_graph_test_pg_res <- monocle3::graph_test(partition_cds, neighbor_graph="principal_graph", cores=8)
partition_genes <- row.names(subset(partition_graph_test_res, q_value < 0.05))

partition_graph_test_res[order(-partition_graph_test_res$morans_I), ]

partition_graph_test_res


# TODO identify top deg ids then plot with
plot_cells(cds, genes=c("RPS18", "RPL13"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)
plot_cells.adjusted(cds, genes = c("RPS18", "SETD9"),
                    show_trajectory_graph=T,
                    label_cell_groups=T,
                    label_leaves=FALSE)





monocle3::plot_genes_in_pseudotime(cds[c("RPS18", "SETD9", "MKI67", "MALAT1", "APOE", "SPARCL1"),],
                         min_expr=0.5, label_by_short_name = F)

# TODO look at specific branch point by select cells/clusters, then learning_graph and identifying interesting genes (markers of interest, and top degs by Moran_I value







# Misc --------------------------------------------------------------------

monocle3::graph_test() # find genes that vary over a trajectory
monocle3::find_gene_modules()
# investigate branch points
monocle3::plot_cells(genes = MARKERS)
