
run_label_transfer <- TRUE

# run label transfer for integrated analysis after plotting regular UMAP
if (run_label_transfer) {
  dir.create('Label_transfer/')
  # TODO if data already loaded at individual analysis, then no need to load it here
  ### get Nowakowski as reference (train) data
  nowakowski_data_path <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Nowakowski/exprMatrix.tsv.gz'
  nowakowski_metadata_path <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Nowakowski/meta.tsv'
  nowakowski_project <- 'adultPancreas'
  # get Nowakowski matrix
  mat <- data.table::fread(nowakowski_data_path)
  # add metadata
  meta <- utils::read.table(nowakowski_metadata_path, header=T, sep="\t", as.is=T, row.names=1)
  # get genes
  genes = mat[,1][[1]]
  genes = gsub(".+[|]", "", genes)
  # get overlapping genes between test and train dataset (to reduce noise and run-time)
  genes.overlap <- rownames(integrated@assays$SCT)[rownames(integrated@assays$SCT) %in% genes]
  # save non-overlap genes for retrospection
  ## most genes contain a: ['-', '_', '.', 'orf']
  genes.non.overlap <- sort(genes[!genes %in% rownames(integrated@assays$SCT)])
  utils::write.csv2(genes.non.overlap, file = paste0("Label_transfer/non_overlapping_ref_genes_", sample_name, ".csv"))

  # set data.frame
  mat = data.frame(mat[,-1], row.names=genes)
  # set data.frame to matrix to dgCMatrix for memory issues
  mat <- as(object = as(mat, "matrix"), Class = "dgCMatrix")
  # get SO from matrix object
  so <- SeuratObject::CreateSeuratObject(counts = mat, project = nowakowski_project, meta.data=meta)
  # subset train data by overlapping genes
  so <- subset(so, cells = names(so$WGCNAcluster)[so$WGCNAcluster != ''], features = genes.overlap)
  # set SO to SCE object
  so.sce <- Seurat::as.SingleCellExperiment(so)
  # subset test data by overlapping genes
  integrated.sub <- subset(integrated, features = genes.overlap)
  # set SO to SCE object
  integrated.sce <- Seurat::as.SingleCellExperiment(integrated.sub, assay = 'SCT')

  # perform SingleR label transfer
  integrated.sce.labels <- SingleR::SingleR(
    test=integrated.sce,
    ref=so.sce,
    labels=so.sce$WGCNAcluster,
    de.method='wilcox',
    aggr.ref=TRUE)
  # TODO or check add aggr.ref = T here too?
  integrated.sce.clusters <- SingleR::SingleR(
    test=integrated.sce,
    ref=so.sce,
    labels=so.sce$WGCNAcluster,
    clusters=integrated.sce$seurat_clusters,
    de.method='wilcox')

  # get a quick a idea of label distribution
  table(integrated.sce.labels$labels)
  table(integrated.sce.clusters$labels)
  # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
  png(paste0("Label_transfer/cell_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  plot(SingleR::plotScoreHeatmap(integrated.sce.labels))
  dev.off()
  png(paste0("Label_transfer/cluster_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  png(SingleR::plotScoreHeatmap(integrated.sce.clusters))
  dev.off()
  # ambiguity 'predicted' by 'deltas', i.e. difference from mean scores for all labels (higher=better resolution)
  # png(paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"))
  p <- SingleR::plotDeltaDistribution(integrated.sce.labels, ncol = 6)
  ggplot2::ggsave(file=paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # png(paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"))
  p <- SingleR::plotDeltaDistribution(integrated.sce.clusters, ncol = 6)
  ggplot2::ggsave(file=paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # print warning if too much pruning happens off of default settings, could mean bad quality of label transfer
  if (table(SingleR::pruneScores(integrated.sce.labels))[1] < length(integrated.sce.labels$labels)*0.99) {
    print("WARNING: SingleR::pruneScores function showed more cells pruned than 99% margin of Gauss dist, look at quality of label transfer")
  }
  # save the cells label transfer results in the Seurat object
  integrated$SingleR.cell.labels <- integrated.sce.labels$labels
  integrated$SingleR.cell.first.labels <- integrated.sce.labels$first.labels

  # add seurat cluster to transfer label for later identification with DE analysis
  integrated.sce.clusters$labels <- paste(as.character(integrated.sce.clusters$labels), levels(integrated), sep='.')
  integrated.sce.clusters$first.labels <- paste(as.character(integrated.sce.clusters$first.labels), levels(integrated), sep='.')

  # original Seurat clusters are needed for RenameIdents function
  names(integrated.sce.clusters$first.labels) <- levels(integrated)
  # set idents/names of clusters to labels of transfer
  integrated <- SeuratObject::RenameIdents(integrated, integrated.sce.clusters$first.labels)
  p1 <- Seurat::DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.first.labels', repel = TRUE) +
    ggplot2::labs(title = "Cells first labels") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) +
    ggplot2::labs(title = "Clusters first labels") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # original Seurat clusters are needed for RenameIdents function
  names(integrated.sce.clusters$labels) <- levels(integrated)
  # set idents/names of clusters to labels of transfer
  integrated <- SeuratObject::RenameIdents(integrated, integrated.sce.clusters$labels)
  # create UMAP on cell and cluster level with transferred labels!
  p3 <- Seurat::DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.labels', repel = TRUE) +
    ggplot2::labs(title = "Cells tuned labels") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p4 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) +
    ggplot2::labs(title = "Clusters tuned labels") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  # create and save total visualizations
  plot_list <- list(p2, p4, p1, p3)
  p <- do.call("grid.arrange", c(plot_list, ncol=2))
  ggplot2::ggsave(paste0("Label_transfer/labels_UMAP_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))

}
