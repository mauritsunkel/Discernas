# run SingleR label transfer before visualizing UMAP, if run take into account labels
if (run_label_transfer) {
  dir.create('Label_transfer/')
  ### get Nowakowski as reference (train) data
  nowakowski_data_path <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Nowakowski/exprMatrix.tsv.gz'
  nowakowski_metadata_path <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Nowakowski/meta.tsv'
  nowakowski_project <- 'adultPancreas'
  # get Nowakowski matrix
  mat <- data.table::fread(nowakowski_data_path)
  # add metadata
  meta <- utils::read.table(nowakowski_metadata_path, header=T, sep="\t", as.is=T, row.names=1)
  # get genes
  genes = mat[,1][[1]]
  genes = gsub(".+[|]", "", genes)

  # get overlapping genes between test and train dataset (to reduce noise and run-time)
  genes.overlap <- rownames(data@assays$SCT)[rownames(data@assays$SCT) %in% genes]
  # save non-overlap genes for retrospection
  ## most genes contain a: ['-', '_', '.', 'orf']
  genes.non.overlap <- sort(genes[!genes %in% rownames(data@assays$SCT)])
  utils::write.csv2(genes.non.overlap, file = paste0("Label_transfer/non_overlapping_ref_genes_", sample_name, ".csv"))

  # set matrix
  mat = data.frame(mat[,-1], row.names=genes)
  # get SO from matrix object
  so <- SeuratObject::CreateSeuratObject(counts = mat, project = nowakowski_project, meta.data=meta)
  # subset train data by overlapping genes
  ## TODO base::subset or SeuratObject::subset
  so <- subset(so, cells = names(so$WGCNAcluster)[so$WGCNAcluster != ''], features = genes.overlap)
  # set SO to SCE object
  so.sce <- Seurat::as.SingleCellExperiment(so)

  # subset test data by overlapping genes
  data.sub <- subset(data, features = genes.overlap)
  # set SO to SCE object
  data.sce <- Seurat::as.SingleCellExperiment(data.sub, assay = 'SCT')

  # perform SingleR label transfer
  data.sce.labels <- SingleR(test=data.sce, ref=so.sce, labels=so.sce$WGCNAcluster, de.method='wilcox')
  data.sce.clusters <- SingleR(test=data.sce, ref=so.sce, labels=so.sce$WGCNAcluster,
                               clusters=data.sce$seurat_clusters, de.method='wilcox')
  # get a quick a idea of label distribution
  table(data.sce.labels$labels)
  table(data.sce.clusters$labels)
  # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
  png(paste0("Label_transfer/cell_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  plot(plotScoreHeatmap(data.sce.labels))
  dev.off()
  png(paste0("Label_transfer/cluster_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  plot(plotScoreHeatmap(data.sce.clusters))
  dev.off()
  # ambiguity 'predicted' by 'deltas', i.e. difference from mean scores for all labels (higher=better resolution)
  # png(paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"))
  p <- plotDeltaDistribution(data.sce.labels, ncol = 6)
  ggsave(file=paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # png(paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"))
  p <- plotDeltaDistribution(data.sce.clusters, ncol = 6)
  ggsave(file=paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # print warning if too much pruning happens off of default settings, could mean bad quality of label transfer
  if (table(pruneScores(data.sce.labels))[1] < length(data.sce.labels$labels)*0.99) {
    print("WARNING: SingleR::pruneScores function showed more cells pruned than 99% margin of Gauss dist, look at quality of label transfer")
  }
  # save the label transfer results in the Seurat object
  data$SingleR.cell.labels <- data.sce.labels$labels
  data$SingleR.cell.first.labels <- data.sce.labels$first.labels

  # add seurat cluster to transfer label for later identification with DE analysis
  data.sce.clusters$labels <- paste(as.character(data.sce.clusters$labels), levels(data), sep='.')
  data.sce.clusters$first.labels <- paste(as.character(data.sce.clusters$first.labels), levels(data), sep='.')

  # original Seurat clusters are needed for RenameIdents function
  names(data.sce.clusters$first.labels) <- levels(data)
  # set idents/names of clusters to labels of transfer
  data <- RenameIdents(data, data.sce.clusters$first.labels)
  p1 <- DimPlot(data, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.first.labels', repel = TRUE)
  p2 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) # + NoLegend()

  # original Seurat clusters are needed for RenameIdents function
  names(data.sce.clusters$labels) <- levels(data)
  # set idents/names of clusters to labels of transfer
  data <- RenameIdents(data, data.sce.clusters$labels)
  # create UMAP on cell and cluster level with transferred labels!
  p3 <- DimPlot(data, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.labels', repel = TRUE)
  p4 <- DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) # + NoLegend()

  # create and save total visualizations
  plot_list <- list(p3, p1)
  p <- do.call("grid.arrange", c(plot_list, ncol=2))
  ggsave(paste0("Label_transfer/cell.labels_UMAP_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))
  plot_list <- list(p4, p2)
  p <- do.call("grid.arrange", c(plot_list, ncol=2))
  ggsave(paste0("Label_transfer/cluster.labels_UMAP_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))
}
