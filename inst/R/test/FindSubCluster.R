data <- qs::qread("C:/SynologyDrive/Projects/scRNAseqR/results/SB_V1/SB/SB.qs")
selection_markers = c("MAP2", "DCX", "NEUROG2")
# selection_markers = c("VIM", "S100B", "SOX9")
# selection_markers = c("AIF1", "CSF1R", "SPI1")

Seurat::Idents(data) <- data$seurat_clusters
data <- Seurat::FindSubCluster(
  object = data,
  cluster = levels(data$seurat_clusters),
  graph.name = "SCT_nn",  # Graphs(data)
  subcluster.name = "seurat_subclusters"
)
Seurat::Idents(data) <- data$seurat_subclusters
# Seurat::Idents(data) <- data$seurat_clusters
p <- Seurat::DotPlot(data, features = selection_markers)
# png("test.png")
hist(p$data$pct.exp, breaks = seq(0, 100, 5), main = paste0("Cells percentage expressed: ", paste(selection_markers, collapse = ", ")))
# dev.off()

# 10/15/20 astrocytes
# 10 neurons
# microglia 10 > 15 > 20

percent_expressed <- 10
cluster_selection <- names(which(table(p$data[p$data$pct.exp > percent_expressed,]$id) == length(unique(p$data$features.plot))))
cell_selection <- data$seurat_subclusters %in% cluster_selection
# cell_selection <- data$seurat_clusters %in% cluster_selection
Seurat::Idents(data) <- cell_selection
Seurat::DefaultAssay(data) <- "RNA"
Seurat::DimPlot(data, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cell_selection)["TRUE"]))

# ---- CELL LEVEL SELECTION
# layer_data <- SeuratObject::LayerData(data)
# select each cell that has expresses each gene from selection_markers
# cell_selection <- sapply(as.data.frame(layer_data[selection_markers, ] > 0), sum) == length(rownames(layer_data[selection_markers, ]))
# Seurat::Idents(data) <- cell_selection
# Seurat::DimPlot(data, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cell_selection)["TRUE"]))

