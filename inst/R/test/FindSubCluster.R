data <- readRDS("C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV8/NSM-NS-NC-M/NSM-NS-NC-M.rds")
Seurat::Idents(data) <- data$seurat_clusters
data <- Seurat::FindSubCluster(
  object = data,
  cluster = levels(data$seurat_clusters),
  graph.name = "SCT_nn",  # Graphs(data)
  subcluster.name = "seurat_subclusters"
)

# selection_markers = c("MAP2", "DCX", "NEUROG2") # 30
# selection_markers = c("VIM", "S100B", "SOX9") # 30
# selection_markers = c("AIF1", "CSF1R", "CX3CR1") # old: 20
selection_markers = c("AIF1", "CSF1R", "SPI1") # 30
percent_expressed <- 30

Seurat::Idents(data) <- data$seurat_subclusters
# Seurat::Idents(data) <- data$seurat_clusters
p <- Seurat::DotPlot(data, features = selection_markers)
# hist(p$data$pct.exp, breaks = seq(0, 100, 5))
png("test.png")
hist(p$data$pct.exp, breaks = seq(0, 100, 5), main = paste0("Cells percentage expressed: ", paste(selection_markers, collapse = " ")))
dev.off()

cluster_selection <- names(which(table(p$data[p$data$pct.exp > percent_expressed,]$id) == length(unique(p$data$features.plot))))
cell_selection <- data$seurat_subclusters %in% cluster_selection
# cell_selection <- data$seurat_clusters %in% cluster_selection
Seurat::Idents(data) <- cell_selection
Seurat::DefaultAssay(data) <- "RNA"
Seurat::DimPlot(data, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cell_selection)["TRUE"]))

# ---- CELL LEVEL SELECTION
layer_data <- SeuratObject::LayerData(data)
# select each cell that has expresses each gene from selection_markers
cell_selection <- sapply(as.data.frame(layer_data[selection_markers, ] > 0), sum) == length(rownames(layer_data[selection_markers, ]))
Seurat::Idents(data) <- cell_selection
Seurat::DimPlot(data, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cell_selection)["TRUE"]))
