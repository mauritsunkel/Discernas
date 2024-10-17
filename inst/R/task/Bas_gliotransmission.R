library(Seurat)
library(SeuratObject)


# file <- "C:/Users/mauri/Desktop/scRNAseqR/results/Bas test/BL_C.qs"
# data <- qs::qread(file)
# FeaturePlot(data, features = markers)
# DotPlot(data, features = markers)


markers <- c("SLC17A7", "SNAP25", "SYT1", "SLC1A2", "SLC1A3")
file <- "C:/Users/mauri/Desktop/scRNAseqR/results/Bas test/BL_A + BL_C.qs"
data <- qs::qread(file)
FeaturePlot(data, features = markers)
DotPlot(data, features = markers)


markers_count <- c("SLC17A7", "SNAP25", "SYT1")
length(colnames(data@assays$RNA@data[markers_count,])) # total n cells # 8048
sum(colSums(data@assays$RNA@data[markers_count,] > 0) >= length(markers_count)) # 129 any cells with all markers

co_culture_astrocytes <- data@assays$RNA@data[,data$orig.ident == "BL_C"] # total BL_C cells 2890
BL_C_gliotransmission <- names(which(colSums(co_culture_astrocytes[markers_count,] > 0) >= length(markers_count))) # 129 cells
BL_C_non_gliotransmission <- names(which(!colSums(co_culture_astrocytes[markers_count,] > 0) >= length(markers_count))) # 2761 cells

sample_markers <- Seurat::FindMarkers(
  data, assay = "SCT", cells.1 = BL_C_gliotransmission, cells.2 = BL_C_non_gliotransmission, only.pos = FALSE, verbose = T,
                                      logfc.threshold = 0, min.pct = 0, test.use = "DESeq2")


sample_markers <- Seurat::FindMarkers(
  data, assay = "SCT", ident.1 = BL_C_gliotransmission, ident.2 = BL_C_non_gliotransmission, only.pos = FALSE, verbose = T,
  logfc.threshold = 0, min.pct = 0)
library(dplyr)
sample_markers2 <- sample_markers
sample_markers <- sample_markers2
# filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
sample_markers <- sample_markers[!sample_markers$p_val_adj > 0.05,]
# order by avg_log2FC
sample_markers_pval_adj <- sample_markers %>% dplyr::arrange(dplyr::desc(avg_log2FC)) # DEPRECATED: filter(pct > 0.1)
write.csv2(sample_markers_pval_adj, file = "C:/Users/mauri/Desktop/scRNAseqR/results/Bas test/gliotransmission_DEG.csv", row.names = TRUE)
