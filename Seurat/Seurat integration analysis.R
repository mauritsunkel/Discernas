### start initialization
library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)

### USER PARAMETERS
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
### END USER PARAMETERS

### READ RDS METHOD
## BL_N - BL_C
rds.files <- c("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden/BL_N/neuronal-subset.rds",
               "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden/BL_C/mixed-neuronal-subset.rds")
sample_name <- "BL_N + BL_C"
## BL_A + BL_C
rds.files <- c("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden/BL_A/astrocytical-subset.rds",
               "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden/BL_C/mixed-astrocytical-subset.rds")
sample_name <- "BL_A + BL_C"
### END READ RDS







work_dir <- paste0(work_dir, 'results/')
dir.create(work_dir)
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- paste0(work_dir, start_time, '/')
dir.create(work_dir)
work_dir <- paste0(work_dir, 'integrated/')
dir.create(work_dir)
work_dir <- paste0(work_dir, sample_name)
dir.create(work_dir)

setwd(work_dir)
### end initialization



###### NOTES
# need start from: raw, latest and user list methods - if latest not available then perform raw for that sample!





data.list <- lapply(X = rds.files, FUN = function(x) {
  readRDS(file = x)
})
### END






### LATEST TIMEPOINT WITH RDS FILE EXISTING METHOD
# sample_name <- "BL_N"
# datetimes <- list.dirs(paste0("../../", sample_name), full.names = FALSE, recursive = FALSE)
# max(datetimes)
### END





# select features that are repeatedly variable across datasets for integration
## nfeatures = 3000 & PrepSCTIntegration, because using SCTransform
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)


# use Canonical COrrelation Analysis (CCA) to find 'anchors' between datasets
## with reference is faster and finds the same amount of anchors
# define combined sample index in data.list position
c <- 1
for (obj in data.list) {
  if (levels(obj$orig.ident) == "BL_C") {
    ind <- c
  }
  c <- c + 1
}
# use combined sample as reference dataset
anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, reference = c(ind))

# this command creates an 'integrated' data assay
## note that the original unmodified data still resides in the 'RNA' assay
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

## start integrated analysis

# specify that we will perform downstream analysis on the corrected data
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
## DEPRECATED: DO NOT run ScaleData here when using SCTransform
# integrated <- ScaleData(integrated, features = rownames(integrated), verbose = FALSE)
# integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated), npcs = 50, verbose = FALSE)

# ## DEPRECATED: we use SCTransform now
# # determine dimensionality of the dataset by the Jackstraw procedure (if takees too long, try something else)
# integrated <- JackStraw(integrated, num.replicate = 100, dims = length(integrated[["pca"]]))
# integrated <- ScoreJackStraw(integrated, dims = 1:length(integrated[["pca"]]))
# # determine amount of PCs based on p-value 0.05 of the jackstraw based method
# jackstraw_p_value <- 0.05
# PC_p_values <- integrated[["pca"]]@jackstraw@overall.p.values[,'Score'] < jackstraw_p_value
# choose_n_PC <- which(PC_p_values==FALSE)[1]-1


choose_N_PCs <- 20
integrated <- FindNeighbors(integrated, dims = 1:choose_N_PCs) # default 20, DEPRECATED length(integrated[["pca"]])
# # algorithm = 4 (= Leiden algorithm) & use: method = "igraph" (for large datasets when using Leiden algorithm)
integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 4)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:choose_N_PCs) # default 20, DEPRECATED length(integrated[["pca"]])
# Visualization
p <- DimPlot(integrated, reduction = "umap", group.by = 'orig.ident')
p <- p + DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
ggplot2::ggsave(file = paste0("UMAP_", sample_name, "_grouped.png"), width = 30, height = 20, units = "cm")
p <- DimPlot(integrated, reduction = "umap", label = TRUE, split.by = "orig.ident")
ggplot2::ggsave(file = paste0("UMAP_", sample_name, "_split.png"), width = 30, height = 20, units = "cm")


# For performing differential expression after integration, we switch back to the original data (RNA, not SCT)
## SCT are now the normalized (and scaled?) Pearson residuals for each data set, prior to integration
DefaultAssay(integrated) <- "RNA"
# based on the test used with any of the FindMarkers or derived Seurat functions the RNA counts or normalized data will be used, which are both in different data slots
## because of using SCTransform the RNA assay data is not yet normalized, the data need not be scaled as the scale.data slot is never used for DE
### proofs by Seurat responses in Github issue numbers: 1836, 2023, 3839, 4032
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)


all_markers <- FindAllMarkers(integrated, min.pct = 0.1)
write.csv2(all_markers, file = paste0("DEG-analysis_all-markers.csv"))
# TODO use this or regular FindMarkers function for integration DEG analysis?
conserved_markers <- FindConservedMarkers(integrated, ident.1 = 1, ident.2 = NULL,
                                grouping.var = "orig.ident", meta.method = metap::minimump, verbose = TRUE)
write.csv2(conserved_markers, file = paste0("DEG-analysis_conserved-markers.csv"))

# select topn genes per cluster for quick plots
topn <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  ungroup() %>%
  pull(gene)
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
sloan_2017_interest <- c("AQP4", "ALDH1L1", "AGXT2L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                         "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                         "HOPX", "FAM107A", "AGT")

plot_DEF <- function(data, features, name) {
  dir.create(paste0(work_dir, "/Feature_split_plots_", name))
  for (i in seq_along(features)) {
    p <- FeaturePlot(data, features = features[i], split.by = "orig.ident", cols = c("grey", "red"))
    ggplot2::ggsave(file=paste0("Feature_split_plots_", name ,"/DEG-analysis_", name, "_feature-plot-split_", features[i], ".png"), width = 30, height = 20, units = "cm")
  }
  p <- VlnPlot(data, features = features, split.by = "orig.ident")
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_violin-plot-split.png"), width = 30, height = 20, units = "cm")

  p <- FeaturePlot(data, features = features)
  ggplot2::ggsave(file=paste0("DEG-analysis_", name, "_feature-plot.png"), width = 30, height = 20, units = "cm")
  p <- VlnPlot(data, features = features)
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_violin-plot.png"), width = 30, height = 20, units = "cm")
  p <- RidgePlot(data, features = features, ncol = 3)
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_ridge-plot.png"), width = 30, height = 20, units = "cm")

  # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
  cell.num <- table(data$seurat_clusters)
  cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(Idents(data)) <- cluster.labels
  p <- DotPlot(data, features = features) + RotatedAxis() + WhiteBackground()
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_dot-plot.png"), width = 30, height = 20, units = "cm")
  p <- DotPlot(data, features = features, split.by = "orig.ident") + RotatedAxis() + WhiteBackground()
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_dot-plot-split.png"), width = 30, height = 20, units = "cm")
  levels(Idents(data)) <- c(0:(length(levels(Idents(data)))-1))

  DefaultAssay(data) <- "SCT"
  p <- DoHeatmap(data, features = features) + NoLegend()
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_heatmap.png"), width = 30, height = 20, units = "cm")
}
plot_DEF(data = integrated, features = unique(topn), name = "topn-features")
plot_DEF(data = integrated, features = astrocyte_interest, name = "astrocyte")
plot_DEF(data = integrated, features = neuron_interest, name = "neuron")
plot_DEF(data = integrated, features = schema_psych_interest, name = "SCHEMA")
plot_DEF(data = integrated, features = sloan_2017_interest, name = "Sloan2017")

# save integrated Seurat object
saveRDS(integrated, file = "integrated.rds")
