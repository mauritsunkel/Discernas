### start initialization
library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)
library(future)
library(SingleR)

beep <- function(n = 5) {
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
}

### USER PARAMETERS
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"

# load future library and set plan to run certain functions with multiprocessing
plan("multisession", workers = 1) # DEVNOTE: n_workers > 1 for parallelization (for me, 5 is max, 4 is safe)
# TODO check if paralellization holds up, otherwise maybe try to apply it per function? (also in individual)
### END USER PARAMETERS

### READ RDS METHOD
## BL_N - BL_C
# rds.files <- c("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C_neuronal_subset.rds",
#                "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_N/BL_N_neuronal_subset.rds")
# sample_name <- "BL_N + BL_C"
## BL_A + BL_C
rds.files <- c("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_A/BL_A_astrocytical_subset.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C_astrocytical_subset.rds")
sample_name <- "BL_A + BL_C"
### END READ RDS

start_time <- format(Sys.time(), "%F %H-%M-%S")
dir.create(paste0(work_dir, 'results/'))
dir.create(paste0(work_dir, 'results/', start_time, '/'))
dir.create(paste0(work_dir, 'results/', start_time, '/integrated/'))
dir.create(paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/"))
setwd(paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/"))
dir.create("DE_analysis/")
dir.create("DE_analysis/markers/")
dir.create("DE_analysis/sample_markers/")
dir.create("DE_analysis/condition_markers/")
dir.create("DE_analysis/conserved_markers/")
dir.create("Plots/")
dir.create("GSEA_analysis/")
### end initialization
### END USER PARAMETERS




data.list <- lapply(X = rds.files, FUN = function(x) {
  readRDS(file = x)
})

# select features that are repeatedly variable across datasets for integration
## nfeatures = 3000 & PrepSCTIntegration, because using SCTransform
features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)

# use Canonical COrrelation Analysis (CCA) to find 'anchors' between datasets
## with reference is faster and finds the same amount of anchors
# define combined sample index in data.list position
ind <- NULL
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
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated), npcs = 50, verbose = FALSE)

choose_N_PCs <- 20
integrated <- FindNeighbors(integrated, dims = 1:choose_N_PCs) # default 20
# # algorithm = 4 (= Leiden algorithm) & use: method = "igraph" (for large datasets when using Leiden algorithm)
integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 4)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:choose_N_PCs) # default 20
# Visualization
p <- DimPlot(integrated, reduction = "umap", group.by = 'orig.ident')
p <- p + DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
ggplot2::ggsave(file = paste0("UMAP_grouped_", sample_name, ".png"), width = 30, height = 20, units = "cm")
p <- DimPlot(integrated, reduction = "umap", label = TRUE, split.by = "orig.ident")
ggplot2::ggsave(file = paste0("UMAP_split_", sample_name, ".png"), width = 30, height = 20, units = "cm")


# run label transfer for integrated analysis after plotting regular UMAP
run_label_transfer <- TRUE
if (run_label_transfer) {
  dir.create('Label_transfer/')
  # TODO if data already loaded at individual analysis, then no need to load it here
  ### get Nowakowski as reference (train) data
  nowakowski_data_path <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Nowakowski/exprMatrix.tsv.gz'
  nowakowski_metadata_path <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Nowakowski/meta.tsv'
  nowakowski_project <- 'adultPancreas'
  # get Nowakowski matrix
  mat <- data.table::fread(nowakowski_data_path)
  # add metadata
  meta <- read.table(nowakowski_metadata_path, header=T, sep="\t", as.is=T, row.names=1)
  # get genes
  genes = mat[,1][[1]]
  genes = gsub(".+[|]", "", genes)

  # get overlapping genes between test and train dataset (to reduce noise and run-time)
  genes.overlap <- rownames(integrated@assays$SCT)[rownames(integrated@assays$SCT) %in% genes]
  # save non-overlap genes for retrospection
  ## most genes contain a: ['-', '_', '.', 'orf']
  genes.non.overlap <- sort(genes[!genes %in% rownames(integrated@assays$SCT)])
  write.csv2(genes.non.overlap, file = paste0("Label_transfer/non_overlapping_ref_genes_", sample_name, ".csv"))

  # set matrix
  mat = data.frame(mat[,-1], row.names=genes)
  # get SO from matrix object
  so <- CreateSeuratObject(counts = mat, project = nowakowski_project, meta.data=meta)
  # subset train data by overlapping genes
  so <- subset(so, cells = names(so$WGCNAcluster)[so$WGCNAcluster != ''], features = genes.overlap)
  # set SO to SCE object
  so.sce <- as.SingleCellExperiment(so)

  # subset test data by overlapping genes
  integrated.sub <- subset(integrated, features = genes.overlap)
  # set SO to SCE object
  integrated.sce <- as.SingleCellExperiment(integrated.sub, assay = 'SCT')

  # perform SingleR label transfer
  integrated.sce.labels <- SingleR(test=integrated.sce, ref=so.sce, labels=so.sce$WGCNAcluster, de.method='wilcox')
  integrated.sce.clusters <- SingleR(test=integrated.sce, ref=so.sce, labels=so.sce$WGCNAcluster,
                               clusters=integrated.sce$seurat_clusters, de.method='wilcox')
  # get a quick a idea of label distribution
  table(integrated.sce.labels$labels)
  table(integrated.sce.clusters$labels)
  # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
  png(paste0("Label_transfer/cell_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  plotScoreHeatmap(integrated.sce.labels)
  dev.off()
  png(paste0("Label_transfer/cluster_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
  plotScoreHeatmap(integrated.sce.clusters)
  dev.off()
  # ambiguity 'predicted' by 'deltas', i.e. difference from mean scores for all labels (higher=better resolution)
  # png(paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"))
  p <- plotDeltaDistribution(integrated.sce.labels, ncol = 6)
  ggsave(file=paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # png(paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"))
  p <- plotDeltaDistribution(integrated.sce.clusters, ncol = 6)
  ggsave(file=paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dev.off()
  # print warning if too much pruning happens off of default settings, could mean bad quality of label transfer
  if (table(pruneScores(integrated.sce.labels))[1] < length(integrated.sce.labels$labels)*0.99) {
    print("WARNING: SingleR::pruneScores function showed more cells pruned than 99% margin of Gauss dist, look at quality of label transfer")
  }
  # save the label transfer results in the Seurat object
  integrated$SingleR.cell.labels <- integrated.sce.labels$labels
  integrated$SingleR.cell.first.labels <- integrated.sce.labels$first.labels

  # add seurat cluster to transfer label for later identification with DE analysis
  integrated.sce.clusters$labels <- paste(as.character(integrated.sce.clusters$labels), levels(integrated), sep='.')
  integrated.sce.clusters$first.labels <- paste(as.character(integrated.sce.clusters$first.labels), levels(integrated), sep='.')

  # original Seurat clusters are needed for RenameIdents function
  names(integrated.sce.clusters$first.labels) <- levels(integrated)
  # set idents/names of clusters to labels of transfer
  integrated <- RenameIdents(integrated, integrated.sce.clusters$first.labels)
  p <- DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.first.labels', repel = TRUE)
  ggsave(file=paste0("Label_transfer/cell_first.labels_UMAP_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) # + NoLegend()
  ggsave(file=paste0("Label_transfer/cluster_first.labels_UMAP_", sample_name, ".png"), width = 30, height = 20, units = "cm")

  # original Seurat clusters are needed for RenameIdents function
  names(integrated.sce.clusters$labels) <- levels(integrated)
  # set idents/names of clusters to labels of transfer
  integrated <- RenameIdents(integrated, integrated.sce.clusters$labels)
  # create UMAP on cell and cluster level with transferred labels!
  p <- DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.labels', repel = TRUE)
  ggsave(file=paste0("Label_transfer/cell_tuned.labels_UMAP_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) # + NoLegend()
  ggsave(file=paste0("Label_transfer/cluster_tuned.labels_UMAP_", sample_name, ".png"), width = 30, height = 20, units = "cm")
}


# For performing differential expression after integration, we switch back to the original data (RNA, not SCT)
## SCT are now the normalized (and scaled?) Pearson residuals for each data set, prior to integration
DefaultAssay(integrated) <- "RNA"
# based on the test used with any of the FindMarkers or derived Seurat functions the RNA counts or normalized data will be used, which are both in different data slots
## because of using SCTransform the RNA assay data is not yet normalized, the data need not be scaled as the scale.data slot is never used for DE
### proofs by Seurat responses in Github issue numbers: 1836, 2023, 3839, 4032
integrated <- NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)


## perform visualization
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                         "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                         "HOPX", "FAM107A", "AGT") # AGXT2L1 not in our data
plot_DEF <- function(data, features, name) {
  dir.create(paste0("Plots/", name, "/"))
  dir.create(paste0("Plots/", name, "/Feature_split/"))

  for (i in seq_along(features)) {
    p <- FeaturePlot(data, features = features[i], split.by = "orig.ident", cols = c("grey", "red"))
    ggplot2::ggsave(file=paste0("Plots/", name ,"/Feature_split/", features[i], ".png"), width = 30, height = 20, units = "cm")
  }
  p <- VlnPlot(data, features = features, split.by = "orig.ident")
  ggplot2::ggsave(file = paste0("Plots/", name, "/violin-split.png"), width = 30, height = 20, units = "cm")

  p <- FeaturePlot(data, features = features)
  ggplot2::ggsave(file=paste0("Plots/", name, "/features.png"), width = 30, height = 20, units = "cm")
  p <- VlnPlot(data, features = features)
  ggplot2::ggsave(file = paste0("Plots/", name, "/violins.png"), width = 30, height = 20, units = "cm")
  p <- RidgePlot(data, features = features, ncol = 3)
  ggplot2::ggsave(file = paste0("Plots/", name, "/ridges.png"), width = 30, height = 20, units = "cm")

  # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
  cell.num <- table(data$seurat_clusters)
  cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(Idents(data)) <- cluster.labels
  p <- DotPlot(data, features = features) + RotatedAxis() + WhiteBackground()
  ggplot2::ggsave(file = paste0("Plots/", name, "/dots.png"), width = 30, height = 20, units = "cm")
  p <- DotPlot(data, features = features, split.by = "orig.ident") + RotatedAxis() + WhiteBackground()
  ggplot2::ggsave(file = paste0("Plots/", name, "/dots-split.png"), width = 30, height = 20, units = "cm")
  levels(Idents(data)) <- c(0:(length(levels(Idents(data)))-1))

  DefaultAssay(data) <- "SCT"
  p <- DoHeatmap(data, features = features) + NoLegend()
  ggplot2::ggsave(file = paste0("Plots/", name, "/heatmap.png"), width = 30, height = 20, units = "cm")
}

plot_DEF(data = integrated, features = astrocyte_interest, name = "astrocyte")
plot_DEF(data = integrated, features = neuron_interest, name = "neuron")
plot_DEF(data = integrated, features = schema_psych_interest, name = "SCHEMA")
plot_DEF(data = integrated, features = sloan_2017_interest, name = "Sloan2017")

# save integrated Seurat object
saveRDS(integrated, file = paste0("integrated_", sample_name, ".rds"))



# TODO remove
beep()





### Perform custom Differential Expression (DE) & Gene Set Enrichment Analysis (GSEA)







# TODO uncomment when done testing
# plot_DEF(data = integrated, features = unique(topn), name = "topn-features")

## TODO build in the custom differential expression (also folder/file_name.ext)
# all_markers <- FindAllMarkers(integrated, min.pct = 0.1)
# write.csv2(all_markers, file = paste0("DE-analysis_all-markers.csv"))
# # TODO use this or regular FindMarkers function for integration DEG analysis?
# conserved_markers <- FindConservedMarkers(integrated, ident.1 = 1, ident.2 = NULL,
#                                 grouping.var = "orig.ident", meta.method = metap::minimump, verbose = TRUE)
# write.csv2(conserved_markers, file = paste0("DE-analysis_conserved-markers.csv"))


## TODO fix this with building in custom DE functions
# select topn genes per cluster for quick plots
# topn <- all_markers %>%
#   group_by(cluster) %>%
#   top_n(n = 1, wt = avg_log2FC) %>%
#   ungroup() %>%
#   pull(gene)
