library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)

## TODO test if these libraries are needed here based on a test set
library(patchwork)
library(data.table)
library(celldex)
library(scRNAseq)
library(scuttle)
library(scater)
library(scran)



### INITIALIZATION
## USER PARAMETERS
sample_name <- 'BL_C'

# set to run cell cycle regression (based on: https://satijalab.org/seurat/articles/cell_cycle_vignette.html)
run_cell_cycle_regression <- FALSE
## END USER PARAMETERS

# create output directories based on start time and sample name and set working directory
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
dir.create(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"), recursive = T)
setwd(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"))
dir.create('Quality_Control/')
dir.create('Principal_Component_Analysis/')
dir.create('DE_analysis/')
### END INITIALIZATION



### START PIPELINE
# read 10X data (preprocessed by 10X Cellranger pipeline) and convert to Seurat object
data.data <- Seurat::Read10X(data.dir = paste0("../../../data/samples/", sample_name, "/filtered_feature_bc_matrix"), strip.suffix = TRUE)
data <- Seurat::CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)

## QUALITY CONTROL
# calculate percentage of all counts belonging to mitochondrial (^MT-) DNA, for filtering
data <- Seurat::PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
# Visualize quality control metrics
png(paste0("Quality_Control/QC_nFeat_nCount_percent.mt_", sample_name, ".png"))
plot(Seurat::VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("#85d0f5", "#2b2f70")))
dev.off()
plot1 <- Seurat::FeatureScatter(data, feature1 = "percent.mt", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
plot2 <- Seurat::FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
png(paste0("Quality_Control/QC_feature-scatter_", sample_name, ".png"))
plot(plot1 + plot2)
dev.off()
# subset cells with less than 200 unique genes (NFeature_RNA)
## TODO check if base::subset or SeuratObject::subset
data <- subset(data, subset = nFeature_RNA > 200)

## SCTransform-v2 replaces NormalizeData + FindVariableFeatures + ScaleData & sets default assay to "SCT"
## https://satijalab.org/seurat/articles/sctransform_vignette.html & https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
# normalize gene expression counts per cell by the total expression and applying a
# scaling factor (default: 10.000) and adding a pseudocount before log-transforming the result
# this global linear scaling on the data sets mean expression across cells is 0 and variance across cells is 1 as to
# provide equal weight in downstream processing such that highly variable genes do not dominate results
# SCTransform-v2 excludes the need for this heuristic pseudocount addition, log-transformation and optimizes variation
# the top 3000 (default) variable genes are kept for improving downstream processing efficiency
# vars.to.regress = regress out variability originating from reads mapped to mitochondrial DNA
# return.only.var.genes = TRUE, as non-sparse matrix is returned and used in PCA
# set transformed data as default data assay for downstream processing
data <- Seurat::SCTransform(data, vst.flavor = "v2", vars.to.regress = "percent.mt", return.only.var.genes = TRUE)

# plot variable features, label top 10
plot1 <- Seurat::VariableFeaturePlot(data, cols = c("#85d0f5", "#2b2f70"))
plot2 <- Seurat::LabelPoints(plot = plot1, points = head(SeuratObject::VariableFeatures(data), 10), repel = TRUE)
png(paste0("Quality_Control/Feature-selection_variable-genes_", sample_name, ".png"))
plot(plot2)
dev.off()

if (run_cell_cycle_regression) {
  dir.create('Cell_Cycle/')

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
  # We can segregate this list into markers of G2/M phase and markers of S phase.
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  # First, we assign each cell a score, based on its expression of G2/M and S phase markers.
  # These marker sets should be anticorrelated in their expression levels,
  # cells expressing neither are likely not cycling and in G1 phase.
  data <- Seurat::CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # Visualize the distribution of cell cycle markers across
  png(paste0("Cell_Cycle/Cell_cycle_markers_ridgeplot_", sample_name, ".png"))
  plot(Seurat::RidgePlot(data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2))
  dev.off()
  data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(object = data), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot-all-features_", sample_name, ".png"))
  plot(DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
  data <- Seurat::RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot-s-and-g2m-features_", sample_name, ".png"))
  plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  # When running a PCA on only cell cycle genes after regression, cells no longer separate by cell-cycle phase
  data <- Seurat::ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data))
  data <- Seurat::RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot_after-regression-s-and-g2m-features_", sample_name, ".png"))
  plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(data), npcs = 20, nfeatures.print = 10, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot_after-regression-variable-features_", sample_name, ".png"))
  plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  # set sample identity to data, instead of cell cycle identity
  data@active.ident <- data@meta.data$old.ident
}









# run Principal Component Analysis as linear dimension reduction
data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(object = data), npcs = 50, verbose = FALSE)
png(paste0("Principal_Component_Analysis/PCA-scores_", sample_name, ".png"))
plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
dev.off()
png(paste0("Principal_Component_Analysis/PCA-loadings_", sample_name, ".png"))
plot(Seurat::VizDimLoadings(data, dims = 1:2, reduction = "pca"))
dev.off()
png(paste0("Principal_Component_Analysis/PCA-genes-heatmap_", sample_name, ".png"))
plot(Seurat::DimHeatmap(data, dims = 1:2, cells = 2000, balanced = TRUE, fast = FALSE))
dev.off()
# custom Elbow (or Scree) plot -> Variance explained
varExplained <- (data[["pca"]]@stdev)^2 / data[["pca"]]@misc$total.variance # Eigenvalues (current subset) / total_variance (Whole dataset)
plotdf <- data.frame('Cumulative' = round(cumsum(varExplained / sum(varExplained)), 3),
                     'Individual' = varExplained / sum(varExplained))
# for geom_bar stacking effect with using stat="identity"
plotdf$diff <- plotdf$Cumulative - plotdf$Individual
plotdf$Cumulative <- NULL
plotdf <- plotdf[, c(2, 1)]
colnames(plotdf) <- c('Cumulative', 'Individual')
longdf <- reshape2::melt(plotdf)
png(paste0("Principal_Component_Analysis/PCA-variance_", sample_name, ".png"))
p <- ggplot2::ggplot(data = longdf, ggplot2::aes(x=rep(1:length(varExplained), times=2), y = value*100, fill = variable, color = variable)) +
  ggplot2::geom_bar(stat="identity", width = .7) +
  # ggplot2::geom_point(stat="identity") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggplot2::labs(fill = "Variance type", color = "Variance type") +
  ggplot2::ylab("Variance explained (%)") +
  ggplot2::xlab("Principal component (#)") +
  ggplot2::scale_fill_manual(values = c('darkgreen', 'darkred')) +
  ggplot2::scale_color_manual(values = c('darkgreen', 'darkred'))
plot(p)
dev.off()
png(paste0("Principal_Component_Analysis/PCA_elbow-plot_", sample_name, ".png"))
plot(Seurat::ElbowPlot(data))
dev.off()

# cell clustering: Levine2015 - Xu & Su2015
choose_N_PCs <- 20 # default: 20 (out of default 50 generated with Seurat::RunPCA)
# construct nearest neighbor graph for clustering
data <- Seurat::FindNeighbors(data, dims = 1:choose_N_PCs)
# use Leiden algorithm for clustering (https://www.nature.com/articles/s41598-019-41695-z/)
## method = "igraph" (for large datasets when using Leiden algorithm)
data <- Seurat::FindClusters(data, resolution = 0.5, algorithm = 4)
# visualize clustering with Uniform Manifold Projection Approximation (UMAP) as non-linear dimension reduction
data <- Seurat::RunUMAP(data, reduction = "pca", dims = 1:choose_N_PCs)
png(paste0("UMAP_unsupervised_", sample_name, ".png"))
plot(Seurat::DimPlot(data, reduction = "umap", label = TRUE))
dev.off()





# FindAllMarkers for indiviual sample --> later DEG for integrated analysis (then explain test.use: Wilcox rst)

# run FindAllMarkers for Differential Gene Expression analysis
## finds markers for every cluster compared to all remaining cells
### report only up-regulated genes as down-regulated genes represent all other cells/clusters here
data.markers <- Seurat::FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.1)
utils::write.csv2(data.markers, file = paste0("DE_analysis/marker-list_", sample_name, ".csv"))

# select top gene per cluster for exploration
topn <- data.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 1, wt = avg_log2FC) %>%
  dplyr::ungroup() %>%
  dplyr::pull(gene)
# marker panels of interest
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE")
astrocyte_maturity <- c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                        "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
neuron_maturity <- c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                     "SYN1", "SNAP25", "SYT1", "APOE")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                         "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                         "HOPX", "FAM107A", "AGT")

plot_DEG <- function(data, features, name) {
  dir.create(paste0("DE_analysis/", name, "/"))
  dir.create(paste0("DE_analysis/", name, "/Feature/"))

  # plot feature expression, if available in Seurat
  for (i in seq_along(features)) {
    tryCatch({
      p <- Seurat::FeaturePlot(data, features = features[i])
      ggplot2::ggsave(file=paste0("DE_analysis/", name ,"/Feature/", features[i], ".png"), width = 30, height = 20, units = "cm")
    },
    error=function(e) {
      message(features[i], ' plot is skipped, as gene was not found with FetchData')
    })
  }

  # expression plots
  p <- Seurat::FeaturePlot(data, features = features)
  ggplot2::ggsave(file=paste0("DE_analysis/", name, "/feature-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- Seurat::VlnPlot(data, features = features)
  ggplot2::ggsave(file = paste0("DE_analysis/", name, "/violin-plot_ ", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- Seurat::DoHeatmap(data, features = features) + NoLegend()
  ggplot2::ggsave(file = paste0("DE_analysis/", name, "/heatmap_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- Seurat::RidgePlot(data, features = features, ncol = 3)
  ggplot2::ggsave(file = paste0("DE_analysis/", name, "/ridge-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  # dotplot with custom labels
  cell.num <- table(SeuratObject::Idents(data))
  cluster.labels = paste(names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(SeuratObject::Idents(data)) <- cluster.labels
  p <- Seurat::DotPlot(data, features = features) + Seurat::RotatedAxis() + Seurat::WhiteBackground()
  ggplot2::ggsave(file = paste0("DE_analysis/", name, "/dot-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  levels(SeuratObject::Idents(data)) <- sapply(stringr::str_split(levels(SeuratObject::Idents(data)), " "), "[[", 1)
}
plot_DEG(data = data, features = unique(topn), name = "topn-features")
plot_DEG(data = data, features = astrocyte_interest, name = "astrocyte")
plot_DEG(data = data, features = astrocyte_maturity, name = "astrocyte_maturity")
plot_DEG(data = data, features = neuron_interest, name = "neuron")
plot_DEG(data = data, features = neuron_maturity, name = "neuron_maturity")
plot_DEG(data = data, features = schema_psych_interest, name = "SCHEMA")
plot_DEG(data = data, features = sloan_2017_interest, name = "Sloan2017")

# plot heatmap for topn genes per cluster
heatmap_features <- data.markers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 8, wt = avg_log2FC) %>%
  dplyr::ungroup() %>%
  dplyr::pull(gene)
p <- Seurat::DoHeatmap(data, features = heatmap_features) + Seurat::NoLegend()
ggplot2::ggsave(file = paste0("DEG-analysis_big-heatmap_", sample_name, ".png"), width = 30, height = 20, units = "cm")

# save data
saveRDS(data, file = paste0(sample_name, ".rds"))
