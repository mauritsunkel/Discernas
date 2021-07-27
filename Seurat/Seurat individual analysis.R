
# initialize
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(data.table)
# library(leiden)

### USER PARAMETERS
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
sample_name <- 'BL_C'
### END USER PARAMETERS



work_dir <- paste0(work_dir, 'results/')
dir.create(work_dir)
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- paste0(work_dir, start_time, '/')
dir.create(work_dir)
work_dir <- paste0(work_dir, sample_name, '/')
dir.create(work_dir)
setwd(work_dir)

### end initialization



# load/read dataset
data.data <- Read10X(data.dir = paste0("../../../data/samples/", sample_name, "/filtered_feature_bc_matrix"), strip.suffix = TRUE)

# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
# head(data@meta.data, 5)

# Visualize QC metrics as a violin plot
png("QC_nFeat_nCount_percent.mt.png")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(data, feature1 = "percent.mt", feature2 = "nCount_RNA")
plot2 <- FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA")
png("QC_feature-scatter.png")
plot1 + plot2
dev.off()

# filter cells with: less then 20% mitochondrial, lower than 200 unique features (genes)
data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 20)

# normalization
## global scaling method: normalizes feature expression measurements for each cell by total expression
## multiply by scale factor (default: 10.000) and log transforms result, stored in -> pbmc[["RNA"]]@data
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data[["RNA"]]@data

# Identification of highly variable features (feature selection)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data), 10)




# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png("Feature-selection_variable-genes.png")
plot2
dev.off()

# data scaling: linear transformation prior to dimension reduction (e.g. PCA)
## reommended new data scaling method SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
## Shifts the expression of each gene, so that the mean expression across cells is 0
## scales the expression of each gene, so that the variance across cells is 1
### This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
## The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(data)
# omitting features = all.genes here would take only the 2000 variable genes into account, speed-up etc
data <- ScaleData(data, features = all.genes)

# perform linear dimension reduction: PCA
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 20, verbose = FALSE)
png(sample_name, "PCA-scores.png")
DimPlot(data, reduction = "pca", label = TRUE) # plot dimensions
dev.off()
png("PCA-loadings.png")
VizDimLoadings(data, dims = 1:2, reduction = "pca") # plot loadings
dev.off()
png("PCA-genes-heatmap.png")
DimHeatmap(data, dims = 1:2, cells = 2000, balanced = TRUE, fast = FALSE) # plot heatmap per PC
dev.off()
# custom Elbow (Scree) plot -> Variance explained
varExplained <- (data[["pca"]]@stdev)^2 / data[["pca"]]@misc$total.variance # Eigenvalues(current subset) / total_variance (Whole dataset)
plotdf <- data.frame('Cumulative' = round(cumsum(varExplained / sum(varExplained)), 3),
                     'Individual' = varExplained / sum(varExplained))
## for geom_bar stacking effect with using stat="identity"
plotdf$diff <- plotdf$Cumulative - plotdf$Individual
plotdf$Cumulative <- NULL
plotdf <- plotdf[, c(2, 1)]
colnames(plotdf) <- c('Cumulative', 'Individual')
##
longdf <- reshape2::melt(plotdf)
png("PCA-variance.png")
p <- ggplot(data = longdf, aes(x=rep(1:length(varExplained), times=2), y = value*100, fill = variable, color = variable)) +
  geom_bar(stat="identity", width = .7) +
  # geom_point(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Variance type", color = "Variance type") +
  ylab("Variance explained (%)") +
  xlab("Principal component (#)") +
  scale_fill_manual(values = c('darkgreen', 'darkred')) +
  scale_color_manual(values = c('darkgreen', 'darkred'))
p
dev.off()
png("PCA_elbow-plot.png")
ElbowPlot(data)
dev.off()

# determine dimensionality of the dataset by the Jackstraw procedure (if takees too long, try something else)
data <- JackStraw(data, num.replicate = 100, dims = length(data[["pca"]]))
data <- ScoreJackStraw(data, dims = 1:length(data[["pca"]]))
# determine amount of PCs based on p-value 0.05 of the jackstraw based method
jackstraw_p_value <- 0.05
PC_p_values <- data[["pca"]]@jackstraw@overall.p.values[,'Score'] < jackstraw_p_value
choose_n_PC <- which(PC_p_values==FALSE)[1]-1
png("PCA_jackstraw.png")
JackStrawPlot(data, dims = 1:length(data[["pca"]]))
dev.off()

# cell clustering: Levine2015 - Xu & Su2015
## +data size = +data sparsity
### primary similarity measures like Euclidean distance less effective than secondary like SNN (w/ Jaccard distance)
data <- FindNeighbors(data, dims = 1:length(data[["pca"]]))
data <- FindClusters(data, resolution = 0.5)

### TRY TO IMPLEMENT LEIDEN ALGORITHM ###
# data <- FindClusters(data, resolution = 0.5, algorithm = "leiden")
## implement Leiden algorithm for clustering instead of standard Louvain or LSM!
# membership <- leiden(pbmc@graphs$RNA_snn)
# table(membership)
# https://cran.r-project.org/web/packages/leiden/vignettes/run_leiden.html
# how to install leiden, leidenalg (pip), igrpah (pip) and reticulate (devtools, via leiden git)
# check https://rstudio.github.io/reticulate/
# redo steps and figure it out
## would this work instead of pip install: reticulate::py_install(packages = "igraph" etc)

# run non-linear dimension reduction (T-SNE/UMAP)
data <- RunUMAP(data, reduction = "pca", dims = 1:length(data[["pca"]]))
png("UMAP_dimensionality-reduction.png")
DimPlot(data, reduction = "umap", label = TRUE)
dev.off()


##  save object without having to run above steps again
# saveRDS(data, file = paste0(sample_name, ".rds"))







## tryout loadRDS and then follow diff gene exp analysis and gene annotaton
# sample_name <- "BL_C"
# work_dir <- "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/BL_C/2021-06-03 16-46-41"
# setwd(work_dir)
# data <- readRDS(file = paste0(sample_name, ".rds"))


# differentially expressed features (genes) analysis
## with the FindAllMarkers function, for each cluster against all other clusters (cells) the DEG will be calculated
### this is useful for finding markers globally to reference to our genes list of interest and then to perform cluster annotation
## with the FindMarkers functions using ident.1 and ident.2 the cluster #/name can be specified solo or in vector c() to define what to compare against each other
### this will be most useful after annotation
## parameters
### min.pct: in either group, minimum percentage of existence of the feature within the cells of either group
#### default: 10%/0.1 - debatable
### min.diff.pct: in either group, how much min.pct is allowed to differ
#### default: -Inf, i.e. no difference - discuss, why even consider changing this parameter, makes no sense to me
### logfc.threshold: minimally X-fold difference between the groups
#### default: 0.25 - higher = faster, debatable
### test.use: see paper https://www.nature.com/articles/nmeth.4612
## proper (log) fold-change explanation: https://www.biostars.org/p/342756/


# find markers for every cluster compared to all remaining cells
## report only the positive ones because otherwise unclear in what specific cluster the gene is more expressed!
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.1)
write.csv2(data.markers, file = paste0("DEG-analysis_genes-list.csv"))


# to select genes from DE list
# selected_genes <- data.markers %>%
#   filter(gene %in% c("NEUROG2"))


# select top X genes per cluster for quick plots
topn <- data.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  ungroup() %>%
  pull(gene)
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")


plot_DEF <- function(data, features, name) {
  p <- FeaturePlot(data, features = features)
  ggsave(file=paste0("DEG-analysis_", name, "_feature-plot.png"), width = 30, height = 20, units = "cm")
  p <- VlnPlot(data, features = features)
  ggsave(file = paste0("DEG-analysis_", name, "_violin-plot.png"), width = 30, height = 20, units = "cm")
  p <- DoHeatmap(data, features = features) + NoLegend()
  ggsave(file = paste0("DEG-analysis_", name, "_heatmap.png"), width = 30, height = 20, units = "cm")
  p <- RidgePlot(data, features = features, ncol = 3)
  ggsave(file = paste0("DEG-analysis_", name, "_ridge-plot.png"), width = 30, height = 20, units = "cm")

  # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
  cell.num <- table(data$seurat_clusters)
  cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(Idents(data)) <- cluster.labels
  p <- DotPlot(data, features = features) + RotatedAxis() + WhiteBackground()
  ggsave(file = paste0("DEG-analysis_", name, "_dot-plot.png"), width = 30, height = 20, units = "cm")
  levels(Idents(data)) <- c(0:(length(levels(Idents(data)))-1))
}
plot_DEF(data = data, features = unique(topn), name = "top-features")
plot_DEF(data = data, features = astrocyte_interest, name = "astrocyte")
plot_DEF(data = data, features = neuron_interest, name = "neuron")
plot_DEF(data = data, features = schema_psych_interest, name = "SCHEMA")

heatmap_features <- data.markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) %>%
  ungroup() %>%
  pull(gene)

p <- DoHeatmap(data, features = heatmap_features) + NoLegend()
ggsave(file = "DEG-analysis_big-heatmap.png", width = 30, height = 20, units = "cm")








### visually inspect/manually check selection/annotation (broad) gene panels on data
## Neuronal panels: NEUROG2, SYN1, RBFOX3 - & MAP2 <--> SYN1
## Astrocyte panel: VIM, S100B, SOX9

# save preselection data
saveRDS(data, file = "mixed.rds") # neuronal, astrocytical or mixed
# subset data
# `%notin%` <- Negate(`%in%`)
subset <- subset(data, seurat_clusters %in% c(0, 4, 5, 6, 7, 8, 10, 12)) # or use idents = c() instead of seurat_clusters
saveRDS(subset, file = "mixed-neuronal-subset.rds") # neuronal, astrocytical, mixed-neuronal, mixed-astrocytical








# read the paper on which test method to use for DEG analysis (which paper is that again?)

# check Slingshot R package for RNA velocity, maybe see development arrows on UMAP
# Pascal papers
## https://bioconductor.org/packages/devel/bioc/vignettes/slingshot/inst/doc/vignette.html
##  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6582955/
## Korsunsky2019 Fast, sensitive and accurate integration of single-cell data with Harmony - on desktop


# READ THIS https://www.nature.com/articles/nmeth.4612
# AND THIS https://www.embopress.org/doi/full/10.15252/msb.20188746
