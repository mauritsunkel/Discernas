# initialize
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(data.table)
library(future)
library(celldex)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(scater)
library(scran)
library(stringr)
# library(leiden)

beep <- function(n = 5) {
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
}

### USER PARAMETERS
sample_name <- 'BL_C'

# work dir should contain forward slashes (/) on Windows
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
dir.create(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"), recursive = T)
setwd(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"))
dir.create('Quality_Control/')
dir.create('Principal_Component_Analysis/')
dir.create('DE_analysis/')



run_label_transfer <- FALSE
run_cell_cycle_regression <- FALSE

## MULTIPROCESSING
# load future library and set plan to run certain functions with multiprocessing
plan("multisession", workers = 1) # DEVNOTE: n_workers > 1 for parallelization (for me, 5 is max, 4 is safe)
### END USER PARAMETERS

## DEPRECATED: use Seurat default for n.var.features instead of median from measurements, as this will be more reproducible and comparable to other work
# # read in amount of variable features to use to from the sample quality report Median Genes Per Cell to make pipeline more data driven
# quality_report <- read.csv(paste0("data/samples/", sample_name, "/metrics_summary.csv"))
# n_var_features <- strtoi(stringr::str_replace_all(quality_report$Median.Genes.per.Cell, ',', ''))


### end initialization



### START PIPELINE
# load/read dataset
data.data <- Read10X(data.dir = paste0("../../../data/samples/", sample_name, "/filtered_feature_bc_matrix"), strip.suffix = TRUE)
# Initialize the Seurat object with the raw (non-normalized data).
data <- CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
data <- PercentageFeatureSet(data, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
png(paste0("Quality_Control/QC_nFeat_nCount_percent.mt_", sample_name, ".png"))
plot(VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("#85d0f5", "#2b2f70")))
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(data, feature1 = "percent.mt", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
plot2 <- FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
png(paste0("Quality_Control/QC_feature-scatter_", sample_name, ".png"))
plot(plot1 + plot2)
dev.off()
# filter cells with: less then 20% mitochondrial, lower than 200 unique features (genes)
data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 20)


# # normalization
# ## global scaling method: normalizes feature expression measurements for each cell by total expression
# ## multiply by scale factor (default: 10.000) and log transforms result, stored in -> pbmc[["RNA"]]@data
# data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
# # Identification of highly variable features (feature selection)
# data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
#
# # data scaling: linear transformation prior to dimension reduction (e.g. PCA)
# ## reommended new data scaling method SCTransform: https://satijalab.org/seurat/articles/sctransform_vignette.html
# ## Shifts the expression of each gene, so that the mean expression across cells is 0
# ## scales the expression of each gene, so that the variance across cells is 1
# ### This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# ## The results of this are stored in pbmc[["RNA"]]@scale.data
# all.genes <- rownames(data)
# # omitting features = all.genes here would take only the 2000 variable genes into account, speed-up etc
# data <- ScaleData(data, features = all.genes)

# normalization + FindVariableFeatures + ScaleData (also sets default assay to SCT)
data <- SCTransform(data, vst.flavor = "v2", vars.to.regress = "percent.mt")

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data, cols = c("#85d0f5", "#2b2f70")) # EMC colors
plot2 <- LabelPoints(plot = plot1, points = head(VariableFeatures(data), 10), repel = TRUE)
png(paste0("Quality_Control/Feature-selection_variable-genes_", sample_name, ".png"))
plot(plot2)
dev.off()


if (run_cell_cycle_regression) {
  dir.create('Cell_Cycle/')

  ### TODO for pipeline, regressing out cell cycle effect should be optional
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  # First, we assign each cell a score, based on its expression of G2/M and S phase markers.
  # These marker sets should be anticorrelated in their expression levels, and cells expressing
  # neither are likely not cycling and in G1 phase.
  data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # Visualize the distribution of cell cycle markers across
  png(paste0("Cell_Cycle/Cell_cycle_markers_ridgeplot_", sample_name, ".png"))
  plot(RidgePlot(data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2))
  dev.off()
  data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot-all-features_", sample_name, ".png"))
  plot(DimPlot(data, reduction = "pca", label = TRUE)) # plot dimensions
  dev.off()
  # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
  data <- RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot-s-and-g2m-features_", sample_name, ".png"))
  plot(DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  data <- ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data))
  # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
  data <- RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot_after-regression-s-and-g2m-features_", sample_name, ".png"))
  plot(DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  data <- RunPCA(data, features = VariableFeatures(data), npcs = 20, nfeatures.print = 10, verbose = FALSE)
  png(paste0("Cell_Cycle/Cell_cycle_PCA_dimplot_after-regression-all-features_", sample_name, ".png"))
  plot(DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  # reset sample identity to data, instead of cell cycle identity
  data@active.ident <- data@meta.data$old.ident
}

# perform linear dimension reduction: PCA
data <- RunPCA(data, features = VariableFeatures(object = data), npcs = 50, verbose = FALSE)
png(paste0("Principal_Component_Analysis/PCA-scores_", sample_name, ".png"))
plot(DimPlot(data, reduction = "pca", label = TRUE)) # plot dimensions
dev.off()
png(paste0("Principal_Component_Analysis/PCA-loadings_", sample_name, ".png"))
plot(VizDimLoadings(data, dims = 1:2, reduction = "pca")) # plot loadings
dev.off()
png(paste0("Principal_Component_Analysis/PCA-genes-heatmap_", sample_name, ".png"))
plot(DimHeatmap(data, dims = 1:2, cells = 2000, balanced = TRUE, fast = FALSE)) # plot heatmap per PC
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
png(paste0("Principal_Component_Analysis/PCA-variance_", sample_name, ".png"))
p <- ggplot(data = longdf, aes(x=rep(1:length(varExplained), times=2), y = value*100, fill = variable, color = variable)) +
  geom_bar(stat="identity", width = .7) +
  # geom_point(stat="identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(fill = "Variance type", color = "Variance type") +
  ylab("Variance explained (%)") +
  xlab("Principal component (#)") +
  scale_fill_manual(values = c('darkgreen', 'darkred')) +
  scale_color_manual(values = c('darkgreen', 'darkred'))
plot(p)
dev.off()
png(paste0("Principal_Component_Analysis/PCA_elbow-plot_", sample_name, ".png"))
plot(ElbowPlot(data))
dev.off()

### DEPRECATED: Seurat Jackstraw assumes equal gene variance which is false with using SCTransform normalization (where gene variance is weighted by biological heterogenetiy)
# # determine dimensionality of the dataset by the Jackstraw procedure (if takes too long, try something else)
# data <- JackStraw(data, num.replicate = 100, dims = length(data[["pca"]]))
# data <- ScoreJackStraw(data, dims = 1:length(data[["pca"]]))
# # determine amount of PCs based on p-value 0.05 of the jackstraw based method
# jackstraw_p_value <- 0.05
# PC_p_values <- data[["pca"]]@jackstraw@overall.p.values[,'Score'] < jackstraw_p_value
# choose_n_PC <- which(PC_p_values==FALSE)[1]-1
# png("PCA_jackstraw.png")
# JackStrawPlot(data, dims = 1:length(data[["pca"]]))
# dev.off()

# cell clustering: Levine2015 - Xu & Su2015
## +data size = +data sparsity
### primary similarity measures like Euclidean distance less effective than secondary like SNN (w/ Jaccard distance)
choose_N_PCs <- 20 # default: 20 (out of default 50 generated at RunPCA), DEPRECATED: length(data[["pca"]])
data <- FindNeighbors(data, dims = 1:choose_N_PCs)
# # algorithm = 4 (= Leiden algorithm) & use: method = "igraph" (for large datasets when using Leiden algorithm)
data <- FindClusters(data, resolution = 0.5, algorithm = 4)

# run non-linear dimension reduction (T-SNE/UMAP)
data <- RunUMAP(data, reduction = "pca", dims = 1:choose_N_PCs)

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
  meta <- read.table(nowakowski_metadata_path, header=T, sep="\t", as.is=T, row.names=1)
  # get genes
  genes = mat[,1][[1]]
  genes = gsub(".+[|]", "", genes)

  # get overlapping genes between test and train dataset (to reduce noise and run-time)
  genes.overlap <- rownames(data@assays$SCT)[rownames(data@assays$SCT) %in% genes]
  # save non-overlap genes for retrospection
  ## most genes contain a: ['-', '_', '.', 'orf']
  genes.non.overlap <- sort(genes[!genes %in% rownames(data@assays$SCT)])
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
  data.sub <- subset(data, features = genes.overlap)
  # set SO to SCE object
  data.sce <- as.SingleCellExperiment(data.sub, assay = 'SCT')

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

# visualize UMAP (if SingleR was run, labels are taken into account)
png(paste0("UMAP_unsupervised_", sample_name, ".png"))
plot(DimPlot(data, reduction = "umap", label = TRUE))
dev.off()



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
data.markers <- FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.1)
write.csv2(data.markers, file = paste0("DE_analysis/marker-list_", sample_name, ".csv"))

# select top X genes per cluster for quick plots
topn <- data.markers %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = avg_log2FC) %>%
  ungroup() %>%
  pull(gene)
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                         "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                         "HOPX", "FAM107A", "AGT")
eva_hnRNPC_interest <- c("HNRNPR", "HNRNPU-AS1", "HNRNPU", "HNRNPLL", "HNRNPA3", "HNRNPD",
                         "HNRNPAB", "HNRNPH1", "HNRNPA2B1", "HNRNPH2", "HNRNPDL")
eva_hnRNPC_interest_2 <- c("HNRNPK", "HNRNPUL2", "HNRNPF", "HNRNPH3", "HNRNPA1", "HNRNPA0",
                         "HNRNPA1L2", "HNRNPC", "HNRNPM", "HNRNPL", "HNRNPUL1")
astrocyte_maturity <- c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                          "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE")
neuron_maturity <- c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                       "SYN1", "SNAP25", "SYT1", "APOE")
# juliette paper: https://www.biorxiv.org/content/10.1101/2022.03.04.482992v1
## BCOR genecard: https://www.genecards.org/cgi-bin/carddisp.pl?gene=BCOR
juliette_paper <- c("BCOR", "BCORL1", "BCL6", "MLLT3", "TFAP2A", "KDM2B", "H3K4me3", "H3K36me2")


plot_DEF <- function(data, features, name) {
  dir.create(paste0("DE_analysis/", name, "/"))
  dir.create(paste0("DE_analysis/", name, "/Feature/"))

  for (i in seq_along(features)) {
    tryCatch({
      p <- Seurat::FeaturePlot(data, features = features[i])
      ggplot2::ggsave(file=paste0("DE_analysis/", name ,"/Feature/", features[i], ".png"), width = 30, height = 20, units = "cm")
    },
    error=function(e) {
      message(features[i], ' plot is skipped, as it was not found with FetchData')
    })
  }

  p <- FeaturePlot(data, features = features)
  ggsave(file=paste0("DE_analysis/", name, "/feature-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- VlnPlot(data, features = features)
  ggsave(file = paste0("DE_analysis/", name, "/violin-plot_ ", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- DoHeatmap(data, features = features) + NoLegend()
  ggsave(file = paste0("DE_analysis/", name, "/heatmap_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  p <- RidgePlot(data, features = features, ncol = 3)
  ggsave(file = paste0("DE_analysis/", name, "/ridge-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")

  # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
  cell.num <- table(Idents(data))
  cluster.labels = paste(names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(Idents(data)) <- cluster.labels
  p <- DotPlot(data, features = features) + RotatedAxis() + WhiteBackground()
  ggsave(file = paste0("DE_analysis/", name, "/dot-plot_", name, "_", sample_name, ".png"), width = 30, height = 20, units = "cm")
  levels(Idents(data)) <- sapply(str_split(levels(Idents(data)), " "), "[[", 1)
}
plot_DEF(data = data, features = unique(topn), name = "topn-features")
plot_DEF(data = data, features = astrocyte_interest, name = "astrocyte")
plot_DEF(data = data, features = neuron_interest, name = "neuron")
plot_DEF(data = data, features = schema_psych_interest, name = "SCHEMA")
plot_DEF(data = data, features = sloan_2017_interest, name = "Sloan2017")
plot_DEF(data = data, features = eva_hnRNPC_interest, name = "Eva_hnRNPC")
plot_DEF(data = data, features = eva_hnRNPC_interest_2, name = "Eva_hnRNPC_2")
plot_DEF(data = data, features = astrocyte_maturity, name = "astrocyte_maturity")
plot_DEF(data = data, features = neuron_maturity, name = "neuron_maturity")
plot_DEF(data = data, features = juliette_paper, name = "Juliette paper BCOR")

heatmap_features <- data.markers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) %>%
  ungroup() %>%
  pull(gene)

p <- DoHeatmap(data, features = heatmap_features) + NoLegend()
ggsave(file = paste0("DEG-analysis_big-heatmap_", sample_name, ".png"), width = 30, height = 20, units = "cm")

### visually inspect/manually check selection/annotation (broad) gene panels on data
## Neuronal panels: NEUROG2, SYN1, RBFOX3 - & MAP2 <--> SYN1
## Astrocyte panel: VIM, S100B, SOX9

# save procesed data
saveRDS(data, file = paste0(sample_name, ".rds"))




# TODO remove
beep()



### NOTE: for new selection, data has been slightly adjusted with re-run, so differing amount of cells possible per cluster (shuffled) and therefore different cluster selection
### NOTE: in old selection a manual visual inspection was done over the dot/feature/ridge/violin plots, this is biased, unreproducible and unautomatable
# BL_N old selection: 2-8-10-11-12-13 (MAP2, NEUROG2, RBFOX3, NOT SYN1 as then 0 clusters would be selected)
## new selection: "nEN-early1.2"  "Glyc.6"        "nEN-early1.9"  "nEN-early1.10" "Glyc.12"       "Glyc.13"
# BL_C neuronal old selection: 1-5-7-9-11-12-14 (MAP2, NEUROG2, RBFOX3, SYN1)
## (should've also not used SYN1 here, then cluster 6 would've been included in the selection)
## new selection: "Glyc.1"  "Glyc.4"  "Glyc.6"  "Glyc.7"  "Glyc.8"  "Glyc.9"  "Glyc.11" "Glyc.13"
# BL_C astrocytical old selection: 3-4 (VIM, S100B, SOX9)
## new selection: "Astrocyte.3" "MGE-RG1.5"
# BL_A old selection: 1-6-7 (VIM, S100B, SOX9)
## new_selection: "tRG.1"      "MGE-RG2.5"  "MGE-IPC1.6" "U1.7"








### subset data manually - uncomment code
# setwd("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/")
# data <- readRDS(file = "BL_C.rds")

## `%notin%` <- Negate(`%in%`)


## CORRELATION PLOT
# x = data@assays$SCT@data['XPO7',]
# y = data@assays$SCT@data['SETD1A',]
# plot(x, y, pch = 19, col = "lightblue", xlab = 'XPO7', ylab = 'SETD1A', main = 'SCT normalised count data')
# abline(lm(y ~ x), col = "red", lwd = 3)
# text(paste("Correlation:", round(cor(x, y), 2)), x = .3, y = .2)
# table(y, x)
# table(data@assays$SCT@counts['XPO7',]) # seemingly categorised
# table(data@assays$SCT@counts['SETD1A',])
# table(data@assays$SCT@counts['XPO7',] > 0, data@assays$SCT@counts['SETD1A',] > 0)
# table(data@assays$SCT@counts['ONECUT3',])
##
