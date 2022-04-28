library(celldex)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(scater)
library(scran)
library(Seurat)
library(pheatmap)

# SingleR Book: http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html

### get UCSC browser data as reference (train) set
## load in via their method
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Nowakowski'
setwd(work_dir)
# get matrix
mat <- data.table::fread("exprMatrix.tsv.gz")
## nrows = 0 (for dry run / maybe use for subsetting)
## verbose = T (see if it shows something useful for testing)
## drop = c(col, names, to, keep)
# showProgress = T (for testing)
# TODO try to load in data in chunks or subset before/after load to not store everything in memory

# add metadata
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# get genes
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)

### get our data as test data
file <- paste0("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds")
data <- readRDS(file = file)

# get overlapping genes between test and train dataset (to reduce noise and runtime)
genes.overlap <- rownames(data@assays$SCT)[rownames(data@assays$SCT) %in% genes]
# save non-overlap genes for retrospection
## most genes contain a: ['-', '_', '.', 'orf']
genes.non.overlap <- sort(genes[!genes %in% rownames(data@assays$SCT)])
write.csv2(genes.non.overlap, file = paste0("Label_transfer_non_overlapping_ref_genes.csv"))

# set gene column as rownames and exclude from count matrix
mat = data.frame(mat[,-1], row.names=genes)
# get SO from matrix object
so <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data=meta)
# subset train data by overlapping genes
so <- subset(so, cells = names(so$WGCNAcluster)[so$WGCNAcluster != ''], features = genes.overlap)



# Bas meeting - Nowakowski label grouping
# TODO make this based on a dataframe and custom function such that it can be
## loaded in from a csv or config file
# TODO discuss if better to build new labels into object and save that
## with a script that does that specifically for that reference dataset
### functioning as it's own config
#### only need to do it once, and in the running/plotting check if exists, else default behavior plotting
so$WGCNAcluster <- sapply(so$WGCNAcluster, function(x) {
  switch(x,
         "RG-div1" = "Radial Glia",
         "RG-div2" = "Radial Glia",
         "vRG" = "Radial Glia",
         "tRG" = "Radial Glia",
         "oRG" = "Radial Glia",
         "RG-early" = "Radial Glia (early)",
         "MGE-RG1" = "Radial Glia (MGE)",
         "MGE-RG2" = "Radial Glia (MGE)",
         "IPC-div1" = "IPC",
         "IPC-div2" = "IPC",
         "MGE-IPC1" = "IPC (MGE)",
         "MGE-IPC2" = "IPC (MGE)",
         "MGE-IPC3" = "IPC (MGE)",
         "MGE-div" = "IPC (MGE)",
         "nIN1" = "Newborn Interneuron",
         "nIN2" = "Newborn Interneuron",
         "nIN3" = "Newborn Interneuron",
         "nIN4" = "Newborn Interneuron",
         "nIN5" = "Newborn Interneuron",
         "IN-CTX-MGE1" = "Interneuron (MGE/Cortex)",
         "IN-CTX-MGE2" = "Interneuron (MGE/Cortex)",
         "IN-STR" = "Interneuron (Striatum)",
         "IN-CTX-CGE1" = "Interneuron (CGE/Cortex)",
         "IN-CTX-CGE2" = "Interneuron (CGE/Cortex)",
         "IPC-nEN1" = "IPC-nEN",
         "IPC-nEN2" = "IPC-nEN",
         "IPC-nEN3" = "IPC-nEN",
         "nEN-early1" = "nEN-early",
         "nEN-early2" = "nEN-early",
         "EN-PFC1" = "EN (Prefrontal Cortex)",
         "EN-PFC2" = "EN (Prefrontal Cortex)",
         "EN-PFC3" = "EN (Prefrontal Cortex)",
         "EN-V1-1" = "EN (Visual Cortex/V1)",
         "EN-V1-2" = "EN (Visual Cortex/V1)",
         "EN-V1-3" = "EN (Visual Cortex/V1)",
         x)
})



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

# # get a quick a idea of label distribution
table(data.sce.labels$labels)
table(data.sce.clusters$first.labels)



# source my custom color palettes from utils
source(file="../../R/my_utils/color_palettes.R")
# get custom colors
custom_colors <- my.color.palettes(type = 'mixed')
# set annotation column for transferred labels from reference data
annotation_col <- data.frame(ref.cluster = data.sce.clusters$labels)
# get ordered and unique label names from reference data
ref.colors <- unique(so.sce$WGCNAcluster)[order(unique(so.sce$WGCNAcluster))]
# set colors to each unique label from reference data
names(ref.colors) <- my.color.palettes(type = 'mixed', n = length(ref.colors))
# swap names and values of named vector to proper formatting
ref.colors <- setNames(names(ref.colors), ref.colors)
# set annotation colors for all transferred labels
annotation_colors <- list(ref.cluster = ref.colors[unique(annotation_col$ref.cluster)])
# set rownames for identification of rows during plotting
rownames(data.sce.clusters$scores) <- levels(data.sce$seurat_clusters)
# plot pretty heatmap
p <- pheatmap(t(data.sce.clusters$scores),
              fontsize = 9,
              color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PRGn"))(100),
              labels_col = paste0(levels(data.sce$seurat_clusters), " (n=", table(data.sce$seurat_clusters), ")"),
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_cols=T,
              main="Scores",
              filename="cluster_scores_heatmap_BL_A_v2.png")





plotScoreHeatmap(results = data.sce.clusters)











# # get a quick a idea of label distribution
table(data.sce.labels$labels)
table(data.sce.clusters$first.labels)
# # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
# plotScoreHeatmap(data.sce.labels)
plotScoreHeatmap(data.sce.clusters)
# # ambiguity 'predicted' by 'deltas', i.e. difference from mean scores for all labels (higher=better resolution)
# plotDeltaDistribution(data.sce.labels, ncol = 6)
# plotDeltaDistribution(data.sce.clusters, ncol = 6)
# # print warning if too much pruning happens off of default settings, could mean bad quality of label transfer
# if (table(pruneScores(data.sce.labels))[1] < length(data.sce.labels$labels)*0.99) {
#   print("WARNING: SingleR::pruneScores function showed more cells pruned than 99% margin of Gauss dist, look at quality of label transfer")
# }




# save the label transfer results in the Seurat object
data$SingleR.cell.labels <- data.sce.labels$labels
# add seurat cluster to transfer label for later identification with DE analysis
data.sce.clusters$labels <- paste(as.character(data.sce.clusters$labels), levels(data), sep='.')
data.sce.clusters$first.labels <- paste(as.character(data.sce.clusters$first.labels), levels(data), sep='.')
# original Seurat clusters are needed for RenameIdents function
names(data.sce.clusters$labels) <- levels(data)
names(data.sce.clusters$first.labels) <- levels(data)
# set idents/names of clusters to labels of transfer
data <- RenameIdents(data, data.sce.clusters$labels)
# create UMAP on cell and cluster level with transferred labels!
# DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# DimPlot(data, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.labels')








colSums(data.sce.clusters$scores)[order(colSums(data.sce.clusters$scores), decreasing=T)]


# TODO could create these heatmaps to further look into specific cell label transfer
## TODO could do the same on cluster level (?)
# all.markers <- metadata(data.sce.labels)$de.genes
# data.sce$labels <- data.sce.labels$labels
# # Beta cell-related markers
# plotHeatmap(data.sce, order_columns_by="labels",
#             features=unique(unlist(all.markers$Astrocyte)))
