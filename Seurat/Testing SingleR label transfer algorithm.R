library(celldex)
library(scRNAseq)
library(SingleR)
library(scuttle)
library(scater)
library(scran)
library(Seurat)

# SingleR Book: http://bioconductor.org/books/devel/SingleRBook/using-the-classic-mode.html
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html

### get Nowakowski as reference (train) data
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Nowakowski'
setwd(work_dir)
# get matrix
mat <- data.table::fread("exprMatrix.tsv.gz")
# add metadata
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# get genes
genes = mat[,1][[1]]
genes = gsub(".+[|]", "", genes)

### get our data as test data
file <- "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/integrated/BL_A + BL_C/integrated_BL_A + BL_C.rds"
data <- readRDS(file = file)

# get overlapping genes between test and train dataset (to reduce noise and runtime)
genes.overlap <- rownames(data@assays$SCT)[rownames(data@assays$SCT) %in% genes]
# save non-overlap genes for retrospection
## most genes contain a: ['-', '_', '.', 'orf']
genes.non.overlap <- sort(genes[!genes %in% rownames(data@assays$SCT)])
write.csv2(genes.non.overlap, file = paste0("Label_transfer_non_overlapping_ref_genes.csv"))

# set matrix
mat = data.frame(mat[,-1], row.names=genes)
# get SO from matrix object
so <- CreateSeuratObject(counts = mat, project = "adultPancreas", meta.data=meta)
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


levels(data)
Idents(data)



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
