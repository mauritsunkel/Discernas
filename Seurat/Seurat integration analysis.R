### start initialization
library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)

### USER PARAMETERS
# combined, neuronal and astrocyte samples
samples.list <- c("BL_C", "BL_N", "BL_A")
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
### END USER PARAMETERS

### READ RDS METHOD
## BL_N - BL_C
rds.files <- c("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/BL_N/neuronal-subset.rds",
               "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/BL_C/mixed-neuronal-subset.rds")

sample_name <- "BL_N + BL_C"
## BL_A + BL_C
rds.files <- c("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/BL_A/astrocytical-subset.rds",
               "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/BL_C/mixed-astrocytical-subset.rds")
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
sample_name <- "BL_N"
datetimes <- list.dirs(paste0("../../", sample_name), full.names = FALSE, recursive = FALSE)
max(datetimes)
### END


data.list <- lapply(X = samples.list, FUN = function(x) {

  ### BETTER TO MAKE LIST OF .rds files to read in after per sample analysis
  ## make data.list of those
  ## perform steps at Perform integration

  ### COMPARE this and LoadRDS method to see if outcome is similar

  # load/read dataset
  data.data <- Read10X(data.dir = paste0("../../../data/samples/", x, "/filtered_feature_bc_matrix"), strip.suffix = TRUE)
  # Initialize the Seurat object with the raw (non-normalized data).
  data <- CreateSeuratObject(counts = data.data, project = x, min.cells = 3, min.features = 700)
  # calculate percentage mitochondrial DNA
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  # filter cells with: less then 20% mitochondrial, lower than 200 unique features (genes)
  data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 20)

  # normalization
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

  # Identification of highly variable features (feature selection)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
})



# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)

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
anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features, reference = c(ind))

# this command creates an 'integrated' data assay
## note that the original unmodified data still resides in the 'RNA' assay
integrated <- IntegrateData(anchorset = anchors)

## start integrated analysis

# specify that we will perform downstream analysis on the corrected data
DefaultAssay(integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, features = all.genes, verbose = FALSE)
# integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunPCA(integrated, features = VariableFeatures(object = integrated), npcs = 20, verbose = FALSE)
# determine dimensionality of the dataset by the Jackstraw procedure (if takees too long, try something else)
integrated <- JackStraw(integrated, num.replicate = 100, dims = length(integrated[["pca"]]))
integrated <- ScoreJackStraw(integrated, dims = 1:length(integrated[["pca"]]))
# determine amount of PCs based on p-value 0.05 of the jackstraw based method
jackstraw_p_value <- 0.05
PC_p_values <- integrated[["pca"]]@jackstraw@overall.p.values[,'Score'] < jackstraw_p_value
choose_n_PC <- which(PC_p_values==FALSE)[1]-1


# integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:length(integrated[["pca"]]))
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:length(integrated[["pca"]]))


# Visualization
p <- DimPlot(integrated, reduction = "umap", group.by = 'orig.ident')
p <- p + DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE)
ggplot2::ggsave(file = paste0("UMAP_", sample_name, "_grouped.png"), width = 30, height = 20, units = "cm")


p <- DimPlot(integrated, reduction = "umap", label = TRUE, split.by = "orig.ident")
ggplot2::ggsave(file = paste0("UMAP_", sample_name, "_split.png"), width = 30, height = 20, units = "cm")






# For performing differential expression after integration, we switch back to the original data (RNA)
DefaultAssay(integrated) <- "RNA"
markers <- FindConservedMarkers(integrated, ident.1 = 0, ident.2 = NULL,
                                grouping.var = "orig.ident", meta.method = metap::minimump, verbose = TRUE)
write.csv2(markers, file = paste0("DEG-analysis_genes-list.csv"))

astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")

plot_DEF <- function(data, features, name) {
  p <- FeaturePlot(data, features = features)
  ggplot2::ggsave(file=paste0("DEG-analysis_", name, "_feature-plot.png"), width = 30, height = 20, units = "cm")
  p <- VlnPlot(data, features = features)
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_violin-plot.png"), width = 30, height = 20, units = "cm")
  # p <- DoHeatmap(data, features = features) + NoLegend()
  # ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_heatmap.png"), width = 30, height = 20, units = "cm")
  p <- RidgePlot(data, features = features, ncol = 3)
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_ridge-plot.png"), width = 30, height = 20, units = "cm")

  # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
  cell.num <- table(data$seurat_clusters)
  cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
  levels(Idents(data)) <- cluster.labels
  p <- DotPlot(data, features = features) + RotatedAxis() + WhiteBackground()
  ggplot2::ggsave(file = paste0("DEG-analysis_", name, "_dot-plot.png"), width = 30, height = 20, units = "cm")
  levels(Idents(data)) <- c(0:(length(levels(Idents(data)))-1))
}
# plot_DEF(data = integrated, features = unique(topn), name = "top-features")
plot_DEF(data = integrated, features = astrocyte_interest, name = "astrocyte")
plot_DEF(data = integrated, features = neuron_interest, name = "neuron")
plot_DEF(data = integrated, features = schema_psych_interest, name = "SCHEMA")

# save integrated Seurat object
saveRDS(data, file = "integrated.rds")













### THOUGHTS and ISSUES
## annotation can happen at 3 points, giving more options again for what the ideal sequence of analysis steps might be
### 1: after individual processing, and using the same panel on both samples
#### difference in pure sample, developing vs mature and mix sample, developing vs mature vs other cell type (neurons/astrocytes), might introduce more bias (highlight differences more?)
### 2/3: after RDS+integration/independent+integration
#### first creates structure by processing and then annotation, might align cells more (highlight similarities more?)
## issue: selecting clusters (by annotation), or rather the cells belonging to specific clusters, is not built into Seurat
### they simply don't need it as they would normally only need specific clusters at DEA level and that is possible
### we on the other hand have a design with mixed cell populations that we want to separate before analysis
#### it's possible, it just takes me more time to figure out proper selection and then perform reanalysis
#### because the Seurat object is complicated, using slots for data types and set options on the background
## automating annotation will only be possible after we have defined which genes and their expressions relate to which celltypes
## automating cluster/cell selection requires predefined gene sets combined with relative expression levels as well
### I don't know if automation or manual work here is the norm or proper approach
## understanding the proper theoretical design approach and/or how to objectively compare results from designs
head(integrated@meta.data)
table(integrated$orig.ident)
table(integrated$seurat_clusters)

integrated@active.assay
integrated@assays
integrated@assays$integrated

integrated@assays$RNA@meta.features
integrated@assays$integrated@meta.features
integrated@assays$RNA@var.features
length(integrated@assays$integrated@var.features)

integrated@reductions$pca@jackstraw
integrated







# subset data
`%notin%` <- Negate(`%in%`)
subset = subset(integrated, seurat_clusters %notin% c(3,4,5,11)) # or use idents = c() instead of seurat_clusters
head(subset@meta.data)
table(subset@meta.data$seurat_clusters)
# reprocess subset data
# specify that we will perform downstream analysis on the corrected data
DefaultAssay(subset) <- "integrated"

# Run the standard workflow for visualization and clustering
subset <- ScaleData(subset, features = all.genes, verbose = FALSE)
# subset <- RunPCA(subset, npcs = 30, verbose = FALSE)
subset <- RunPCA(subset, features = VariableFeatures(object = subset), npcs = 20, verbose = FALSE)
# determine dimensionality of the dataset by the Jackstraw procedure (if takees too long, try something else)
subset <- JackStraw(subset, num.replicate = 100, dims = length(subset[["pca"]]))
subset <- ScoreJackStraw(subset, dims = 1:length(subset[["pca"]]))
# determine amount of PCs based on p-value 0.05 of the jackstraw based method
jackstraw_p_value <- 0.05
PC_p_values <- subset[["pca"]]@jackstraw@overall.p.values[,'Score'] < jackstraw_p_value
choose_n_PC <- which(PC_p_values==FALSE)[1]-1


# subset <- FindNeighbors(subset, reduction = "pca", dims = 1:30)
subset <- FindNeighbors(subset, dims = 1:length(subset[["pca"]]))
subset <- FindClusters(subset, resolution = 0.5)
subset <- RunUMAP(subset, reduction = "pca", dims = 1:length(subset[["pca"]]))


# Visualization
p <- DimPlot(subset, reduction = "umap", group.by = 'orig.ident')
p <- p + DimPlot(subset, reduction = "umap", label = TRUE, repel = TRUE)
ggplot2::ggsave(file = paste0("Integration subset alignment UMAP - subset", sample_name, " - ref r1 - equal parameters - grouped and clustered.png"), width = 30, height = 20, units = "cm")


p <- DimPlot(subset, reduction = "umap", split.by = "orig.ident")
ggplot2::ggsave(file = paste0("Integration subset alignment UMAP - subset", sample_name, "- ref r1 - equal parameters - split.png"), width = 30, height = 20, units = "cm")


# For performing differential expression after integration, we switch back to the original data (RNA)
DefaultAssay(subset) <- "RNA"
markers <- FindConservedMarkers(subset, ident.1 = 0, ident.2 = NULL,
                                grouping.var = "orig.ident", meta.method = metap::minimump, verbose = TRUE)

# plotting
# plot_DEF(data = subset, features = unique(topn), name = "subset-top-features")
plot_DEF(data = subset, features = astrocyte_interest, name = "subset-astrocyte")
plot_DEF(data = subset, features = neuron_interest, name = "subset-neuron")
plot_DEF(data = subset, features = schema_psych_interest, name = "subset-SCHEMA")





### TODO use discussed gene panels for cluster selection (high level annotation)
## automate it in pipeline with an option
## see if annotation selection changes when initiated from RDS/individually processed -> integrated
## create functions that reduce redundancy in the pipeline
# - Neurons: NEUROG2, SYN1, RBFOX3
# - Astrocytes: VIM, S100B, SOX9

### TODO check why heatmap doesnt run from RDS + integrated ?
## error: No requested features found in the scale.data slot for the RNA assay.

### TODO
## create topn: loop door alle clusters (verander ident.1) (maak soortvan findallclusters na van de individuele DE analyse)
## zoek topN voor elk cluster en gebruik deze als genpanel voor plotting
### annotatie wordt hier belangrijk, vooral als we specifiek naar subceltypes willen kijken/vergelijken


### TODO
## kan "Identify differential expressed genes across conditions" (vignet) nadoen zonder annotatie om te kijken hoe
## de cellen in het algemeen veranderen
### met annotatie duidelijker en specifieker natuurlijk









### Seurat --> Monocle 3 for pseudo time analysis
# main site: https://cole-trapnell-lab.github.io/monocle3/
# main paper (cite): https://www.nature.com/articles/s41586-019-0969-x
# Seurat -> Monocle vignette: https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
# Monocle tutorial: http://cole-trapnell-lab.github.io/monocle-release/monocle3/#tutorial-1-learning-trajectories-with-monocle-3
# Monocle -> TradeSeq: https://bioconductor.org/packages/release/bioc/vignettes/tradeSeq/inst/doc/Monocle.html
## TradeSeq: An R package that allows analysis of gene expression along trajectories


#### possible pseudotime analysis --- compare to official Seurat vignette approach: https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
# Hi everyone,
# I've just done the import and pseudotime from seurat v3 to monocle v3 using a seurat integrated object. I used the integrated assay because I wanted monocle to map a trajectory onto the merged samples and onto essentially the same UMAP as the seurat object had. I think if I'd used the RNA assay, monocle wouldn't be able to batch-correct and the samples wouldn't even be in the same clusters. Here's what I did. Please let me know if you think I did something wrong lol!
# `
# seurat.object <- readRDS("Integrated.rds")
# data <- as(as.matrix(GetAssayData(seurat.object, assay = "integrated", slot = "scale.data")), 'sparseMatrix')
#
# pd <- data.frame(seurat.object@meta.data)
# #keep only the columns that are relevant
# pData <- pd %>% select(orig.ident, nCount_RNA, nFeature_RNA)
# fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#
# #Construct monocle cds
# monocle.object <- new_cell_data_set(expression_data = data, cell_metadata = pData, gene_metadata = fData)
# #preprocess
# monocle.object = preprocess_cds(monocle.object, num_dim = 100, norm_method = "size_only", pseudo_count = 0)
# monocle.object = reduce_dimension(monocle.object)
# #map pseudotime
# monocle.object = order_cells(monocle.object, reduction_method = "UMAP")
# monocle.object = learn_graph(monocle.object)
# plot_cells(monocle.object, color_cells_by = "pseudotime")
#
# `


