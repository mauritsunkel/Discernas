library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
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
  print(paste0("beeped at: ", format(Sys.time(), "%F %H-%M-%S")))
}

### USER PARAMETERS
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"

# load future library and set plan to run certain functions with multiprocessing
future::plan("multisession", workers = 1) # DEVNOTE: n_workers > 1 for parallelization (for me, 5 is max, 4 is safe)
# TODO check if paralellization holds up, otherwise maybe try to apply it per function? (also in individual)

# set to run Nowakowski label transfer
## TODO should become Kriegstein label transfer
run_label_transfer <- FALSE
# set to perform selection after integration pre-selection
perform_selection <- TRUE
perform_cell_level_selection <- TRUE
perform_cluster_level_selection <- FALSE

## SELECTION PARAMETERS FOR NEURONS
# define selection panel (varies per selection type)
selection_panel <- c("MAP2", "NEUROG2", "RBFOX3") # RBFOX3 <-> DCX
# define panel type
selection_panel_type <- "neuronal"
## SELECTION PARAMETERS FOR ASTROCYTES
# # define selection panel (varies per selection type)
# selection_panel <- c("VIM", "S100B", "SOX9") # SOX9 <-> FABP7
# # define panel type
# selection_panel_type <- "astrocytical"
## SELECTION PARAMETERS IN GENERAL
# define minimally percent expression in cluster for each feature from selection_panel
selection_percent_expressed <- 20

# files and sample names
rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_29-06/BL_N/BL_N.rds",
               "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_29-06/BL_C/BL_C.rds")
sample_name <- "BL_N + BL_C"
ref_sample <- "BL_C"

# initialize start time and directories
start_time <- format(Sys.time(), "%F %H-%M-%S")
dir.create(paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/"), recursive = T)
setwd(paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/"))
dir.create("DE_analysis/")
dir.create("DE_analysis/markers/")
dir.create("DE_analysis/sample_markers/")
dir.create("DE_analysis/condition_markers/")
dir.create("DE_analysis/conserved_markers/")
dir.create("Plots/")
dir.create("GSEA_analysis/")
### END USER PARAMETERS




data.list <- lapply(X = rds.files, FUN = function(x) {
  readRDS(file = x)
})

# select features that are repeatedly variable across datasets for integration
## nfeatures = 3000 & PrepSCTIntegration, because using SCTransform
features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- Seurat::PrepSCTIntegration(object.list = data.list, anchor.features = features)

# use Canonical Correlation Analysis (CCA) to find 'anchors' between datasets
## with reference is faster and finds the same amount of anchors
# define combined sample index in data.list position
for (i in seq_along(1:length(data.list))) if (levels(data.list[[i]]$orig.ident) == ref_sample) ind <- i

# use combined sample as reference dataset
anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, reference = c(ind))
# clear data.list to save RAM
rm(data.list)

# this command creates an 'integrated' data assay
## note that the original unmodified data still resides in the 'RNA' assay
integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")
# clear anchors to save RAM
rm(anchors)

# define integrated analysis
integration_analysis <- function(integrated, selection_performed = FALSE) {
  if (selection_performed) {
    dir.create(paste0(getwd(), '/after_selection/'))
    dir.create(paste0(getwd(), '/after_selection/Plots/'))
    setwd(paste0(getwd(), '/after_selection/'))

  }

  # specify that we will perform downstream analysis on the corrected data
  SeuratObject::DefaultAssay(integrated) <- "integrated"

  # Run the standard workflow for visualization and clustering
  integrated <- Seurat::RunPCA(integrated, features = SeuratObject::VariableFeatures(object = integrated), npcs = 50, verbose = FALSE)
  choose_N_PCs <- 20 # default 20
  integrated <- Seurat::FindNeighbors(integrated, dims = 1:choose_N_PCs)
  # algorithm = 4 (= Leiden algorithm) & use: method = "igraph" (for large datasets when using Leiden algorithm)
  integrated <- Seurat::FindClusters(integrated, resolution = 0.5, algorithm = 4, method = "igraph")
  integrated <- Seurat::RunUMAP(integrated, reduction = "pca", dims = 1:choose_N_PCs)

  # prep data (recorrect counts) for SCT DEG and visualization
  integrated <- PrepSCTFindMarkers(integrated, assay = "SCT")
  DefaultAssay(integrated) <- "SCT"

  # run label transfer for integrated analysis after plotting regular UMAP
  if (run_label_transfer) {
    dir.create('Label_transfer/')
    # TODO if data already loaded at individual analysis, then no need to load it here
    ### get Nowakowski as reference (train) data
    nowakowski_data_path <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Nowakowski/exprMatrix.tsv.gz'
    nowakowski_metadata_path <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Nowakowski/meta.tsv'
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

    # set data.frame
    mat = data.frame(mat[,-1], row.names=genes)
    # set data.frame to matrix to dgCMatrix for memory issues
    mat <- as(object = as(mat, "matrix"), Class = "dgCMatrix")
    # get SO from matrix object
    so <- SeuratObject::CreateSeuratObject(counts = mat, project = nowakowski_project, meta.data=meta)
    # subset train data by overlapping genes
    so <- subset(so, cells = names(so$WGCNAcluster)[so$WGCNAcluster != ''], features = genes.overlap)
    # set SO to SCE object
    so.sce <- Seurat::as.SingleCellExperiment(so)
    # subset test data by overlapping genes
    integrated.sub <- subset(integrated, features = genes.overlap)
    # set SO to SCE object
    integrated.sce <- Seurat::as.SingleCellExperiment(integrated.sub, assay = 'SCT')

    # perform SingleR label transfer
    integrated.sce.labels <- SingleR::SingleR(
      test=integrated.sce,
      ref=so.sce,
      labels=so.sce$WGCNAcluster,
      de.method='wilcox',
      aggr.ref=TRUE)
    # TODO or check add aggr.ref = T here too?
    integrated.sce.clusters <- SingleR::SingleR(
      test=integrated.sce,
      ref=so.sce,
      labels=so.sce$WGCNAcluster,
      clusters=integrated.sce$seurat_clusters,
      de.method='wilcox')

    # get a quick a idea of label distribution
    table(integrated.sce.labels$labels)
    table(integrated.sce.clusters$labels)
    # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
    png(paste0("Label_transfer/cell_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
    plot(SingleR::plotScoreHeatmap(integrated.sce.labels))
    dev.off()
    png(paste0("Label_transfer/cluster_heatmap_scores=first_labels=tuned_", sample_name, ".png"))
    png(SingleR::plotScoreHeatmap(integrated.sce.clusters))
    dev.off()
    # ambiguity 'predicted' by 'deltas', i.e. difference from mean scores for all labels (higher=better resolution)
    # png(paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"))
    p <- SingleR::plotDeltaDistribution(integrated.sce.labels, ncol = 6)
    ggplot2::ggsave(file=paste0("Label_transfer/cell_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
    # dev.off()
    # png(paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"))
    p <- SingleR::plotDeltaDistribution(integrated.sce.clusters, ncol = 6)
    ggplot2::ggsave(file=paste0("Label_transfer/cluster_deltas_distribution_", sample_name, ".png"), width = 30, height = 20, units = "cm")
    # dev.off()
    # print warning if too much pruning happens off of default settings, could mean bad quality of label transfer
    if (table(pruneScores(integrated.sce.labels))[1] < length(integrated.sce.labels$labels)*0.99) {
      print("WARNING: SingleR::pruneScores function showed more cells pruned than 99% margin of Gauss dist, look at quality of label transfer")
    }
    # save the cells label transfer results in the Seurat object
    integrated$SingleR.cell.labels <- integrated.sce.labels$labels
    integrated$SingleR.cell.first.labels <- integrated.sce.labels$first.labels

    # add seurat cluster to transfer label for later identification with DE analysis
    integrated.sce.clusters$labels <- paste(as.character(integrated.sce.clusters$labels), levels(integrated), sep='.')
    integrated.sce.clusters$first.labels <- paste(as.character(integrated.sce.clusters$first.labels), levels(integrated), sep='.')

    # original Seurat clusters are needed for RenameIdents function
    names(integrated.sce.clusters$first.labels) <- levels(integrated)
    # set idents/names of clusters to labels of transfer
    integrated <- SeuratObject::RenameIdents(integrated, integrated.sce.clusters$first.labels)
    p1 <- Seurat::DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.first.labels', repel = TRUE) +
      ggplot2::labs(title = "Cells first labels") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) +
      ggplot2::labs(title = "Clusters first labels") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # original Seurat clusters are needed for RenameIdents function
    names(integrated.sce.clusters$labels) <- levels(integrated)
    # set idents/names of clusters to labels of transfer
    integrated <- SeuratObject::RenameIdents(integrated, integrated.sce.clusters$labels)
    # create UMAP on cell and cluster level with transferred labels!
    p3 <- Seurat::DimPlot(integrated, reduction = "umap", label = FALSE, pt.size = 0.5, group.by = 'SingleR.cell.labels', repel = TRUE) +
      ggplot2::labs(title = "Cells tuned labels") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p4 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, pt.size = 0.5, repel = TRUE) +
      ggplot2::labs(title = "Clusters tuned labels") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    # create and save total visualizations
    plot_list <- list(p2, p4, p1, p3)
    p <- do.call("grid.arrange", c(plot_list, ncol=2))
    ggplot2::ggsave(paste0("Label_transfer/labels_UMAP_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))

  }

  # Visualization
  p1 <- Seurat::DimPlot(integrated, reduction = "umap", group.by = 'orig.ident') +
    ggplot2::labs(title = "Original sample identity") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE) +
    ggplot2::labs(title = "Integrated") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # initiate plot_list for arranging ggplots in final visualization
  plot_list <- list(p1, p2)
  for (sample in stringr::str_split(sample_name, " \\+ ")[[1]]) {
    p <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, cells = names(integrated$orig.ident[integrated$orig.ident == sample])) +
      ggplot2::labs(title = sample) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    # add individual identity plot to plot_list
    plot_list[[length(plot_list)+1]] <- p
  }
  # create total visualization
  p <- do.call("grid.arrange", c(plot_list, ncol=2))
  # save total visualization
  ggplot2::ggsave(paste0("UMAPs_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))
  dev.off() # TODO test if this works to get heatmap instead of old plot in label transfer after selection

  ### DEPRECATED: was for Seurat SCT v1, v2 now performs DEG and visualization on the SCT assay data slot
  # # For performing differential expression after integration, we switch back to the original data (RNA, not SCT)
  # ## SCT are now the normalized (and scaled?) Pearson residuals for each data set, prior to integration
  # SeuratObject::DefaultAssay(integrated) <- "RNA"
  # # based on the test used with any of the FindMarkers or derived Seurat functions the RNA counts or normalized data will be used, which are both in different data slots
  # ## because of using SCTransform the RNA assay data is not yet normalized, the data need not be scaled as the scale.data slot is never used for DE
  # ### proofs by Seurat responses in Github issue numbers: 1836, 2023, 3839, 4032
  # integrated <- Seurat::NormalizeData(integrated, normalization.method = "LogNormalize", scale.factor = 10000)

  ## perform visualization
  astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                          "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE")
  neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1", "DCX")
  schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                             "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
  sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                           "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                           "HOPX", "FAM107A", "AGT") # AGXT2L1 not in our data
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
  # juliette_paper <- c("BCOR", "BCORL1", "BCL6", "MLLT3", "TFAP2A", "KDM2B")
  plot_DEF <- function(data, features, name) {
    dir.create(paste0("Plots/", name, "/"))
    dir.create(paste0("Plots/", name, "/Feature_split/"))

    for (i in seq_along(features)) {
      p <- Seurat::FeaturePlot(data, features = features[i], split.by = "orig.ident", cols = c("grey", "red"))
      ggplot2::ggsave(file=paste0("Plots/", name ,"/Feature_split/", features[i], ".png"), width = 30, height = 20, units = "cm")
    }
    p <- Seurat::VlnPlot(data, features = features, split.by = "orig.ident")
    ggplot2::ggsave(file = paste0("Plots/", name, "/violin-split.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::FeaturePlot(data, features = features)
    ggplot2::ggsave(file=paste0("Plots/", name, "/features.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::VlnPlot(data, features = features)
    ggplot2::ggsave(file = paste0("Plots/", name, "/violins.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::RidgePlot(data, features = features, ncol = 3)
    ggplot2::ggsave(file = paste0("Plots/", name, "/ridges.png"), width = 30, height = 20, units = "cm")
    # Change cluster labels from 0, 1, 2 etc to labels for DotPlot only
    cell.num <- table(data$seurat_clusters)
    cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
    levels(SeuratObject::Idents(data)) <- cluster.labels
    p <- Seurat::DotPlot(data, features = features) + Seurat::RotatedAxis() + Seurat::WhiteBackground()
    ggplot2::ggsave(file = paste0("Plots/", name, "/dots.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::DotPlot(data, features = features, split.by = "orig.ident") + Seurat::RotatedAxis() + Seurat::WhiteBackground()
    ggplot2::ggsave(file = paste0("Plots/", name, "/dots-split.png"), width = 30, height = 20, units = "cm")
    levels(SeuratObject::Idents(data)) <- c(0:(length(levels(SeuratObject::Idents(data)))-1))
    Seurat::DefaultAssay(data) <- "SCT"
    p <- Seurat::DoHeatmap(data, features = features) + Seurat::NoLegend()
    ggplot2::ggsave(file = paste0("Plots/", name, "/heatmap.png"), width = 30, height = 20, units = "cm")
  }
  print("plot astrocyte panel")
  plot_DEF(data = integrated, features = astrocyte_interest, name = "astrocyte")
  print("plot neuron panel")
  plot_DEF(data = integrated, features = neuron_interest, name = "neuron")
  print("plot SCHEMA panel")
  plot_DEF(data = integrated, features = schema_psych_interest, name = "SCHEMA")
  print("plot Sloan2017 panel")
  plot_DEF(data = integrated, features = sloan_2017_interest, name = "Sloan2017")
  # print("plot Juliette paper BCOR panel")
  # plot_DEF(data = integrated, features = juliette_paper, name = "Juliette paper BCOR")

  # save integrated Seurat object
  saveRDS(integrated, file = paste0(sample_name, ".rds"))

  if (!selection_performed) {
    return(integrated)
  }
}

# perform integrated analysis
integrated <- integration_analysis(integrated, selection_performed = FALSE)

# perform marker selection & rerun integration analysis?
if (perform_selection) {
  print("start performing marker selection and then integration on selected data")

  if (perform_cell_level_selection) {
    assay_data <- GetAssayData(integrated, slot = "data", assay = "SCT")
    cellsToSelect <- sapply(as.data.frame(assay_data[selection_panel, ] > 0), sum) == length(selection_panel)
    integrated <- integrated[, cellsToSelect]
  }
  if (perform_cluster_level_selection) {
    # create dotplot to extract percent expressed (and average expression) data
    p <- Seurat::DotPlot(integrated, features = selection_panel)
    # get cluster names where percent expressed is above (25)% for each feature (gene) of selection_panel
    cluster_selection <- names(which(table(p$data[p$data$pct.exp > selection_percent_expressed,]$id) == length(selection_panel)))
    # if no clusters selected, set selection to NULL
    if (length(cluster_selection) == 0) {
      cluster_selection <- NULL
    }
    # subset/keep only clusters based on cluster_selection
    integrated <- subset(integrated, idents = cluster_selection)
  }

  # remove empty clusters from original seurat_clusters
  integrated$seurat_clusters <- factor(integrated$seurat_clusters)
  # add selection panel and type to subset file as metadata (for later reference)
  integrated@misc$selection_panel <- as.data.frame(selection_panel, replicate(length(selection_panel), selection_panel_type))
  # perform integration_analysis after selection
  integration_analysis(integrated, selection_performed = TRUE)
}

# TODO remove after automating the pipeline
beep()



### TODO Perform custom Differential Expression (DE) & Gene Set Enrichment Analysis (GSEA)

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


