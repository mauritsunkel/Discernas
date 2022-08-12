library(Seurat)
library(ggplot2)
library(SingleR)
library(gridExtra)



### INITIALIZATION
## USER PARAMETERS
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"

# files and sample names
rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/2022-08-12 13-50-22/BL_N/BL_N.rds",
               "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/2022-08-12 14-42-57/BL_C/BL_C.rds")
sample_name <- "BL_N + BL_C"
ref_sample <- "BL_C"

# set to perform selection after integration and re-run integration
perform_cell_level_selection <- FALSE
perform_cluster_level_selection <- TRUE
selection_panel <- c("MAP2", "DCX", "NEUROG2") # RBFOX3 <-> DCX
selection_panel_type <- "neuronal"
# selection_panel <- c("VIM", "S100B", "SOX9") # SOX9 <-> FABP7
# selection_panel_type <- "astrocytical"
# define minimal expression percentage in each cluster for each feature from selection_panel
selection_percent_expressed <- 20
## END USER PARAMETERS

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
### END INITIALIZATION



# read/load all sample data
data.list <- lapply(X = rds.files, FUN = function(x) {
  readRDS(file = x)
})

# select repeatedly variable features across data sets
features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
data.list <- Seurat::PrepSCTIntegration(object.list = data.list, anchor.features = features)

# set ref sample index
for (i in seq_along(1:length(data.list))) if (levels(data.list[[i]]$orig.ident) == ref_sample) ind <- i

# run Canonical Correlation Analysis (CCA) to find 'anchors' between data sets
anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, reference = c(ind))
rm(data.list)

# integrate data by overlapping data spaces by integration anchors, creating 'integrated' data assay
integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors)

integration_analysis <- function(integrated, selection_performed = FALSE) {
  if (selection_performed) {
    dir.create(paste0(getwd(), '/after_selection/'))
    dir.create(paste0(getwd(), '/after_selection/Plots/'))
    setwd(paste0(getwd(), '/after_selection/'))
  }

  # run the workflow for visualization and clustering on integrated assay
  SeuratObject::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::RunPCA(integrated, features = SeuratObject::VariableFeatures(object = integrated), npcs = 50, verbose = FALSE)
  choose_N_PCs <- 20
  integrated <- Seurat::FindNeighbors(integrated, dims = 1:choose_N_PCs)
  integrated <- Seurat::FindClusters(integrated, resolution = 0.5, algorithm = 4, method = "igraph")
  integrated <- Seurat::RunUMAP(integrated, reduction = "pca", dims = 1:choose_N_PCs)

  # prepare data (recorrect counts) for SCT assay DEG and visualization
  integrated <- Seurat::PrepSCTFindMarkers(integrated, assay = "SCT")
  SeuratObject::DefaultAssay(integrated) <- "SCT"

  # Visualization
  p1 <- Seurat::DimPlot(integrated, reduction = "umap", group.by = 'orig.ident') +
    ggplot2::labs(title = "Original sample identity") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE) +
    ggplot2::labs(title = "Integrated") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  # initiate plot_list for arranging ggplot objects in final visualization
  plot_list <- list(p1, p2)
  for (sample in stringr::str_split(sample_name, " \\+ ")[[1]]) {
    p <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, cells = names(integrated$orig.ident[integrated$orig.ident == sample])) +
      ggplot2::labs(title = sample) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    plot_list[[length(plot_list)+1]] <- p
  }
  # create arranged visualization
  p <- do.call("grid.arrange", c(plot_list, ncol=2))
  ggplot2::ggsave(paste0("UMAPs_", sample_name, ".png"), plot = p, width = c(12,12), height = c(12,12))
  dev.off()

  # define marker panels of interest
  astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                          "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE")
  astrocyte_maturity <- c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                          "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE")
  neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1", "DCX")
  neuron_maturity <- c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                       "SYN1", "SNAP25", "SYT1", "APOE")
  schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                             "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
  sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                           "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                           "HOPX", "FAM107A", "AGT")

  # define expression visualization function
  plot_DEG <- function(data, features, name) {
    dir.create(paste0("Plots/", name, "/"))
    dir.create(paste0("Plots/", name, "/Feature/"))
    dir.create(paste0("Plots/", name, "/Feature_split/"))

    # plot feature expression, if available in Seurat
    for (i in seq_along(features)) {
      tryCatch({
        p <- Seurat::FeaturePlot(data, features = features[i])
        ggplot2::ggsave(file=paste0("Plots/", name ,"/Feature/", features[i], ".png"), width = 30, height = 20, units = "cm")

        p <- Seurat::FeaturePlot(data, features = features[i], split.by = "orig.ident", cols = c("grey", "red"))
        ggplot2::ggsave(file=paste0("Plots/", name ,"/Feature_split/", features[i], ".png"), width = 30, height = 20, units = "cm")
      },
      error=function(e) {
        message(features[i], ' plot is skipped, as feature was not found with FetchData')
      })
    }

    # expression plots
    p <- Seurat::VlnPlot(data, features = features, split.by = "orig.ident")
    ggplot2::ggsave(file = paste0("Plots/", name, "/violin-split.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::FeaturePlot(data, features = features)
    ggplot2::ggsave(file=paste0("Plots/", name, "/features.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::VlnPlot(data, features = features)
    ggplot2::ggsave(file = paste0("Plots/", name, "/violins.png"), width = 30, height = 20, units = "cm")
    p <- Seurat::RidgePlot(data, features = features, ncol = 3)
    ggplot2::ggsave(file = paste0("Plots/", name, "/ridges.png"), width = 30, height = 20, units = "cm")
    # dotplot with custom labels
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
  plot_DEG(data = integrated, features = astrocyte_interest, name = "astrocyte")
  plot_DEG(data = integrated, features = astrocyte_maturity, name = "astrocyte_maturity")
  plot_DEG(data = integrated, features = neuron_interest, name = "neuron")
  plot_DEG(data = integrated, features = neuron_maturity, name = "neuron_maturity")
  plot_DEG(data = integrated, features = schema_psych_interest, name = "SCHEMA")
  plot_DEG(data = integrated, features = sloan_2017_interest, name = "Sloan2017")

  saveRDS(integrated, file = paste0(sample_name, ".rds"))

  # if selection not yet and to be ran, return object
  if (!selection_performed && perform_cluster_level_selection || perform_cell_level_selection) {
    return(integrated)
  }
}

# run integrated analysis
integrated <- integration_analysis(integrated, selection_performed = FALSE)

# run marker selection & rerun integration analysis
if (perform_cluster_level_selection || perform_cell_level_selection) {
  message("start performing marker selection and then rerun integration on selected data")

  if (perform_cell_level_selection) {
    assay_data <- SeuratObject::GetAssayData(integrated, slot = "data", assay = "SCT")
    # select each cell that has expresses each gene from selection_panel
    cellsToSelect <- sapply(as.data.frame(assay_data[selection_panel, ] > 0), sum) == length(selection_panel)
    # perform cell selection
    integrated <- integrated[, cellsToSelect]
  }
  if (perform_cluster_level_selection) {
    # create dotplot to extract percent expressed information
    p <- Seurat::DotPlot(integrated, features = selection_panel)
    # get cluster names where percent expressed is above %threshold for each gene of selection_panel
    cluster_selection <- names(which(table(p$data[p$data$pct.exp > selection_percent_expressed,]$id) == length(selection_panel)))
    # if no clusters selected, set selection to NULL
    if (length(cluster_selection) == 0) {
      cluster_selection <- NULL
    }
    # perform cluster selection
    integrated <- subset(integrated, idents = cluster_selection)
  }

  # remove empty clusters from original seurat_clusters
  integrated$seurat_clusters <- factor(integrated$seurat_clusters)
  # add selection panel and type as metadata
  integrated@misc$selection_panel <- as.data.frame(selection_panel, replicate(length(selection_panel), selection_panel_type))
  # rerun integration_analysis post selection
  integration_analysis(integrated, selection_performed = TRUE)
}
