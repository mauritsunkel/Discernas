#' Perform Seurat sample integration.
#'
#' Seurat sample integration by Canonical Correlation Analysis.
#'
#' @param sample_files character vector with paths to .rds data files to be integrated
#' @param sample_names character vector with sample names of .rds data files
#' @param output_dir Package home directory, used to create output directory for results.
#' @param ref_sample default: NULL. Otherwise, a sample within sample_names to be used as CCA reference anchor dataset.
#' @param perform_cluster_level_selection perform marker selection based on selection panel and selection percent expressed for each cluster, then reintegrate (note: advise to use either cell or cluster based selection, not both)
#' @param perform_cell_level_selection perform marker selection based on selection panel and selection percent expressed for each cell, then reintegrate (note: advise to use either cell or cluster based selection, not both)
#' @param selection_panel default: c(). Gene/feature marker panel for selection and reintegration
#' @param selection_percent_expressed default: 20. Minimal percentage threshold for each marker in panel to be expressed in order to take whole cluster/cell into account
#'
#' @export
#'
#' @examplesIf FALSE
#' output_dir <- file.path("EMC-SKlab-scRNAseq", "results")
#'
#' # files and sample names
#' sample_files <- c(file.path("EMC-SKlab-scRNAseq", "results", "t1.rds"),
#'                file.path("EMC-SKlab-scRNAseq", "results", "t2.rds"),
#'                file.path("EMC-SKlab-scRNAseq", "results", "t3.rds"))
#' sample_names <- c("t1", "t2", "t3")
#' selection_panel <- c("GENES", "OF", "INTEREST")
#' samples_integration(sample_files, sample_names, output_dir,
#'                     selection_panel = selection_panel)
#'
#' @note During development notes
#' Explored different selectionp panels
#' - Neurons: c("MAP2", "DCX", "NEUROG2") # RBFOX3 <-> DCX
#' - Astrocytes: c("VIM", "S100B", "SOX9") # SOX9 <-> FABP7
#'
#' Calling leidenalg via reticulate to run Leiden algorithm instead of Louvain
#' algorithm with Seurat::FindClusters() used to work, now it has stopped
#' working. Still trying to find the cause of this, issued on their Github.
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import RcppAnnoy
samples_integration <- function(sample_files, sample_names, output_dir,
                                ref_sample = NULL,
                                perform_cluster_level_selection = TRUE,
                                perform_cell_level_selection = FALSE,
                                selection_panel = c(),
                                selection_percent_expressed = 20,
                                features_of_interest = features_of_interest) {
  # read/load all sample data
  data.list <- lapply(X = sample_files, FUN = function(x) {
    readRDS(file = x)
  })

  # set sample name for integration
  sample_name <- paste(sample_names, collapse = "-")
  # set reference sample for integration
  ref_sample <- if (is.null(ref_sample)) stringr::str_split(sample_name, "-")[[1]][1] else ref_sample
  # set ref sample index
  for (i in seq_along(1:length(data.list))) if (levels(data.list[[i]]$orig.ident) == ref_sample) ind <- i

  # initialize start time and directories
  output_dir <- file.path(output_dir, 'integrated', sample_name)
  dir.create(file.path(output_dir, 'plots'), recursive = T)
  ### END INITIALIZATION



  # select repeatedly variable features across data sets
  features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- Seurat::PrepSCTIntegration(object.list = data.list, anchor.features = features)


  # run Canonical Correlation Analysis (CCA) to find 'anchors' between data sets
  anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, reference = c(ind))
  rm(data.list)

  # integrate data by overlapping data spaces by integration anchors, creating 'integrated' data assay
  integrated <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")
  rm(anchors)

  integration_analysis <- function(integrated, selection_performed = FALSE) {
    if (selection_performed) {
      output_dir <- file.path(output_dir, 'postSelect')
      dir.create(file.path(output_dir, 'plots'), recursive = T)
    }

    # run the workflow for visualization and clustering on integrated assay
    SeuratObject::DefaultAssay(integrated) <- "integrated"
    integrated <- Seurat::RunPCA(integrated, features = SeuratObject::VariableFeatures(object = integrated), npcs = 50, verbose = FALSE)
    choose_N_PCs <- 20
    integrated <- Seurat::FindNeighbors(integrated, dims = 1:choose_N_PCs)
    # could give warning: "NAs introduced by coercion" as '.' in data will be coerced to NA
    integrated <- Seurat::FindClusters(integrated, resolution = 0.5, algorithm = 1)
    # TODO use Leiden > Louvain algorithm again (not working because of unable to call leidenalg via reticulate anymore) integrated <- Seurat::FindClusters(integrated, resolution = 0.5, algorithm = 4, method = "igraph")
    integrated <- Seurat::RunUMAP(integrated, reduction = "pca", dims = 1:choose_N_PCs)

    # prepare data (recorrect counts) for SCT assay DEG and visualization
    integrated <- Seurat::PrepSCTFindMarkers(integrated, assay = "SCT")
    SeuratObject::DefaultAssay(integrated) <- "SCT"

    # save Seurat object in .RDS data file
    saveRDS(integrated, file = file.path(output_dir, paste0(sample_name, ".rds")))



    ### VISUALIZATION
    p1 <- Seurat::DimPlot(integrated, reduction = "umap", group.by = 'orig.ident') +
      ggplot2::labs(title = "Original sample identity") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE) +
      ggplot2::labs(title = "Integrated") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    # initiate plot_list for arranging ggplot objects in final visualization
    plot_list <- list(p1, p2)
    for (sample in sample_names) {
      p <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, cells = names(integrated$orig.ident[integrated$orig.ident == sample])) +
        ggplot2::labs(title = sample) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
      plot_list[[length(plot_list)+1]] <- p
    }
    # create arranged visualization
    p <- do.call("grid.arrange", c(plot_list, ncol=2))
    ggplot2::ggsave(file.path(output_dir, paste0("UMAPs_", sample_name, ".png")), plot = p, width = c(12,12), height = c(12,12))
    dev.off()

    # define expression visualization function
    plot_DEG <- function(data, features, name, sample_order = NULL) {
      dir.create(file.path(output_dir, 'plots' , name, 'feature'), recursive = T)
      dir.create(file.path(output_dir, 'plots', name, 'feature_split'))

      # set plot sample order
      if(!is.null(sample_order)) {
        data$orig.ident <- factor(data$orig.ident, levels = sample_order)
      }

      # plot feature expression, if available in Seurat
      for (i in seq_along(features)) {
        tryCatch({
          p <- Seurat::FeaturePlot(data, features = features[i])
          ggplot2::ggsave(file = file.path(output_dir, 'plots', name , 'feature', paste0(features[i], ".png")), width = 30, height = 20, units = "cm")

          p <- Seurat::FeaturePlot(data, features = features[i], split.by = "orig.ident", by.col = FALSE, order = TRUE, cols = c("grey", "red"))
          ggplot2::ggsave(file = file.path(output_dir, 'plots', name , 'feature_split', paste0(features[i], ".png")), width = 30, height = 20, units = "cm")
        },
        error=function(e) {
          message(features[i], ' plot is skipped, as feature was not found with FetchData')
        })
      }

      # expression plots
      p <- Seurat::VlnPlot(data, features = features, split.by = "orig.ident") + Seurat::RestoreLegend()
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'violin-split.png'), width = 30, height = 20, units = "cm")
      p <- Seurat::FeaturePlot(data, features = features, order = TRUE)
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'features.png'), width = 30, height = 20, units = "cm")
      p <- Seurat::VlnPlot(data, features = features)
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'violins.png'), width = 30, height = 20, units = "cm")
      p <- Seurat::RidgePlot(data, features = features, ncol = 3)
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'ridges.png'), width = 30, height = 20, units = "cm")

      # dotplot with custom labels
      cell.num <- table(data$seurat_clusters)
      cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
      levels(SeuratObject::Idents(data)) <- cluster.labels
      p <- Seurat::DotPlot(data, features = features) + Seurat::RotatedAxis() + Seurat::WhiteBackground()
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'dots.png'), width = 30, height = 20, units = "cm")

      p <- Seurat::DotPlot(data, features = features, split.by = "orig.ident", cols="RdYlGn") + Seurat::RotatedAxis() + Seurat::WhiteBackground()
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'dots-split.png'), width = 30, height = 20, units = "cm")

      levels(SeuratObject::Idents(data)) <- c(0:(length(levels(SeuratObject::Idents(data)))-1))
      Seurat::DefaultAssay(data) <- "SCT"

      p <- Seurat::DoHeatmap(data, features = features) + Seurat::NoLegend()
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'heatmap.png'), width = 30, height = 20, units = "cm")
    }

    for (feat_name in names(features_of_interest)) {
      plot_DEG(data = data, features = features_of_interest[[feat_name]], name = feat_name, sample_order = sample_names)
    }

    saveRDS(integrated, file = file.path(output_dir, paste0(sample_name, ".rds")))

    # if selection not yet and to be ran, return object
    if (!selection_performed && perform_cluster_level_selection || perform_cell_level_selection) {
      return(integrated)
    }
  }

  # run integrated analysis
  integrated <- integration_analysis(integrated, selection_performed = FALSE)

  # run marker selection & rerun integration analysis
  if (perform_cluster_level_selection || perform_cell_level_selection) {
    if (perform_cluster_level_selection && perform_cell_level_selection) {
      message("WARNING: be advised to rather use cell or cluster based selection, not both simultaneously")
    }
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
    integrated@misc$selection_panel <- selection_panel
    # rerun integration_analysis post selection
    integration_analysis(integrated, selection_performed = TRUE)
  }
}
