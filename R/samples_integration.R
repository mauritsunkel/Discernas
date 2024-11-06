#' Perform Seurat sample integration.
#'
#' Seurat sample integration by Canonical Correlation Analysis.
#'
#' @importFrom gridExtra grid.arrange arrangeGrob
#'
#' @param sample_files character vector with paths to .qs data files to be integrated
#' @param sample_names character vector with sample names of .qs data files
#' @param output_dir Package home directory, used to create output directory for results.
#' @param features_of_interest list of gene marker panels used for plotting
#' @param integration_method default: "RPCA", one of c("RPCA", "CCA", "harmony") "harmony", harmony is run from harmony::RunHarmony(), other are Seurat::IntegrateLayers()
#'
#' @export
#'
#' @note During development notes
#' Explored different selection panels
#' - Neurons: c("MAP2", "DCX", "NEUROG2") # RBFOX3 <-> DCX
#' - Astrocytes: c("VIM", "S100B", "SOX9") # SOX9 <-> FABP7
#'
#' Calling leidenalg via reticulate to run Leiden algorithm instead of Louvain
#' algorithm with Seurat::FindClusters() used to work, now it has stopped
#' working. Still trying to find the cause of this, issued on their Github.
#'
#' Seurat integration notes
#' RPCA > CCA when: "a substantial fraction of the cells in one dataset have no matching type in the other" and/or "datasets originate from the same platform (i.e. multiple lanes of 10x genomics)"
#'  the major difference with Seurat V4 is that in Seurat v5, we perform the integration in low-dimensional space. What that means is that prior to performing any integration - we run a PCA on the full dataset (V5: IntegrateLayers() vs V4: IntegrateData(), no more 'integrated' assay)
#'  RPCA struggles with varying cell amounts between samples because of mutual neighborhood constraints, in that case harmony could be better
#' Using harmony::RunHarmony instead of Seurat::IntegrateLayers to have more parameter control.
#' Harmony is fast, memory efficient, removes batch effects, scalable, and ranks among the top in many different type of benchmarking studies
#' - only for the most complex use cases, consider scVI/scANVI/Scanorama/scGen (all Python based)
#' - https://www.nature.com/articles/s41592-021-01336-8
#' - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6964114/
#' - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8053088/
#' - https://www.nature.com/articles/s41467-023-41855-w
#' - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9#Sec9
#' IntegrateLayers(normalization.method = "SCT") replaces: SelectIntegrationFeatures, PrepSCTIntegration, FindIntegrationAnchors, IntegrateData
#' - For IntegrateLayers, specifically using method = SeuratWrappers::FastMNNIntegration, got error  - following issue: https://github.com/satijalab/seurat/issues/8631
#' - For IntegrateLayers, specifically using method = SeuratWrappers::scVIIntegration, need reticulate/conda setup for scVI
samples_integration <- function(
    sample_files, sample_names, output_dir, features_of_interest,
    integration_method = "RPCA", integrated_sample_name = NULL) {
  set.seed(42)
  library(Seurat) # added because of error
  # Error: package or namespace load failed for ‘Seurat’ in .doLoadActions(where, attach):
  #   error in load action .__A__.1 for package RcppAnnoy: loadModule(module = "AnnoyAngular", what = TRUE, env = ns, loadNow = TRUE): Unable to load module "AnnoyAngular": attempt to apply non-function
  # Error in .requirePackage(package) :
  #   unable to find required package ‘Seurat’

  # read/load all sample data
  data.list <- lapply(X = sample_files, FUN = function(x) {
    qs::qread(x)
  })
  ### END INITIALIZATION



  # select repeatedly variable features across data sets
  data.features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  ## DEVNOTE: deprecated, used with Seurat::FindIntegrationAnchors and Seurat::IntegrateData
  ## data.list <- Seurat::PrepSCTIntegration(object.list = data.list, anchor.features = data.features)

  ## DEVNOTE: using base::merge instead of SeuratObject::merge.Seurat as the latter lost cells
  data.merged <- base::merge(data.list[[1]], data.list[2:length(data.list)])
  Seurat::VariableFeatures(data.merged) <- data.features

  data.merged <- Seurat::RunPCA(object = data.merged, assay = "SCT", features = data.features, npcs = min(c(dim(data.merged)[2], 50)))

  data.merged <- run_integration(so = data.merged, integration_method = integration_method)

  if (is.null(integrated_sample_name)) integrated_sample_name <- paste(sample_names, collapse = "-")
  # run integrated analysis
  integration_analysis(
    integrated = data.merged, output_dir = output_dir,
    sample_name = integrated_sample_name, features_of_interest = features_of_interest,
    sample_names = sample_names)
}

#' Run sample layers integration
#'
#' @param so seurat object with layers to be integrated
#' @param integration_method default: "RPCA", one of c("RPCA", "CCA", "harmony") "harmony", harmony is run from harmony::RunHarmony(), other are Seurat::IntegrateLayers()
#'
#' @return so
run_integration <- function(so, integration_method) {
  message("\n RUN integration: ", integration_method, "\n")
  if (integration_method == "harmony") {
    integrated_so <- harmony::RunHarmony(
      object = so,
      group.by.vars = "orig.ident",
      assay.use = "SCT",
      reduction.use = "pca",
      reduction.save = "integrated.dr",
      dims.use = 1:50,
      plot_convergence = FALSE,
      epsilon.cluster=-Inf,
      epsilon.harmony=-Inf,
      verbose = TRUE)
  } else if (integration_method == "RPCA") {
    integrated_so <- tryCatch({
      integrated_so <- Seurat::IntegrateLayers(
        object = so,
        method = Seurat::RPCAIntegration,
        normalization.method = "SCT",
        orig.reduction = "pca",
        new.reduction = "integrated.dr",
        dims = 1:(min(c(table(so$orig.ident), 50))-1),
        k.weight = min(c(table(so$orig.ident), 100))-1,
        verbose = TRUE)
    },
    error=function(e) {
      message(e)
      message("\n Removing smallest sample as data cannot be integrated, rerunning without")
      Seurat::Idents(so) <- so@meta.data[, "orig.ident"]
      sizeSorted_sample_names <- names(sort(table(so$orig.ident)))
      so <- subset(so, idents = sizeSorted_sample_names[2:length(sizeSorted_sample_names)])
      if (length(unique(so$orig.ident)) == 1) {
        stop('Only a single sample leftover, no integration needed')
      }
      integrated_so <- run_integration(so = so, integration_method = integration_method)
      return(integrated_so)
    })
  } else if (integration_method == "CCA") {
    integrated_so <- tryCatch({
      integrated_so <- Seurat::IntegrateLayers(
        object = so,
        method = Seurat::CCAIntegration,
        normalization.method = "SCT",
        orig.reduction = "pca",
        new.reduction = "integrated.dr",
        dims = 1:(min(c(table(so$orig.ident), 50))-1),
        k.weight = min(c(table(so$orig.ident), 100))-1,
        verbose = TRUE)
    },
    error=function(e) {
      message(e)
      message("\n Removing smallest sample as data cannot be integrated, rerunning without")
      Seurat::Idents(so) <- so@meta.data[, "orig.ident"]
      sizeSorted_sample_names <- names(sort(table(so$orig.ident)))
      so <- subset(so, idents = sizeSorted_sample_names[2:length(sizeSorted_sample_names)])
      if (length(unique(so$orig.ident)) == 1) {
        stop('Only a single sample leftover, no integration needed')
      }
      integrated_so <- run_integration(so = so, integration_method = integration_method)
      return(integrated_so)
    })
  }
  return(integrated_so)
}

#' Analysis of integrated samples
#'
#' @param integrated Integrated Seurat object
#' @param output_dir Package home directory, used to create output directory for results.
#' @param sample_names names of samples to be integrated
#' @param sample_name name of integrated sample, combined of samples_names
#' @param features_of_interest gene of interest
#'
#' @export
integration_analysis <- function(integrated, output_dir, sample_name, features_of_interest, sample_names = NULL) {
  message("\n RUNNING integration_analysis: ", output_dir, "\n")
  if (is.null(sample_names)) sample_names <- sample_name
  # initialize start time and directories
  if (!grepl(paste0("/", sample_name, "/"), output_dir)) output_dir <- file.path(output_dir, sample_name)
  dir.create(file.path(output_dir, 'plots'), recursive = TRUE)
  dir.create(file.path(output_dir, 'UMAPs'), recursive = TRUE)

  # prepare data (recorrect counts) for SCT assay DEG: https://satijalab.org/seurat/articles/integration_introduction
  ## Seurat recommends to use recorrected counts for visualization: https://github.com/satijalab/seurat/issues/6675
  integrated <- Seurat::PrepSCTFindMarkers(integrated, assay = "SCT")
  # run the workflow for visualization and clustering
  integrated <- Seurat::FindNeighbors(integrated, assay = "SCT", reduction = "integrated.dr", dims = 1:50)
  # could give warning: "NAs introduced by coercion" as '.' in data will be coerced to NA
  message("\n RUN FindClusters: ", output_dir, "\n")
  integrated <- Seurat::FindClusters(integrated, resolution = 0.8, algorithm = 1)
  message("\n RUN RunUMAP: ", output_dir, "\n")
  integrated <- Seurat::RunUMAP(integrated, assay = "SCT", reduction = "integrated.dr", dims = 1:50)

  ### VISUALIZATION
  p1 <- Seurat::DimPlot(integrated, reduction = "umap", group.by = 'orig.ident') +
    ggplot2::labs(title = "Original sample identity") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p1 <- p1 + ggplot2::xlim(min(integrated@reductions$umap@cell.embeddings[,1]),max(integrated@reductions$umap@cell.embeddings[,1]))
  p1 <- p1 + ggplot2::ylim(min(integrated@reductions$umap@cell.embeddings[,2]),max(integrated@reductions$umap@cell.embeddings[,2]))
  ggplot2::ggsave(file.path(output_dir, "UMAPs", "original-identity.png"), plot = p1, width = c(12,12), height = c(12,12))
  p2 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE) +
    ggplot2::labs(title = "Integrated") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  p2 <- p2 + ggplot2::xlim(min(integrated@reductions$umap@cell.embeddings[,1]),max(integrated@reductions$umap@cell.embeddings[,1]))
  p2 <- p2 + ggplot2::ylim(min(integrated@reductions$umap@cell.embeddings[,2]),max(integrated@reductions$umap@cell.embeddings[,2]))
  ggplot2::ggsave(file.path(output_dir, "UMAPs", "integrated.png"), plot = p2, width = c(12,12), height = c(12,12))
  # initiate plot_list for arranging ggplot objects in final visualization
  plot_list <- list(p1, p2)

  for (sample in unique(integrated$orig.ident)) {
    p3 <- Seurat::DimPlot(integrated, reduction = "umap", label = TRUE, repel = TRUE, cells = names(integrated$orig.ident[integrated$orig.ident == sample])) +
      ggplot2::labs(title = sample) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    p3 <- p3 + ggplot2::xlim(min(integrated@reductions$umap@cell.embeddings[,1]),max(integrated@reductions$umap@cell.embeddings[,1]))
    p3 <- p3 + ggplot2::ylim(min(integrated@reductions$umap@cell.embeddings[,2]),max(integrated@reductions$umap@cell.embeddings[,2]))

    plot_list[[length(plot_list)+1]] <- p3
    ggplot2::ggsave(file.path(output_dir, "UMAPs", paste0(sample, ".png")), plot = p3, width = c(12,12), height = c(12,12))
  }
  # create arranged visualization
  p4 <- do.call(gridExtra::grid.arrange, c(plot_list, ncol=2))
  ggplot2::ggsave(file.path(output_dir, paste0("UMAPs_", sample_name, ".png")), plot = p4, width = c(12,12), height = c(12,12))
  dev.off()

  # define expression visualization function
  plot_DEG <- function(data, data.features, name, sample_order = NULL, output_dir) {
    dir.create(file.path(output_dir, 'plots' , name, 'feature'), recursive = T)
    dir.create(file.path(output_dir, 'plots', name, 'feature_split'))

    # set plot sample order
    if(!is.null(sample_order)) {
      data$orig.ident <- factor(data$orig.ident, levels = sample_order)
    }

    # plot feature expression, if available in Seurat
    for (i in seq_along(data.features)) {
      tryCatch({
        p <- Seurat::FeaturePlot(data, features = data.features[i])
        ggplot2::ggsave(file = file.path(output_dir, 'plots', name , 'feature', paste0(data.features[i], ".png")), width = 30, height = 20, units = "cm")

        p <- Seurat::FeaturePlot(data, features = data.features[i], split.by = "orig.ident", by.col = FALSE, order = TRUE, cols = c("grey", "red"))
        ggplot2::ggsave(file = file.path(output_dir, 'plots', name , 'feature_split', paste0(data.features[i], ".png")), width = 30, height = 20, units = "cm")
      },
      error=function(e) {
        message(data.features[i], ' plot is skipped, as feature was not found with FetchData')
      })
    }

    # expression plots
    p <- Seurat::VlnPlot(data, features = data.features, split.by = "orig.ident") + Seurat::RestoreLegend()
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'violin-split.png'), width = 30, height = 20, units = "cm")
    p <- Seurat::FeaturePlot(data, features = data.features, order = TRUE)
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'features.png'), width = 30, height = 20, units = "cm")
    p <- Seurat::VlnPlot(data, features = data.features)
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'violins.png'), width = 30, height = 20, units = "cm")
    p <- Seurat::RidgePlot(data, features = data.features, ncol = 3)
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'ridges.png'), width = 30, height = 20, units = "cm")

    # dotplot with custom labels
    cell.num <- table(data$seurat_clusters)
    cluster.labels = paste("Cluster", names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
    levels(SeuratObject::Idents(data)) <- cluster.labels
    p <- Seurat::DotPlot(data, features = data.features) + Seurat::RotatedAxis() + Seurat::WhiteBackground()
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'dots.png'), width = 30, height = 20, units = "cm")

    p <- Seurat::DotPlot(data, features = data.features, split.by = "orig.ident", cols="RdYlGn") + Seurat::RotatedAxis() + Seurat::WhiteBackground()
    ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'dots-split.png'), width = 30, height = 20, units = "cm")

    levels(SeuratObject::Idents(data)) <- c(0:(length(levels(SeuratObject::Idents(data)))-1))
    Seurat::DefaultAssay(data) <- "SCT"


    tryCatch({
      p <- Seurat::DoHeatmap(data, features = data.features) + Seurat::NoLegend()
      ggplot2::ggsave(file = file.path(output_dir, 'plots', name, 'heatmap.png'), width = 30, height = 20, units = "cm")
    },
    error=function(e) {
      message(data.features, ' features were not found for DoHeatmap')
    })

  }

  if (!is.null(features_of_interest)) {
    for (feat_name in names(features_of_interest)) {
      plot_DEG(data = integrated, data.features = features_of_interest[[feat_name]], name = feat_name, sample_order = sample_names, output_dir = output_dir)
    }
  }

  # save Seurat object in .qs data file
  qs::qsave(integrated, file = file.path(output_dir, paste0(sample_name, ".qs")), preset = 'custom', algorithm = "zstd_stream", compress_level = 4, shuffle_control = 15, nthreads = 1)
}

#' Perform selection of cells/clusters/annotation
#'
#' Perform selection of cells/clusters/annotation and then reintegrate and reanalyze
#'
#' @param so_filename filename of seurat object to perform selection on
#' @param selection_markers genes of interest to be used for selection
#' @param percent_expressed threshold for percentage of cells that have to express all selection_markers
#' @param reference_annotations reference database and annotation label of it to perform selection on
#' @param output_dir Package home directory, used to create output directory for results
#' @param sample_name name of integrated sample, combined of samples_names
#' @param features_of_interest genes of interest
#' @param exclude_samples default NULL, else vector of sample names to exclude
#' @param integration method default = "RPCA", else "CCA" or "harmony"
#'
#' @export
#'
#' @note
#' SCT + integration workflow: https://github.com/satijalab/seurat/issues/4896
#' subset + reintegration workflow: https://github.com/satijalab/seurat/issues/1883
#' RNA assay set --> selection --> DietSeurat to clean up https://github.com/satijalab/seurat/issues/8073
#' - reanalyze using same workflow you used for the entire dataset as a blank slate
#' no need to do MT_features as they are already regressed out, no need to redo doublet and ambient RNA removal
selection_reintegration <- function(
    so_filename,
    output_dir, sample_name, features_of_interest = NULL, exclude_samples = NULL, integration_method = 'RPCA',
    selection_markers = NULL, percent_expressed = NULL, reference_annotations = NULL, subset_reintegration = TRUE) {

  if (is.null(selection_markers) && is.null(percent_expressed) && is.null(reference_annotations)) {
    stop("To perform selection and reintegration, pass in the parameters")
  }

  dir.create(output_dir, recursive = TRUE)

  message(paste0("\n reading .qs... --> ", output_dir, "\n"))
  so <- qs::qread(file = so_filename)

  if (!is.null(exclude_samples)) {
    Seurat::DefaultAssay(so) <- "RNA"
    Seurat::Idents(so) <- so$orig.ident
    so <- subset(so, idents = setdiff(unique(so$orig.ident), exclude_samples))
  }
  Seurat::DefaultAssay(so) <- "SCT"

  if (!is.null(reference_annotations)) {
    if (length(unique(names(reference_annotations))) != 1) stop("Pass a single reference")
    message("\n Selecting specified annotation from reference as idents \n")
    # set idents to reference name
    Seurat::Idents(so) <- so@meta.data[, names(reference_annotations)]
    to_select <- unique(Seurat::Idents(so)[grepl(reference_annotations[[names(reference_annotations)]], Seurat::Idents(so), ignore.case = TRUE)])
    # select idents by annotation (must be in reference)
    Seurat::DefaultAssay(so) <- "RNA"
    so <- subset(so, idents = to_select)
  } else if (!is.null(percent_expressed) && !is.null(selection_markers)) {
    message("\n Selecting subclusters based on markers and percent expressed \n")
    ## perform subclustering for increased resolution
    Seurat::Idents(so) <- so$seurat_clusters
    so <- Seurat::FindSubCluster(
      object = so,
      cluster = factor(so$seurat_clusters),
      graph.name = "SCT_nn",
      subcluster.name = "seurat_subclusters"
    )
    # set idents to subclustering
    Seurat::Idents(so) <- so$seurat_subclusters
    # create dotplot to extract percent expressed information
    p <- Seurat::DotPlot(so, features = selection_markers)
    png(file.path(output_dir, "subclusters_percentage_expressed.png"))
    hist(p$data$pct.exp, breaks = seq(0, 100, 5), main = paste0("Cells percentage expressed: ", paste(selection_markers, collapse = ", ")))
    dev.off()
    # get cluster names where percent expressed is above %threshold for each gene of selection_markers
    subcluster_selection <- names(which(table(p$data[p$data$pct.exp > percent_expressed,]$id) == length(unique(p$data$features.plot))))
    # if no clusters selected, set selection to NULL
    if (length(subcluster_selection) == 0) subcluster_selection <- NULL
    # set RNA assay
    Seurat::DefaultAssay(so) <- "RNA"
    # plot selected cells
    cell_selection <- so$seurat_subclusters %in% subcluster_selection
    Seurat::Idents(so) <- cell_selection
    p <- Seurat::DimPlot(so, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cell_selection)["TRUE"]))
    ggplot2::ggsave(file.path(output_dir, "subclusters_selected_cells.png"), plot = p, width = c(12,12), height = c(12,12))
    # perform subcluster selection
    Seurat::Idents(so) <- so$seurat_subclusters
    so <- subset(so, idents = subcluster_selection)
    # add selection panel and type as metadata
    so@misc$selection_markers <- selection_markers
  } else if (!is.null(selection_markers)) {
    message("\n Selecting cells based on expressing all markers \n")
    layer_data <- SeuratObject::LayerData(so)
    # select each cell that has expresses each gene from selection_markers
    if (length(levels(so$orig.ident)) == 1) {
      cells_to_select <- names(which(layer_data[selection_markers, ] > 0))
      cells_to_plot <- layer_data[selection_markers, ] > 0
    } else {
      cells_to_select <- sapply(as.data.frame(layer_data[selection_markers, ] > 0), sum) == length(rownames(layer_data[selection_markers, ]))
      # TODO
      cells_to_plot <- NULL
    }

    # set RNA assay
    Seurat::DefaultAssay(so) <- "RNA"
    # plot selected cells
    Seurat::Idents(so) <- cells_to_plot
    p <- Seurat::DimPlot(so, reduction = "umap", label = F, repel = TRUE) + ggplot2::ggtitle(paste0("selected cells: ", table(cells_to_plot)["TRUE"]))
    ggplot2::ggsave(file.path(output_dir, "selected_cells.png"), plot = p, width = c(12,12), height = c(12,12))
    # perform cell selection
    so <- so[, cells_to_select]
    # add selection panel and type as metadata
    so@misc$selection_markers <- selection_markers
  }

  # remove empty clusters from original seurat_clusters
  so$seurat_clusters <- factor(so$seurat_clusters)
  Seurat::Idents(so) <- so$seurat_clusters

  # cleanup filtered Seurat object
  so <- Seurat::DietSeurat(so, assays = c("RNA"))
  # split (recommended by Seurat V5: https://github.com/satijalab/seurat/issues/8406) -> layers
  message("\ splitting seurat object into layers \n")
  split_so <- function(so) {
    so <- tryCatch({
      so[["RNA"]] <- split(so[["RNA"]], f = so$orig.ident)
      return(so)
    },
    error=function(e) {
      message(e)
      message("\n Removing smallest sample as data cannot be split, rerunning without")
      Seurat::Idents(so) <- so@meta.data[, "orig.ident"]
      sizeSorted_sample_names <- names(sort(table(so$orig.ident)))
      so <- subset(so, idents = sizeSorted_sample_names[2:length(sizeSorted_sample_names)])
      if(length(names(table(so$orig.ident))) == 0) stop("All samples to small to split in layers")
      so <- split_so(so)
      return(so)
    })
    return(so)
  }
  so <- split_so(so)
  message("\n RUN SCTransform \n")
  options(future.globals.maxSize = 8000 * 1024^2)
  SCTransform_subset <- function(so) {
    so <- tryCatch({
      so <- Seurat::SCTransform(so, vst.flavor = "v2", method = "glmGamPoi", return.only.var.genes = FALSE)
    },
    error=function(e) {
      message(e)
      message("\n Removing smallest sample as data cannot be processed by SCTransform, rerunning without")
      Seurat::Idents(so) <- so@meta.data[, "orig.ident"]
      sizeSorted_sample_names <- names(sort(table(so$orig.ident)))
      if (length(sizeSorted_sample_names) > 1) {
        so <- subset(so, idents = sizeSorted_sample_names[2:length(sizeSorted_sample_names)])
        so <- SCTransform_subset(so)
        return(so)
      } else {
        stop("All samples to small to run SCTransform on subset, keeping old SCTransformed data")
      }
    })
  }
  so <- SCTransform_subset(so)
  message("\n RUN PCA \n")
  so <- Seurat::RunPCA(so, features = SeuratObject::VariableFeatures(object = so), npcs = min(c(dim(so)[2], 50)), verbose = TRUE)

  if (subset_reintegration) {
    # do integration and integration_analysis post selection
    so <- run_integration(so = so, integration_method = integration_method)
    integration_analysis(
      integrated = so, output_dir = output_dir,
      sample_name = sample_name, features_of_interest = features_of_interest)
  } else {
    return (so)
  }
}
