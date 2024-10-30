#' Analyse individual samples
#'
#' Performs Seurat SCTv2 analysis workflow on individual samples.
#'
#' @param samples_dir Root directory of samples.
#' @param sample_name Current sample name, used to match in samples_dir.
#' @param output_dir Package home directory, used to create output directory for results.
#' @param features_of_interest list of gene marker panels used for plotting
#' @param run_cell_cycle_regression True/False to regress out genes to do with cell cycle, based on Tirosh et al, 2015.
#' @param run_doublet_removal default: TRUE, to run scDblFinder and call and remove doublets, plotted in QC
#' @param doublet_removal_rate default: NULL, automatic doublet removal detection rate, NULL for 10x data okay, otherwise base on sequencing protocol and amount of cells
#' @param run_ambient_RNA_removal default: TRUE, to run SoupX with automatic estimation of contamination rate and correcting counts matrix
#'
#' @importFrom dplyr .data
#'
#' @export
#'
#' @examplesIf FALSE
#' # create output directories based on start time and sample name and set working directory
#' start_time <- format(Sys.time(), "%F %H-%M-%S")
#' output_dir <- file.path("EMC-SKlab-scRNAseq", "results", start_time)
#' # directory where samples are located
#' samples_dir <- file.path("EMC-SKlab-scRNAseq", "data", "samples", "project")
#' # selected sample names from sample dir
#' sample_names <- c('t1', 't2', 't3')
#'
#' for (sample_name in sample_names) {
#'   individual_analysis(samples_dir, sample_name, output_dir)
#' }
#'
#' @note Paper writing & during development notes
#'
#' This SCT->Harmony workflow is inspired by: https://github.com/immunogenomics/harmony/issues/41
#'
#' Cell cycle regression based on: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
#'
#' SCTransform-v2 replaces NormalizeData + FindVariableFeatures + ScaleData & sets default assay to SCT
#' https://satijalab.org/seurat/articles/sctransform_vignette.html & https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
#' normalize gene expression counts per cell by the total expression and applying a
#' scaling factor (default: 10.000) and adding a pseudocount before log-transforming the result
#' this global linear scaling on the data sets mean expression across cells is 0 and variance across cells is 1 as to
#' provide equal weight in downstream processing such that highly variable genes do not dominate results
#' SCTransform-v2 excludes the need for this heuristic pseudocount addition, log-transformation and optimizes variation
#' the top 3000 (default) variable genes are kept for improving downstream processing efficiency
#' vars.to.regress = regress out variability originating from reads mapped to mitochondrial DNA
#' return.only.var.genes = TRUE, as non-sparse matrix is returned and used in PCA
#' set transformed data as default data assay for downstream processing
#'
#' Choosing cell-free droplet detection algorithm, EmptyNN based on comparison with CellRanger embedded EmptyDrops: https://www.cell.com/patterns/fulltext/S2666-3899(21)00154-9
#' Run before removing doublets/multiplets, possibly even before base QC of reming min.cells and min.features
#' Error, issued at: https://github.com/lkmklsmn/empty_nn/issues/4
#'
#' Choosing doublet detection algorithm, DoubletFinder based on benchmerk: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7897250/
#' scDblFinder based on comparison with doubletFinder: https://f1000research.com/articles/10-979/v2
#' run scDblFinder between fundamental QC and SCTransform: https://github.com/plger/scDblFinder/issues/86
#' # for 10x data, you can leave scDblFinder(dbr) parameter NULL: https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
#'
#' Choosing ambient RNA removal algorithm, some are heavily dependent on pre-clustering etc, so SoupX: https://academic.oup.com/gigascience/article/9/12/giaa151/6049831?login=true
#' Good to run doublet removal first and supply only singlets to SoupX.
#' 2-5-10-20% contamination rate low-usual-moderate-high
#'
#' order of doublet, empty, ambient RNA removal unclear: https://github.com/plger/scDblFinder/issues/93
#' order by 10x Genomics: https://www.10xgenomics.com/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data
#'
#' discussion point: if we were to run ambient RNA removal first, then doublet/multiplet detection could be easier and more accurate
sample_analysis <- function(
    samples_dir, sample_name, output_dir, features_of_interest,
    run_cell_cycle_regression = F,
    run_doublet_removal = T, doublet_removal_rate = NULL,
    run_ambient_RNA_removal = T) {
  set.seed(42)

  sample_path <- file.path(output_dir, sample_name)
  if (dir.exists(sample_path)) {
    stop("Sample already exists in output directory, please choose another to avoid overwriting results...")
  }
  dir.create(sample_path, recursive = T)
  dir.create(file.path(sample_path, 'quality_control'), recursive = T)
  dir.create(file.path(sample_path, 'Principal_Component_Analysis'), recursive = T)
  dir.create(file.path(sample_path, 'DE_analysis'), recursive = T)

  ## find samples dir automatically, based on unique sample names from workflow
  pattern <- paste0("/", sample_name, "/")
  samples_dirs_list <- list.dirs(samples_dir, recursive = TRUE)
  samples_dir <- samples_dirs_list[grep(pattern, samples_dirs_list)][1]
  samples_dir <- substr(samples_dir, 1, stringr::str_locate(samples_dir, pattern)[2]-1)

  # read 10X data (preprocessed by 10X Cellranger pipeline) and convert to Seurat object
  data.data <- Seurat::Read10X(data.dir = file.path(samples_dir, "filtered_feature_bc_matrix"))
  ## DEVNOTE: min.features = 500 (CellRanger default)
  data <- Seurat::CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 500)
  rm(data.data)

  if (run_doublet_removal) {
    plot_and_remove_doublets <- function(data, sample_path, sample_name, doublet_removal_rate = NULL) {
      temp_QC_data <- data

      sce <- scDblFinder::scDblFinder(temp_QC_data@assays$RNA$counts, clusters = TRUE, dbr = doublet_removal_rate)
      temp_QC_data$scDblFinder.score <- sce$scDblFinder.score
      temp_QC_data$scDblFinder.class <- sce$scDblFinder.class

      temp_QC_data <- Seurat::NormalizeData(temp_QC_data)
      temp_QC_data <- Seurat::ScaleData(temp_QC_data)
      temp_QC_data <- Seurat::FindVariableFeatures(temp_QC_data)
      temp_QC_data <- Seurat::RunPCA(temp_QC_data, features = SeuratObject::VariableFeatures(object = temp_QC_data), npcs = 50, verbose = FALSE)
      temp_QC_data <- Seurat::RunUMAP(temp_QC_data, reduction = "pca", dims = 1:30)

      p <- Seurat::FeaturePlot(temp_QC_data, features = "scDblFinder.score") +
        ggplot2::labs(subtitle = paste0(
          "singlets: ",
          table(temp_QC_data$scDblFinder.class)[1],
          " - doublets: ",
          paste0(table(temp_QC_data$scDblFinder.class)[2],
                 " - doublet rate: ",
                 round(table(temp_QC_data$scDblFinder.class)[2]/table(temp_QC_data$scDblFinder.class)[1], digits = 3))))
      ggplot2::ggsave(file=file.path(sample_path, "quality_control", paste0("scDblFinder_scores_", sample_name, ".png")), width = 30, height = 20, units = "cm")

      # keep only singlets for downstream processing
      data <- data[, temp_QC_data$scDblFinder.class == "singlet"]
      return(data)
    }
    data <- plot_and_remove_doublets(data, sample_path, sample_name, doublet_removal_rate)
  }

  if (run_ambient_RNA_removal) {
    clean_ambient_RNA <- function(data, samples_dir, sample_name) {
      table_of_droplets = Seurat::Read10X(data.dir = file.path(samples_dir, "raw_feature_bc_matrix"))

      # add non-overlapping genes tot table_of_droplets with 0 counts to satisfy having same amount of genes
      non_overlapping_genes <- rownames(data@assays$RNA$counts)[which(!rownames(data@assays$RNA$counts) %in% rownames(table_of_droplets))]
      l <- list()
      l <- t(sapply(non_overlapping_genes, function(gene) {
        l[[gene]] <- rep(0, length(colnames(table_of_droplets)))
      }))
      colnames(l) <- colnames(table_of_droplets)
      table_of_droplets <- rbind(table_of_droplets, l)

      overlapping_genes <- rownames(data@assays$RNA$counts)[which(rownames(data@assays$RNA$counts) %in% rownames(table_of_droplets))]
      table_of_droplets <- table_of_droplets[overlapping_genes,]
      data@assays$RNA$counts <- data@assays$RNA$counts[overlapping_genes,]

      sc <- SoupX::SoupChannel(tod = table_of_droplets, toc = data@assays$RNA$counts) # estimateSoup()

      graphclust_dir <- file.path(samples_dir, "analysis", "clustering")
      graphclust_dir <- list.files(graphclust_dir, full.names = T)[grep(pattern = "graphclust$", list.files(graphclust_dir))]
      tenx_graphclust <- read.csv(file.path(graphclust_dir, "clusters.csv"))
      tenx_clusters <- tenx_graphclust$Cluster[tenx_graphclust$Barcode %in% colnames(data@assays$RNA$counts)]
      sc <- SoupX::setClusters(sc, tenx_clusters)

      png(filename = file.path(sample_path, "quality_control", paste0("SoupX_contamination_", sample_name, ".png")))
      sc <- SoupX::autoEstCont(sc, doPlot = TRUE) # if fails, contamination rate default: 10% -> 0.1, or determine gene (sets) to estimate fraction
      dev.off()


      out <- SoupX::adjustCounts(sc, roundToInt = FALSE) # TODO set roundToInt = T is downstream processing algortihms need integers
      colnames(out) <- colnames(data@assays$RNA$counts)

      data@misc[["tenx_counts"]] <- data@assays$RNA$counts
      data@assays$RNA$counts <- out
      data@misc[["SoupX_contamination_percentage"]] <- sc$fit$rhoEst * 100
      return(data)
    }
    data <- clean_ambient_RNA(data, samples_dir, sample_name)
  }

  ## QUALITY CONTROL
  # calculate percentage of all counts belonging to mitochondrial (^MT-) DNA, for filtering
  mt_features <- grep("^MT-", rownames(data@assays$RNA), value = TRUE)
  data <- Seurat::PercentageFeatureSet(data, features = mt_features, col.name = "percent.mt")
  # Visualize quality control metrics
  png(file.path(sample_path, "quality_control", paste0("QC_nFeat_nCount_percent.mt_", sample_name, ".png")))
  plot(Seurat::VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("#85d0f5", "#2b2f70")))
  dev.off()
  plot1 <- Seurat::FeatureScatter(data, feature1 = "percent.mt", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
  plot2 <- Seurat::FeatureScatter(data, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", cols = c("#85d0f5", "#2b2f70"))
  png(file.path(sample_path, "quality_control", paste0("QC_feature-scatter_", sample_name, ".png")))
  plot(plot1 + plot2)
  dev.off()

  ## SCTransform-v2 (default in Seurat V5) replaces NormalizeData + FindVariableFeatures + ScaleData
  ### return.only.var.genes = FALSE: https://github.com/satijalab/seurat/issues/4896
  ## default assay set to "SCT"
  options(future.globals.maxSize = 8000 * 1024^2)
  data <- Seurat::SCTransform(data, vst.flavor = "v2", vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes = FALSE)


  # plot variable features, label top 10
  plot1 <- Seurat::VariableFeaturePlot(data, cols = c("#85d0f5", "#2b2f70"), selection.method = 'sctransform')
  plot2 <- Seurat::LabelPoints(plot = plot1, points = head(SeuratObject::VariableFeatures(data), 10), repel = TRUE)
  png(file.path(sample_path, "quality_control", paste0("Feature-selection_variable-genes_", sample_name, ".png")))
  plot(plot2)
  dev.off()

  # TODO check: https://github.com/satijalab/seurat/discussions/4259
  if (run_cell_cycle_regression) {
    dir.create(file.path(sample_path, 'Cell_Cycle'), recursive = T)

    # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.
    # We can segregate this list into markers of G2/M phase and markers of S phase.
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    # First, we assign each cell a score, based on its expression of G2/M and S phase markers.
    # These marker sets should be anticorrelated in their expression levels,
    # cells expressing neither are likely not cycling and in G1 phase.
    data <- Seurat::CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    # Visualize the distribution of cell cycle markers across
    png(file.path("Cell_Cycle", paste0("Cell_cycle_markers_ridgeplot_", sample_name, ".png")))
    plot(Seurat::RidgePlot(data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2))
    dev.off()
    data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(object = data), npcs = 20, verbose = FALSE)
    png(file.path("Cell_Cycle", paste0("Cell_cycle_PCA_dimplot-all-features_", sample_name, ".png")))
    plot(DimPlot(data, reduction = "pca", label = TRUE))
    dev.off()
    # Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
    data <- Seurat::RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
    png(file.path("Cell_Cycle", paste0("Cell_cycle_PCA_dimplot-s-and-g2m-features_", sample_name, ".png")))
    plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
    dev.off()
    # When running a PCA on only cell cycle genes after regression, cells no longer separate by cell-cycle phase
    data <- Seurat::ScaleData(data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(data))
    data <- Seurat::RunPCA(data, features = c(s.genes, g2m.genes), npcs = 20, verbose = FALSE)
    png(file.path("Cell_Cycle", paste0("Cell_cycle_PCA_dimplot_after-regression-s-and-g2m-features_", sample_name, ".png")))
    plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
    dev.off()
    data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(data), npcs = 20, nfeatures.print = 10, verbose = FALSE)
    png(file.path("Cell_Cycle", paste0("Cell_cycle_PCA_dimplot_after-regression-variable-features_", sample_name, ".png")))
    plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
    dev.off()
    # set sample identity to data, instead of cell cycle identity
    data@active.ident <- data@meta.data$old.ident
  }

  # run Principal Component Analysis as linear dimension reduction
  data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(object = data), npcs = 50, verbose = FALSE)
  png(file.path(sample_path, "Principal_Component_Analysis", paste0("/PCA-scores_", sample_name, ".png")))
  plot(Seurat::DimPlot(data, reduction = "pca", label = TRUE))
  dev.off()
  png(file.path(sample_path, "Principal_Component_Analysis", paste0("PCA-loadings_", sample_name, ".png")))
  plot(Seurat::VizDimLoadings(data, dims = 1:2, reduction = "pca"))
  dev.off()
  png(file.path(sample_path, "Principal_Component_Analysis", paste0("PCA-genes-heatmap_", sample_name, ".png")))
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
  png(file.path(sample_path, "Principal_Component_Analysis", paste0("PCA-variance_", sample_name, ".png")))
  p <- ggplot2::ggplot(data = longdf, ggplot2::aes(x=rep(1:length(varExplained), times=2), y = .data$value*100, fill = .data$variable, color = .data$variable)) +
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
  png(file.path(sample_path, "Principal_Component_Analysis", paste0("PCA_elbow-plot_", sample_name, ".png")))
  plot(Seurat::ElbowPlot(data))
  dev.off()

  # cell clustering: Levine2015 - Xu & Su2015
  # construct nearest neighbor graph for clustering
  data <- Seurat::FindNeighbors(data, dims = 1:30)
  # use Leiden algorithm for clustering (https://www.nature.com/articles/s41598-019-41695-z/)
  ## method = "igraph" (for large datasets when using Leiden algorithm)
  data <- Seurat::FindClusters(data, resolution = 0.8, algorithm = 1)
  # visualize clustering with Uniform Manifold Projection Approximation (UMAP) as non-linear dimension reduction
  data <- Seurat::RunUMAP(data, reduction = "pca", dims = 1:30)
  png(file.path(sample_path, paste0("UMAP_unsupervised_", sample_name, ".png")))
  plot(Seurat::DimPlot(data, reduction = "umap", label = TRUE))
  dev.off()





  # FindAllMarkers for indiviual sample --> later DEG for integrated analysis (then explain test.use: Wilcox rst)

  # run FindAllMarkers for Differential Gene Expression analysis
  ## finds markers for every cluster compared to all remaining cells
  ### report only up-regulated genes as down-regulated genes represent all other cells/clusters here
  data.markers <- Seurat::FindAllMarkers(data, assay = "SCT", only.pos = TRUE, min.pct = 0.1)
  utils::write.csv2(data.markers, file = file.path(sample_path, "DE_analysis", paste0("marker-list_", sample_name, ".csv")))

  # select top gene per cluster for exploration
  topn <- data.markers %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::top_n(n = 1, wt = .data$avg_log2FC) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data$gene)
  features_of_interest[["topn-features"]] <- topn

  plot_DEG <- function(data, features, name) {
    dir.create(file.path(sample_path, "DE_analysis", name), recursive = T)
    dir.create(file.path(sample_path, "DE_analysis", name, "Feature"), recursive = T)

    # plot feature expression, if available in Seurat
    for (i in seq_along(features)) {
      tryCatch({
        p <- Seurat::FeaturePlot(data, features = features[i])
        ggplot2::ggsave(file=file.path(sample_path, "DE_analysis", name ,"Feature", paste0(features[i], ".png")), width = 30, height = 20, units = "cm")
      },
      error=function(e) {
        message(features[i], ' plot is skipped, as gene was not found with FetchData')
      })
    }

    # expression plots
    if (any(features %in% rownames(data@assays$SCT@scale.data))) {
      p <- Seurat::DoHeatmap(data, features = features) + Seurat::NoLegend()
      ggplot2::ggsave(file = file.path(sample_path, "DE_analysis", name, paste0("heatmap_", name, "_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    }
    p <- Seurat::FeaturePlot(data, features = features)
    ggplot2::ggsave(file=file.path(sample_path, "DE_analysis", name, paste0("feature-plot_", name, "_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    p <- Seurat::VlnPlot(data, features = features)
    ggplot2::ggsave(file = file.path(sample_path, "DE_analysis", name, paste0("violin-plot_ ", name, "_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    p <- Seurat::RidgePlot(data, features = features, ncol = 3)
    ggplot2::ggsave(file = file.path(sample_path, "DE_analysis", name, paste0("ridge-plot_", name, "_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    # dotplot with custom labels
    cell.num <- table(SeuratObject::Idents(data))
    cluster.labels = paste(names(cell.num), paste0("(", round(cell.num/sum(cell.num), 2)*100, "%, n = ", cell.num, ")"))
    levels(SeuratObject::Idents(data)) <- cluster.labels
    p <- Seurat::DotPlot(data, features = features) + Seurat::RotatedAxis() + Seurat::WhiteBackground()
    ggplot2::ggsave(file = file.path(sample_path, "DE_analysis", name, paste0("dot-plot_", name, "_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    levels(SeuratObject::Idents(data)) <- sapply(stringr::str_split(levels(SeuratObject::Idents(data)), " "), "[[", 1)
  }

  for (feat_name in names(features_of_interest)) {
    plot_DEG(data = data, features = features_of_interest[[feat_name]], name = feat_name)
  }

  # plot heatmap for topn genes per cluster
  heatmap_features <- data.markers %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::top_n(n = 8, wt = .data$avg_log2FC) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data$gene)
  p <- Seurat::DoHeatmap(data, features = heatmap_features) + Seurat::NoLegend()
  ggplot2::ggsave(file = file.path(sample_path, paste0("DEG-analysis_big-heatmap_", sample_name, ".png")), width = 30, height = 20, units = "cm")

  # save data
  qs::qsave(data, file = file.path(sample_path, paste0(sample_name, ".qs")), preset = 'custom', algorithm = "zstd_stream", compress_level = 4, shuffle_control = 15, nthreads = 1)
}
