#' Run DEA on processed Seurat object
#'
#' Perform differential expression analysis on processed and annotated Seurat object.
#'
#' @param sample_name sample name, string
#' @param qs_file character string file path to .qs file of processed Seurat object
#' @param output_dir output directory for plots, string
#' @param sample_celltype_DEA list of sample_celltype comparisons, as exampled
#' @param features_of_interest marker features to plot as violins and dots per DE comparison
#'
#' @export
#'
#' @description Changed following functions in mauritsunkel/EMC-SKlab-seurat:
#' Seurat::FindMarkers.default
#' Seurat:::WilcoxDETest
#' Seurat::Foldchange.default
#'
#' For benchmarking DE methods see: Junttilla et al 2022 (summarizes both Zimmerman et al 2021/2022 and Squair et al 2021!)
#' If you have replicate samples: pseudobulk (ROTS with sum > mean) > mixed models (best: MAST_RE) > naive models > latent variable models
#' If no replicate samples: use Seurat default Wilcoxon rank-sum test (with Presto installed for speed) on SCT normalized assay (@data slot)
#' -  set: recorrect_umi = FALSE, when running on a subset
#' AUROC suggests that p-values are generally in order for all methods, however large amount of false positives for naive/latent variable models.
#' Don't use ComBat for batch normalization for DE. If taking into account batch effects for DE use as a latent variable.
#' Squair et al 2021: "the central principle underlying valid DE analysis is the ability of statistical methods to account for the intrinsic variability of biological replicates".
#'
#' @examplesIf FALSE
#' differential_expression_analysis(
#'   sample_name = "T1",
#'   qs_file = file.path("EMC-SKlab-scRNAseq", "results", "T1.qs"),
#'   output_dir = file.path("EMC-SKlab-scRNAseq", "results", 'integrated', 'sample_name')
#' )
#'
#' sample_celltype_DEA <- list(
#' c("SampleA_CelltypeA", "SampleB_CelltypeA"),
#' c("SampleA_CelltypeB", "SampleB_CelltypeB"))
differential_expression_analysis <- function(
    sample_name, qs_file, output_dir,
    sample_celltype_DEA = NULL,
    features_of_interest = NULL) {
  library(Seurat) # added because of error
  # Error: package or namespace load failed for ‘Seurat’ in .doLoadActions(where, attach):
  #   error in load action .__A__.1 for package RcppAnnoy: loadModule(module = "AnnoyAngular", what = TRUE, env = ns, loadNow = TRUE): Unable to load module "AnnoyAngular": attempt to apply non-function
  # Error in .requirePackage(package) :
  #   unable to find required package ‘Seurat’

  ## create output directories
  output_dir <- file.path(output_dir, 'DEA')
  dir.create(output_dir, recursive = TRUE)
  ## read data
  integrated <- qs::qread(qs_file)

  if (!length(unique(integrated$orig.ident)) > 1) {
    message("Need multiple samples in Seurat object to perform sample-level DEA")
  } else {
    # plot scatters with correlation for AggregatedExpression of unique sample pairs
    aggregate_data <- Seurat::AggregateExpression(integrated, group.by = "orig.ident", assays = "SCT", return.seurat = TRUE)
    # get unique sample-pairs
    sample_combinations <- combn(names(aggregate_data$orig.ident), 2, simplify = F)
    scatter_plots <- patchwork::wrap_plots(lapply(sample_combinations, function(x) {
      Seurat::CellScatter(aggregate_data, x[1], x[2])
    }))
    ggplot2::ggsave(plot = scatter_plots, file = file.path(output_dir, 'pseudobulk_scatter_plots.png'), width = 30, height = 20, units = "cm")
  }

  ## sample(s)-sample(s) & sample(s)-celltype(s)-level DE
  if (is.null(sample_celltype_DEA)) {
    message("To perform sample(s)-sample(s) and sample(s)-celltype(s) DE comparisons, provide sample_celltype_DEA (named) list with comparisons, see workflow for instructions")
  } else {
    for (meta_name in names(sample_celltype_DEA)) {
      message(meta_name)
      if (meta_name == "orig.ident") {
        SeuratObject::Idents(integrated) <- integrated$orig.ident
      } else {
        idents <- paste(integrated$orig.ident, integrated@meta.data[, meta_name], sep = "_")
        SeuratObject::Idents(integrated) <- idents
      }

      for (i in seq_along(sample_celltype_DEA[[meta_name]])) {
        ref_ident_DEA <- sample_celltype_DEA[[meta_name]][[i]]$ref
        vs_ident_DEA <- sample_celltype_DEA[[meta_name]][[i]]$vs
        ref_ident <- c()
        for (id in ref_ident_DEA) {
          ref_ident <- c(ref_ident, grep(id, levels(SeuratObject::Idents(integrated)), value = TRUE))
        }
        if (vs_ident_DEA == 'rest') {
          vs_ident <- 'rest'
        } else {
          vs_ident <- c()
          for (id in vs_ident_DEA) {
            vs_ident <- c(vs_ident, grep(id, levels(SeuratObject::Idents(integrated)), value = TRUE))
          }
        }

        comp_name <- names(sample_celltype_DEA[[meta_name]])[i]
        if (comp_name == "name") {
          comp_name <- paste0(ref_ident_DEA, "_vs_", vs_ident_DEA)
        }
        meta_dirname <- switch(
          meta_name,
          "kriegstein.seurat.custom.clusters.mean" = "Kriegstein.Seurat",
          "orig.ident" = "Sample-level",
          "mapmycells_supercluster" = "MapMyCells"
        )
        DE_output_dir <- file.path(output_dir, meta_dirname, comp_name)
        dir.create(DE_output_dir, recursive = TRUE)

        message("DE: ", comp_name)
        message("ref_idents: ", ref_ident)
        message("vs_idents: ", vs_ident)
        DE_EnhancedVolcano(integrated, ref_ident, vs_ident, DE_output_dir, comp_name)
        if (!is.null(features_of_interest)) {
          DE_MarkerExpression(integrated, features_of_interest, idents = c(ref_ident, vs_ident), DE_output_dir)
        }
      }
    }
  }
}

#' Plot marker features as violins and dots per DE comparison
#'
#' @param seurat_object SO with integrated samples
#' @param features_of_interest List with named vectors with features (genes)
#' @param idents for which DE idents to plot
#' @param output_dir output directory for plots, string
DE_MarkerExpression <- function(seurat_object, features_of_interest, idents, output_dir) {
  output_dir_split <- strsplit(output_dir, "/")[[1]]
  output_filename <- paste0(output_dir_split[[length(output_dir_split)]])

  for (feat_name in (names(features_of_interest))) {
    message("DE Violin & DotPLot: ", feat_name)
    output_dir_path <- file.path(output_dir, "markers", feat_name)
    dir.create(output_dir_path, recursive = TRUE)

    features <- features_of_interest[[feat_name]]

    Seurat::VlnPlot(seurat_object, features, idents = idents, assay = "SCT", same.y.lims = T)
    ggplot2::ggsave(file.path(output_dir_path, paste0(output_filename, "_ViolinPlot.png")), width = 30, height = 20, units = "cm")
    Seurat::DotPlot(seurat_object, features, idents = idents, assay = "SCT")
    ggplot2::ggsave(file.path(output_dir_path, paste0(output_filename, "_DotPlot.png")), width = 30, height = 20, units = "cm")
  }
}

#' Perform DE with Seurat FindMarkers() and plot EnhancedVolcano
#'
#' @param seurat_object integrated Seurat object. Containing ref_ident and vs_ident in Idents(seurat_object)
#' @param ref_ident reference sample name, pct.1 in Seurat DE result
#' @param vs_ident versus sample name, pct.2 in Seurat DE result
#' @param directory to save plot in, filename is based on reference and versus sample
#' @param comp_name used to handle filenaming and EnhancedVolcano plot title
DE_EnhancedVolcano <- function(seurat_object, ref_ident, vs_ident, DE_output_dir, comp_name) {
  DE_res <- Seurat::FindMarkers(
    seurat_object,
    assay = "SCT",
    ident.1 = ref_ident,
    ident.2 = if (vs_ident != 'rest') vs_ident else NULL,
    only.pos = if (vs_ident == 'rest') TRUE else FALSE,
    verbose = TRUE,
    logfc.threshold = 0,
    min.pct = 0)

  ## sort by average log2 fold-change
  DE_res <- DE_res %>% dplyr::arrange(dplyr::desc(avg_log2FC))
  ## filter by p-val-adj (Bonferroni corrected)
  DE_res_adj <- DE_res[DE_res$p_val_adj < 0.05,]
  ## write raw and p-val-adj filtered sample-level DE
  if (comp_name == "name") {
    filename <- file.path(DE_output_dir, paste0("1=", ref_ident, "_vs_2=", vs_ident, ".xlsx"))
  } else {
    split_name <- strsplit(comp_name, "_vs_")
    ref_ident_name <- split_name[[1]][1]
    vs_ident_name <- split_name[[1]][2]
    filename <- file.path(DE_output_dir, paste0("1=", ref_ident_name, "_vs_2=", vs_ident_name, ".xlsx"))
  }
  openxlsx::write.xlsx(x = DE_res, file = filename, row.names = TRUE)
  openxlsx::write.xlsx(x = DE_res_adj, file = sub(".xlsx$", "_adj.xlsx", filename), row.names = TRUE)

  ## plot EnhancedVolcano per sample DE
  plotEnhancedVolcano(
    seurat_object, DE_res, ref_ident, vs_ident,
    filedir = DE_output_dir, comp_name
  )
}
