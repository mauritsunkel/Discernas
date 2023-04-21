#' Run DEA on processed Seurat object
#'
#' Perform differential expression analysis on processed and annotated Seurat object.
#'
#' @param sample_name sample name, string
#' @param rds_file character string file path to .rds file of processed Seurat object
#' @param output_dir output directory for plots, string
#'
#' @export
#'
#' @description Changed following functions in mauritsunkel/EMC-SKlab-seurat:
#' Seurat::FindMarkers.default
#' Seurat:::WilcoxDETest
#' Seurat::Foldchange.default
#'
#' @examplesIf FALSE
#' differential_expression_analysis(
#'   sample_name = "T1",
#'   rds_file = file.path("EMC-SKlab-scRNAseq", "results", "T1.rds"),
#'   file.path("EMC-SKlab-scRNAseq", "results")
#' )
differential_expression_analysis <- function(sample_name, rds_file, output_dir) {
  # set working directory and create output directories
  dir.create(paste0(output_dir, 'results/integrated/', sample_name, "/DE_analysis/"), recursive = T)
  output_dir <- paste0(output_dir, 'results/integrated/', sample_name, "/DE_analysis/")
  setwd(output_dir)
  dir.create(paste0(output_dir, "../GSE_analysis/"))
  dir.create(paste0(output_dir, "markers/"))
  integrated <- readRDS(rds_file)
  if (length(table(integrated$orig.ident)) != 1) {
    dir.create(paste0(output_dir, "sample_markers/"))
    dir.create(paste0(output_dir, "conserved_markers/"))
    dir.create(paste0(output_dir, "condition_markers/"))
  }

  ## perform sample level comparison for integration
  # check if all samples for comparison have more then 3 cells
  if (!any(table(integrated$orig.ident) < 3)) {
    # check if multiple samples
    if (length(table(integrated$orig.ident)) != 1) {
      # set idents to compare cells at sample level instead of cluster level
      SeuratObject::Idents(integrated) <- integrated$orig.ident

      # create sample marker dataframes to count n cells used in comparisons
      sample_markers_columns <- c(paste0("n_cells_", names(table(integrated$orig.ident))[1]),
                                  paste0("n_cells_", names(table(integrated$orig.ident))[2]))
      sample_markers_df <- data.frame(matrix(nrow = 0, ncol = length(sample_markers_columns)))
      colnames(sample_markers_df) <- sample_markers_columns
      sample_markers_df[nrow(sample_markers_df) + 1,] = c(table(integrated$orig.ident)[1],
                                                          table(integrated$orig.ident)[2])
      # write n cells for comparison to CSV files
      write.csv2(sample_markers_df, file = "sample_markers/n_cells_for_comparison_m.csv", row.names = FALSE)

      # get sample markers (note: p_val_adj based on Bonferroni correction using ALL genes)
      sample_markers <- Seurat::FindMarkers(integrated, assay = "SCT", ident.1 = names(table(integrated$orig.ident))[1], only.pos = FALSE, verbose = T,
                                    logfc.threshold = 0, min.pct = 0)

      # set pct variable based on BL_C orig.identity index
      pct <- if(names(table(integrated$orig.ident))[1] == "BL_C") "pct.1" else "pct.2"

      # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
      sample_markers <- sample_markers[!(sample_markers$p_val_adj > 0.05 & sample_markers$nz_p_val_adj > 0.05),]
      # order by avg_log2FC
      sample_markers_pval_adj <- sample_markers %>% dplyr::arrange(dplyr::desc(avg_log2FC)) # DEPRECATED: filter(pct > 0.1)

      write.csv2(sample_markers_pval_adj, file = paste0("sample_markers/pct1=", names(table(integrated$orig.ident))[1], "-pct2=", names(table(integrated$orig.ident))[2], " - (nz-)p-val st 0.05.csv"), row.names = TRUE)

      # add sample markers and n cells count as miscellaneous data to Seurat object
      SeuratObject::Misc(object = integrated, slot = paste0("DEG.sample_markers")) <- sample_markers
      SeuratObject::Misc(object = integrated, slot = paste0("DEG.sample_markers_n")) <- sample_markers_df
    }
  }

  # set idents back to cluster level for downstream comparisons
  SeuratObject::Idents(integrated) <- integrated$seurat_clusters

  # markers
  markers_columns <- c('cluster_ID', 'n_cells_cluster', 'n_all_other_cells')
  markers_df <- data.frame(matrix(nrow = 0, ncol = length(markers_columns)))
  colnames(markers_df) <- markers_columns
  # conserved markers
  conserved_markers_columns <- c('cluster_ID',
                                 paste0('n_cells_cluster_identity_', names(table(integrated$orig.ident))[1]),
                                 paste0('n_cells_cluster_identity_', names(table(integrated$orig.ident))[2]),
                                 paste0('n_all_other_cells_identity_', names(table(integrated$orig.ident))[1]),
                                 paste0('n_all_other_cells_identity_', names(table(integrated$orig.ident))[2]))
  conserved_markers_df <- data.frame(matrix(nrow = 0, ncol = length(conserved_markers_columns)))
  colnames(conserved_markers_df) <- conserved_markers_columns
  # condition markers
  condition_markers_columns <- c('cluster_id',
                                 paste0('n_cells_', names(table(integrated$orig.ident))[1]),
                                 paste0('n_cells_', names(table(integrated$orig.ident))[2]))
  condition_markers_df <- data.frame(matrix(nrow = 0, ncol = length(condition_markers_columns)))
  colnames(condition_markers_df) <- condition_markers_columns

  # get cluster IDs to loop over
  cluster_ids <- levels(integrated$seurat_clusters)
  for (i in cluster_ids) {
    message('Cluster ID:', i)

    # add amount of cells used for markers comparison to df
    markers_df[nrow(markers_df) + 1,] = c(i,
                                          table(integrated$seurat_clusters)[i],
                                          sum(table(integrated$seurat_clusters)[-as.integer(i)]))
    # write n cells for comparison to CSV files
    write.csv2(markers_df, file = "markers/n_cells_for_comparison_m.csv", row.names = FALSE)
    # check more than 2 cells in ident (cluster) before comparison
    if (table(integrated$seurat_clusters)[i] < 3) {
      message("For markers, skipping ident (cluster) ", i, " comparison because < 3 cells")
    } else {
      # create markers for integrated data for each cluster vs all other clusters
      markers <- Seurat::FindMarkers(integrated, assay = "SCT", ident.1 = i, only.pos = FALSE, verbose = T)
      # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
      markers <- markers[!(markers$p_val_adj > 0.05 & markers$nz_p_val_adj > 0.05),]

      write.csv2(markers, file = paste0("markers/all_cluster", i, "_m.csv"))

      # add as miscellaneous data to Seurat object
      SeuratObject::Misc(object = integrated, slot = paste0("DEG.markers")) <- markers
      SeuratObject::Misc(object = integrated, slot = paste0("DEG.markers_n")) <- markers_df
      message("wrote markers")
    }

    # check if multiple samples
    if (length(table(integrated$orig.ident)) != 1) {
      # add amount of cells used for conserved_markers comparison to df
      df <- data.frame('orig.ident' = integrated$orig.ident, 'seurat_clusters' = integrated$seurat_clusters)
      # condition 1 & match cluster
      cm_val1 <- nrow(df %>% dplyr::filter(orig.ident == names(table(integrated$orig.ident))[1] & seurat_clusters == i))
      # condition 2 & match cluster
      cm_val2 <- nrow(df %>% dplyr::filter(orig.ident == names(table(integrated$orig.ident))[2] & seurat_clusters == i))
      # condition 1 & no match cluster
      cm_val3 <- nrow(df %>% dplyr::filter(orig.ident == names(table(integrated$orig.ident))[1] & seurat_clusters != i))
      # condition 2 & no match cluster
      cm_val4 <- nrow(df %>% dplyr::filter(orig.ident == names(table(integrated$orig.ident))[2] & seurat_clusters != i))
      conserved_markers_df[nrow(conserved_markers_df) + 1,] = c(i, cm_val1, cm_val2, cm_val3, cm_val4)
      # write n cells for comparison to CSV files
      write.csv2(conserved_markers_df, file = "conserved_markers/n_cells_for_comparison_cm.csv", row.names = FALSE)
      # check more than 2 cells in ident (cluster) for each group before comparison
      if (any(c(cm_val1, cm_val2, cm_val3, cm_val4) < 3)) {
        message("For conserved markers, skipping ident (cluster) ", i, " comparison because < 3 cells")
      } else {
        # create markers conserved between groups (conditions) for integrated data for each cluster vs all other clusters
        conserved_markers <- Seurat::FindConservedMarkers(integrated, assay = "SCT", ident.1 = i, only.pos = FALSE,
                                                  grouping.var = "orig.ident", verbose = T)


        # filters rows (genes) for both compared samples if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
        conserved_markers <- conserved_markers[!(conserved_markers[, paste0(names(table(integrated$orig.ident))[1], "_p_val_adj")] > 0.05 & conserved_markers[, paste0(names(table(integrated$orig.ident))[1], "_nz_p_val_adj")] > 0.05),]
        conserved_markers <- conserved_markers[!(conserved_markers[, paste0(names(table(integrated$orig.ident))[2], "_p_val_adj")] > 0.05 & conserved_markers[, paste0(names(table(integrated$orig.ident))[2], "_nz_p_val_adj")] > 0.05),]

        write.csv2(conserved_markers, file = paste0("conserved_markers/all_cluster", i, "_cm.csv"))

        # add as miscellaneous data to Seurat object
        SeuratObject::Misc(object = integrated, slot = paste0("DEG.conserved_markers")) <- conserved_markers
        SeuratObject::Misc(object = integrated, slot = paste0("DEG.conserved_markers_n")) <- conserved_markers_df
        message("wrote conserved markers")
      }

      # create condition markers for integrated data within each cluster between each condition
      subset <- subset(integrated, seurat_clusters == i)
      # change cluster identity to original identity to find markers between conditions
      SeuratObject::Idents(subset) <- subset$orig.ident
      # check if subset contains cells for at least 2 conditions/samples for comparison
      if (length(names(table(subset$orig.ident))) == 1) {
        message('In cluster ', i, ' only cells for condition/sample ',
                names(table(subset$orig.ident)), ' were found, cannot create condition markers for this cluster.')
        next
      }
      # add amount of cells used for condition_markers comparison to df
      condition_markers_df[nrow(condition_markers_df) + 1,] = c(i, table(subset$orig.ident)[1], table(subset$orig.ident)[2])
      # write n cells for comparison to CSV files
      write.csv2(condition_markers_df, file = "condition_markers/n_cells_for_comparison.csv", row.names = FALSE)
      # check more than 2 cells in subset ident (cluster) before comparison
      if (any(c(table(subset$orig.ident)[1], table(subset$orig.ident)[2]) < 3)) {
        message("For condition markers, skipping ident (cluster) ", i, " comparison because < 3 cells")
      } else {
        # create condition_markers for subset data for within each cluster to compare conditions
        condition_markers <- Seurat::FindMarkers(subset, assay = "SCT", recorrect_umi = FALSE, ident.1 = "BL_C", verbose = T, only.pos = FALSE)
        # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
        condition_markers <- condition_markers[!(condition_markers$p_val_adj > 0.05 & condition_markers$nz_p_val_adj > 0.05),]

        write.csv2(condition_markers, file = paste0("condition_markers/all_cluster", i, ".csv"))

        # add as miscellaneous data to Seurat object
        SeuratObject::Misc(object = integrated, slot = paste0("DEG.condition_markers")) <- condition_markers
        SeuratObject::Misc(object = integrated, slot = paste0("DEG.condition_markers_n")) <- condition_markers_df
        message("wrote condition markers")
      }
    }
  }

  # save data for possible adjustments
  saveRDS(integrated, file = rds_file)
}
