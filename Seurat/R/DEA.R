library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)
library(DESeq2)
library(MAST)

# set sample and set and get data
DE_method = 'DESeq2'

sample_name <- "BL_A + BL_C"
rds.file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/BL_A + BL_C.rds"
integrated <- readRDS(rds.file)

# set working directory and create ouput directories
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/DE_analysis/")
dir.create(work_dir, recursive = T)
setwd(work_dir)
dir.create(paste0(work_dir, "../GSE_analysis/"))
dir.create(paste0(work_dir, "markers/"))
if (length(table(integrated$orig.ident)) != 1) {
  dir.create(paste0(work_dir, "sample_markers/"))
}

message(start_time)
### END USER PARAMETERS







## perform sample level comparison for integration
# check if multiple samples AND all samples for comparison have more then 3 cells
if (length(table(integrated$orig.ident)) != 1 && !any(table(integrated$orig.ident) < 3)) {
  # set idents to compare cells at sample level instead of cluster level
  Idents(integrated) <- integrated$orig.ident

  # TODO if 2+ samples then this hard code will error
  # create sample marker dataframes to count n cells used in comparisons
  sample_markers_columns <- c(paste0("n_cells_", names(table(integrated$orig.ident))[1]),
                              paste0("n_cells_", names(table(integrated$orig.ident))[2]))
  sample_markers_df <- data.frame(matrix(nrow = 0, ncol = length(sample_markers_columns)))
  colnames(sample_markers_df) <- sample_markers_columns
  sample_markers_df[nrow(sample_markers_df) + 1,] = c(table(integrated$orig.ident)[1],
                                                      table(integrated$orig.ident)[2])

  # write n cells for comparison to CSV files
  write.csv2(sample_markers_df, file = "sample_markers/n_cells_for_comparison_m.csv", row.names = FALSE)

  # TODO need a user put name instead of: names(table(integrated$orig.ident))[1]
  # get sample markers (note: p_val_adj based on Bonferroni correction using ALL genes)
  sample_markers <- FindMarkers(integrated, assay = "SCT", ident.1 = names(table(integrated$orig.ident))[1], only.pos = FALSE, verbose = T,
                                logfc.threshold = 0, min.pct = 0, test.use = DE_method)

  # TODO check if filter on feature not in dataframe filters all rows
  # # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
  # sample_markers <- sample_markers[!(sample_markers$p_val_adj > 0.05 & sample_markers$nz_p_val_adj > 0.05),]

  # order by avg_log2FC
  sample_markers_pval_adj <- sample_markers %>% arrange(desc(avg_log2FC)) # DEPRECATED: filter(pct > 0.1)

  write.csv2(sample_markers_pval_adj, file = paste0("sample_markers/method=", DE_method, "-pct1=", names(table(integrated$orig.ident))[1], "-pct2=", names(table(integrated$orig.ident))[2], " - (nz-)p-val st 0.05.csv"), row.names = TRUE)
}
