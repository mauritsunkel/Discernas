### start initialization
library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)

### USER PARAMETERS
# read an integrated saved RDS file
sample_name <- "BL_N + BL_C"
integrated <- readRDS(paste0("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden - Cellcycle/integrated/", sample_name, "/integrated.rds"))

# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
### END USER PARAMETERS

work_dir <- paste0(work_dir, 'results/')
dir.create(work_dir)
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- paste0(work_dir, start_time, '/')
dir.create(work_dir)
work_dir <- paste0(work_dir, 'integrated/')
dir.create(work_dir)
work_dir <- paste0(work_dir, sample_name, "/")
dir.create(work_dir)
work_dir <- paste0(work_dir, 'DE_analysis/')
dir.create(work_dir)
setwd(work_dir)
### end initialization ###








dir.create(paste0(work_dir, "markers/"))
dir.create(paste0(work_dir, "conserved_markers/"))
dir.create(paste0(work_dir, "condition_markers/"))

cluster_ids <- levels(integrated$seurat_clusters)
for (i in cluster_ids) {
  # create markers for integrated data for each cluster vs all other clusters
  markers <- FindMarkers(integrated, ident.1 = i, only.pos = TRUE, verbose = T)
  markers_top10 <- markers %>% slice_head(n = 10)
  markers_top100 <- markers %>% slice_head(n = 100)
  write.csv2(markers, file = paste0("markers/all_cluster", i, "_m.csv"))
  write.csv2(markers_top10, file = paste0("markers/top10_cluster", i, "_m.csv"))
  write.csv2(markers_top100, file = paste0("markers/top100_cluster", i, "_m.csv"))

  # create markers conserved between groups (conditions) for integrated data for each cluster vs all other clusters
  conserved_markers <- FindConservedMarkers(integrated, ident.1 = i, only.pos = TRUE,
                                            grouping.var = "orig.ident", verbose = T)
  conserved_markers_top10 <- conserved_markers %>% slice_head(n = 10)
  conserved_markers_top100 <- conserved_markers %>% slice_head(n = 100)
  write.csv2(conserved_markers, file = paste0("conserved_markers/all_cluster", i, "_cm.csv"))
  write.csv2(conserved_markers_top10, file = paste0("conserved_markers/top10_cluster", i, "_cm.csv"))
  write.csv2(conserved_markers_top100, file = paste0("conserved_markers/top100_cluster", i, "_cm.csv"))

  # create condition markers for integrated data within each cluster between each condition
  ## DEV NOTE: this is not pairwise if more than 2 conditions are integrated at the same time
  subset <- subset(integrated, seurat_clusters == i)
  # change cluster identity to original identity to find markers between conditions
  Idents(subset) <- subset$orig.ident
  condition_markers <- FindMarkers(subset, ident.1 = "BL_C", verbose = T)
  condition_markers_top10 <- condition_markers %>% slice_head(n = 10)
  condition_markers_top100 <- condition_markers %>% slice_head(n = 100)
  write.csv2(condition_markers, file = paste0("condition_markers/all_cluster", i, "_cm.csv"))
  write.csv2(condition_markers_top10, file = paste0("condition_markers/top10_cluster", i, "_cm.csv"))
  write.csv2(condition_markers_top100, file = paste0("condition_markers/top100_cluster", i, "_cm.csv"))

  # if want to assign each table to its own variable, use assign() and get()
  # assign(paste0("cluster", i, "_markers"), markers)
  ## get(paste0("cluster", i, "_markers"))
}
#cleanup environment
rm("markers", "markers_top10", "markers_top100",
   "conserved_markers", "conserved_markers_top10", "conserved_markers_top100",
   "condition_markers", "condition_markers_top10", "condition_markers_top100",
   "cluster_ids")
