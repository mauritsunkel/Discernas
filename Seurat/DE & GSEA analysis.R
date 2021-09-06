### start initialization
library(Seurat)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)

# FGSEA package: vignette: http://127.0.0.1:31440/library/fgsea/doc/fgsea-tutorial.html
FGSEA_analysis <- function(markers, working_directory, marker_type, cluster) {
  library(biomaRt)
  library(fgsea)
  library(data.table)
  library(ggplot2)
  ## library(reactome.db) # not needed in test run

  dir.create(paste0(work_dir, "../GSEA_analysis/", marker_type, "/"))
  dir.create(paste0(work_dir, "../GSEA_analysis/", marker_type, "/", cluster, "/"))

  ## fix infinite values later by applying -log10 function
  markers$p_val[markers$p_val == 0] <- min(markers$p_val[markers$p_val != 0])

  ## calculate metric by FoldChangeSign and -LogPvalue
  markers$fcsign <- sign(markers$avg_log2FC)
  markers$logPval <- -log10(markers$p_val)



  ## create ranked vector
  fgsea_ranks <- markers$logPval/markers$fcsign

  ## get Entrez IDs by HGNC symbol to match gene names and provide a translation map
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  mapping <- getBM(
    attributes = c('entrezgene_id', 'hgnc_symbol'),
    filters = 'hgnc_symbol',
    values = rownames(markers),
    mart = hsmart
  )
  names(fgsea_ranks) <- match(rownames(markers), mapping$hgnc_symbol)

  ## get Reactome pathways by Entrez IDs
  pathways <- reactomePathways(names(fgsea_ranks))

  fgsea_results <- fgsea(pathways = pathways,
                         stats    = fgsea_ranks,
                         eps      = 0.0,
                         minSize  = 15,
                         maxSize  = 500)

  topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  png(filename=paste0(working_directory, "../GSEA_analysis/", marker_type, "/", cluster, "/overview_table.png"), width = 1600)
  plotGseaTable(pathways[topPathways], fgsea_ranks, fgsea_results,
                     gseaParam=0.5)
  dev.off()

  ## can try to collapse pathways if there are many seemlingly alike in the plot above
  collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < 0.01],
                                        pathways, fgsea_ranks)
  mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]

  ## check if mainPathways is empty (likely collapsedPathways is empty too)
  if (length(mainPathways) > 0) {
    png(filename=paste0(working_directory, "../GSEA_analysis/", marker_type, "/", cluster, "/collapsed_table.png"), width = 1600)
    p <- plotGseaTable(pathways[mainPathways], fgsea_ranks, fgsea_results,
                       gseaParam = 0.5)
    dev.off()
  }

  # TODO check why NES column is wrong and check what columns mean
  # TODO uncomment DE function code
  print(fgsea_results$NES)
  fwrite(fgsea_results, file=paste0(working_directory, "../GSEA_analysis/", marker_type, "/", cluster, "/overview.xls"), sep="\t", sep2=c("", ",", ""))

  for (i in seq_along(topPathways)) {
    # png(filename=paste0(working_directory, "GSEA/cluster_", cluster, "/enriched_", i, ".png"), width = 1600)
    p <- plotEnrichment(pathways[[topPathways[i]]],
                        fgsea_ranks) + labs(title=topPathways[[i]])
    ggsave(file = paste0(working_directory, "../GSEA_analysis/", marker_type, "/", cluster, "/enriched_", i, ".png"), width = 30, height = 20, units = "cm")
  }
}


### USER PARAMETERS
# read an integrated saved RDS file
sample_name <- "BL_A + BL_C"
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
dir.create(paste0(work_dir, "markers/positive/"))
dir.create(paste0(work_dir, "markers/negative/"))
dir.create(paste0(work_dir, "conserved_markers/"))
dir.create(paste0(work_dir, "conserved_markers/positive/"))
dir.create(paste0(work_dir, "conserved_markers/negative/"))
dir.create(paste0(work_dir, "condition_markers/"))
dir.create(paste0(work_dir, "condition_markers/positive/"))
dir.create(paste0(work_dir, "condition_markers/negative/"))
dir.create(paste0(work_dir, "../GSEA_analysis/"))

cluster_ids <- levels(integrated$seurat_clusters)
for (i in cluster_ids) {
  print(paste('Cluster ID:', i))
  ## create markers for integrated data for each cluster vs all other clusters
  markers <- FindMarkers(integrated, ident.1 = i, only.pos = FALSE, verbose = T)
  # pos_markers <- markers %>% filter(avg_log2FC > 0)
  # pos_markers_top10 <- markers %>% filter(avg_log2FC > 0) %>% slice_head(n = 10)
  # pos_markers_top100 <- markers %>% filter(avg_log2FC > 0) %>% slice_head(n = 100)
  # neg_markers <- markers %>% filter(avg_log2FC < 0)
  # neg_markers_top10 <- markers %>% filter(avg_log2FC < 0) %>% slice_head(n = 10)
  # neg_markers_top100 <- markers %>% filter(avg_log2FC < 0) %>% slice_head(n = 100)
  # write.csv2(pos_markers, file = paste0("markers/positive/all_cluster", i, "_m.csv"))
  # write.csv2(pos_markers_top10, file = paste0("markers/positive/top10_cluster", i, "_m.csv"))
  # write.csv2(pos_markers_top100, file = paste0("markers/positive/top100_cluster", i, "_m.csv"))
  # write.csv2(neg_markers, file = paste0("markers/negative/all_cluster", i, "_m.csv"))
  # write.csv2(neg_markers_top10, file = paste0("markers/negative/top10_cluster", i, "_m.csv"))
  # write.csv2(neg_markers_top100, file = paste0("markers/negative/top100_cluster", i, "_m.csv"))
  FGSEA_analysis(markers = markers, working_directory = work_dir, marker_type = 'markers', cluster = i)

  ## create markers conserved between groups (conditions) for integrated data for each cluster vs all other clusters
  # conserved_markers <- FindConservedMarkers(integrated, ident.1 = i, only.pos = FALSE,
  #                                           grouping.var = "orig.ident", verbose = T)
  # pos_conserved_markers <- conserved_markers %>% filter((.[[2]] > 0) & (.[[7]] > 0))
  # pos_conserved_markers_top10 <- conserved_markers %>% filter((.[[2]] > 0) & (.[[7]] > 0)) %>% slice_head(n = 10)
  # pos_conserved_markers_top100 <- conserved_markers %>% filter((.[[2]] > 0) & (.[[7]] > 0)) %>% slice_head(n = 100)
  # neg_conserved_markers <- conserved_markers %>% filter((.[[2]] < 0) | (.[[7]] < 0))
  # neg_conserved_markers_top10 <- conserved_markers %>% filter((.[[2]] < 0) | (.[[7]] < 0)) %>% slice_head(n = 10)
  # neg_conserved_markers_top100 <- conserved_markers %>% filter((.[[2]] < 0) | (.[[7]] < 0)) %>% slice_head(n = 100)
  # write.csv2(pos_conserved_markers, file = paste0("conserved_markers/positive/all_cluster", i, "_cm.csv"))
  # write.csv2(pos_conserved_markers_top10, file = paste0("conserved_markers/positive/top10_cluster", i, "_cm.csv"))
  # write.csv2(pos_conserved_markers_top100, file = paste0("conserved_markers/positive/top100_cluster", i, "_cm.csv"))
  # write.csv2(neg_conserved_markers, file = paste0("conserved_markers/negative/all_cluster", i, "_cm.csv"))
  # write.csv2(neg_conserved_markers_top10, file = paste0("conserved_markers/negative/top10_cluster", i, "_cm.csv"))
  # write.csv2(neg_conserved_markers_top100, file = paste0("conserved_markers/negative/top100_cluster", i, "_cm.csv"))
  #



  # Error in ValidateCellGroups(object = object, cells.1 = cells.1, cells.2 = cells.2,  :
  # Cell group 2 is empty - no cells with identity class
  # TODO check if ident.2 has cells before comparison is done, otherwise get above error message

  # ## create condition markers for integrated data within each cluster between each condition
  # ## DEV NOTE: this is not pairwise if more than 2 conditions are integrated at the same time
  subset <- subset(integrated, seurat_clusters == i)
  # ##  change cluster identity to original identity to find markers between conditions
  Idents(subset) <- subset$orig.ident
  condition_markers <- FindMarkers(subset, ident.1 = "BL_C", verbose = T, only.pos = FALSE)
  # pos_condition_markers <- condition_markers %>% filter(avg_log2FC > 0)
  # pos_condition_markers_top10 <- condition_markers %>% filter(avg_log2FC > 0) %>% slice_head(n = 10)
  # pos_condition_markers_top100 <- condition_markers %>% filter(avg_log2FC > 0) %>% slice_head(n = 100)
  # neg_condition_markers <- condition_markers %>% filter(avg_log2FC < 0)
  # neg_condition_markers_top10 <- condition_markers %>% filter(avg_log2FC < 0) %>% slice_head(n = 10)
  # neg_condition_markers_top100 <- condition_markers %>% filter(avg_log2FC < 0) %>% slice_head(n = 100)
  # write.csv2(pos_condition_markers, file = paste0("condition_markers/positive/all_cluster", i, "_cm.csv"))
  # write.csv2(pos_condition_markers_top10, file = paste0("condition_markers/positive/top10_cluster", i, "_cm.csv"))
  # write.csv2(pos_condition_markers_top100, file = paste0("condition_markers/positive/top100_cluster", i, "_cm.csv"))
  # write.csv2(neg_condition_markers, file = paste0("condition_markers/negative/all_cluster", i, "_cm.csv"))
  # write.csv2(neg_condition_markers_top10, file = paste0("condition_markers/negative/top10_cluster", i, "_cm.csv"))
  # write.csv2(neg_condition_markers_top100, file = paste0("condition_markers/negative/top100_cluster", i, "_cm.csv"))
  print(paste('Cluster ID:', i, ' before condition_markers call'))
  FGSEA_analysis(markers = condition_markers, working_directory = work_dir, marker_type = 'condition_markers', cluster = i)

  # DEVNOTE if want to assign each table to its own variable, use assign() and get()
  # assign(paste0("cluster", i, "_markers"), markers)
  ## get(paste0("cluster", i, "_markers"))
}
## cleanup environment
# rm("markers", "pos_markers", "pos_markers_top10", "pos_markers_top100", "neg_markers", "neg_markers_top10", "neg_markers_top100",
#    "conserved_markers", "pos_conserved_markers", "pos_conserved_markers_top10", "pos_conserved_markers_top100", "neg_conserved_markers", "neg_conserved_markers_top10", "neg_conserved_markers_top100",
#    "condition_markers", "pos_condition_markers", "pos_condition_markers_top10", "pos_condition_markers_top100", "neg_condition_markers", "neg_condition_markers_top10", "neg_condition_markers_top100",
#    "cluster_ids", "subset")




