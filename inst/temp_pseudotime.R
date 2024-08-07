data <- readRDS("C:/Users/mauri/Desktop/scRNAseqR/results/sakshi_4/NSM-NS-NC-M/NSM-NS-NC-M.rds")
cds <- SeuratWrappers::as.cell_data_set(data, assay = "SCT")

# load match_clustering function
cds <- match_clustering(data, cds, match_seurat_clustering = FALSE)

cds <- monocle3::learn_graph(cds, use_partition = TRUE)



for (partition in seq_along(levels(monocle3::partitions(cds)))) {
  message("partition: ", partition)





  # # create df of cell names and gene(s) of interest values
  # df <- SeuratObject::FetchData(data, genes_of_interest)
  # # get cell name with highest summed expression for gene(s) of interest
  # if (length(genes_of_interest) > 1) {
  #   cell_name <- names(which.max(apply(df[monocle3::partitions(cds) == partition, ], 1, sum)))
  # } else if (length(genes_of_interest) == 1) {
  #   cell_name <- rownames(df)[which.max(df[monocle3::partitions(cds) == partition, ])]
  # }
  # # calculate pseudo-time from the cell with highest expression
  # cds <- monocle3::order_cells(cds, root_cells = cell_name)
  # # get cell name at furthest point to invert pseudo-time calculations
  # cell_name <- names(which.max(monocle3::pseudotime(cds)[monocle3::partitions(cds) == partition]))
  # # recalculate pseudo-time from root cell furthest from cell with highest expression - as trajectory starting point
  # cds <- monocle3::order_cells(cds, root_cells = cell_name)



  ### PLOTTING
}


pseudotime_root_markers <- list(
  "Microglia" = c("AIF1"),
  "Astrocyte" = c("VIM", "S100B", "SOX9"),
  "Neuron"    = c("MAP2", "DCX", "NEUROG2"),
  "Dividing"  = c("MKI67"),
  "other"     = c("FOXJ1")
)
