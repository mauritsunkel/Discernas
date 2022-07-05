library(Seurat)

### USER SETTINGS
file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated - old selection/BL_A + BL_C/BL_A + BL_C.rds"
data <- readRDS(file)

# selection_panel <- c("MAP2", "NEUROG2", "DCX") # RBFOX3 <-> DCX
selection_panel <- c("VIM", "S100B", "FABP7") # SOX9 <-> FABP7
### END USER SETTINGS



# initialize matrixS
rownames <- c(levels(data$seurat_clusters), "all")
colnames <- c("total", "match 3", "match 3%", "match 2 additionally", "match 2 additionally%")
m <- matrix(NA, nrow = length(rownames), ncol = length(colnames), dimnames = list(rownames, colnames))

# initialize update matrix function
update_matrix <- function(matrix, data, selection_panel, cluster = 0) {
  assay_data <- GetAssayData(data, slot = "data", assay = "SCT")

  n_cells_2 <- sum(sapply(as.data.frame(assay_data[selection_panel, ] > 0), sum) == length(selection_panel)-1)
  n_cells_3 <- sum(sapply(as.data.frame(assay_data[selection_panel, ] > 0), sum) == length(selection_panel))
  n_cells_total <- dim(data)[2]

  matrix[cluster, ] <- c(n_cells_total, n_cells_3, n_cells_3/n_cells_total*100, n_cells_2, (n_cells_2)/n_cells_total*100)

  return(matrix)
}

# count selection matching n markers from panel
for (cluster in seq_along(levels(data$seurat_clusters))) {
  message(cluster)
  subset <- subset(data, seurat_clusters == cluster)
  m <- update_matrix(m, subset, selection_panel, cluster)
}
m <- update_matrix(m, data, selection_panel, length(levels(data$seurat_clusters))+1)

m
message(selection_panel)
