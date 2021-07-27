library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(magrittr)
library(ggplot2)

# This strategy works by constructing a graph W on all cells using the principal graph as a guide,
# and then computing the pseudotime of each cell as its geodesic distance back to one or more user-selected root nodes
# in the trajectory.

# A geodesic curve commonly represents the shortest distance (in some sense) on a surface
## a Riemannian manifold (like a UMAP) in our case.




### TODO try CCA instead of RPCA?

### TODO try monoculture samples for trajectory analysis as well


# read an integrated saved RDS file
## BL_A + BLC_C
integrated <- readRDS("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/integrated/BL_A + BL_C/integrated.rds")
## BL_N + BLC_C
integrated <- readRDS("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-07-20 13-23-16/integrated/BL_N + BL_C/integrated.rds")



cds <- SeuratWrappers::as.cell_data_set(integrated)
cds <- monocle3::cluster_cells(cds, cluster_method = "leiden") # it's a requirement that Monocle 3 runs its own clustering!
p1 <- monocle3::plot_cells(cds, show_trajectory_graph = FALSE, group_label_size = 4)
p2 <- monocle3::plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, group_label_size = 4)
patchwork::wrap_plots(p1, p2)

integrated.sub <- subset(Seurat::as.Seurat(cds, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment"), monocle3_partitions == 1)
cds <- SeuratWrappers::as.cell_data_set(integrated.sub)
cds <- monocle3::learn_graph(cds)
monocle3::plot_cells(cds, label_groups_by_cluster = T, label_leaves = T, label_branch_points = F, group_label_size = 6, graph_label_size = 4)

max.avp <- which.max(unlist(Seurat::FetchData(integrated.sub, "AVP")))
max.avp <- colnames(integrated.sub)[max.avp]
cds <- monocle3::order_cells(cds, root_cells = max.avp)
monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T,
                     label_branch_points = F, group_label_size = 6, graph_label_size = 4)

# Set the assay back as 'integrated'
integrated.sub <- Seurat::as.Seurat(cds, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment")
FeaturePlot(integrated.sub, "monocle3_pseudotime")


