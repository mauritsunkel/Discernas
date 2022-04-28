library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(magrittr)
library(ggplot2)

# Monocle 3 docs: https://cole-trapnell-lab.github.io/monocle3/docs/differential/
## Seurat to Monocle 3 vignette (OLD): https://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
# Monocle 3 Google group issues: https://groups.google.com/g/monocle-3-users

# From Seurat: when to use which assay (intuitive):
## Use integrated assay when 'aligning' cell states shared across datasets (i.e. clustering, plotting, pseudotime).
### Use RNA assay when exploring genes that change either across clusters, trajectories, or conditions.
## In this case: assay and NORMALIZED data used for pseudo-time analysis (in this case these 2 options are analogous for pseudo-time trajectory)
### RNA assay --> data slot (NormalizeData(object, normalization.method = "LogNormalize") = treated as log-normalized, corrected data
### integrated assay --> scale.data slot (IntegrateData(object, normalization.method = "SCT") = treated as centered, corrected Pearson residuals

# Seurat --> Monocle 3 for pseudo time analysis - main site: https://cole-trapnell-lab.github.io/monocle3/
# Monocle theory - from paper: https://www.nature.com/articles/s41586-019-0969-x
## reduce dimensionality with UMAP (> t-SNE, also preserves global distance AND complexity O(N) versus O(Nlog[N]))
## organize into partitioned approximate graph abstraction (PAGA) (construct k-nearest neighbor graph)
## identify 'communities' (clusters) with Louvain (or Leiden) algorithm
## PAGA -> graph where clusters = nodes, edged when more neighborly vs expected binomial (coarse-grained trajectory)
## Monocle3 -> principal graph (SimplePPT-like, +faster, +large datasets, +trajectory loops, +prune branches)
### SimplePPT paper coined reverse graph embedding (heavy math, I do not understand most of it)
### uses 'landmark' cells (locally dense k-mean representation), more cells -> higher resolution but also runtime
### pruning for smoothing and adding loops is performed at the end
## compute pseudotime: geodesic distance (shortest-path edge distances) from node to root node (principal node)
### map each cell to closest principal point with Euclidean distance in UMAP space
### for each principal graph edge: map similarly to it's endpoints orthogonally (ref 22)
### then order is defined and geodesic distance is computed as pseudotime
## identify genes that vary in expression over a trajectory (spatial pseudotime)
### Use Moran's I statistic: multi-directional and multi-dimensional spatial autocorrelation
#### intuitive explanation: for each cell pair and therein each gene
##### sum(multiply differences of gene expression vs gene mean expression for cell pairs & multiply by cell connection weight)
###### multiply sum by total cell pairs, normalize by total weight of cell interactions and static gene expression for each cell
####### 0 = no change over trajectory, higher = relatively more change over trajectory









### practical questions to test in code
# TODO how to get Seurat object (meta)features --> Monocle?
## and the inverse
## TODO test Kriegstein labels --> integrated RDS
## TODO use plot_cells then with argument color_cells_by="meta_feature_name" (name of any col in colData(cds))
### TODO add feature manually by colData(cds)$feature_name <- data?

# TODO plot monoculture vs coculture coloring (also with both astrocytes and neurons simultaneously)
## TODO create 'pseudobulk' samples from these 3
### c: https://www.nature.com/articles/s41586-019-0969-x/figures/1
### r: https://www.nature.com/articles/s41586-019-0969-x/figures/7
### TODO create 'pseudobulk' samples from their clusters?
### TODO create kinetics plots for these
#### h: https://www.nature.com/articles/s41586-019-0969-x/figures/3


# TODO after testing practical questions and applications, finish docs @DEA applied to trajectories especially
## TODO select subset of cells based on principal edge/node
### TODO Monocle 3 selecting cells/segments manually
#### choose_cells()
#### choose_graph_segments()




### INITIALIZE ###
# read an integrated saved RDS file
sample_name <- "BL_A + BL_C"
file <- paste0("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection/integrated/", sample_name, "/", sample_name, ".rds")
integrated <- readRDS(file)
SeuratObject::DefaultAssay(integrated) <- "integrated"

# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
start_time <- format(Sys.time(), "%F %H-%M-%S")
dir.create(paste0(work_dir, 'results/'))
dir.create(paste0(work_dir, 'results/', start_time, '/'))
dir.create(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/'))
dir.create(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/', sample_name, "/"))
setwd(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/', sample_name, "/"))













### end initialization ###



## convert from Seurat to Monocle3 object
cds <- SeuratWrappers::as.cell_data_set(integrated)







## Monocle 3 requires to run it's own clustering (flag in custom function allows to match Seurat_n_clusters)
match_clustering <- function(seurat_obj, monocle_obj, match_seurat_clustering) {
  cds <- monocle3::cluster_cells(monocle_obj, cluster_method = "leiden", resolution = NULL, num_iter = 10, verbose = F)
  if (match_seurat_clustering) {
    ## in wanting to match clustering, iterate with different resolutions, converging to Seurat_n_clusters
    resolution = 1e-03
    n_clusters_seurat <- length(levels(Seurat::Idents(seurat_obj)))
    n_clusters_monocle <- length(levels(clusters(cds)))
    while (n_clusters_seurat != n_clusters_monocle) {
      if (n_clusters_seurat > n_clusters_monocle) {
        resolution <- resolution * 2
      } else {
        resolution <- resolution / 5
      }
      cds <- monocle3::cluster_cells(monocle_obj, cluster_method = "leiden", resolution = resolution, num_iter = 10, verbose = F)
      n_clusters_seurat <- length(levels(Seurat::Idents(seurat_obj)))
      n_clusters_monocle <- length(levels(clusters(cds)))
      print(paste0("resolution: ", resolution, " - n_seurat_clusters: ", n_clusters_seurat, " - n_monocle_clusters: ", n_clusters_monocle))
    }
    colData(cds)$monocle3_clustering_resolution <- resolution
    return(cds)
  } else {
    return(cds)
  }
}
cds <- match_clustering(integrated, cds, match_seurat_clustering = TRUE)
p1 <- monocle3::plot_cells(cds, show_trajectory_graph = FALSE, group_label_size = 4)
p2 <- monocle3::plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE, group_label_size = 4)









# TODO is this step even necessarily needed? it just selects a monocle partition, can also plot without selection
## cds <- cds[, partitions(cds) == 1] # for subsetting on monocle object itself
integrated.sub <- subset(Seurat::as.Seurat(cds, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment"), monocle3_partitions == 1)
## select single cell with highest expression for certain (marker) gene to use as trajectory starting point
## TODO currently use PAX6 as stem cell marker, is that okay?
### TODO or select whole cluster of cells to use as starting point (how to select with marker genes?)
#### # cluster2.cells <- WhichCells(object = integrated.sub, idents = 2)
cells.max.pax6 <- which.max(unlist(SeuratObject::FetchData(integrated.sub, "PAX6")))
cells.max.pax6 <- colnames(integrated.sub)[cells.max.pax6]



## learn principal graph and plot the trajectory
cds <- monocle3::learn_graph(cds)
p3 <- monocle3::plot_cells(cds, label_groups_by_cluster = T, label_leaves = T, label_branch_points = F, group_label_size = 6, graph_label_size = 4)

## set root_cells to select trajectory starting point
cds <- monocle3::order_cells(cds, root_cells = cells.max.pax6)
p4 <- monocle3::plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T,
                     label_branch_points = F, group_label_size = 6, graph_label_size = 4)

## add titles to figures and plot together
p1 <- p1 + labs(title="Clustering")
p2 <- p2 + labs(title="Partition(s)")
p3 <- p3 + labs(title="Clustering + trajectory")
p4 <- p4 + labs(title="Pseudotime + trajectory")
pw <- patchwork::wrap_plots(p1, p2, p3, p4)
ggplot2::ggsave(file = paste0("Overview_", sample_name, ".png"), width = 30, height = 20, units = "cm")









### for doing marker gene (differential?) expression analysis in Monocle 3
marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", cores=1)

top_specific_markers <- marker_test_res %>%
  dplyr::filter(fraction_expressing >= 0.10) %>%
  dplyr::group_by(cell_group) %>%
  dplyr::top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% dplyr::pull(gene_id))

rowData(cds)$gene_short_name <- row.names(rowData(cds))
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="none",
                    max.size=3)








## get data back to Seurat object
# integrated.sub <- Seurat::as.Seurat(cds, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment")


## TODO check Monocle3 3d plotting - Github issued: https://github.com/cole-trapnell-lab/monocle3/issues/590
## inspiration: https://www.nature.com/articles/s41586-019-0969-x/figures/4
# cds_3d <- reduce_dimension(cds, max_components = 3, preprocess_method = "PCA")
# match_clustering <- function(seurat_obj, monocle_obj, match_seurat_clustering) {
#   cds <- monocle3::cluster_cells(monocle_obj, cluster_method = "leiden", resolution = NULL, num_iter = 10, verbose = F)
#   if (match_seurat_clustering) {
#     ## in wanting to match clustering, iterate with different resolutions, converging to Seurat_n_clusters
#     resolution = 1e-03
#     n_clusters_seurat <- length(levels(Seurat::Idents(seurat_obj)))
#     n_clusters_monocle <- length(levels(clusters(cds)))
#     while (n_clusters_seurat != n_clusters_monocle) {
#       if (n_clusters_seurat > n_clusters_monocle) {
#         resolution <- resolution * 2
#       } else {
#         resolution <- resolution / 5
#       }
#       cds <- monocle3::cluster_cells(monocle_obj, cluster_method = "leiden", resolution = resolution, num_iter = 10, verbose = F)
#       n_clusters_seurat <- length(levels(Seurat::Idents(seurat_obj)))
#       n_clusters_monocle <- length(levels(clusters(cds)))
#       print(paste0("resolution: ", resolution, " - n_seurat_clusters: ", n_clusters_seurat, " - n_monocle_clusters: ", n_clusters_monocle))
#     }
#     colData(cds)$monocle3_clustering_resolution <- resolution
#     return(cds)
#   } else {
#     return(cds)
#   }
# }
# cds_3d <- match_clustering(integrated, cds_3d, match_seurat_clustering = TRUE)
# cds_3d <- learn_graph(cds_3d)
# cds_3d <- order_cells(cds_3d)
# source(file="C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/R/my_utils/color_palettes.R")
# my.color.palettes <- my.color.palettes # TODO check if works: explicit mention because of scoping error
# # get custom colors
# custom_colors <- my.color.palettes(type = 'mixed')
# plot_cells_3d(cds_3d,
#               color_cells_by="cluster",
#               color_palette = custom_colors[2:2+length(levels(clusters(cds_3d)))],
#               show_trajectory_graph = T,
#               trajectory_graph_segment_size = 5,
#               norm_method = "log",
#               alpha=0.5,
#               min_expr=0)
