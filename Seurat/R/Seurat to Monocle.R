library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(monocle3)
library(patchwork)
library(magrittr)
library(ggplot2)



plot_cells.adjusted <- function(cds, x = 1, y = 2,
                                 reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                                 color_cells_by = "cluster", group_cells_by = c("cluster", "partition"),
                                 genes = NULL, show_trajectory_graph = TRUE,
                                 trajectory_graph_color = "grey28", trajectory_graph_segment_size = 0.75,
                                 norm_method = c("log", "size_only"), label_cell_groups = TRUE,
                                 label_groups_by_cluster = TRUE, group_label_size = 2, labels_per_group = 1,
                                 label_branch_points = TRUE, label_roots = TRUE, label_leaves = TRUE,
                                 graph_label_size = 2, cell_size = 0.35, cell_stroke = I(cell_size/2),
                                 alpha = 1, min_expr = 0.1, rasterize = FALSE, scale_to_range = FALSE,
                                 label_principal_points = FALSE)
{
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for", reduction_method,
                                      "calculated.", "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method, "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >= max(x, y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension", "space."))
  if (!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster",
                                                  "partition", "pseudotime") | color_cells_by %in%
                              names(colData(cds)), msg = paste("color_cells_by must one of",
                                                               "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."))
    if (color_cells_by == "pseudotime") {
      tryCatch({
        pseudotime(cds, reduction_method = reduction_method)
      }, error = function(x) {
        stop(paste("No pseudotime for", reduction_method,
                   "calculated. Please run order_cells with",
                   "reduction_method =", reduction_method, "before attempting to color by pseudotime."))
      })
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))
  norm_method = match.arg(norm_method)
  ### MONOCLE3 ORIGINAL CODE
  # group_cells_by = match.arg(group_cells_by)
  ### MY INJECTED CODE
  if (!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% c("cluster", "partition") | group_cells_by %in%
                              names(colData(cds)), msg = paste("group_cells_by must be one of",
                                                               "'cluster', 'partition',", "or a column in the colData table."))}




  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))
  if (show_trajectory_graph && is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  if (label_principal_points && is.null(principal_graph(cds)[[reduction_method]])) {
    message("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
    label_principal_points = FALSE
  }
  if (label_principal_points) {
    label_branch_points <- FALSE
    label_leaves <- FALSE
    label_roots <- FALSE
  }
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize) {
    plotting_func <- ggrastr::geom_point_rast
  }
  else {
    plotting_func <- ggplot2::geom_point
  }
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    data_df$cell_group <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    ### MY INJECTED CODE
    data_df$cell_group <- colData(cds)[data_df$sample_name,
                                       group_cells_by]
    ### MONOCLE3 ORIGINAL CODE
    # stop("Error: unrecognized way of grouping cells.")
  }
  if (color_cells_by == "cluster") {
    data_df$cell_color <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    data_df$cell_color <- colData(cds)[data_df$sample_name,
                                       color_cells_by]
  }
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = x,
                                         prin_graph_dim_2 = y) %>% dplyr::mutate(sample_name = rownames(.),
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[[reduction_method]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from",
                                                                     target = "to") %>% dplyr::left_join(ica_space_df %>%
                                                                                                           dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
                                                                                                                          source_prin_graph_dim_2 = "prin_graph_dim_2"),
                                                                                                         by = "source") %>% dplyr::left_join(ica_space_df %>%
                                                                                                                                               dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1",
                                                                                                                                                              target_prin_graph_dim_2 = "prin_graph_dim_2"),
                                                                                                                                             by = "target")
  }
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2) {
      markers = unlist(genes[, 1], use.names = FALSE)
    }
    else {
      markers = genes
    }
    markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in%
                                       markers) | (row.names(rowData(cds)) %in% markers),
                                    , drop = FALSE]
    markers_rowData <- as.data.frame(markers_rowData)
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData),
                                                     , drop = FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
      if (!is.null(dim(genes)) && dim(genes) >= 2) {
        genes = as.data.frame(genes)
        row.names(genes) = genes[, 1]
        genes = genes[row.names(cds_exprs), ]
        agg_mat = as.matrix(aggregate_gene_expression(cds,
                                                      genes, norm_method = norm_method, scale_agg_values = FALSE))
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id",
                                          "cell_id")
        if (is.factor(genes[, 2]))
          markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                            levels = levels(genes[, 2]))
        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      }
      else {
        cds_exprs@x = round(10000 * cds_exprs@x)/10000
        markers_exprs = matrix(cds_exprs, nrow = nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id",
                                          "cell_id")
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y = "row.names")
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <- as.character(markers_exprs$feature_id)
        }
        else {
          markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        }
        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) |
                                                !as.character(markers_exprs$feature_label) %in%
                                                markers, as.character(markers_exprs$feature_id),
                                              as.character(markers_exprs$feature_label))
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
        if (norm_method == "size_only")
          expression_legend_label = "Expression"
        else expression_legend_label = "log10(Expression)"
      }
      if (scale_to_range) {
        markers_exprs = dplyr::group_by(markers_exprs,
                                        feature_label) %>% dplyr::mutate(max_val_for_feature = max(value),
                                                                         min_val_for_feature = min(value)) %>% dplyr::mutate(value = 100 *
                                                                                                                               (value - min_val_for_feature)/(max_val_for_feature -
                                                                                                                                                                min_val_for_feature))
        expression_legend_label = "% Max"
      }
    }
  }
  if (label_cell_groups && is.null(color_cells_by) == FALSE) {
    if (is.null(data_df$cell_color)) {
      if (is.null(genes)) {
        message(paste(color_cells_by, "not found in colData(cds), cells will",
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }
    else {
      if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        if (label_groups_by_cluster && is.null(data_df$cell_group) ==
            FALSE) {
          text_df = data_df %>% dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster = dplyr::n()) %>%
            dplyr::group_by(cell_color, .add = TRUE) %>%
            dplyr::mutate(per = dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(),
                                                         text_x = stats::median(x = data_dim_1),
                                                         text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        }
        else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>%
            dplyr::mutate(per = 1)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(),
                                                         text_x = stats::median(x = data_dim_1),
                                                         text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }
        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
      }
      else {
        message(paste("Cells aren't colored in a way that allows them to",
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }
  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name",
                     by.y = "cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr,
                                          value, NA))
    ya_sub <- data_df[!is.na(data_df$value), ]
    na_sub <- data_df[is.na(data_df$value), ]
    if (norm_method == "size_only") {
      g <- ggplot(data = data_df, aes(x = data_dim_1,
                                      y = data_dim_2)) + plotting_func(aes(data_dim_1,
                                                                           data_dim_2), size = I(cell_size), stroke = I(cell_stroke),
                                                                       color = "grey80", alpha = alpha, data = na_sub) +
        plotting_func(aes(color = value), size = I(cell_size),
                      stroke = I(cell_stroke), data = ya_sub[order(ya_sub$value),
                      ]) + viridis::scale_color_viridis(option = "viridis",
                                                        name = expression_legend_label, na.value = NA,
                                                        end = 0.8, alpha = alpha) + guides(alpha = FALSE) +
        facet_wrap(~feature_label)
    }
    else {
      g <- ggplot(data = data_df, aes(x = data_dim_1,
                                      y = data_dim_2)) + plotting_func(aes(data_dim_1,
                                                                           data_dim_2), size = I(cell_size), stroke = I(cell_stroke),
                                                                       color = "grey80", data = na_sub, alpha = alpha) +
        plotting_func(aes(color = log10(value + min_expr)),
                      size = I(cell_size), stroke = I(cell_stroke),
                      data = ya_sub[order(ya_sub$value), ], alpha = alpha) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label, na.value = NA,
                                     end = 0.8, alpha = alpha) + guides(alpha = FALSE) +
        facet_wrap(~feature_label)
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    if (color_cells_by %in% c("cluster", "partition")) {
      if (is.null(data_df$cell_color)) {
        g <- g + geom_point(color = I("gray"), size = I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      }
      else {
        g <- g + geom_point(aes(color = cell_color),
                            size = I(cell_size), stroke = I(cell_stroke),
                            na.rm = TRUE, alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
    else if (class(data_df$cell_color) == "numeric") {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by,
                                            option = "C")
    }
    else {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
  }
  if (show_trajectory_graph) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1",
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                                     yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size,
                          color = I(trajectory_graph_color), linetype = "solid",
                          na.rm = TRUE, data = edge_df)
    if (label_principal_points) {
      mst_branch_nodes <- branch_nodes(cds, reduction_method)
      mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
      mst_root_nodes <- root_nodes(cds, reduction_method)
      pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
      princ_point_df <- ica_space_df %>% dplyr::slice(match(names(pps),
                                                            sample_name))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "white", fill = "black", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, princ_point_df) + ggrepel::geom_text_repel(aes_string(x = "prin_graph_dim_1",
                                                                                                                                                y = "prin_graph_dim_2", label = "sample_name"),
                                                                                                                                     size = I(graph_label_size * 1.5), color = "Black",
                                                                                                                                     na.rm = TRUE, princ_point_df)
    }
    if (label_branch_points) {
      mst_branch_nodes <- branch_nodes(cds, reduction_method)
      branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes),
                                                             sample_name)) %>% dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "white", fill = "black", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1",
                                                                                                                                  y = "prin_graph_dim_2", label = "branch_point_idx"),
                                                                                                                       size = I(graph_label_size), color = "white",
                                                                                                                       na.rm = TRUE, branch_point_df)
    }
    if (label_leaves) {
      mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
      leaf_df <- ica_space_df %>% dplyr::slice(match(names(mst_leaf_nodes),
                                                     sample_name)) %>% dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "black", fill = "lightgray", size = I(graph_label_size *
                                                                          1.5), na.rm = TRUE, leaf_df) + geom_text(aes_string(x = "prin_graph_dim_1",
                                                                                                                              y = "prin_graph_dim_2", label = "leaf_idx"),
                                                                                                                   size = I(graph_label_size), color = "black",
                                                                                                                   na.rm = TRUE, leaf_df)
    }
    if (label_roots) {
      mst_root_nodes <- root_nodes(cds, reduction_method)
      root_df <- ica_space_df %>% dplyr::slice(match(names(mst_root_nodes),
                                                     sample_name)) %>% dplyr::mutate(root_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "black", fill = "white", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, root_df) + geom_text(aes_string(x = "prin_graph_dim_1",
                                                                                                                          y = "prin_graph_dim_2", label = "root_idx"),
                                                                                                               size = I(graph_label_size), color = "black",
                                                                                                               na.rm = TRUE, root_df)
    }
  }
  if (label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df, mapping = aes_string(x = "text_x",
                                                                           y = "text_y", label = "label"), size = I(group_label_size))
    if (is.null(markers_exprs))
      g <- g + theme(legend.position = "none")
  }
  g <- g + monocle_theme_opts() + xlab(paste(reduction_method,
                                             x)) + ylab(paste(reduction_method, y)) + theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill = "white"))
  g
}
# DEVNOTE: use if want to overwrite original function namespace
# namespace of customFunction is R_globalenv, where it is defined, bur should be Seurat as that is ns of my targeted function
environment(plot_cells.adjusted) <- asNamespace("monocle3")
# # then with assignInNameSpace I can basically inject my code their copied function and then substitute it back in their environment
# assignInNamespace("plot_cells", plot_cells.adjusted, ns = "monocle3")

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
# TODO after testing practical questions and applications, finish docs @DEA applied to trajectories especially
## TODO create kinetics plots for these
### h: https://www.nature.com/articles/s41586-019-0969-x/figures/3
## TODO select subset of cells based on principal edge/node
### TODO Monocle 3 selecting cells/segments manually
#### choose_cells()
#### choose_graph_segments()




### INITIALIZE ###
# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
start_time <- format(Sys.time(), "%F %H-%M-%S")
dir.create(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/'), recursive = TRUE)




## pre-selection rds files
### TODO testing
sample_names <- c("BL_A + BL_C")
rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06/integrated - old selection/BL_A + BL_C/after_selection/BL_A + BL_C.rds")

# sample_names <- c("BL_C", "BL_A", "BL_N", "BL_A + BL_C", "BL_N + BL_C")
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_C/BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_A/BL_A.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_N/BL_N.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_N + BL_C/BL_N + BL_C.rds"
# )
## after selection rds files
# sample_names <- c("BL_A + BL_C", "BL_N + BL_C")
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/after_selection/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_N + BL_C/after_selection/BL_N + BL_C.rds"
# )
names(rds.files) <- sample_names

# set genes of interest for trajectory starting point selection
astrocytical_interest_markers <- c("VIM", "S100B", "SOX9", "SLC1A3")
neuronal_interest_markers <- c("MAP2", "RBFOX3", "NEUROG2")
### end initialization ###



# initialize plots to wrap
plots <- list()
for (sample_name in sample_names) {
  # create sample specific directory
  dir.create(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/', sample_name, "/"))
  # set working directory to sample specific directory
  setwd(paste0(work_dir, 'results/', start_time, '/monocle-pseudotime/', sample_name, "/"))

  # set genes of interest
  if ("BL_A" %in% sample_name) {
    genes_of_interest <- astrocytical_interest_markers
  } else if ("BL_N" %in% sample_name) {
    genes_of_interest <- neuronal_interest_markers
    # TODO hoogste SOX2 (voorlopen cellen marker) proberen

    # TODO kleur som marker panels in Seurat (visuele validatie)
  } else {
    genes_of_interest <- c(astrocytical_interest_markers, neuronal_interest_markers)
  }

  # read an integrated saved RDS file
  integrated <- readRDS(rds.files[[sample_name]])

  ## convert from Seurat to Monocle3 object
  cds <- SeuratWrappers::as.cell_data_set(integrated, assay = "SCT")

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
  cds <- match_clustering(integrated, cds, match_seurat_clustering = FALSE)
  if (grepl("\\+", sample_name)) {
    p1 <- plot_cells.adjusted(cds,
                              color_cells_by = "orig.ident",
                              group_cells_by = "orig.ident",
                              show_trajectory_graph = FALSE,
                              label_cell_groups = FALSE, # if false, show legend
                              group_label_size = 3,
                              graph_label_size = 3)
    p1 <- p1 + labs(title="Seurat sample identity")
    plots[[length(plots)+1]] <- p1
  }

  p2 <- monocle3::plot_cells(cds,
                             color_cells_by = "partition",
                             show_trajectory_graph = FALSE,
                             label_cell_groups = FALSE, # if false, show legend
                             group_label_size = 3,
                             graph_label_size = 3)
  p2 <- p2 + labs(title="Monocle partition(s)")
  plots[[length(plots)+1]] <- p2

  ## learn principal graph and plot the trajectory
  cds <- monocle3::learn_graph(cds, use_partition = TRUE)

  for (partition in seq_along(table(partitions(cds)))) {
    # create df of cell names and marker of interest values
    df <- SeuratObject::FetchData(integrated, genes_of_interest)
    # get cell name with highest summed expression for marker(s) of interest
    cell_name <- names(which.max(apply(df[partitions(cds) == partition, ], 1, sum)))
    # calculate pseudo-time from cell with highest expression
    cds <- monocle3::order_cells(cds, root_cells = cell_name)
    # get cell name at furthest point to invert pseudo-time calculations
    cell_name <- names(which.max(pseudotime(cds)[partitions(cds) == partition]))
    # recalculate pseudo-time from root cell furthest from cell with highest expression - as trajectory starting point
    cds <- monocle3::order_cells(cds, root_cells = cell_name)

    # plot once
    if (partition == "1") {
      p3 <- monocle3::plot_cells(cds,
                                 label_cell_groups = FALSE, # if false, show legend
                                 label_groups_by_cluster = TRUE,
                                 label_leaves = TRUE,
                                 label_branch_points = FALSE,
                                 group_label_size = 3,
                                 graph_label_size = 3)
      p3 <- p3 + guides(color = guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
      p3 <- p3 + labs(title="Monocle clusters + trajectory")
      plots[[length(plots)+1]] <- p3
      p4 <- plot_cells.adjusted(cds,
                                color_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                group_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE, # if false, show legend
                                label_groups_by_cluster = TRUE,
                                label_roots = FALSE,
                                label_leaves = TRUE,
                                label_branch_points = FALSE,
                                group_label_size = 3,
                                graph_label_size = 3)
      p4 <- p4 + guides(color = guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
      p4 <- p4 + labs(title="Kriegstein clusters UMAP")
      plots[[length(plots)+1]] <- p4
      p5 <- plot_cells.adjusted(cds,
                                color_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                group_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                show_trajectory_graph = TRUE,
                                label_cell_groups = FALSE, # if false, show legend
                                label_groups_by_cluster = TRUE,
                                label_roots = FALSE,
                                label_leaves = TRUE,
                                label_branch_points = FALSE,
                                group_label_size = 3,
                                graph_label_size = 3)
      p5 <- p5 + guides(color = guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
      p5 <- p5 + labs(title="Kriegstein clusters + trajectory")
      plots[[length(plots)+1]] <- p5
    }
    p6 <- monocle3::plot_cells(cds,
                               color_cells_by = "pseudotime",
                               label_cell_groups = TRUE,
                               label_roots = FALSE,
                               label_leaves = TRUE,
                               label_branch_points = FALSE,
                               group_label_size = 3,
                               graph_label_size = 3)
    p6 <- p6 + labs(title=paste0("Pseudotime + trajectory: partition ", partition))
    plots[[length(plots)+1]] <- p6
  }

  # wrap plots
  pw <- patchwork::wrap_plots(plots, ncol = 2)
  ## save plots
  ggplot2::ggsave(file = paste0("Overview_", sample_name, ".png"), width = 30, height = 20, units = "cm")

  # reset plots list
  plots <- list()
}














## add a feature to Monocle3/cell_data_set object manually
# colData(cds)$monoculture_coculture <- colData(cds)$orig.ident

## get data back to Seurat object
# integrated.sub <- Seurat::as.Seurat(cds, counts = "counts", data = "logcounts", assay = NULL, project = "SingleCellExperiment")

### TODO (?) for doing marker gene (differential?) expression analysis in Monocle 3
# marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", cores=1)
#
# top_specific_markers <- marker_test_res %>%
#   dplyr::filter(fraction_expressing >= 0.10) %>%
#   dplyr::group_by(cell_group) %>%
#   dplyr::top_n(1, pseudo_R2)
#
# top_specific_marker_ids <- unique(top_specific_markers %>% dplyr::pull(gene_id))
#
# rowData(cds)$gene_short_name <- row.names(rowData(cds))
# plot_genes_by_group(cds,
#                     top_specific_marker_ids,
#                     group_cells_by="cluster",
#                     ordering_type="none",
#                     max.size=3)

## TODO check aggregating gene expression to create cluster-level pseudobulk samples to create schematic trajectory
# cell_group_df <- data.frame("cell_ids" = names(clusters(cds)),
#                             "seurat_clusters" = clusters(cds))
# monocle3::aggregate_gene_expression(cds, cell_group_df = cell_group_df)

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
