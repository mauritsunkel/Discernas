#' Run Monocle3 from Seurat object.
#'
#' Run monocle3 from Seurat object, performing its own clustering, partitioning
#' and then creating a pseudotime trajectory and plots.
#'
#' @param input_files Character vector with paths to Seurat .rds files.
#' @param input_names Character vector with sample names of input files.
#' @param output_dir String containing output directory.
#' @param genes_of_interest Character vector with gene name(s) of interest, used to select pseudotime trajectory starting point.
#'
#' @export
#'
#' @note Notes during development
#'
#' Seurat --> Monocle3 for pseudo time analysis - main site: https://cole-trapnell-lab.github.io/monocle3/
#'
#' Monocle theory - from paper: https://www.nature.com/articles/s41586-019-0969-x
#' reduce dimensionality with UMAP (> t-SNE, also preserves global distance AND complexity O(N) versus O(Nlog[N]))
#' organize into Partitioned Approximate Graph Abstraction (PAGA) (construct k-nearest neighbor graph)
#' identify 'communities' (clusters) with Louvain (or Leiden) algorithm
#' PAGA -> graph where clusters = nodes, edged when more neighborly vs expected binomial (coarse-grained trajectory)
#'
#' Monocle3 -> principal graph (SimplePPT-like, +faster, +large datasets, +trajectory loops, +prune branches)
#' SimplePPT paper coined reverse graph embedding (heavy math, I do not understand most of it)
#' uses 'landmark' cells (locally dense k-mean representation), more cells -> higher resolution but also runtime
#' pruning for smoothing and adding loops is performed at the end
#'
#' compute pseudotime: geodesic distance (shortest-path edge distances) from node to root node (principal node)
#' map each cell to closest principal point with Euclidean distance in UMAP space
#' for each principal graph edge: map similarly to it's endpoints orthogonally (ref 22)
#' then order is defined and geodesic distance is computed as pseudotime
#'
#' identify genes that vary in expression over a trajectory (spatial pseudotime)
#' Use Moran's I statistic: multi-directional and multi-dimensional spatial autocorrelation
#' intuitive explanation: for each cell pair and therein each gene
#' sum(multiply differences of gene expression vs gene mean expression for cell pairs & multiply by cell connection weight)
#' multiply sum by total cell pairs, normalize by total weight of cell interactions and static gene expression for each cell
#' 0 = no change over trajectory, higher = relatively more change over trajectory
#'
#' @examplesIf FALSE
#' # EMC.SKlab.scRNAseq::pseudotime(
#' #  input_files = c(file.path("path", "to", "seurat.rds")),
#' #  input_names = c("sample_name"),
#' #  output_dir = file.path("path", "to", results"),
#' #  genes_of_interest = c("GENES", "OF", "INTEREST"))
pseudotime <- function(input_files, input_names, output_dir, genes_of_interest) {
  # set input names on files
  names(input_files) <- input_names

  base_output_dir <- output_dir
  for (input_name in input_names) {
    # initialize and reset plots to wrap per sample
    plots <- list()

    # create sample specific directory
    output_dir <- file.path(base_output_dir, input_name, 'pseudotime')
    dir.create(output_dir, recursive = TRUE)

    # get data
    data <- readRDS(input_files[[input_name]])

    # convert from Seurat to cell data set object
    cds <- SeuratWrappers::as.cell_data_set(data, assay = "SCT")

    # Monocle 3 requires to run its own clustering
    match_clustering <- function(seurat_obj, monocle_obj, match_seurat_clustering) {
      cds <- monocle3::cluster_cells(monocle_obj, cluster_method = "leiden", resolution = NULL, num_iter = 10, verbose = F)
      if (match_seurat_clustering) {
        # in wanting to match clustering, iterate with different resolutions, converging to Seurat_n_clusters
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
          message(paste0("resolution: ", resolution, " - n_seurat_clusters: ", n_clusters_seurat, " - n_monocle_clusters: ", n_clusters_monocle))
        }
        SummarizedExperiment::colData(cds)$monocle3_clustering_resolution <- resolution
        return(cds)
      } else {
        return(cds)
      }
    }
    cds <- match_clustering(data, cds, match_seurat_clustering = FALSE)

    # if integrated sample, plot sample identity
    if (grepl("\\-", input_name)) {
      p1 <- plot_cells.adjusted(cds,
                                color_cells_by = "orig.ident",
                                group_cells_by = "orig.ident",
                                show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE, # if false, show legend
                                group_label_size = 3,
                                graph_label_size = 3,
                                cell_size = 1,
                                trajectory_graph_segment_size = 2)
      p1 <- p1 + ggplot2::labs(title="Seurat sample identity")
      ggplot2::ggsave(file = file.path(output_dir, 'seurat_sample_identity.png'), width = 30, height = 20, units = "cm")
      plots[[length(plots)+1]] <- p1
    }

    # plot Monocle3 partitions
    p2 <- plot_cells.adjusted(cds,
                               color_cells_by = "partition",
                               group_cells_by = "cluster",
                               show_trajectory_graph = FALSE,
                               label_cell_groups = FALSE, # if false, show legend
                               group_label_size = 3,
                               graph_label_size = 3,
                               cell_size = 1,
                               trajectory_graph_segment_size = 2)
    p2 <- p2 + ggplot2::labs(title="Monocle partition(s)")
    ggplot2::ggsave(file = file.path(output_dir, 'monocle_partitions.png'), width = 30, height = 20, units = "cm")
    plots[[length(plots)+1]] <- p2

    # learn principal graph
    cds <- monocle3::learn_graph(cds, use_partition = TRUE)

    for (partition in seq_along(table(monocle3::partitions(cds)))) {
      message("partition: ", partition)
      # create df of cell names and gene(s) of interest values
      df <- SeuratObject::FetchData(data, genes_of_interest)
      # get cell name with highest summed expression for gene(s) of interest
      if (length(genes_of_interest) > 1) {
        cell_name <- names(which.max(apply(df[monocle3::partitions(cds) == partition, ], 1, sum)))
      } else if (length(genes_of_interest) == 1) {
        cell_name <- rownames(df)[which.max(df[monocle3::partitions(cds) == partition, ])]
      }
      # calculate pseudo-time from the cell with highest expression
      cds <- monocle3::order_cells(cds, root_cells = cell_name)
      # get cell name at furthest point to invert pseudo-time calculations
      cell_name <- names(which.max(monocle3::pseudotime(cds)[monocle3::partitions(cds) == partition]))
      # recalculate pseudo-time from root cell furthest from cell with highest expression - as trajectory starting point
      cds <- monocle3::order_cells(cds, root_cells = cell_name)

      # plot once
      if (partition == "1") {
        # plot Monocle3 clusters and trajectory
        p3 <- plot_cells.adjusted(cds,
                                   group_cells_by = "cluster",
                                   label_cell_groups = FALSE, # if false, show legend
                                   label_groups_by_cluster = TRUE,
                                   label_roots = TRUE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE,
                                   group_label_size = 3,
                                   graph_label_size = 5,
                                   cell_size = 1,
                                   trajectory_graph_segment_size = 2)
        p3 <- p3 + ggplot2::guides(color = ggplot2::guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
        p3 <- p3 + ggplot2::labs(title="Monocle clusters + trajectory")
        ggplot2::ggsave(file = file.path(output_dir, 'monocle_clusters_trajectory.png'), width = 30, height = 20, units = "cm")
        plots[[length(plots)+1]] <- p3

        # plot reference clusters UMAP
        p4 <- plot_cells.adjusted(cds,
                                  color_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                  group_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                  show_trajectory_graph = F,
                                  label_cell_groups = T, # if false, show legend
                                  label_groups_by_cluster = TRUE,
                                  label_roots = TRUE,
                                  label_leaves = FALSE,
                                  label_branch_points = FALSE,
                                  group_label_size = 3,
                                  graph_label_size = 5,
                                  cell_size = 1,
                                  trajectory_graph_segment_size = 2)
        p4 <- p4 + ggplot2::guides(color = ggplot2::guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
        p4 <- p4 + ggplot2::labs(title="Kriegstein clusters UMAP")
        ggplot2::ggsave(file = file.path(output_dir, 'kriegstein_clusters_UMAP.png'), width = 30, height = 20, units = "cm")
        plots[[length(plots)+1]] <- p4

        # plot reference data trajectory
        p5 <- plot_cells.adjusted(cds,
                                  color_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                  group_cells_by = "kriegstein.seurat.custom.clusters.mean",
                                  show_trajectory_graph = TRUE,
                                  label_cell_groups = FALSE, # if false, show legend
                                  label_groups_by_cluster = TRUE,
                                  label_roots = TRUE,
                                  label_leaves = FALSE,
                                  label_branch_points = FALSE,
                                  group_label_size = 3,
                                  graph_label_size = 5,
                                  cell_size = 1,
                                  trajectory_graph_segment_size = 2)
        p5 <- p5 + ggplot2::guides(color = ggplot2::guide_legend(title = "", ncol=2, override.aes = list(size = 4)))
        p5 <- p5 + ggplot2::labs(title="Kriegstein clusters + trajectory")
        ggplot2::ggsave(file = file.path(output_dir, 'kriegstein_clusters_trajectory.png'), width = 30, height = 20, units = "cm")
        plots[[length(plots)+1]] <- p5
      }

      # plot pseudotime and trajectory per Monocle3 partition
      p6 <- plot_cells.adjusted(cds,
                                 color_cells_by = "pseudotime",
                                 group_cells_by = "cluster",
                                 label_cell_groups = TRUE,
                                 label_roots = TRUE,
                                 label_leaves = FALSE,
                                 label_branch_points = FALSE,
                                 group_label_size = 3,
                                 graph_label_size = 5,
                                 cell_size = 1,
                                 trajectory_graph_segment_size = 2)
      p6 <- p6 + ggplot2::labs(title=paste0("Pseudotime + trajectory: partition ", partition))
      ggplot2::ggsave(file = file.path(output_dir, paste0('pseudotime_trajectory_partition_', partition, '.png')), width = 30, height = 20, units = "cm")
      plots[[length(plots)+1]] <- p6
    }

    # wrap and save plots
    pw <- patchwork::wrap_plots(plots, ncol = 2)
    ggplot2::ggsave(file = file.path(output_dir, paste0('Overview_', input_name, '.png')), width = 30, height = 20, units = "cm")
  }
}


#' Run adjusted monocle3::plot_cells.
#'
#' Fix monocle3::plot_cells for coloring plot by groups in colData.
#'
#' @description see monocle3::plot_cells for documentation
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
  assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for", reduction_method,
                                      "calculated.", "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method, "before attempting to plot."))
  low_dim_coords <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >= max(x, y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension", "space."))
  if (!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster",
                                                  "partition", "pseudotime") | color_cells_by %in%
                              names(SummarizedExperiment::colData(cds)), msg = paste("color_cells_by must one of",
                                                               "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."))
    if (color_cells_by == "pseudotime") {
      tryCatch({
        monocle3::pseudotime(cds, reduction_method = reduction_method)
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
                              names(SummarizedExperiment::colData(cds)), msg = paste("group_cells_by must be one of",
                                                               "'cluster', 'partition', or a column in the colData table."))}



  ### MONOCLE3 CODE
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))
  if (show_trajectory_graph && is.null(monocle3::principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  if (label_principal_points && is.null(monocle3::principal_graph(cds)[[reduction_method]])) {
    message("Cannot label principal points when no trajectory to plot. Has learn_graph() been called yet?")
    label_principal_points = FALSE
  }
  if (label_principal_points) {
    label_branch_points <- FALSE
    label_leaves <- FALSE
    label_roots <- FALSE
  }
  gene_short_name <- NA
  input_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize) {
    plotting_func <- ggrastr::geom_point_rast
  }
  else {
    plotting_func <- ggplot2::geom_point
  }
  S_matrix <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$input_name <- row.names(data_df)
  data_df <- as.data.frame(cbind(data_df, SummarizedExperiment::colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$input_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    data_df$cell_group <- tryCatch({
      monocle3::partitions(cds, reduction_method = reduction_method)[data_df$input_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    ### MY INJECTED CODE
    data_df$cell_group <- SummarizedExperiment::colData(cds)[data_df$input_name,
                                       group_cells_by]
    ### MONOCLE3 ORIGINAL CODE
    # stop("Error: unrecognized way of grouping cells.")
  }
  if (color_cells_by == "cluster") {
    data_df$cell_color <- tryCatch({
      monocle3::clusters(cds, reduction_method = reduction_method)[data_df$input_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      monocle3::partitions(cds, reduction_method = reduction_method)[data_df$input_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      monocle3::pseudotime(cds, reduction_method = reduction_method)[data_df$input_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    data_df$cell_color <- SummarizedExperiment::colData(cds)[data_df$input_name,
                                       color_cells_by]
  }
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = x,
                                         prin_graph_dim_2 = y) %>% dplyr::mutate(input_name = rownames(.),
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[[reduction_method]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from",
                                                                     target = "to") %>% dplyr::left_join(ica_space_df %>%
                                                                                                           dplyr::select_(source = "input_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
                                                                                                                          source_prin_graph_dim_2 = "prin_graph_dim_2"),
                                                                                                         by = "source") %>% dplyr::left_join(ica_space_df %>%
                                                                                                                                               dplyr::select_(target = "input_name", target_prin_graph_dim_1 = "prin_graph_dim_1",
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
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/monocle3::size_factors(cds))
      if (!is.null(dim(genes)) && dim(genes) >= 2) {
        genes = as.data.frame(genes)
        row.names(genes) = genes[, 1]
        genes = genes[row.names(cds_exprs), ]
        agg_mat = as.matrix(monocle3::aggregate_gene_expression(cds,
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
    data_df <- merge(data_df, markers_exprs, by.x = "input_name",
                     by.y = "cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr,
                                          value, NA))
    ya_sub <- data_df[!is.na(data_df$value), ]
    na_sub <- data_df[is.na(data_df$value), ]
    if (norm_method == "size_only") {
      g <- ggplot2::ggplot(data = data_df, ggplot2::aes(x = data_dim_1,
                                      y = data_dim_2)) + plotting_func(ggplot2::aes(data_dim_1,
                                                                           data_dim_2), size = I(cell_size), stroke = I(cell_stroke),
                                                                       color = "grey80", alpha = alpha, data = na_sub) +
        plotting_func(ggplot2::aes(color = value), size = I(cell_size),
                      stroke = I(cell_stroke), data = ya_sub[order(ya_sub$value),
                      ]) + viridis::scale_color_viridis(option = "viridis",
                                                        name = expression_legend_label, na.value = NA,
                                                        end = 0.8, alpha = alpha) + ggplot2::guides(alpha = FALSE) +
        ggplot2::facet_wrap(~feature_label)
    }
    else {
      g <- ggplot2::ggplot(data = data_df, ggplot2::aes(x = data_dim_1,
                                      y = data_dim_2)) + plotting_func(ggplot2::aes(data_dim_1,
                                                                           data_dim_2), size = I(cell_size), stroke = I(cell_stroke),
                                                                       color = "grey80", data = na_sub, alpha = alpha) +
        plotting_func(ggplot2::aes(color = log10(value + min_expr)),
                      size = I(cell_size), stroke = I(cell_stroke),
                      data = ya_sub[order(ya_sub$value), ], alpha = alpha) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label, na.value = NA,
                                     end = 0.8, alpha = alpha) + guides(alpha = FALSE) +
        ggplot2::facet_wrap(~feature_label)
    }
  }
  else {
    g <- ggplot2::ggplot(data = data_df, ggplot2::aes(x = data_dim_1, y = data_dim_2))
    if (color_cells_by %in% c("cluster", "partition")) {
      if (is.null(data_df$cell_color)) {
        g <- g + ggplot2::geom_point(color = I("gray"), size = I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      }
      else {
        g <- g + ggplot2::geom_point(ggplot2::aes(color = cell_color),
                            size = I(cell_size), stroke = I(cell_stroke),
                            na.rm = TRUE, alpha = alpha)
      }
      g <- g + ggplot2::guides(color = ggplot2::guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
    else if (class(data_df$cell_color) == "numeric") {
      g <- g + ggplot2::geom_point(ggplot2::aes(color = cell_color), size = I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by,
                                            option = "C")
    }
    else {
      g <- g + ggplot2::geom_point(ggplot2::aes(color = cell_color), size = I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + ggplot2::guides(color = ggplot2::guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }
  }
  if (show_trajectory_graph) {
    g <- g + ggplot2::geom_segment(ggplot2::aes_string(x = "source_prin_graph_dim_1",
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
                                     yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size,
                          color = I(trajectory_graph_color), linetype = "solid",
                          na.rm = TRUE, data = edge_df)
    if (label_principal_points) {
      mst_branch_nodes <- monocle3:::branch_nodes(cds, reduction_method)
      mst_leaf_nodes <- monocle3:::leaf_nodes(cds, reduction_method)
      mst_root_nodes <- monocle3:::root_nodes(cds, reduction_method)
      pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
      princ_point_df <- ica_space_df %>% dplyr::slice(match(names(pps),
                                                            input_name))
      g <- g + ggplot2::geom_point(ggplot2::aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "white", fill = "black", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, princ_point_df) + ggrepel::geom_text_repel(ggplot2::aes_string(x = "prin_graph_dim_1",
                                                                                                                                                y = "prin_graph_dim_2", label = "input_name"),
                                                                                                                                     size = I(graph_label_size * 1.5), color = "Black",
                                                                                                                                     na.rm = TRUE, princ_point_df)
    }
    if (label_branch_points) {
      mst_branch_nodes <- monocle3:::branch_nodes(cds, reduction_method)
      branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes),
                                                             input_name)) %>% dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      g <- g + ggplot2::geom_point(ggplot2::aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "white", fill = "black", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, branch_point_df) + ggplot2::geom_text(ggplot2::aes_string(x = "prin_graph_dim_1",
                                                                                                                                  y = "prin_graph_dim_2", label = "branch_point_idx"),
                                                                                                                       size = I(graph_label_size), color = "white",
                                                                                                                       na.rm = TRUE, branch_point_df)
    }
    if (label_leaves) {
      mst_leaf_nodes <- monocle3:::leaf_nodes(cds, reduction_method)
      leaf_df <- ica_space_df %>% dplyr::slice(match(names(mst_leaf_nodes),
                                                     input_name)) %>% dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      g <- g + ggplot2::geom_point(ggplot2::aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "black", fill = "lightgray", size = I(graph_label_size *
                                                                          1.5), na.rm = TRUE, leaf_df) + ggplot2::geom_text(ggplot2::aes_string(x = "prin_graph_dim_1",
                                                                                                                              y = "prin_graph_dim_2", label = "leaf_idx"),
                                                                                                                   size = I(graph_label_size), color = "black",
                                                                                                                   na.rm = TRUE, leaf_df)
    }
    if (label_roots) {
      mst_root_nodes <- monocle3:::root_nodes(cds, reduction_method)
      root_df <- ica_space_df %>% dplyr::slice(match(names(mst_root_nodes),
                                                     input_name)) %>% dplyr::mutate(root_idx = seq_len(dplyr::n()))
      g <- g + ggplot2::geom_point(ggplot2::aes_string(x = "prin_graph_dim_1",
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size),
                          color = "black", fill = "white", size = I(graph_label_size *
                                                                      1.5), na.rm = TRUE, root_df) + ggplot2::geom_text(ggplot2::aes_string(x = "prin_graph_dim_1",
                                                                                                                          y = "prin_graph_dim_2", label = "root_idx"),
                                                                                                               size = I(graph_label_size), color = "black",
                                                                                                               na.rm = TRUE, root_df)
    }
  }
  if (label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df, mapping = ggplot2::aes_string(x = "text_x",
                                                                           y = "text_y", label = "label"), size = I(group_label_size))
    if (is.null(markers_exprs))
      g <- g + ggplot2::theme(legend.position = "none")
  }
  g <- g + monocle3:::monocle_theme_opts() + ggplot2::xlab(paste(reduction_method,
                                             x)) + ggplot2::ylab(paste(reduction_method, y)) + ggplot2::theme(legend.key = ggplot2::element_blank()) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))
  g
}
# environment(plot_cells.adjusted) <- asNamespace("monocle3") # TODO no need for this line if no need Monoclo3 namespace, used package::functions where possible
# assignInNamespace("plot_cells", plot_cells.adjusted, ns = "monocle3") # TODO keep this line?
