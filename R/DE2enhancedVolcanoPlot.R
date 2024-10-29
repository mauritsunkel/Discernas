#' Create enhanced Volcano plot from Seurat DE
#'
#' @param seurat_object integrated Seurat object. Containing ref_ident and vs_ident in Idents(seurat_object)
#' @param seurat_DE Seurat DE result, from a FindMarkers() function
#' @param ref_ident reference sample name, pct.1 in Seurat DE result
#' @param vs_ident versus sample name, pct.2 in Seurat DE result
#' @param filedir directory to save plot in, filename is based on reference and versus sample
#' @param comp_name used to handle filenaming and EnhancedVolcano plot title
#' @param genes_to_label default: NULL, if gene vector given, only those are labeled
#'
#' @export
plotEnhancedVolcano <- function(
    seurat_object, seurat_DE,
    ref_ident, vs_ident,
    filedir, comp_name,
    genes_to_label = NULL, label_size = 4) {

  if (vs_ident != 'rest') {
    if(!all(c(ref_ident, vs_ident) %in% levels(SeuratObject::Idents(seurat_object)))) {
      stop("Either ", ref_ident, " or ", vs_ident, " not in Idents(seurat_object")
    }
  } else {
    if(!all(c(ref_ident) %in% levels(SeuratObject::Idents(seurat_object)))) {
      stop("Any of ", ref_ident, " not in Idents(seurat_object")
    }
  }




  message("EnhancedVolcano:", comp_name)
  if (comp_name == "name") {
    title <- paste0(vs_ident, " vs ", ref_ident)
  } else {
    title_split <- strsplit(comp_name, "_vs_")
    title <- paste0(title_split[[1]][2], " vs ", title_split[[1]][1])
  }

  rescale <- function(x, from, to) {(x - min(x))/(max(x)-min(x)) * (to - from) + from}
  n_cells_ref <- sum(SeuratObject::Idents(seurat_object) %in% ref_ident)
  if (vs_ident != 'rest') {
    n_cells_vs <- sum(SeuratObject::Idents(seurat_object) %in% vs_ident)
  } else {
    n_cells_vs <- sum(!SeuratObject::Idents(seurat_object) %in% ref_ident)
  }

  n_cells.1 <- seurat_DE$pct.1 * n_cells_ref
  n_cells.2 <- seurat_DE$pct.2 * n_cells_vs
  absolute_diff_scaled <- rescale(abs(n_cells.1 - n_cells.2), 0, 1)
  relative_diff_scaled <- rescale(abs(seurat_DE$pct.1-seurat_DE$pct.2), 0, 1)
  point_size <- rescale(absolute_diff_scaled + relative_diff_scaled, 1, 4)

  x_axis_min <- floor(min(seurat_DE$avg_log2FC))
  x_axis_max <- ceiling(max(seurat_DE$avg_log2FC))

  p <- EnhancedVolcano::EnhancedVolcano(
    seurat_DE,
    lab = rownames(seurat_DE),
    x = "avg_log2FC",
    y = "p_val_adj",
    selectLab = genes_to_label,
    boxedLabels = TRUE,
    labFace = "bold",
    drawConnectors = TRUE,
    arrowheads = FALSE,
    widthConnectors = 0.25,
    pCutoff = 0.005,
    FCcutoff = 1,
    pointSize = point_size,
    colAlpha = .35,
    labSize = label_size,
    borderWidth = 1.5,
    legendPosition = "right",
    legendIconSize = 4.5,
    legendLabels = c("NS", expression(log[2]~FC), bquote(italic(P)), bquote(italic(P)~+~log[2]~FC)),
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    xlab = expression(log[2]~FC),
    ylab = bquote(~-log[10] ~ italic(P)),
    title = title,
    subtitle = paste0(n_cells_vs, " vs ", n_cells_ref),
    caption = paste0("N genes: ", nrow(seurat_DE))
  ) + ggplot2::coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
    ggplot2::scale_x_continuous(
      breaks=seq(x_axis_min, x_axis_max, 1))

  filename <- file.path(filedir, paste0(sub(" vs ", "_vs_", title), '_EVP.png'))
  ggplot2::ggsave(plot = p, file = filename, width = 30, height = 20, units = "cm")
}
