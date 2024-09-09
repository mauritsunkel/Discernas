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
    genes_to_label = NULL) {

  if(!all(c(ref_ident, vs_ident) %in% levels(SeuratObject::Idents(seurat_object)))) {
    stop("Either ", ref_ident, " or ", vs_ident, " not in Idents(seurat_object")
  }

  if (comp_name == "name") { # is.null(comp_name || comp_name == ""
    title <- paste0(ref_ident, " vs ", vs_ident)
  } else {
    title <- sub("_vs_", " vs ", comp_name)
  }
  ## to plot avg_log2FC in more intuitive sense of left for ident.1 and right for ident.2
  seurat_DE$avg_log2FC <- seurat_DE$avg_log2FC * -1

  ## set plot point size to represent relative difference in amount of cells
  rescale <- function(val, from, to) {
    from + to * ((val - min(val)) / (max(val) - min(val)))
  }
  pseudocount <- 0.00001
  pct_ratio <- abs(log((seurat_DE$pct.1 + pseudocount)/(seurat_DE$pct.2 + pseudocount)))
  point_size <- rescale(pct_ratio, from = 1, to = 5)

  n_cells_1 <- sum(SeuratObject::Idents(seurat_object) %in% ref_ident)
  n_cells_2 <- sum(SeuratObject::Idents(seurat_object) %in% vs_ident)
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
    pCutoff = 1e-05,
    FCcutoff = 1,
    pointSize = point_size, # DEVNOTE: colAlpha can also be a vector
    labSize = 4,
    borderWidth = 1.5,
    legendPosition = "right",
    legendIconSize = 4.5,
    legendLabels = c("NS", expression(log[2]~FC), bquote(italic(P)), bquote(italic(P)~+~log[2]~FC)),
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    xlab = expression(log[2]~FC),
    ylab = bquote(~-log[10] ~ italic(P)),
    title = title,
    subtitle = paste0(n_cells_1, " vs ", n_cells_2),
    caption = paste0("N genes: ", nrow(seurat_DE))
  ) + ggplot2::coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
    ggplot2::scale_x_continuous(
      breaks=seq(x_axis_min, x_axis_max, 1))

  filename <- file.path(filedir, paste0(sub(" vs ", "_vs_", title), '_EVP.png'))
  ggplot2::ggsave(plot = p, file = filename, width = 30, height = 20, units = "cm")
}
