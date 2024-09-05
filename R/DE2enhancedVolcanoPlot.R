#' Create enhanced Volcano plot from Seurat DE
#'
#' @param seurat_object integrated Seurat object. Containing ref_sample and vs_sample in Idents(seurat_object)
#' @param seurat_DE Seurat DE result, from a FindMarkers() function
#' @param ref_sample reference sample name, pct.1 in Seurat DE result
#' @param vs_sample versus sample name, pct.2 in Seurat DE result
#' @param filedir directory to save plot in, filename is based on reference and versus sample
#' @param genes_to_label default: NULL, if gene vector given, only those are labeled
#'
#' @export
plotEnhancedVolcano <- function(
    seurat_object, seurat_DE,
    ref_sample, vs_sample,
    filedir, genes_to_label = NULL) {

  if(!all(c(ref_sample, vs_sample) %in% levels(SeuratObject::Idents(integrated)))) {
    stop("Either ", ref_sample, " or ", vs_sample, " not in Idents(seurat_object")
  }

  n_cells_1 <- sum(SeuratObject::Idents(seurat_object) == ref_sample)
  n_cells_2 <- sum(SeuratObject::Idents(seurat_object) == vs_sample)
  x_axis_min <- floor(min(seurat_DE$avg_log2FC))
  x_axis_max <- ceiling(max(seurat_DE$avg_log2FC))
  p <- EnhancedVolcano::EnhancedVolcano(
    temp_DE,
    lab = rownames(temp_DE),
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
    pointSize = 1.5,
    labSize = 4,
    borderWidth = 1.5,
    legendPosition = "right",
    legendIconSize = 4.5,
    legendLabels = c("NS", expression(log[2]~FC), bquote(italic(P)), bquote(italic(P)~+~log[2]~FC)),
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    xlab = expression(log[2]~FC),
    ylab = bquote(~-log[10] ~ italic(P)),
    title = paste0(ref_sample, " vs ", vs_sample),
    subtitle = paste0(n_cells_1, " vs ", n_cells_2)
  ) + ggplot2::coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
    ggplot2::scale_x_continuous(
      breaks=seq(x_axis_min, x_axis_max, 1))
  ggplot2::ggsave(plot = p, file = file.path(filedir, paste0(ref_sample, "_vs_", vs_sample, '_EVP.png')), width = 30, height = 20, units = "cm")
}
