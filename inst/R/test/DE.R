library(dplyr)

integrated <- qs::qread("C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_4/NSM-NS-NC-M/NSM-NS-NC-M.qs")
integrated$kriegstein.seurat.custom.clusters.mean
integrated$mapmycells_supercluster
idents <- paste(integrated$orig.ident, integrated@meta.data[, "kriegstein.seurat.custom.clusters.mean"], sep = "_")
SeuratObject::Idents(integrated) <- idents

# SeuratObject::Idents(integrated) <- integrated$orig.ident

ref_ident <- c("NSM_Microglia.4","NSM_Microglia.5","NSM_Microglia.6")
vs_ident <- c("M_Microglia.4","M_Microglia.5","M_Microglia.6")
DE_res <- Seurat::FindMarkers(
  integrated,
  assay = "SCT",
  ident.1 = ref_ident,
  ident.2 = vs_ident,
  only.pos = FALSE,
  verbose = TRUE,
  logfc.threshold = 0,
  min.pct = 0)


## sort by average log2 fold-change
DE_res <- DE_res %>% dplyr::arrange(dplyr::desc(avg_log2FC))
head(DE_res)

idents <- c(ref_ident, vs_ident)
features <- c("CCL4", "TYROBP", "TJP1")
Seurat::VlnPlot(integrated, features, idents = idents, assay = "SCT", same.y.lims = T)
Seurat::DotPlot(integrated, features, idents = idents, assay = "SCT")



rescale <- function(x, from, to) {(x - min(x))/(max(x)-min(x)) * (to - from) + from}
n_cells_ref <- sum(SeuratObject::Idents(integrated) %in% ref_ident)
n_cells_vs <- sum(SeuratObject::Idents(integrated) %in% vs_ident)
n_cells.1 <- DE_res$pct.1 * n_cells_ref
n_cells.2 <- DE_res$pct.2 * n_cells_vs
absolute_diff_scaled <- rescale(abs(n_cells.1 - n_cells.2), 0, 1)
relative_diff_scaled <- rescale(abs(DE_res$pct.1-DE_res$pct.2), 0, 1)
point_size <- rescale(absolute_diff_scaled + relative_diff_scaled, 2, 4)
EnhancedVolcano::EnhancedVolcano(
  DE_res,
  lab = rownames(DE_res),
  x = "avg_log2FC",
  y = "p_val_adj",
  # selectLab = genes_to_label,
  boxedLabels = TRUE,
  labFace = "bold",
  drawConnectors = TRUE,
  arrowheads = FALSE,
  widthConnectors = 0.25,
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = point_size, # DEVNOTE: colAlpha can also be a vector
  colAlpha = .35,
  labSize = 4,
  borderWidth = 1.5,
  legendPosition = "right",
  legendIconSize = 4.5,
  legendLabels = c("NS", expression(log[2]~FC), bquote(italic(P)), bquote(italic(P)~+~log[2]~FC)),
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  xlab = expression(log[2]~FC),
  ylab = bquote(~-log[10] ~ italic(P)),
  # title = title,
  # subtitle = paste0(n_cells_1, " vs ", n_cells_2),
  caption = paste0("N genes: ", nrow(DE_res))
)
