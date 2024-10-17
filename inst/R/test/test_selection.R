## NOTE: integrated becomes integrated.test during selection testing
integrated <- qs::qread("C:/Users/Maurits/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV8/NSM-NS-NC-M/NSM-NS-NC-M.qs")
# perform_cluster_level_selection <- T
# "CX3CR1" %in% rownames(integrated)
selection_panel <- c("VIM", "S100B", "SOX9") # 19: 14.7, 24: 9.6, 26: 40
# selection_panel <- c("MAP2", "DCX", "NEUROG2") # 8: 16.5, 12: 16.4
# selection_panel <- c("AIF1", "CSF1R", "CX3CR1") # CX3CR1/GPR34
selection_percent_expressed <- 15



if (perform_cluster_level_selection) {
  # create dotplot to extract percent expressed information
  p <- Seurat::DotPlot(integrated, features = selection_panel)
  p
  # get cluster names where percent expressed is above %threshold for each gene of selection_panel
  cluster_selection <- names(which(table(p$data[p$data$pct.exp > selection_percent_expressed,]$id) == length(unique(p$data$features.plot))))
  # if no clusters selected, set selection to NULL
  if (length(cluster_selection) == 0) {
    cluster_selection <- NULL
  }
  cluster_selection
  # perform cluster selection
  Seurat::DefaultAssay(integrated) <- "RNA"
  integrated_test <- subset(integrated, idents = cluster_selection)
}
Seurat::FeaturePlot(integrated, features = selection_panel)
if (perform_cell_level_selection) {
  Seurat::DefaultAssay(integrated) <- "SCT"
  layer_data <- SeuratObject::LayerData(integrated)
  # select each cell that has expresses each gene from selection_panel
  cellsToSelect <- sapply(as.data.frame(layer_data[selection_panel, ] > 0), sum) == length(rownames(layer_data[selection_panel, ]))
  # perform cell selection
  Seurat::DefaultAssay(integrated) <- "RNA"
  integrated_test <- integrated[, cellsToSelect]
}
reference_annotations <- list(mapmycells_supercluster = c("Microglia"))
if (!is.null(reference_annotations)) {
  message("Selecting specified annotation from reference as idents")
  # set idents to reference name
  Seurat::Idents(integrated) <- integrated@meta.data[, names(reference_annotations)]
  # select idents by annotation (must be in reference)
  Seurat::DefaultAssay(integrated) <- "RNA"
  integrated_test <- subset(integrated, idents = reference_annotations[[names(reference_annotations)]])
}

# cleanup filtered Seurat object
reset <- Seurat::DietSeurat(integrated_test, assays = c("RNA"))
# split (recommended by Seurat V5: https://github.com/satijalab/seurat/issues/8406) -> layers
reset[["RNA"]] <- split(reset[["RNA"]], f = reset$orig.ident)

# default assay set to "SCT"
options(future.globals.maxSize = 8000 * 1024^2)
reset <- Seurat::SCTransform(reset, vst.flavor = "v2", method = "glmGamPoi", return.only.var.genes = FALSE)
reset <- Seurat::RunPCA(reset, features = SeuratObject::VariableFeatures(object = reset), npcs = 50, verbose = TRUE)
reset.integrated <- Seurat::IntegrateLayers(
  object = reset,
  method = Seurat::RPCAIntegration,
  orig.reduction = "pca",
  normalization.method = "SCT",
  new.reduction = "integrated.dr",
  verbose = TRUE,
  dims = 1:50)
# post integration processing
reset.integrated <- Seurat::PrepSCTFindMarkers(reset.integrated, assay = "SCT")
reset.integrated <- Seurat::RunUMAP(reset.integrated, dims = 1:50, reduction = "integrated.dr")
reset.integrated <- Seurat::FindNeighbors(reset.integrated, reduction = "integrated.dr", dims = 1:50)
reset.integrated <- Seurat::FindClusters(reset.integrated, resolution = 0.8)
Seurat::DimPlot(reset.integrated, reduction = "umap", group.by = c("orig.ident"))
