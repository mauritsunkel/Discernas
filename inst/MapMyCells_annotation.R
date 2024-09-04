
# reticulate::use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python311")

### USER INPUT ###
seurat_file <- "C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_4/NSM-NS-NC-M/NSM-NS-NC-M.rds"
##################



data <- readRDS(seurat_file)
if ("integrated" %in% names(data@assays)) Seurat::DefaultAssay(data) <- "integrated"

ensemblIDs_geneSymbols <- read.csv(
  file = "C:/Users/mauri/Desktop/scRNAseqR/inst/extdata/ensembl_genes.tsv",
  header = F,
  sep = "\t")
colnames(ensemblIDs_geneSymbols) <- c("ENSEMBL_ID", "GENE_SYMBOL", "GENE_TYPE")

# anndata: rows are cells, columns are genes (Ensembl IDs), data expected to be raw counts
datat <- Matrix::t(data@assays$RNA$counts)
colnames(datat) <- plyr::mapvalues(
  x = colnames(datat),
  from = ensemblIDs_geneSymbols$GENE_SYMBOL,
  to = ensemblIDs_geneSymbols$ENSEMBL_ID,
  warn_missing = FALSE
)

ad <- anndata::AnnData(
  X = datat,
  obs = data.frame(group = rownames(datat), row.names = rownames(datat)),
  var = data.frame(type = colnames(datat), row.names = colnames(datat))
)

dir.create(file.path(dirname(seurat_file), "annotation_MapMyCells"), recursive = T)
anndata_outfile <- file.path(dirname(seurat_file), "annotation_MapMyCells", sub(".rds", ".h5ad", basename(seurat_file)))
anndata::write_h5ad(
  anndata = ad,
  filename = anndata_outfile,
  compression = "gzip")

# RUN MapMyCells: https://knowledge.brain-map.org/mapmycells/process/
## Human Whole Brain
## Hierarchical Mapping
### unpack downloaded results.zip to: dirname(anndata_outfile)

# read MapMyCells result and save to Seurat object (.rds)
mmc_files <- list.files(dirname(anndata_outfile), full.names = T)
mmc_csv <- mmc_files[grepl(".csv", mmc_files)]
mmc <- read.csv(
  file = mmc_csv,
  skip = 4)
mmc[mmc == "Splatter"] <- "Neuron"
data$mapmycells_supercluster <- mmc$supercluster_name
data$mapmycells_cluster <- mmc$cluster_name
data$mapmycells_subcluster <- mmc$subcluster_name
data$mapmycells_supercluster_bootstrapping_probability <- mmc$supercluster_bootstrapping_probability
data$mapmycells_cluster_bootstrapping_probability <- mmc$cluster_bootstrapping_probability
data$mapmycells_subcluster_bootstrapping_probability <- mmc$subcluster_bootstrapping_probability

write.csv2(table(data$mapmycells_supercluster), file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_supercluster_table.csv")))
write.csv2(table(data$mapmycells_cluster), file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_cluster_table.csv")))
write.csv2(table(data$mapmycells_subcluster), file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_subcluster_table.csv")))

p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_supercluster", pt.size = .75)
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_supercluster.png")), width = 30, height = 20, units = "cm")
p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_supercluster", label = T) + Seurat::NoLegend()
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_supercluster_noLegend.png")), width = 30, height = 20, units = "cm")

p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_cluster", pt.size = .75)
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_cluster.png")), width = 30, height = 20, units = "cm")
p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_cluster", label = T) + Seurat::NoLegend()
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_cluster_noLegend.png")), width = 30, height = 20, units = "cm")

p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_subcluster", pt.size = .75)
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_subcluster.png")), width = 30, height = 20, units = "cm")
p <- Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_subcluster", label = T) + Seurat::NoLegend()
ggplot2::ggsave(file = file.path(dirname(seurat_file), "annotation_MapMyCells", paste0("UMAP_subcluster_noLegend.png")), width = 30, height = 20, units = "cm")

# plot composition
composition_df <- as.data.frame(table(data@meta.data[,c("orig.ident", "mapmycells_supercluster")]))
composition_df <- composition_df %>%
  dplyr::arrange(orig.ident, mapmycells_supercluster) %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::mutate(csum = cumsum(Freq))
composition_df$mapmycells_supercluster <- factor(composition_df$mapmycells_supercluster, levels = rev(levels(composition_df$mapmycells_supercluster)))
# remove labels with less than 10 cells
composition_df$Freq[composition_df$Freq < 10] <- NA
p <- ggplot2::ggplot(composition_df, ggplot2::aes(x = orig.ident, y = Freq, fill = mapmycells_supercluster)) +
  ggplot2::geom_bar(stat = "identity", color = "black") +
  # ggrepel::geom_text_repel(ggplot2::aes(x = orig.ident, y = csum, label=csum), color="white", size=3, max.overlaps = Inf) +
  ggplot2::geom_text(ggplot2::aes(x = orig.ident, y = csum, label=Freq), nudge_y = -75, color="white", size=3) + #
  # ggplot2::geom_text() +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "MapMyCells supercluster")) +
  ggplot2::labs(x = "", y = "# cells") +
  ggplot2::theme_minimal()
ggplot2::ggsave(plot = p, file = file.path(dirname(seurat_file), "annotation_MapMyCells", 'samples_composition.png'), width = 30, height = 20, units = "cm")

saveRDS(data, file = seurat_file)
