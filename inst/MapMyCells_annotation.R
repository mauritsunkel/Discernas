
reticulate::use_python("C:/Users/mauri/AppData/Local/Programs/Python/Python311")

# seurat_file <- "C:/Users/mauri/Desktop/scRNAseqR/results/sakshi_3/integrated/NSM-NS-NC-M/NSM-NS-NC-M.rds"
# anndata_outfile <- "C:/Users/mauri/Desktop/scRNAseqR/results/sakshi_3/integrated/NSM-NS-NC-M/MapMyCells_AIBS.h5ad"

seurat_file <- "K:/Maurits/Drive/Kushnerlab/Projects/scRNAseq/results/Pipe_29-06 cell_level_selection/BL_C.rds"
anndata_outfile <- "K:/Maurits/Drive/Kushnerlab/Projects/scRNAseq/results/Pipe_29-06 cell_level_selection/BL_A+BL_C_MapMyCells_AIBS.h5ad"

data <- readRDS(seurat_file)
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

anndata::write_h5ad(
  anndata = ad,
  filename = anndata_outfile,
  compression = "gzip")

# RUN MapMyCells: https://knowledge.brain-map.org/mapmycells/process/
## Hierarchical Mapping
## Human Whole Brain

# read MapMyCells result and save to Seurat object (.rds)
mmc <- read.csv(
  file = "K:/Maurits/Drive/Kushnerlab/Projects/scRNAseq/results/Pipe_29-06 cell_level_selection/BL_C_MapMyCells_AIBS_10xWholeHumanBrain(CCN202210140)_HierarchicalMapping_UTC_1721819541296.csv",
  skip = 4)
mmc[mmc == "Splatter"] <- "Neuron"
data$mapmycells_supercluster <- mmc$supercluster_name
data$mapmycells_cluster <- mmc$cluster_name
data$mapmycells_subcluster <- mmc$subcluster_name
data$mapmycells_supercluster_bootstrapping_probability <- mmc$supercluster_bootstrapping_probability
data$mapmycells_cluster_bootstrapping_probability <- mmc$cluster_bootstrapping_probability
data$mapmycells_subcluster_bootstrapping_probability <- mmc$subcluster_bootstrapping_probability
saveRDS(data, file = seurat_file)

# manually export figure
Seurat::DefaultAssay(data) <- "integrated"
Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_supercluster", pt.size = .75)
Seurat::DimPlot(data, reduction = "umap", group.by = "mapmycells_supercluster", label = T) + Seurat::NoLegend()

table(data$mapmycells_supercluster)
# table(mmc$supercluster_bootstrapping_probability > 0.6)
