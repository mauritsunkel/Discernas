library(SoupX)

sample_name <- "NS"
samples_dir <- "C:/Users/mauri/Desktop/scRNAseqR/data/samples/sakshi NGN2 microglia"

data.data <- Seurat::Read10X(data.dir = file.path(samples_dir, sample_name, "filtered_feature_bc_matrix"))
data <- Seurat::CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)
rm(data.data)



sample_path <- "C:/Users/mauri/Desktop/scRNAseqR/results/test"
doublet_removal_rate <- NULL
run_doublet_removal <- T
run_ambient_RNA_removal <- T
#TODO functionalize, remove manual rm()s
if (run_doublet_removal) {
  temp_QC_data <- data

  set.seed(1)
  sce <- scDblFinder(temp_QC_data@assays$RNA$counts, clusters = TRUE, dbr = doublet_removal_rate)
  temp_QC_data$scDblFinder.score <- sce$scDblFinder.score
  temp_QC_data$scDblFinder.class <- sce$scDblFinder.class
  rm(sce)

  temp_QC_data <- Seurat::NormalizeData(temp_QC_data)
  temp_QC_data <- Seurat::ScaleData(temp_QC_data)
  temp_QC_data <- Seurat::FindVariableFeatures(temp_QC_data)
  temp_QC_data <- Seurat::RunPCA(temp_QC_data, features = SeuratObject::VariableFeatures(object = temp_QC_data), npcs = 50, verbose = FALSE)
  temp_QC_data <- Seurat::RunUMAP(temp_QC_data, reduction = "pca", dims = 1:30)

  png(file.path(sample_path, "quality_control", paste0("scDblFinder_scores_", sample_name, ".png")))
  Seurat::FeaturePlot(temp_QC_data, features = "scDblFinder.score") +
    ggplot2::labs(subtitle = paste0(
      "singlets: ",
      table(temp_QC_data$scDblFinder.class)[1],
      " - doublets: ",
      paste0(table(temp_QC_data$scDblFinder.class)[2],
             " - doublet rate: ",
             round(table(temp_QC_data$scDblFinder.class)[2]/table(temp_QC_data$scDblFinder.class)[1], digits = 3))))
  dev.off()

  # keep only singlets for downstream processing
  data <- data[, temp_QC_data$scDblFinder.class == "singlet"]
  rm(temp_QC_data)
}

if (run_ambient_RNA_removal) {
  table_of_droplets = Seurat::Read10X(data.dir = file.path(samples_dir, sample_name, "raw_feature_bc_matrix"))
  overlapping_genes <- rownames(data@assays$RNA$counts)[which(rownames(data@assays$RNA$counts) %in% rownames(table_of_droplets))]
  table_of_droplets <- table_of_droplets[overlapping_genes,]
  data@assays$RNA$counts <- data@assays$RNA$counts[overlapping_genes,]
  rm(overlapping_genes)

  tenx_graphclust <- read.csv(file.path(samples_dir, sample_name, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"))
  tenx_clusters <- tenx_graphclust$Cluster[tenx_graphclust$Barcode %in% colnames(data@assays$RNA$counts)]
  rm(tenx_graphclust)
  sc = SoupX::SoupChannel(tod = table_of_droplets, toc = data@assays$RNA$counts) # estimateSoup()
  rm(table_of_droplets)
  sc = SoupX::setClusters(sc, tenx_clusters)
  rm(tenx_clusters)
  sc = SoupX::autoEstCont(sc) # if fails, contamination rate default: 10% -> 0.1, or determine gene (sets) to estimate fraction
  out = SoupX::adjustCounts(sc, roundToInt = F) # TODO roundToInt = T if downstream process needs int counts, instead of float
  colnames(out) <- colnames(data@assays$RNA$counts)
  data@assays$RNA$tenx_counts <- data@assays$RNA$counts
  data@assays$RNA$counts <- out
  data@misc[["SoupX_contamination_percentage"]] <- sc$fit$rhoEst * 100
  rm(sc, out)
}


## with SoupX::load10x()
# test <- SoupX::load10X(file.path(samples_dir, sample_name))
# test <- SoupX::autoEstCont(test)
# test_out <- SoupX::adjustCounts(test)

## Manual sanity checks
# SoupX::plotMarkerMap(test, c("VIM", "SOX9"))
# SoupX::plotChangeMap(test, test_out, c("VIM", "SOX9"))
