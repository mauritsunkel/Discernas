library(scDblFinder)

samples_dir <- "C:/Users/mauri/Desktop/scRNAseqR/data/samples/sakshi NGN2 microglia"
sample_name <- "NSM"

data.data <- Seurat::Read10X(data.dir = file.path(samples_dir, sample_name, "filtered_feature_bc_matrix"), strip.suffix = TRUE)
data <- Seurat::CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)

plot_and_remove_doublets <- function(data, sample_path, sample_name) {
  orig_data <- data

  set.seed(1)
  sce <- scDblFinder(data@assays$RNA$counts, clusters = TRUE, dbr = NULL)
  data$scDblFinder.score <- sce$scDblFinder.score
  data$scDblFinder.class <- sce$scDblFinder.class

  data <- Seurat::NormalizeData(data)
  data <- Seurat::ScaleData(data)
  data <- Seurat::FindVariableFeatures(data)
  data <- Seurat::RunPCA(data, features = SeuratObject::VariableFeatures(object = data), npcs = 50, verbose = FALSE)
  data <- Seurat::RunUMAP(data, reduction = "pca", dims = 1:30)

  png(file.path(sample_path, "quality_control", paste0("scDblFinder_scores_", sample_name, ".png")))
  Seurat::FeaturePlot(data, features = "scDblFinder.score") +
    ggplot2::labs(subtitle = paste0(
      "singlets: ",
      table(orig_data$scDblFinder.class)[1],
      " - doublets: ",
      paste0(table(orig_data$scDblFinder.class)[2],
      " - doublet rate: ",
      round(table(orig_data$scDblFinder.class)[2]/table(orig_data$scDblFinder.class)[1], digits = 3))))
  dev.off()

  # keep only singlets for downstream processing
  orig_data <- orig_data[,data$scDblFinder.class == "singlet"]
  return(orig_data)
}
plot_and_remove_doublets(data, sample_path, sample_name)
