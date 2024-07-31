library(Seurat)
library(EmptyNN)
library(tensorflow)

tensorflow::install_tensorflow(python_version="4.4")

# data.data <- Seurat::Read10X(data.dir = file.path(samples_dir, sample_name, "filtered_feature_bc_matrix"), strip.suffix = TRUE)
# data <- Seurat::CreateSeuratObject(counts = data.data, project = sample_name, min.cells = 3, min.features = 700)


# download .hdf5 dataset from: https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neurons_900?
counts <- Seurat::Read10X_h5("C:/Users/mauri/Desktop/neurons_900_raw_gene_bc_matrices_h5.h5", use.names = TRUE, unique.features = TRUE)

# Transpose the count matrix, so rows are cells and columns are genes
counts <- t(counts)

# Run emptynn()
nn.res <- EmptyNN::emptynn(counts, threshold = 100, k = 10, iteration = 10, verbose = TRUE)

# Downstream analysis
retained <- runSeurat(counts = counts[, nn.res$nn.keep], resolution = 0.2)
DimPlot(retained,reduction = 'tsne') + ggtitle("EmptyNN") + NoLegend()
