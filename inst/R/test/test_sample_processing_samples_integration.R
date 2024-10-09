

## implemented and merged with sample_analysis
NSM.data <- Seurat::Read10X(data.dir = "C:/SynologyDrive/Projects/scRNAseqR/data/samples/sakshi NGN2 microglia/NSM/filtered_feature_bc_matrix")
M.data <- Seurat::Read10X(data.dir = "C:/SynologyDrive/Projects/scRNAseqR/data/samples/sakshi NGN2 microglia/M/filtered_feature_bc_matrix")
NS.data <- Seurat::Read10X(data.dir = "C:/SynologyDrive/Projects/scRNAseqR/data/samples/sakshi NGN2 microglia/NS/filtered_feature_bc_matrix")
NC.data <- Seurat::Read10X(data.dir = "C:/SynologyDrive/Projects/scRNAseqR/data/samples/sakshi NGN2 microglia/NC/filtered_feature_bc_matrix")
NSM.data <- Seurat::CreateSeuratObject(counts = NSM.data, project = "NSM", min.cells = 3, min.features = 500)
M.data <- Seurat::CreateSeuratObject(counts = M.data, project = "M", min.cells = 3, min.features = 500)
NS.data <- Seurat::CreateSeuratObject(counts = NS.data, project = "NS", min.cells = 3, min.features = 500)
NC.data <- Seurat::CreateSeuratObject(counts = NC.data, project = "NC", min.cells = 3, min.features = 500)
data.list <- list(NSM.data, M.data, NS.data, NC.data)


# TODO works, need?
# merged.test <- CreateSeuratObject(counts = c("NSM" = NSM.data@assays$RNA$counts, "NS" = NS.data@assays$RNA$counts))

# doublet removal
if (T) {
  plot_and_remove_doublets <- function(data, sample_path, sample_name, doublet_removal_rate = NULL) {
    temp_QC_data <- data
    
    set.seed(1)
    sce <- scDblFinder::scDblFinder(temp_QC_data@assays$RNA$counts, clusters = TRUE, dbr = doublet_removal_rate)
    temp_QC_data$scDblFinder.score <- sce$scDblFinder.score
    temp_QC_data$scDblFinder.class <- sce$scDblFinder.class
    
    temp_QC_data <- Seurat::NormalizeData(temp_QC_data)
    temp_QC_data <- Seurat::ScaleData(temp_QC_data)
    temp_QC_data <- Seurat::FindVariableFeatures(temp_QC_data)
    temp_QC_data <- Seurat::RunPCA(temp_QC_data, features = SeuratObject::VariableFeatures(object = temp_QC_data), npcs = 50, verbose = FALSE)
    temp_QC_data <- Seurat::RunUMAP(temp_QC_data, reduction = "pca", dims = 1:30)
    
    p <- Seurat::FeaturePlot(temp_QC_data, features = "scDblFinder.score") +
      ggplot2::labs(subtitle = paste0(
        "singlets: ",
        table(temp_QC_data$scDblFinder.class)[1],
        " - doublets: ",
        paste0(table(temp_QC_data$scDblFinder.class)[2],
               " - doublet rate: ",
               round(table(temp_QC_data$scDblFinder.class)[2]/table(temp_QC_data$scDblFinder.class)[1], digits = 3))))
    ggplot2::ggsave(file=file.path(sample_path, "quality_control", paste0("scDblFinder_scores_", sample_name, ".png")), width = 30, height = 20, units = "cm")
    
    # keep only singlets for downstream processing
    data <- data[, temp_QC_data$scDblFinder.class == "singlet"]
    return(data)
  }
  data <- plot_and_remove_doublets(data, sample_path, sample_name, doublet_removal_rate)
}

# ambient RNA removal
if (T) {
  clean_ambient_RNA <- function(data, samples_dir, sample_name) {
    table_of_droplets = Seurat::Read10X(data.dir = file.path(samples_dir, sample_name, "raw_feature_bc_matrix"))
    
    # add non-overlapping genes tot table_of_droplets with 0 counts to satisfy having same amount of genes
    non_overlapping_genes <- rownames(data@assays$RNA$counts)[which(!rownames(data@assays$RNA$counts) %in% rownames(table_of_droplets))]
    l <- list()
    l <- t(sapply(non_overlapping_genes, function(gene) {
      l[[gene]] <- rep(0, length(colnames(table_of_droplets)))
    }))
    colnames(l) <- colnames(table_of_droplets)
    table_of_droplets <- rbind(table_of_droplets, l)
    
    overlapping_genes <- rownames(data@assays$RNA$counts)[which(rownames(data@assays$RNA$counts) %in% rownames(table_of_droplets))]
    table_of_droplets <- table_of_droplets[overlapping_genes,]
    data@assays$RNA$counts <- data@assays$RNA$counts[overlapping_genes,]
    
    sc <- SoupX::SoupChannel(tod = table_of_droplets, toc = data@assays$RNA$counts) # estimateSoup()
    
    graphclust_dir <- file.path(samples_dir, sample_name, "analysis", "clustering")
    graphclust_dir <- list.files(graphclust_dir, full.names = T)[grep(pattern = "graphclust$", list.files(graphclust_dir))]
    tenx_graphclust <- read.csv(file.path(graphclust_dir, "clusters.csv"))
    tenx_clusters <- tenx_graphclust$Cluster[tenx_graphclust$Barcode %in% colnames(data@assays$RNA$counts)]
    sc <- SoupX::setClusters(sc, tenx_clusters)
    
    png(filename = file.path(sample_path, "quality_control", paste0("SoupX_contamination_", sample_name, ".png")))
    sc <- SoupX::autoEstCont(sc, doPlot = TRUE) # if fails, contamination rate default: 10% -> 0.1, or determine gene (sets) to estimate fraction
    dev.off()
    
    
    out <- SoupX::adjustCounts(sc, roundToInt = FALSE) # TODO set roundToInt = T is downstream processing algortihms need integers
    colnames(out) <- colnames(data@assays$RNA$counts)
    
    data@misc[["tenx_counts"]] <- data@assays$RNA$counts
    data@assays$RNA$counts <- out
    data@misc[["SoupX_contamination_percentage"]] <- sc$fit$rhoEst * 100
    return(data)
  }
  data <- clean_ambient_RNA(data, samples_dir, sample_name)
}

options(future.globals.maxSize = 8000 * 1024^2)
for (i in seq_along(data.list)) {
  so <- data.list[[i]]
  so_mt_features <- grep("^MT-", rownames(so@assays$RNA), value = TRUE)
  so <- Seurat::PercentageFeatureSet(so, features = so_mt_features, col.name = "percent.mt")
  ## return.only.var.genes = FALSE (https://github.com/satijalab/seurat/issues/4896)
  so <- Seurat::SCTransform(so, vars.to.regress = "percent.mt", method = "glmGamPoi", return.only.var.genes = FALSE)
  data.list[[i]] <- so
}

## list .rds files and make data.list with loading them
sample_files <- c(
  "C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV6/NSM/NSM.rds",
  "C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV6/NS/NS.rds",
  "C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV6/NC/NC.rds",
  "C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_pipeV6/M/M.rds"
)
data.list <- lapply(X = sample_files, FUN = function(x) {
  readRDS(file = x)
})

## implemented and merged with samples_integration
data.features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
# devnote: using base::merge instead of SeuratObject::merge.Seurat as the latter lost cells
data.merged <- base::merge(data.list[[1]], data.list[2:length(data.list)])
Seurat::VariableFeatures(data.merged) <- data.features
data.merged <- Seurat::RunPCA(object = data.merged, assay = "SCT", features = data.features, npcs = 50)


data.merged <- Seurat::IntegrateLayers(
  object = data.merged,
  method = Seurat::RPCAIntegration,
  orig.reduction = "pca",
  normalization.method = "SCT",
  new.reduction = "integrated.dr",
  verbose = TRUE,
  dims = 1:50)

# data.merged <- harmony::RunHarmony(
#   object = data.merged,
#   group.by.vars = "orig.ident",
#   assay.use = "SCT",
#   reduction.use = "pca",
#   reduction.save = "integrated.dr",
#   dims.use = 1:50,
#   plot_convergence = FALSE,
#   epsilon.cluster=-Inf,
#   epsilon.harmony=-Inf)

data.merged <- Seurat::FindNeighbors(object = data.merged, assay = "SCT", reduction = "integrated.dr", dims = 1:50)
data.merged <- Seurat::FindClusters(object = data.merged, resolution = 0.8)
data.merged <- Seurat::PrepSCTFindMarkers(data.merged, assay = "SCT")
data.merged <- Seurat::RunUMAP(object = data.merged, assay = "SCT", reduction = "integrated.dr", dims = 1:50)
Seurat::DimPlot(data.merged, reduction = "umap", group.by = c("orig.ident"))
# Seurat::DimPlot(data.merged, reduction = "umap", group.by = c("seurat_clusters"))
# saveRDS(data.merged, file = "harmony.rds")