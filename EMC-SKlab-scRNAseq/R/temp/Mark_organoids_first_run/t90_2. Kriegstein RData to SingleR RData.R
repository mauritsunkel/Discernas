library(Seurat)
library(SingleR)



### USER PARAMETERS ###
# set annotations for training
annotations <- c("age", "structure", "custom.clusterv2")

# set sample names and files
sample_names <- c("t90")
rds.files <- c(
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Mark organoids first run/t90/t90.rds")

# get start time
startTime <- format(Sys.time(), "%F %H-%M-%S")

# set work and output dirs
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein'
setwd(work_dir)
RData.folder <- paste0(getwd(), "/RData/chunks_25/")
RData.output.folder <- paste0("RData/SingleR_RData_", startTime, "/")
dir.create(RData.output.folder)
### END USER PARAMETERS ###



# initialize and get genes and meta functions and data
getGenes <- function() {
  genesFile <- "genes.csv"
  if (!file.exists(genesFile)) {
    stop(genesFile, " does not exist, generate this file during reference data chunking.")
  } else {
    genes <- utils::read.table(genesFile, sep = "\t", col.names = 'gene')
    return(genes$gene[2:length(genes$gene)])
  }
}
getMeta <- function() {
  metaFile <- "custom.meta.tsv"
  if (!file.exists(metaFile)) {
    stop(metaFile, " does not exist, generate this file during reference data chunking.")
  } else {
    return(utils::read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
  }
}
genes <- getGenes()
meta <- getMeta()

# iterate files and perform SingleR for annotation with Pearson correlation
for (j in 1:length(rds.files)) {
  # initialize iteration timer
  start_iter_time <- Sys.time()

  sample_name <- sample_names[j]
  message("start iteration of ", sample_name)

  sample_data <- readRDS(file = rds.files[j])
  # set RNA assay to have genes in same feature space as the reference data
  DefaultAssay(sample_data) <- "RNA"

  # get overlapping genes between data and reference
  sample_genes <- rownames(sample_data)
  overlapping_genes <- intersect(sample_genes, genes)

  # add overlapping genes to .rds and set back to SCT assay
  SeuratObject::Misc(object = sample_data, slot = "Kriegstein.gene.overlap.assayRNA") <- overlapping_genes
  DefaultAssay(sample_data) <- "SCT"
  saveRDS(sample_data, file = rds.files[j])

  message("create Seurat sample -> SCE")
  sample_data <- Seurat::as.SingleCellExperiment(sample_data, assay = 'RNA')
  sample_data <- sample_data[overlapping_genes,]
  # perform transformation to have genes in the same feature space as the reference data
  sample_data <- scuttle::logNormCounts(sample_data)

  for (i in 1:length(list.files(path = RData.folder))) {
    message("start iteration of ", sample_name, " at cell RData chunk:  ", RData.folder, "iter.", i, ".RData")

    load(paste0(RData.folder, "chunk.", i, ".RData")) # cell_data (R object name)
    cell_data <- cell_data[overlapping_genes,]

    # set filename base for RData saving
    filename_base <- paste0(RData.output.folder, sample_name, ".iter.", i)
    # run SingleR for each annotation and save RData
    for (annotation in annotations) {
      message("sample=", sample_name, " iter=", i, " annotation=", annotation)
      result <- SingleR::SingleR(
        test = sample_data,
        ref = cell_data,
        labels = SummarizedExperiment::colData(cell_data)[, annotation],
        clusters = SummarizedExperiment::colData(sample_data)[, "seurat_clusters"],
        de.method = 'wilcox',
        aggr.ref = FALSE)
      save(result, file = paste0(filename_base, ".", annotation, ".RData"))
      # remove result before next iteration to save memory
      rm(result)
    }
  }
  message("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes - memory in use:", as.data.frame(gc())$'(Mb)'[2])
  message('point 1 mem ', memory.size(), ' ', memory.size(max=TRUE))
}
