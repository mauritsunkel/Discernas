library(Seurat)
library(SingleR)
library(pheatmap)
library(data.table)



### USER PARAMETERS ###
# set annotations for training
annotations <- c("age", "structure", "custom.clusterv2")

# load our data (as test data)
sample_names <- c("BL_N + BL_C")
rds.files <- c(
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06/integrated - old selection/BL_N + BL_C/BL_N + BL_C.rds"
)
## pre-selection data
# sample_names <- c("BL_C", "BL_A", "BL_N", "BL_A + BL_C", "BL_N + BL_C")
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_C/BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_A/BL_A.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_N/BL_N.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_N + BL_C/BL_N + BL_C.rds"
# )
## post-selection data
# sample_names <- c("BL_A + BL_C", "BL_N + BL_C")
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/after_selection/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_N + BL_C/after_selection/BL_N + BL_C.rds"
# )


# get start time
startTime <- format(Sys.time(), "%F %H-%M-%S")

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein'
setwd(work_dir)
RData.folder <- paste0(getwd(), "/RData/chunks_25/")
RData.output.folder <- paste0("RData/SingleR_RData_", startTime, "/")
dir.create(RData.output.folder)
### END USER PARAMETERS ###



# initialize get genes function
getGenes <- function() {
  genesFile <- "genes.csv"
  if (!file.exists(genesFile)) {
    stop(genesFile, " does not exist, this file is generated during data chunking.")
  } else {
    genes <- read.table(genesFile, sep = "\t", col.names = 'gene')
    return(genes$gene)
  }
}
# initialize get meta function
getMeta <- function() {
  metaFile <- "custom.meta.tsv"
  if (!file.exists(metaFile)) {
    stop(metaFile, " does not exist, this file is generated during data chunking.")
  } else {
    return(read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
  }
}
# get genes and meta
genes <- getGenes()
meta <- getMeta()

# start iterating file processing and perform SingleR for correlation annotations
for (j in 1:length(rds.files)) {
  # initialize iteration timer
  start_iter_time <- Sys.time()
  # setup filename for saving RData
  sample_name <- sample_names[j]
  # print user iteration status
  message("start iteration of ", sample_name)

  # read in sample data
  sample_data <- readRDS(file = rds.files[j])

  # get sample genes
  sample_genes <- rownames(sample_data)
  # get genes overlapping with Kriegstein genes
  overlapping_genes <- intersect(sample_genes, genes)

  # add overlapping genes as misc(elleneous) annotation to Seurat object (sample_data@misc$Kriegstein.gene.overlap)
  SeuratObject::Misc(object = sample_data, slot = "Kriegstein.gene.overlap") <- overlapping_genes
  # overwrite rds file with new misc(elleneous) annotation (note: NOT metadata, as that is about cells here)
  saveRDS(sample_data, file = rds.files[j])

  # set SO to SCE object while getting raw counts
  message("create Seurat sample -> sce")
  sample_data <- Seurat::as.SingleCellExperiment(sample_data, assay = 'RNA')

  # subset after convert to sce
  sample_data <- sample_data[intersect(rownames(sample_data), overlapping_genes),]

  for (i in 1:length(list.files(path = RData.folder))) {
    # print user iteration status
    message("start iteration of ", sample_name, " at cell RData chunk:  ", RData.folder, "iter.", i, ".RData")

    # load cell data
    load(paste0(RData.folder, "chunk.", i, ".RData")) # cell_data (R object name)
    cell_data <- cell_data
    # subset only overlapping genes from cell data
    cell_data <- cell_data[overlapping_genes,]

    # set filename base for RData saving
    filename_base <- paste0(RData.output.folder, sample_name, ".iter.", i)
    for (annotation in annotations) {
      message("sample=", sample_name, " iter=", i, " annotation=", annotation)
      # run SingleR (train + classify) to generate RData
      result <- SingleR::SingleR(
        test = sample_data,
        ref = cell_data,
        labels = colData(cell_data)[, annotation],
        clusters = colData(sample_data)[, "seurat_clusters"],
        de.method = 'wilcox',
        aggr.ref = FALSE)
      save(result, file = paste0(filename_base, ".", annotation, ".RData"))

      # remove result before next iteration to save memory
      rm(result)
    }

  }
  for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

  print(paste("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes - memory in use:", as.data.frame(gc())$'(Mb)'[2]))
  cat('point 1 mem', memory.size(), memory.size(max=TRUE), '\n')
}
