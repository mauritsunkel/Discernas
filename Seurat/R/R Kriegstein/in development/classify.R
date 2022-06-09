library(Seurat)
library(SingleR)
library(pheatmap)
library(data.table)


### USER PARAMETERS ###

# set work, output and data directories
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
start_time <- format(Sys.time(), "%F %H-%M-%S")
out_dir <- paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/')
dir.create(out_dir, recursive = TRUE)
setwd(out_dir)
data_dir <- paste0(work_dir, 'data/Kriegstein/RData/')
SingleR.trained.input.folder <- paste0(data_dir, "SingleR_trainedDataChunks/")
SingleR.tested.output.folder <- paste0(data_dir, "SingleR_testedDataChunks/")
# TODO testing - keep only non-aggregated or aggregated later
dir.create(SingleR.tested.output.folder)
SingleR.tested.aggregated.output.folder <- paste0(data_dir, "SingleR_testedAggregatedDataChunks/")
dir.create(SingleR.tested.aggregated.output.folder)

# source beep from my_utils
source(file=paste0(work_dir, "R/my_utils/beep.R"))

# set annotations
annotations <- c("celltype", "age", "individual", "structure", "area", "cell.type.v2", "custom.clusterv2")

# TODO TESTING
# # set sample names (NOTE: same order as rds.files)
# sample_names <- c("BL_A", "BL_C", "BL_N", "BL_A + BL_C", "BL_N + BL_C")
# # set rds files ((NOTE: same order as samples))
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_A/BL_A.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_C/BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_N/BL_N.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_N + BL_C/BL_N + BL_C.rds"
# )
# set sample names (NOTE: same order as rds.files)
sample_names <- c("BL_A", "BL_C", "BL_N")
# set rds files ((NOTE: same order as samples))
rds.files <- c(
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_A/BL_A.rds",
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_C/BL_C.rds",
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_N/BL_N.rds"
)
# sample_names <- c("BL_N")
# rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/BL_N/BL_N.rds")
# set sample names on rds.files
names(rds.files) <- sample_names
### END USER PARAMETERS ###






# TODO testing
i <- 1

for (i in 1:length(rds.files)) {
  sample_data <- readRDS(file = rds.files[i])
  print(length(rownames(sample_data)))
  if (i == 1) {
    genes_superset <- rownames(sample_data)
  } else {
    genes_superset <- intersect(genes_superset, rownames(sample_data))
  }
}
length(genes_superset)
astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                        "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1")
neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
schema_psych_interest <- c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                           "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11")
sloan_2017_interest <- c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                         "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                         "HOPX", "FAM107A", "AGT")
eva_hnRNPC_interest <- c("HNRNPR", "HNRNPU-AS1", "HNRNPU", "HNRNPLL", "HNRNPA3", "HNRNPD",
                         "HNRNPAB", "HNRNPH1", "HNRNPA2B1", "HNRNPH2", "HNRNPDL")
eva_hnRNPC_interest_2 <- c("HNRNPK", "HNRNPUL2", "HNRNPF", "HNRNPH3", "HNRNPA1", "HNRNPA0",
                           "HNRNPA1L2", "HNRNPC", "HNRNPM", "HNRNPL", "HNRNPUL1")
astrocyte_maturity <- c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                        "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE")
neuron_maturity <- c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                     "SYN1", "SNAP25", "SYT1", "APOE")
neuron_interest %in% genes_superset









# initialize get genes function
getGenes <- function() {
  genesFile <- paste0(data_dir, "../genes.csv")
  if (!file.exists(genesFile)) {
    message(genesFile, " does not exist, this file is generated during data chunking.")
    stop()
  } else {
    genes <- read.table(genesFile, sep = "\t", col.names = 'gene')
    return(genes$gene)
  }
}
# initialize get meta function
getMeta <- function() {
  metaFile <- paste0(data_dir, "../custom.meta.tsv")
  if (!file.exists(metaFile)) {
    message(metaFile, " does not exist, this file is generated during data chunking.")
    stop()
  } else {
    return(read.table(metaFile, header=T, sep="\t", as.is=T, row.names=1))
  }
}
# get genes and meta
genes <- getGenes()
meta <- getMeta()



# TODO for testing
i <- 1

# start iterating file processing and perform SingleR for correlation annotations
for (i in 1:length(rds.files)) {
  # initialize iteration timer
  start_iter_time <- Sys.time()
  # setup filename for saving RData
  sample_name <- sample_names[i]
  # print user iteration status
  print(paste0("start iteration of ", sample_name))

  # read in sample data
  sample_data <- readRDS(file = rds.files[i])
  # get sample genes
  sample_genes <- rownames(sample_data@assays$SCT)

  # get genes overlapping with Kriegstein genes
  overlapping_genes <- intersect(sample_genes, genes)


  # add overlapping genes as misc(elleneous) annotation to Seurat object
  SeuratObject::Misc(object = sample_data, slot = "Kriegstein.gene.overlap") <- overlapping_genes
  # DEVNOTE: usage:
  ## Misc(sample_data)
  ## names(sample_data@misc)
  ## sample_data@misc$Kriegstein.gene.overlap
  # overwrite rds file with new misc(elleneous) annotation (note: NOT metadata, as that is about cells here)
  saveRDS(sample_data, file = rds.files[i])

  # subset test data by overlapping genes after saving those in miscellaneous slot
  print("subset sample data")
  sample_data <- subset(sample_data, features = overlapping_genes)
  # set SO to SCE object
  print("Seurat sample -> sce")
  sample_data <- Seurat::as.SingleCellExperiment(sample_data, assay = 'SCT')



  # TODO for loop here on data chunks
  ## TODO next if not (non-)aggregated




  # TODO check that rows of $original.exprs are cells (because 30741) just like amount of genes in genes
  ## TODO so perform overlapping_genes by genes not with trainedData, that comes later??
  # load trained data
  load(paste0(SingleR.trained.input.folder, "chunk.1.age.aggregateRef.RData")) # trainedData (R object name)
  . <- trainedData

  trainedData$original.exprs[["14"]]

  dim(trainedData$nn.indices[[1]])
  identical(test, trainedData$common.genes)
  test <- trainedData$common.genes
  max(as.integer(trainedData$common.genes))


  load("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein/RData/chunks/chunk.1.RData")
  cell_data
  # set genes as row names
  cell_data <- data.frame(cell_data, row.names=genes) # genes = 30741
  rownames(cell_data)
  dim(cell_data)








  ## set overlapping genes?
  # overlapping_genes <- intersect(trainedData$common.genes, rownames(sample_data))







  # subset common genes and indices of trainedData on overlapping genes with testData
  trainedData$nn.indices <- lapply(trainedData$nn.indices, FUN = function(x) {
    x@data <- x@data[trainedData$common.genes %in% overlapping_genes,]
    return(x)
  })
  trainedData$nn.indices <- SimpleList(trainedData$nn.indices)
  trainedData$search$extra <- lapply(trainedData$search$extra, lapply, intersect, y = overlapping_genes)
  trainedData$common.genes <- overlapping_genes

  # TODO fine.tune on and off for comparison (it definitely costs a lot of extra time)
  ## TODO prune on/off? (like fine.tune, or work together?)
  result <- SingleR::classifySingleR(test = sample_data,
                                     trained = trainedData)

  "COPE" %in% genes[as.integer(trainedData$common.genes)]



# TODO capture run time per file (chunk, annotation, trainParam (aggr.ref/genes.all))








  # for (j in 1:length(list.files(path = RData.folder))) {
  #   # print user iteration status
  #   print(paste0("start iteration of ", sample_name, " at cell RData chunk:  ", RData.folder, "iter.", j, ".RData"))
  #
  #   # set filename base for RData saving
  #   filename_base <- paste0(RData.output.folder, sample_name, ".iter.", j)
  #
  #   # load cell data
  #   load(paste0(RData.folder, "iter.", j, ".RData")) # cell_data (R object name)
  #   cell_data <- cell_data
  #   # set genes as row names
  #   cell_data <- data.frame(cell_data, row.names=genes)
  #
  #
  #
  #   print(cell_data)
  #
  #   # subset only overlapping genes from cell data
  #   cell_data <- cell_data[overlapping_genes_ids,]
  #
  # }
}





# TODO
## classify test (each)
## combineCommonResults
## visualize
## test same/similar label transfer with previous version, for both aggregated and un-aggregated
### keep one (in training script)
