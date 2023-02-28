library(Seurat)
library(data.table)
library(scuttle)

# DEVNOTE: can't train first and then intersect to common genes space (sadly)



### USER PARAMETERS ###
# set amount of chunks, higher = more chunks = less memory per iteration
n_chunks <- 25

# set seed for randomness in SingleR aggregation
set.seed(0)

# get start time
startTime <- format(Sys.time(), "%F %H-%M-%S")
# set work dir
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein'
setwd(work_dir)
# set output folder
RData.chunks.output.folder <- paste0("RData/", startTime, "_chunks/")
dir.create(RData.chunks.output.folder)
### END USER PARAMETERS ###



# get start time to measure run time
start_run_time <- Sys.time()

# get custom meta data from file, create from original if doesn't exist yet
if (file.exists("custom.meta.tsv")) {
  meta <- utils::read.table("custom.meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
} else {
  meta <- utils::read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
  # create custom merged cluster annotations for clusterv2 from reference data
  anno_df <- read.csv2("custom_annotation.txt")
  meta$custom.clusterv2 <- plyr::mapvalues(meta$clusterv2, anno_df$from, anno_df$to)

  # save custom annotation for use in processing and visualization
  utils::write.table(meta, file = "custom.meta.tsv", sep="\t", row.names=T)
}

# get genes, if non-existent create from reference data
if (file.exists("genes.csv") & file.exists("n_cols.txt")) {
  genes = utils::read.table("genes.csv", sep = "\t", col.names = 'gene')
  genes = genes$gene
  n_cols <- read.csv2("n_cols.txt", header = F)
  n_cols <- n_cols$V1
} else {
  # read only columns for efficiency
  initial <- data.table::fread("exprMatrix.tsv.gz", nrows=0, colClasses=c("gene"="character"))
  n_cols <- dim(initial)[2]
  # read first row only for efficiency
  initial <- data.table::fread("exprMatrix.tsv.gz", select=c(1), colClasses=c("gene"="character"))
  genes = initial[,1][[1]]
  genes = gsub(".+[|]", "", genes)
  utils::write.csv2(genes, file = paste0("genes.csv"), row.names = FALSE)
  write(n_cols, "n_cols.txt")
}

# set constants for data chunking
n_cells <- n_cols - 1
cell_ids <- 2:n_cols
n_rows <- length(genes)
n_sample <- ceiling(n_cells/n_chunks)

# print initial run time and memory usage
message('pre-chunk run time: ', Sys.time()-start_run_time)
message('point 1 mem ', utils::memory.size(), ' - ', utils::memory.size(max=TRUE))
# reset run time for iterations
start_run_time <- Sys.time()
# initialize sampled_cells (to store sampled cell_ids)
sampled_cells <- vector()

# iterate to chunk Kriegstein data
for (i in 1:n_chunks) {
  message("iter ", i, " start")
  start_iter_time <- Sys.time()

  # for final iteration, set amount of cell_ids to sample to amount of leftover cell_ids
  if (length(cell_ids)-length(sampled_cells) < n_sample) {
    n_sample <- length(cell_ids) - length(sampled_cells)
  }
  # sample cell_ids unsampled before
  sample_cells <- sample(cell_ids[!cell_ids %in% sampled_cells], n_sample, replace = F)
  if (length(sample_cells) < 2) {
    message('WARNING: sample taken of less then 2 cell_ids -> skipping iteration')
    next
  }
  # store sampled cells
  sampled_cells <- c(sampled_cells, sample_cells)

  message("load sampled cells data")
  cell_data <- data.table::fread("exprMatrix.tsv.gz", select = c(unique(sample_cells)), colClasses=c("gene"="character"))
  rownames(cell_data) <- genes
  # convert table to SO and add metadata matching by cell name
  cell_data <- SeuratObject::CreateSeuratObject(counts = cell_data, meta.data = meta[names(cell_data),])
  cell_data <- Seurat::as.SingleCellExperiment(cell_data)
  # SingleR() expects REFERENCE datasets to be normalized and log-transformed
  cell_data <- scuttle::logNormCounts(cell_data)

  save(cell_data, file = paste0(RData.chunks.output.folder, "chunk.", i, ".RData"))
  # remove cell data before next iteration to save memory
  rm(cell_data)

  message("Iteration ", i ," runtime was ", Sys.time()-start_iter_time, " minutes")
  }
message("Total run time: ", Sys.time() - start_run_time)
