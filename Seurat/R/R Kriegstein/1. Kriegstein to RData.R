library(Seurat)
library(SingleR)
library(pheatmap)
library(data.table)

start_run_time <- Sys.time()




# TODO could load in earlier defined overlapping genes instead of redoing that all the time
## TODO build in flag if file (by filename) exists then load, otherwise load rds.file and perform overlap
### TODO otherwise no need for loading RDS files at this step anymore (and all downstream processing of those)
#### TODO also since Seurat genes are identical over samples, need only one RDS file instead one per sample

# TODO there is hardcoding in this script because of testing and memory issues, if published, adjust and redo


# load our data (as test data)
rds.files <- c("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_A/BL_A.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_N/BL_N.rds")
data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Kriegstein'
setwd(work_dir)
# get meta
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# get nrows
# initial <- data.table::fread("exprMatrix.tsv.gz", nrows=0, colClasses=c("gene"="character"))
n_cols <- 691929 # dim(initial)[2]
n_cells <- n_cols - 1
cells <- 2:n_cols
# get genes
# initial <- fread("exprMatrix.tsv.gz", select=c(1), colClasses=c("gene"="character"))
# genes = initial[,1][[1]]
# genes = gsub(".+[|]", "", genes)
genes = read.table("genes.csv", sep = "\t", col.names = 'gene')
genes = genes$gene
# set #rows
n_rows <- length(genes)

# get overlapping genes between test and train dataset (to reduce noise and runtime)
genes.overlap <- rownames(data.list[[1]]@assays$SCT)[rownames(data.list[[1]]@assays$SCT) %in% genes]
which.genes.overlap <- which(rownames(data.list[[1]]@assays$SCT) %in% genes)
# save non-overlap genes for retrospection
## most genes contain a: ['-', '_', '.', 'orf']
genes.non.overlap <- sort(genes[!genes %in% rownames(data.list[[1]]@assays$SCT)])
write.csv2(genes.non.overlap, file = paste0("Label_transfer_non_overlapping_ref_genes.csv"))

data.list <- lapply(X = data.list, FUN = function(x) {
  # subset test data by overlapping genes
  x <- subset(x, features = genes.overlap)
  # set SO to SCE object
  x <- Seurat::as.SingleCellExperiment(x, assay = 'SCT')
  return(x)
})

# total RAM needed is around 96GB (my laptop has 16GB)
n_iter <- 40 # higher = more chunks for less memory per iteration
n_sample <- ceiling(n_cells/n_iter)

print(Sys.time()-start_run_time)
cat('point 1 mem', memory.size(), memory.size(max=TRUE), '\n')
start_total_runtime <- Sys.time()
sampled_cells <- c()
# initialize resulting list structure
results.list <- list()
for (i in 1:length(rds.files)) {
  results.list[[i]] <- list()
}
for (i in 1:n_iter) {
  start_iter_time <- Sys.time()
  print("iter start")

  # for final iteration, set amount of cells to sample to amount of leftover cells
  if (length(cells)-length(sampled_cells) < n_sample) {
    n_sample <- length(cells)-length(sampled_cells)
  }
  # sample cells
  if (i == 1) {
    sample_cells <- sample(cells, n_sample, replace = F)
  } else {
    # TODO needs to be with ! to remove already sampled cells right?
    ## TODO if that is true, then why did it work before?
    sample_cells <- sample(cells[-sampled_cells], n_sample, replace = F)
  }
  if (length(sample_cells) < 2) {
    print('WARNING: sample taken of less then 2 cells -> skipping iteration')
    next
  }
  sampled_cells <- c(sample_cells, sampled_cells)

  print("load cell data")
  # load cell data
  cell_data <- data.table::fread("exprMatrix.tsv.gz", select = c(unique(sample_cells)), colClasses=c("gene"="character"))

  # TODO if keeps crashing before full run, then need to try the non-sampling linear solution for data saving
  ## look into skip= and nrows= !!!
  ## also look at cmd= parameter and think about being able to only take certain genes first
  ### https://gitlab.com/jozefhajnala/fread-benchmarks/-/blob/master/rscripts/04_fread_grep.R
  save(cell_data, file = paste0("RData/Kriegstein_RData/iter.", i, ".RData"))

  # print size of all current objects in the environment
  for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

  # remove from environment to save memory for next iteration
  rm(cell_data)

  print(paste("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes - memory in use:", as.data.frame(gc())$'(Mb)'[2]))
  cat('point 1 mem', memory.size(), memory.size(max=TRUE), '\n')
  }
print(paste("Total run time:", Sys.time() - start_total_runtime))
## cells per iteration
# 100: 4.2-4.6 mins - 9.1 total
# 1000: 4.7-4.8 - 9.6 total      - mem: 992-1000
# 10.000: 5-5 - 10.1 total       - mem: 1274-1392
# 100.000 crash?
# 50.000: Error: cannot allocate vector of size 8.0 Gb
# 20.000 forgot to note
# 30.000 (3-iter): Error: cannot allocate vector of size 4.8 Gb
# 15.000 (3-iter): 6.6-8.0-5.8     1530-1535-1535
# 10.000 (2-iter, all samples): 5.4-5.6 - 2220-2265
# 20.000 (2-iter, all samples): 7.1-5.8 - 2222-2221
