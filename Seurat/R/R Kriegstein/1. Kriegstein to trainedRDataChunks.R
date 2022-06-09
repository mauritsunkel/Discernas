library(Seurat)
library(SingleR)
library(pheatmap)
library(data.table)



### USER PARAMETERS ###
# set seed (SingleR aggregation has randomness)
set.seed(0)

# get start time
startTime <- format(Sys.time(), "%F %H-%M-%S")

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein'
setwd(work_dir)

# set saving RData chunks
saveChunks <- TRUE
# higher = more chunks (less memory per iteration)
n_chunks <- 25
# set training RData chunks
trainData <- FALSE

if (saveChunks) {
  RData.chunks.output.folder <- paste0("RData/", startTime, "_chunks_test/")
  dir.create(RData.chunks.output.folder)
}
if (trainData) {
  RData.output.folder <- "RData/SingleR_trainedData_customReference_20-05-2022/"
  dir.create(RData.output.folder)
}

# set annotations for training
annotations <- c("age", "structure", "custom.clusterv2")
### END USER PARAMETERS ###

# get script start time to print runtime later
start_run_time <- Sys.time()

# get meta
if (file.exists("custom.meta.tsv")) {
  meta <- read.table("custom.meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
} else {
  meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
  # create custom merged cluster annotations for clusterv2
  meta$custom.clusterv2 <- sapply(meta$clusterv2, function(x) {
    switch(x,
           "Neuron_8" = "Excitatory neuron 1",
           "Neuron_28" = "Excitatory neuron 1",
           "Neuron_31" = "Excitatory neuron 1",
           "Neuron_5" = "Excitatory neuron 1",
           "Neuron_6" = "Excitatory neuron 1",
           "Neuron_92Neuron" = "Excitatory neuron 1",
           "Neuron_93Neuron" = "Excitatory neuron 1",
           "Neuron_4Neuron" = "Excitatory neuron 1",
           "Neuron_29" = "Excitatory neuron 2",
           "Neuron_30" = "Excitatory neuron 2",
           "Neuron_35Neuron" = "Excitatory neuron 2",
           "Neuron_13" = "Excitatory neuron 3",
           "Neuron_16" = "Excitatory neuron 3",
           "Neuron_15" = "Excitatory neuron 3",
           "Neuron_14" = "Excitatory neuron 3",
           "Neuron_22" = "Excitatory neuron 3",
           "Neuron_20" = "Excitatory neuron 3",
           "Neuron_21" = "Excitatory neuron 3",
           "Neuron_18" = "Excitatory neuron 3",
           "Neuron_17" = "Excitatory neuron 3",
           "Neuron_23" = "Excitatory neuron 3",
           "Neuron_24" = "Excitatory neuron 3",
           "Neuron_26" = "Excitatory neuron 3",
           "Neuron_25" = "Excitatory neuron 3",
           "Neuron_27" = "Excitatory neuron 3",
           "Neuron_36" = "Excitatory neuron 3",
           "Neuron_37" = "Excitatory neuron 3",
           "Neuron_19" = "Excitatory neuron 3",
           "Neuron_35" = "Excitatory neuron 3",
           "Neuron_9" = "Excitatory neuron 3",
           "Neuron_63Neuron" = "Excitatory neuron 3",
           "Neuron_100Neuron" = "Excitatory neuron 3",
           "Neuron_68Neuron" = "Excitatory neuron 3",
           "Neuron_98Neuron" = "Excitatory neuron 3",
           "Neuron_102Neuron" = "Excitatory neuron 3",
           "Neuron_53Neuron" = "Excitatory neuron 3",
           "Neuron_101Neuron" = "Excitatory neuron 3",
           "Neuron_72Neuron" = "Excitatory neuron 3",
           "Neuron_50Neuron" = "Excitatory neuron 3",
           "Neuron_5Neuron" = "Excitatory neuron 3",
           "IPC_18IPC" = "IPC",
           "IPC_29IPC" = "IPC",
           "IPC_34IPC" = "IPC",
           "Neuron_4" = "Excitatory neuron 4",
           "Neuron_7" = "Excitatory neuron 4",
           "Neuron_32" = "Excitatory neuron 4",
           "Neuron_33" = "Excitatory neuron 4",
           "Neuron_34" = "Excitatory neuron 4",
           "Neuron_33Neuron" = "Excitatory neuron 4",
           "Neuron_39Neuron" = "Excitatory neuron 4",
           "Neuron_2Neuron" = "Excitatory neuron 4",
           "Neuron_66Neuron" = "Excitatory neuron 4",
           "GW19_2_29NeuronNeuron" = "Excitatory neuron 4",
           "GW18_2_42NeuronNeuron" = "Excitatory neuron 4",
           "Neuron_33Neuron" = "Excitatory neuron 4",
           "GW19_2_1OutlierrOutlierr" = "Outlier (in excitatory neuron 4)",
           "GW19_2_45OutlierrOutlierr" = "Outlier (in excitatory neuron 4)",
           "Neuron_1" = "Excitatory neuron 5",
           "Neuron_2" = "Excitatory neuron 5",
           "Neuron_3" = "Excitatory neuron 5",
           "Neuron_12" = "Excitatory neuron 5",
           "Neuron_10Neuron" = "Excitatory neuron 5",
           "Neuron_37Neuron" = "Excitatory neuron 5",
           "Neuron_11Neuron" = "Excitatory neuron 6",
           "Neuron_36Neuron" = "Excitatory neuron 6",
           "Neuron_64Neuron" = "Excitatory neuron 6",
           "Neuron_10" = "Excitatory neuron 7",
           "Neuron_11" = "Excitatory neuron 7",
           "Interneuron_6" = "Interneuron 1",
           "Interneuron_9" = "Interneuron 1",
           "Interneuron_31Interneuron" = "Interneuron 1",
           "Interneuron_1" = "Interneuron 2",
           "Interneuron_2" = "Interneuron 2",
           "Interneuron_4" = "Interneuron 2",
           "Interneuron_3" = "Interneuron 2",
           "Interneuron_5" = "Interneuron 2",
           "Interneuron_7" = "Interneuron 2",
           "Interneuron_10" = "Interneuron 2",
           "Interneuron_8" = "Interneuron 2",
           "Interneuron_24Interneuron" = "Interneuron 2",
           "Interneuron_43Interneuron" = "Interneuron 2",
           "GW20_Interneuron_32InterneuronInterneuron" = "Interneuron 2",
           "GW20_Interneuron_34InterneuronInterneuron" = "Interneuron 2",
           "Endo_1" = "Endo 1",
           "Endo_2" = "Endo 1",
           "Endo_3" = "Endo 1",
           "Endo_4" = "Endo 1",
           "Endo_4Endo" = "Endo 1",
           "Endo_26Endo" = "Endo 1",
           "Endo_5" = "Endo 2",
           "Endo_6" = "Endo 2",
           "Endo_11" = "Endo 2",
           "Endo_34Endo" = "Endo 2",
           "Endo_7" = "Endo 3",
           "Endo_8" = "Endo 3",
           "Endo_9" = "Endo 3",
           "Endo_10" = "Endo 3",
           "Endo_12" = "Endo 3",
           "Endo_13" = "Endo 3",
           "Endo_14" = "Endo 3",
           "Endo_15" = "Endo 3",
           "Endo_16" = "Endo 3",
           "GW17_Vascular_14VascularVascular" = "Vascular",
           "GW17_Vascular_21VascularVascular" = "Vascular",
           "Oligo_2" = "Oligodendrocyte 1",
           "Oligo_3" = "Oligodendrocyte 1",
           "Oligo_4" = "Oligodendrocyte 1",
           "Oligo_6" = "Oligodendrocyte 1",
           "Oligo_1" = "Oligodendrocyte 2",
           "Oligo_5" = "Oligodendrocyte 2",
           "Oligo_7" = "Oligodendrocyte 2",
           "Oligo_9" = "Oligodendrocyte 2",
           "GW18_2_45OligoOligo" = "Oligodendrocyte 2",
           "Oligo_8" = "Oligodendrocyte 2",
           "Oligo_11" = "Oligodendrocyte 2",
           "Oligo_12" = "Oligodendrocyte 2",
           "Microglia_1" = "Microglia",
           "Microglia_2" = "Microglia",
           "Microglia_3" = "Microglia",
           "Microglia_5" = "Microglia",
           "Microglia_6" = "Microglia",
           "Microglia_4" = "Microglia",
           "Microglia_10Microglia" = "Microglia",
           "Microglia_16Microglia" = "Microglia",
           "Microglia_76Microglia" = "Microglia",
           "Microglia_15Microglia" = "Microglia",
           "Microglia_8Microglia" = "Microglia",
           "GW19_2_43MicrogliaMicroglia" = "Microglia",
           "GW14_Microglia_27MicrogliaMicroglia" = "Microglia",
           "Astrocyte_1" = "Astrocyte 1",
           "Astrocyte_2" = "Astrocyte 2",
           "Astrocyte_3" = "Astrocyte 3",
           "Astrocyte_4" = "Astrocyte 3",
           "Astrocyte_5" = "Astrocyte 3",
           "Astrocyte_6" = "Astrocyte 3",
           "GW18_Astrocyte_17AstrocyteAstrocyte" = "Astrocyte 3",
           "GW18_Astrocyte_51AstrocyteAstrocyte" = "Astrocyte 3",
           "GW18_Astrocyte_26AstrocyteAstrocyte" = "Astrocyte 3",
           "GW19_Astrocyte_53AstrocyteAstrocyte" = "Astrocyte 3",
           "GW20_Astrocyte_39AstrocyteAstrocyte" = "Astrocyte 3",
           "GW19_Astrocyte_35AstrocyteAstrocyte" = "Astrocyte 3",
           "RG_1" = "Radial glia 1",
           "RG_17" = "Radial glia 1",
           "RG_12" = "Radial glia 1",
           "RG_13RG" = "Radial glia 1",
           "RG_37RG" = "Radial glia 1",
           "RG_6" = "Radial glia 2",
           "RG_7" = "Radial glia 2",
           "RG_9" = "Radial glia 2",
           "RG_11RG" = "Radial glia 2",
           "RG_59RG" = "Radial glia 2",
           "RG_10" = "Radial glia 2",
           "RG_52RG" = "Radial glia 2",
           "GW19_2_30RGRG" = "Radial glia 2",
           "RG_2RG" = "Radial glia 2",
           "RG_4" = "Radial glia 2",
           "RG_13" = "Radial glia 2",
           "RG_14" = "Radial glia 2",
           "RG_46RG" = "Radial glia 2",
           "GW22both_RG_38RGRG" = "Radial glia 2",
           "GW22both_RG_23RGRG" = "Radial glia 2",
           "RG_5" = "Radial glia 2",
           "RG_47RG" = "Radial glia 2",
           "RG_22RG" = "Radial glia 2",
           "RG_3RG" = "Radial glia 2",
           "RG_4RG" = "Radial glia 2",
           "RG_8" = "Radial glia 2",
           "RG_3" = "Radial glia 2",
           "RG_45RG" = "Radial glia 2",
           "RG_50RG" = "Radial glia 2",
           "RG_11" = "Radial glia 2",
           "RG_16" = "Radial glia 2",
           "RG_15" = "Radial glia 2",
           "RG_15RG" = "Radial glia 2",
           "Dividing_1" = "Dividing",
           "Dividing_21Dividing" = "Dividing",
           "Dividing_2" = "Dividing",
           "Dividing_4" = "Dividing",
           "Dividing_5" = "Dividing",
           "Dividing_6" = "Dividing",
           "Dividing_8" = "Dividing",
           "Dividing_13Dividing" = "Dividing",
           "Dividing_7" = "Dividing",
           "Dividing_41Dividing" = "Dividing",
           "Dividing_3" = "Dividing",
           "Dividing_11" = "Dividing",
           "Dividing_12" = "Dividing",
           "Dividing_13" = "Dividing",
           "Dividing_10" = "Dividing",
           "Dividing_48Dividing" = "Dividing",
           "Dividing_15Dividing" = "Dividing",
           "Dividing_9" = "Dividing",
           "Dividing_14" = "Dividing",
           "Dividing_15" = "Dividing",
           "Dividing_16" = "Dividing",
           "Dividing_42Dividing" = "Dividing",
           "Dividing_18Dividing" = "Dividing",
           "GW19_Dividing_39DividingDividing" = "Dividing",
           x)
  })
  # save custom annotation to be used in processing and visualization
  write.table(meta, file = "custom.meta.tsv", sep="\t", row.names=T)
}

# if "genes.csv" exists, then n_rows, n_cells, and cells are defined in earlier run for this version of data
if (file.exists("genes.csv")) {
  n_cols <- 691929
  genes = read.table("genes.csv", sep = "\t", col.names = 'gene')
  genes = genes$gene
} else {
  # read in only columns for efficiency
  initial <- data.table::fread("exprMatrix.tsv.gz", nrows=0, colClasses=c("gene"="character"))
  n_cols <- dim(initial)[2]

  # read in first row only for efficiency
  initial <- fread("exprMatrix.tsv.gz", select=c(1), colClasses=c("gene"="character"))
  genes = initial[,1][[1]]
  genes = gsub(".+[|]", "", genes)
  write.csv2(genes, file = paste0("genes.csv"))
}
n_cells <- n_cols - 1
cell_ids <- 2:n_cols
n_rows <- length(genes)

n_sample <- ceiling(n_cells/n_chunks)

# print initial runtime
print(Sys.time()-start_run_time)
# print initial memory usage
cat('point 1 mem', memory.size(), memory.size(max=TRUE), '\n')
# initialize total runtime for iterations
start_total_runtime <- Sys.time()
# initialize sampled_cells (for cell_ids)
sampled_cells <- vector()

# start iterating for chunking Kriegstein data
for (i in 1:n_chunks) {
  start_iter_time <- Sys.time()
  print("iter start")

  # for final iteration, set amount of cell_ids to sample to amount of leftover cell_ids
  if (length(cell_ids)-length(sampled_cells) < n_sample) {
    n_sample <- length(cell_ids) - length(sampled_cells)
  }
  # sample cell_ids
  sample_cells <- sample(cell_ids[!cell_ids %in% sampled_cells], n_sample, replace = F)

  if (length(sample_cells) < 2) {
    print('WARNING: sample taken of less then 2 cell_ids -> skipping iteration')
    next
  }
  sampled_cells <- c(sample_cells, sampled_cells)

  print("load cell data")
  # load cell data
  ## TODO if keeps crashing before full run, then need to try the non-sampling linear solution for data saving
  ### look into skip= and nrows= !!!
  ### also look at cmd= parameter and think about being able to only take certain genes first
  #### https://gitlab.com/jozefhajnala/fread-benchmarks/-/blob/master/rscripts/04_fread_grep.R
  # TODO look into scan instead of fread? (and Pascal data format?)
  # read selected cells from all Kriegstein data
  cell_data <- data.table::fread("exprMatrix.tsv.gz", select = c(unique(sample_cells)), colClasses=c("gene"="character"))
  # add genes as rownames
  rownames(cell_data) <- genes
  # convert table to SO and add metadata (use cell names to match metadata stored in meta)
  cell_data <- SeuratObject::CreateSeuratObject(counts = cell_data, meta.data = meta[names(cell_data),])
  # convert SO to SCE object
  cell_data <- Seurat::as.SingleCellExperiment(cell_data)

  if (saveChunks) {
    save(cell_data, file = paste0(RData.chunks.output.folder, "chunk.", i, ".RData"))
  }

  # DEVNOTE: can't train first and then intersect to common genes space
  if (trainData) {
    for (annotation in annotations) {
      print(annotation)
      # TODO test if want to use aggr.ref param, delete other, don't keep both data saves
      print(paste0(annotation, ": ", start_iter_time-Sys.time()))
      trainedData <- SingleR::trainSingleR(ref = cell_data,
                                           labels = cell_data@colData[, annotation],
                                           recompute = TRUE,
                                           genes = "de")
      print(paste0("no aggr: ", start_iter_time-Sys.time()))
      save(trainedData, file = paste0(RData.output.folder, "chunk.", i, ".", annotation, ".RData"))
      trainedData <- SingleR::trainSingleR(ref = cell_data,
                                           labels = cell_data@colData[, annotation],
                                           recompute = TRUE,
                                           genes = "de",
                                           aggr.ref = T)
      print(paste0("aggr: ", start_iter_time-Sys.time()))
      save(trainedData, file = paste0(RData.output.folder, "chunk.", i, ".", annotation, ".aggregateRef.RData"))

      # remove trainedData before next iteration to save memory
      rm(trainedData)
    }
  }

  # remove cell data before next iteration to save memory
  rm(cell_data)

  print(paste("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes"))
  }
print(paste("Total run time:", Sys.time() - start_total_runtime))
