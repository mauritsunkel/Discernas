library(Seurat)
library(SingleR)
library(pheatmap)
library(data.table)

### USER PARAMETERS ###
# load our data (as test data)
rds.files <- c("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_A/BL_A.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_N/BL_N.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/integrated/BL_A + BL_C/integrated_BL_A + BL_C.rds",
               "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/integrated/BL_N + BL_C/integrated_BL_N + BL_C.rds")
sample_names <- c("BL_A", "BL_C", "BL_N", "BL_A + BL_C", "BL_N + BL_C")

data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Kriegstein'
setwd(work_dir)
# get meta
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

RData.folder <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Kriegstein/RData/Kriegstein_RData/'

### END USER PARAMETERS ###



# create custom merged cluster annotations for clusterv2
meta$custom.clusterv2 <- sapply(meta$clusterv2, function(x) {
  switch(x,
         "Neuron_8" = "Neuron 1",
         "Neuron_28" = "Neuron 1",
         "Neuron_31" = "Neuron 1",
         "Neuron_5" = "Neuron 1",
         "Neuron_6" = "Neuron 1",
         "Neuron_92Neuron" = "Neuron 1",
         "Neuron_93Neuron" = "Neuron 1",
         "Neuron_4Neuron" = "Neuron 1",
         "Neuron_29" = "Neuron 2",
         "Neuron_30" = "Neuron 2",
         "Neuron_35Neuron" = "Neuron 2",
         "Neuron_13" = "Neuron 3",
         "Neuron_16" = "Neuron 3",
         "Neuron_15" = "Neuron 3",
         "Neuron_14" = "Neuron 3",
         "Neuron_22" = "Neuron 3",
         "Neuron_20" = "Neuron 3",
         "Neuron_21" = "Neuron 3",
         "Neuron_18" = "Neuron 3",
         "Neuron_17" = "Neuron 3",
         "Neuron_23" = "Neuron 3",
         "Neuron_24" = "Neuron 3",
         "Neuron_26" = "Neuron 3",
         "Neuron_25" = "Neuron 3",
         "Neuron_27" = "Neuron 3",
         "Neuron_36" = "Neuron 3",
         "Neuron_37" = "Neuron 3",
         "Neuron_19" = "Neuron 3",
         "Neuron_35" = "Neuron 3",
         "Neuron_9" = "Neuron 3",
         "Neuron_63Neuron" = "Neuron 3",
         "Neuron_100Neuron" = "Neuron 3",
         "Neuron_68Neuron" = "Neuron 3",
         "Neuron_98Neuron" = "Neuron 3",
         "Neuron_102Neuron" = "Neuron 3",
         "Neuron_53Neuron" = "Neuron 3",
         "Neuron_101Neuron" = "Neuron 3",
         "Neuron_72Neuron" = "Neuron 3",
         "Neuron_50Neuron" = "Neuron 3",
         "Neuron_5Neuron" = "Neuron 3",
         "IPC_18IPC" = "IPC",
         "IPC_29IPC" = "IPC",
         "IPC_34IPC" = "IPC",
         "Neuron_4" = "Neuron 4",
         "Neuron_7" = "Neuron 4",
         "Neuron_32" = "Neuron 4",
         "Neuron_33" = "Neuron 4",
         "Neuron_34" = "Neuron 4",
         "Neuron_33Neuron" = "Neuron 4",
         "Neuron_39Neuron" = "Neuron 4",
         "Neuron_2Neuron" = "Neuron 4",
         "Neuron_66Neuron" = "Neuron 4",
         "GW19_2_29NeuronNeuron" = "Neuron 4",
         "GW18_2_42NeuronNeuron" = "Neuron 4",
         "Neuron_33Neuron" = "Neuron 4",
         "GW19_2_1OutlierrOutlierr" = "Outlier (in Neuron 4)",
         "GW19_2_45OutlierrOutlierr" = "Outlier (in Neuron 4)",
         "Neuron_1" = "Neuron 5",
         "Neuron_2" = "Neuron 5",
         "Neuron_3" = "Neuron 5",
         "Neuron_12" = "Neuron 5",
         "Neuron_10Neuron" = "Neuron 5",
         "Neuron_37Neuron" = "Neuron 5",
         "Neuron_11Neuron" = "Neuron 6",
         "Neuron_36Neuron" = "Neuron 6",
         "Neuron_64Neuron" = "Neuron 6",
         "Neuron_10" = "Neuron 7",
         "Neuron_11" = "Neuron 7",
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



# get genes
genes = read.table("genes.csv", sep = "\t", col.names = 'gene')
genes = genes$gene
# get overlapping genes between test and train dataset (to reduce noise and runtime)
genes.overlap <- rownames(data.list[[1]]@assays$SCT)[rownames(data.list[[1]]@assays$SCT) %in% genes]
which.genes.overlap <- which(rownames(data.list[[1]]@assays$SCT) %in% genes)

# subset data and convert to sce
data.list <- lapply(X = data.list, FUN = function(x) {
  # subset test data by overlapping genes
  x <- subset(x, features = genes.overlap)
  # set SO to SCE object
  x <- Seurat::as.SingleCellExperiment(x, assay = 'SCT')
  return(x)
})

kriegstein.data.list <- list.files(path = RData.folder)
for (i in 1:length(kriegstein.data.list)) {
  start_iter_time <- Sys.time()
  print("iter start")

  print("load cell data")

  # load cell data
  load(paste0(RData.folder, "iter.", i, ".RData")) # cell_data (R object name)
  cell_data <- cell_data # TODO check if works: explicit mention because of scoping error

  # set genes as row names
  cell_data <- data.frame(cell_data, row.names=genes)
  # subset only overlapping genes from cell data
  cell_data <- cell_data[which.genes.overlap,]

  print("create seurat -> sce")
  # get SO from matrix
  cell_data <- SeuratObject::CreateSeuratObject(counts = cell_data, meta.data = meta[names(cell_data),]) # so
  # set SO to SCE object
  cell_data <- Seurat::as.SingleCellExperiment(cell_data) # so.sce

  for (j in 1:length(rds.files)) {
    # setup filename for temporary saving data
    sample_name <- sample_names[j]
    filename <- paste0("RData/SingleR_RData +custom.reference/", sample_name, "_iter.", i)

    ## perform SingleR label transfer
    print("celltype")
    # celltype
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$celltype,
                       clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".celltype.RData"))
    print("age")
    # age
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$age,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".age.RData"))
    print("individual")
    # individual
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$individual,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".individual.RData"))
    print("structure")
    # structure
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$structure,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".structure.RData"))
    print("area")
    # area
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$area,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".area.RData"))
    print("cell.type.v2")
    # cell.type.v2
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$cell.type.v2,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".cell.type.v2.RData"))

    # clusterv2
    print("custom.clusterv2")
    result <- SingleR::SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$custom.clusterv2,
                      clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    save(result, file = paste0(filename, ".custom.clusterv2.RData"))



    # print("cellabbreviation")
    # # cellabbreviation
    # result <- SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$cellabbreviation,
    #                   clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    # save(result, file = paste0(filename, ".cellabbreviation.RData"))

    # clusterv1
    # print("clusterv1")
    # result <- SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$clusterv1,
    #                   clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    # save(result, file = paste0(filename, ".clusterv1.RData"))
    # clusterv2
    # print("clusterv2")
    # result <- SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$clusterv2,
    #                   clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    # save(result, file = paste0(filename, ".clusterv2.RData"))

    # print("celltype.structure")
    # # celltype.structure
    # result <- SingleR(test=data.list[[j]], ref=cell_data, labels=cell_data$celltype.structure,
    #                   clusters=data.list[[j]]$seurat_clusters, de.method='wilcox')
    # save(result, file = paste0(filename, ".celltype.structure.RData"))
  }

  for (thing in ls()) { message(thing); print(object.size(get(thing)), units='auto') }

  print(paste("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes - memory in use:", as.data.frame(gc())$'(Mb)'[2]))
  cat('point 1 mem', memory.size(), memory.size(max=TRUE), '\n')
}

# unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
