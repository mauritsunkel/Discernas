
library(Seurat)
library(dplyr)

renameOrigIdent <- function(data, renameFrom, renameTo) {
  data$orig.ident <- plyr::mapvalues(data$orig.ident, from = renameFrom, to = renameTo)
  return(data)
}

setwd("C:/Users/mauri/Desktop/UCSC Cell Browser scRNAseq  data share")

# "BL_A + BL_C.rds"
# "BL_A.rds"
# "BL_C.rds"
# "BL_N + BL_C.rds"
# "BL_N.rds"
data <- readRDS("BL_N.rds")
data <- renameOrigIdent(data, renameFrom = "BL_N", renameTo = "neurons mono-culture")
data <- renameOrigIdent(data, renameFrom = "BL_C", renameTo = "co-culture")
data <- renameOrigIdent(data, renameFrom = "BL_A", renameTo = "astrocytes mono-culture")


# remove unnecessary metadata
data$SCT_snn_res.0.5 <- NULL
data$seurat_clusters <- NULL
data$clusters <- data$kriegstein.seurat.custom.clusters.mean
data$kriegstein.seurat.custom.clusters.mean <- NULL
data$percent_mitochondrial <- data$percent.mt
data$percent.mt <- NULL
data$original_identity <- data$orig.ident
data$orig.ident <- NULL


# integrated astrocytes
# astrocytes
# co-culture
# integrated neurons
# neurons
name <- "neurons"
data@project.name <- name
saveRDS(data, paste0(name, ".rds"))
