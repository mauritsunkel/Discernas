library(SingleR)
library(pheatmap)
library(Seurat)
library(stringr)
library(grDevices)



### USER PARAMETERS ###
# set RData folder
RData_folder <-"C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein/RData/t90"
# set sample(s)
samples <- c("t90")
# set file(s)
rds.files <- c(
  "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Mark organoids first run/t90/t90.rds"
  )
names(rds.files) <- samples

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein/'
setwd(work_dir)

# "max" as SingleR intended (combineCommonResults), max.scores/max.labels across references
# "mean" custom for averaging scores (and labels) across references
RefAggrStrategy <- "max"

# source and get custom color palette from utils
source(file="C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/R/utils/color_palettes.R")
my.color.palettes <- my.color.palettes
custom_colors <- my.color.palettes(type = 'mixed')

# initialize function and get meta data
getMeta <- function() {
  metaFile <- "custom.meta.tsv"
  if (!file.exists(metaFile)) {
    stop(metaFile, " does not exist, this file is generated during data chunking.")
  } else {
    return(utils::read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
  }
}
meta <- getMeta()

# set annotations (used as meta plot data and/or used for plotting its values)
annotations <- c("age", "structure", "custom.clusterv2")
annotations.to.plot <- c("custom.clusterv2")

# set start time and create output directories
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
dir.create(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'), recursive = TRUE)
setwd(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'))
### END USER PARAMETERS ###



# read and get .rds data
data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})
names(data.list) <- samples

# initialize list to store SingleR results
results.list <- list()
# iterate samples and annotations: grab all corresponding files, then load and get corresponding data
for (sample in samples) {
  for (anno in annotations) {
    # list all files in RData_folder based on sample and annotation
    files <- list.files(path = RData_folder, pattern = paste0(stringr::str_replace(sample, " \\+ ", " .* "), ".*", anno), full.names=T)

    # if '+' not in sample (name), then
    if (!grepl("+", sample, fixed = T)) {
      # if '+' in stems of files, then remove those files from list (separate individual from integrated files)
      files <- files[!grepl("+", sapply(stringr::str_split(files, "/"), tail, 1), fixed = T)]
    }

    # load each file, return its actual object
    results <- sapply(files, function(file) {
      load(file) # result (name of R object)
      return(result)
    })

    # put loaded objects into results.list
    results.list[[paste(sample, anno)]] <- results
  }
}

# CombineCommonResults of corresponding data for each single sample-reference comparison
combined.results <- lapply(X = results.list, FUN = function(x) {
  # get scores from combined results
  combined <- SingleR::combineCommonResults(x)
  df <- as.data.frame(combined$scores)

  ## get max score and label for each highest scoring reference
  # for each unique column name, select all its corresponding columns (scores, result chunks of SingleR)
  combined$max.scores <- sapply(unique(names(df)), function(names) {
    # for each row (cluster)
    sapply(1:nrow(df), function(row) {
      # get the column by highest row-wise value within the subset of selected columns
      col <- max.col(df[names(df) == names])[row]
      # select the highest score (value)
      df[names(df) == names][row, col]
    })
  })
  combined$max.labels <- colnames(combined$max.scores)[max.col(combined$max.scores)]

  # write scores for figure reference
  utils::write.csv2(combined$max.scores, file = paste0("Kriegstein_Pearson.correlation.max_", sample, "_", anno ,".csv"))

  ## get mean score and label of all references
  # for each unique column name select all columns
  combined$mean.scores <- sapply(unique(names(df)), function(names) {
    rowMeans(df[names(df) == names]) # mean per cluster (row) per set of named columns (scores)
  })
  combined$mean.labels <- colnames(combined$mean.scores)[max.col(combined$mean.scores)]

  # write scores for figure reference
  write.csv2(combined$mean.scores, file = paste0("Kriegstein_Pearson.correlation.mean_", sample, "_", anno ,".csv"))

  return(combined)
})


# save Kriegstein cluster labels into Seurat object --> rds
for (sample in samples) {
  for (anno in annotations) {
    # Kriegstein data and labels from sample-annotation correlation overlap
    misc_data <- combined.results[[paste(sample, anno)]]
    # add as miscellaneous data to Seurat object
    SeuratObject::Misc(object = data.list[[sample]], slot = paste0("Kriegstein.SingleR.", anno)) <- misc_data
  }

  SeuratObject::Misc(object = data.list[[sample]], slot = "Kriegstein.SingleR.RefAggrStrategy") <- RefAggrStrategy
  # copy seurat clusters metadata to aggregate with Kriegstein labels
  data.list[[sample]]$kriegstein.seurat.custom.clusters.mean <- data.list[[sample]]$seurat_clusters
  # overwrite copied metadata to an aggregation of Kriegstein.Seurat custom cluster labels
  levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean) <- paste0(combined.results[[paste(sample, "custom.clusterv2")]]$mean.labels, ".", levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean))

  # overwrite rds file with new misc annotation
  saveRDS(data.list[[sample]], file = rds.files[[sample]])
}

# custom visualizations per sample-reference comparison for each annotation
for (sample in samples) {
  annotation_col <- data.frame(row.names = levels(data.list[[sample]]$seurat_clusters))
  annotation_colors <- list()

  for (anno in annotations) {
    # set annotation column for transferred labels from reference data
    annotation_col[ , ncol(annotation_col) + 1] <- data.list[[sample]]@misc[[paste0("Kriegstein.SingleR.", anno)]][[paste0(RefAggrStrategy, ".labels")]]
    colnames(annotation_col)[ncol(annotation_col)] <- paste0("ref.", anno)

    # get ordered and unique label names from reference data
    labels <- unique(meta[which(colnames(meta) == anno)][,1])
    annotation_colors[[length(annotation_colors) + 1]] <- labels[order(labels)]
    names(annotation_colors)[length(annotation_colors)] <- paste0("ref.", anno)

    # set colors to each unique label from reference data
    names(annotation_colors[[paste0("ref.", anno)]]) <- my.color.palettes(type = 'mixed', n = length(annotation_colors[[paste0("ref.", anno)]]))

    # swap names and values of named vector to proper formatting
    annotation_colors[[paste0("ref.", anno)]] <- stats::setNames(names(annotation_colors[[paste0("ref.", anno)]]), annotation_colors[[paste0("ref.", anno)]])

    # set annotation colors for all transferred labels
    annotation_colors[[paste0("ref.", anno)]] <- annotation_colors[[paste0("ref.", anno)]][unique(annotation_col[which(colnames(annotation_col) == paste0("ref.", anno))][,1])]
  }

  # order annotation columns for visualization
  names(annotation_col) <- c("fetal.brain.age", "fetal.brain.structure", "fetal.brain.celltype")
  annotation_col <- annotation_col[, c("fetal.brain.celltype", "fetal.brain.structure", "fetal.brain.age")]

  for (anno in annotations.to.plot) {
    filename <- paste0("Kriegstein_heatmap_", sample, "_", anno ,".png")
    print(paste("plotting:", filename))

    # set rownames for identification of rows during plotting
    rownames(combined.results[[paste(sample, anno)]][[paste0(RefAggrStrategy, ".scores")]]) <- levels((data.list[[sample]]$seurat_clusters))

    # plot pretty heatmap
    p <- pheatmap::pheatmap(t(combined.results[[paste(sample, anno)]][[paste0(RefAggrStrategy, ".scores")]]),
                  fontsize = 9,
                  color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PiYG"))(100),
                  labels_col = paste0(levels(data.list[[sample]]$seurat_clusters), " (n=", table(data.list[[sample]]$seurat_clusters), ")"),
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  cluster_cols = T,
                  main = "Scores",
                  filename = filename
    )
  }
}
