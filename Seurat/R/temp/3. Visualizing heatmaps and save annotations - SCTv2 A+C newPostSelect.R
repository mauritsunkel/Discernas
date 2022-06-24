library(SingleR)
library(pheatmap)
library(Seurat)



### USER PARAMETERS ###
## TODO testing
# set RData folder
RData_folder <-"C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein/RData/Pipe_SCTv2_23-06/A+C newPostSelect"
# set sample(s)
samples <- c("BL_A + BL_C")
# set file(s)
rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06/integrated - new selection/BL_A + BL_C/after_selection/BL_A + BL_C.rds")



## pre-selection data
# # set sample names (NOTE: same order as rds.files)
# samples <- c("BL_A", "BL_C", "BL_N", "BL_A + BL_C", "BL_N + BL_C")
# # set rds files ((NOTE: same order as samples))
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/BL_A/BL_A.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/BL_C/BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/BL_N/BL_N.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/integrated/BL_A + BL_C/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/integrated/BL_N + BL_C/BL_N + BL_C.rds"
#   )
## post-selection data
# samples <- c("BL_A + BL_C", "BL_N + BL_C")
# rds.files <- c(
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/integrated/BL_A + BL_C/after_selection/BL_A + BL_C.rds",
#   "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/integrated/BL_N + BL_C/after_selection/BL_N + BL_C.rds"
#   )
names(rds.files) <- samples







# source my custom color palettes from utils
source(file="C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/R/my_utils/color_palettes.R")
my.color.palettes <- my.color.palettes
# get custom colors
custom_colors <- my.color.palettes(type = 'mixed')

# set work dir
work_dir <- 'C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Kriegstein/'
setwd(work_dir)

# initialize function and get meta data
getMeta <- function() {
  metaFile <- "custom.meta.tsv"
  if (!file.exists(metaFile)) {
    stop(metaFile, " does not exist, this file is generated during data chunking.")
  } else {
    return(read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
  }
}
meta <- getMeta()


# set annotations
annotations <- c("age", "structure", "custom.clusterv2")
annotations.to.plot <- c("custom.clusterv2")

# set start time and create output directories
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"
dir.create(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'), recursive = TRUE)
setwd(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'))
### END USER PARAMETERS ###

# read and return rds data
data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})
# add sample names to rds data
names(data.list) <- samples



# initialize list to store results
results.list <- list()
# iterate samples and annotations: grab all corresponding files, then load and return corresponding data
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
  # use SingleR::combineCommonResults (instead of combineRecomputedResults) because
  ## from combineRecomputedResults docs: It is strongly recommended that the
  ## universe of genes be the same across all references
  ### because if the intersection of all genes over all references is highly
  ### different, the availability between refs may have unpredictable results
  combined <- SingleR::combineCommonResults(x)

  # create dataframe to select values of the highest scoring reference
  df <- as.data.frame(combined$scores)
  ## get max score and label for each highest scoring reference
  # for each unique column name select all columns
  combined$max.scores <- sapply(unique(names(df)), function(cols) {
    # for each row (cluster)
    sapply(1:nrow(df), function(row) {
      # get the column by highest row-wise value within the subset of selected columns
      col <- max.col(df[names(df) == cols])[row]
      # select the highest score (value) based on the defined column and row within the subset of selected columns
      df[names(df) == cols][row, col]
    })
  })
  combined$max.labels <- colnames(combined$max.scores)[max.col(combined$max.scores)]
  ## get mean score and label of all references
  # for each unique column name select all columns
  combined$mean.scores <- sapply(unique(names(df)), function(cols) {
    rowMeans(df[names(df) == cols]) # mean per cluster (row) per set of columns (label)
  })
  combined$mean.labels <- colnames(combined$mean.scores)[max.col(combined$mean.scores)]

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
  # copy seurat clusters metadata to aggregate with Kriegstein labels
  data.list[[sample]]$kriegstein.seurat.custom.clusters.mean <- data.list[[sample]]$seurat_clusters
  # overwrite copied metadata to an aggregation of Kriegstein.Seurat custom cluster labels
  levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean) <- paste0(combined.results[[paste(sample, "custom.clusterv2")]]$mean.labels, ".", levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean))

  # overwrite rds file with new misc(elleneous) annotation (note: NOT metadata, as that is about cells in SO)
  saveRDS(data.list[[sample]], file = rds.files[[sample]])
}



# TODO put this (visualizing) in its own script, just need rds data + meta now as SingleR results are stored in misc
# custom visualizations per sample-reference comparison for each annotation
for (sample in samples) {
  annotation_col <- data.frame(row.names = levels(data.list[[sample]]$seurat_clusters))
  annotation_colors <- list()

  for (anno in annotations) {
    # set annotation column for transferred labels from reference data
    annotation_col[ , ncol(annotation_col) + 1] <- data.list[[sample]]@misc[[paste0("Kriegstein.SingleR.", anno)]][["mean.labels"]]
    colnames(annotation_col)[ncol(annotation_col)] <- paste0("ref.", anno)

    # get ordered and unique label names from reference data
    labels <- unique(meta[which(colnames(meta) == anno)][,1])
    annotation_colors[[length(annotation_colors) + 1]] <- labels[order(labels)]
    names(annotation_colors)[length(annotation_colors)] <- paste0("ref.", anno)

    # set colors to each unique label from reference data
    names(annotation_colors[[paste0("ref.", anno)]]) <- my.color.palettes(type = 'mixed', n = length(annotation_colors[[paste0("ref.", anno)]]))

    # swap names and values of named vector to proper formatting
    annotation_colors[[paste0("ref.", anno)]] <- setNames(names(annotation_colors[[paste0("ref.", anno)]]), annotation_colors[[paste0("ref.", anno)]])

    # set annotation colors for all transferred labels
    annotation_colors[[paste0("ref.", anno)]] <- annotation_colors[[paste0("ref.", anno)]][unique(annotation_col[which(colnames(annotation_col) == paste0("ref.", anno))][,1])]
  }

  # order annotation columns for visualization
  names(annotation_col) <- c("fetal.brain.age", "fetal.brain.structure", "fetal.brain.celltype")
  annotation_col <- annotation_col[, c("fetal.brain.celltype", "fetal.brain.structure", "fetal.brain.age")]

  for (anno in annotations.to.plot) {
    filename <- paste0("heatmap_clusters_scores=first_labels=tuned_", sample, "_", anno ,".png")
    print(paste("plotting:", filename))

    # set rownames for identification of rows during plotting
    rownames(combined.results[[paste(sample, anno)]][["mean.scores"]]) <- levels((data.list[[sample]]$seurat_clusters))

    # plot pretty heatmap
    p <- pheatmap::pheatmap(t(combined.results[[paste(sample, anno)]][["mean.scores"]]),
                  fontsize = 9,
                  color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "PRGn"))(100),
                  labels_col = paste0(levels(data.list[[sample]]$seurat_clusters), " (n=", table(data.list[[sample]]$seurat_clusters), ")"),
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  cluster_cols = T,
                  main = "Scores",
                  filename = filename
    )
  }
}
