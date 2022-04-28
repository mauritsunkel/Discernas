
# create Kriegstein custom EMC poster heatmap

library(SingleR)
library(pheatmap)
library(Seurat)

### USER PARAMETERS ###
# set work dir
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Kriegstein/'
setwd(work_dir)
# get meta
meta <- read.table("custom.meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

folder <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Kriegstein/singler_6+custom.clusterv2/'
kriegstein.singler.data.list <- list.files(path = folder)

samples <- c("BL_C")
annotations <- c("custom.clusterv2")

# load our data (as test data)
rds.files <- c(
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds"
)
### END USER PARAMETERS ###


data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})
names(data.list) <- samples

results.list <- list()
for (sample in samples) {
  for (anno in annotations) {
    # list all files in folder based on sample and annotation
    files <- list.files(path = folder, pattern = paste0(stringr::str_replace(sample, " \\+ ", " .* "), ".*", anno), full.names=T)

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

combined.results <- lapply(X = results.list, FUN = function(x) {
  # use SingleR::combineCommonResults (instead of combineRecomputedResults) because
  ## from combineRecomputedResults docs: It is strongly recommended that the
  ## universe of genes be the same across all references
  ### because if the intersection of all genes over all references is highly
  ### different, the availability between refs may have unpredictable results
  combined <- SingleR::combineCommonResults(x)

  # set scores (cbind of different refs (iterations) to combined.scores)
  combined$combined.scores <- combined$scores
  # take row mean for each ref and label as visualization score
  df <- as.data.frame(combined$scores)
  combined$scores <- as.matrix( # sapply returns a list here, so we convert it to a data.frame
    sapply(unique(names(df)), # for each unique column name
           function(col) rowMeans(df[names(df) == col]) # calculate row means
    )
  )
  return(combined)
})

# simple visualizations per type
for (sample in samples) {
  for (anno in annotations) {
    filename <- paste0("results/EMC_poster/cell_heatmap_scores=first_labels=tuned_", sample, "_", anno, ".png")

    # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
    plotScoreHeatmap(results = combined.results[[paste(sample, anno)]], scores.use=0, calls.use=0,
                     normalize=T, filename = filename)
  }
}






# custom visualizations per sample with all annotations
# source my custom color palettes from utils
source(file="C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/R/my_utils/color_palettes.R")
# get custom colors
custom_colors <- my.color.palettes(type = 'mixed')





for (sample in samples) {
  for (anno in annotations) {
    filename <- paste0("results/EMC_poster/heatmap_clusters_scores=first_labels=tuned_", sample, "_", anno ,".png")

    # set annotation column for transferred labels from reference data
    annotation_col <- data.frame(ref.custom.clusterv2 = combined.results[[paste(sample, "custom.clusterv2")]]$first.labels)
    # get ordered and unique label names from reference data
    colors.custom.clusterv2 <- unique(meta$custom.clusterv2)[order(unique(meta$custom.clusterv2))]
    # set colors to each unique label from reference data
    names(colors.custom.clusterv2) <- my.color.palettes(type = 'mixed', n = length(colors.custom.clusterv2))
    # swap names and values of named vector to proper formatting
    colors.custom.clusterv2 <- setNames(names(colors.custom.clusterv2), colors.custom.clusterv2)
    # set annotation colors for all transferred labels
    annotation_colors <- list(ref.custom.clusterv2 = colors.custom.clusterv2[unique(annotation_col$ref.custom.clusterv2)])
    # order annotation columns for visualization
    # annotation_col <- annotation_col[, c("ref.custom.clusterv2")]
    # set rownames for identification of rows during plotting
    rownames(combined.results[[paste(sample, anno)]]$scores) <- levels((data.list[[sample]]$seurat_clusters))
    # plot pretty heatmap
    pheatmap(t(combined.results[[paste(sample, anno)]]$scores),
                  fontsize = 9,
                  color = colorRampPalette(RColorBrewer::brewer.pal(n=7,"Blues"))(100),
                  labels_col = paste0(levels(data.list[[sample]]$seurat_clusters), " (n=", table(data.list[[sample]]$seurat_clusters), ")"),
                  annotation_col = annotation_col,
                  annotation_colors = annotation_colors,
                  cluster_cols = T,
                  # main = "Scores",
                  # filename = filename,
                  border_color = "#283271",
                  number_col = "yellow",
                  angle_col = 315
    )

  }
}


# number_col (expect black, want EMC blue)
# angle_col (expect 90, try skewed 45 (presets: 0, 45, 90, 270))



# annotation_colors (expect black, try EMC blue)


# combined.results

# BL_C
data <- readRDS("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds")



labels <- combined.results$`BL_C custom.clusterv2`$first.labels
old_labels <- sapply(strsplit(levels(data), "\\."), function(x) x[2])
labels <- paste(labels, old_labels, sep = ".")
names(labels) <- levels(data)
data <- RenameIdents(data, labels)

cols <- c("#2071b5", "#00c4ff", "#c4edfc", "#85d0f5", "#c4edfc", "#85d0f5",
          "#85d0f5", "#85d0f5", "#85d0f5", "#00c4ff", "#2071b5", "#85d0f5",
          "#85d0f5")
DimPlot(data, reduction = "umap", label = TRUE, cols = cols, label.size = 5, label.color = "#283271")



# neuron 5 #2071b5
# neuron 4 #85d0f5
# dividing #00c4ff
# radial glia 2 #c4edfc


# feature plot EMC poster
features <- c("VIM", "S100B", "SOX9", "FABP7", "SLC1A2",
              "MAP2", "NEUROG2", "SYN1", "RBFOX3", "SLC1A3")
FeaturePlot(data, features = features, cols = c("#85d0f5", "#2b2f70"), ncol = 5, label.color = "#283271")
