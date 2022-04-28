library(SingleR)
library(pheatmap)
library(Seurat)

### USER PARAMETERS ###
# set work dir
work_dir <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Kriegstein/'
setwd(work_dir)
# get meta
meta <- read.table("custom.meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

RData_folder <- 'C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/data/Kriegstein/RData/SingleR_RData +custom.reference/'
kriegstein.singler.data.list <- list.files(path = RData_folder)

samples <- c("BL_A", "BL_C", "BL_N", "BL_A + BL_C", "BL_N + BL_C")
annotations <- c("celltype", "age", "individual", "structure", "area", "cell.type.v2", "custom.clusterv2")

# load our data (as test data)
rds.files <- c(
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_A/BL_A.rds",
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_C/BL_C.rds",
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/BL_N/BL_N.rds",
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/integrated/BL_A + BL_C/integrated_BL_A + BL_C.rds",
  "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/results/SCTransform + Leiden - Cellcycle + Annotation/integrated/BL_N + BL_C/integrated_BL_N + BL_C.rds"
)

start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/"
dir.create(paste0(work_dir, 'results/'))
dir.create(paste0(work_dir, 'results/', 'Kriegstein/'))
dir.create(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'))
setwd(paste0(work_dir, 'results/', 'Kriegstein/', start_time, '/'))
### END USER PARAMETERS ###



data.list <- lapply(X = rds.files, FUN = function(x) {
  return(readRDS(file = x))
})
names(data.list) <- samples

results.list <- list()
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
      result <- result # TODO check if works: explicit mention because of scoping error
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
    filename <- paste0("cluster_heatmap_scores=first_labels=tuned_", sample, "_", anno, ".png")

    # create scores heatmap to see what label scores highest for a given cell, ambiguity is possible
    SingleR::plotScoreHeatmap(results = combined.results[[paste(sample, anno)]], scores.use=0, calls.use=0,
                     normalize=T, filename = filename)
  }
}

# custom visualizations per sample with all annotations
# source my custom color palettes from utils
source(file="C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/Seurat/R/my_utils/color_palettes.R")
my.color.palettes <- my.color.palettes # TODO check if works: explicit mention because of scoping error
# get custom colors
custom_colors <- my.color.palettes(type = 'mixed')

for (sample in samples) {
  for (anno in annotations) {
    filename <- paste0("heatmap_clusters_scores=first_labels=tuned_", sample, "_", anno ,".png")

    # set annotation column for transferred labels from reference data
    annotation_col <- data.frame(ref.celltype = combined.results[[paste(sample, "celltype")]]$labels,
                                 ref.age = combined.results[[paste(sample, "age")]]$labels,
                                 ref.individual = combined.results[[paste(sample, "individual")]]$labels,
                                 ref.structure = combined.results[[paste(sample, "structure")]]$labels,
                                 ref.area = combined.results[[paste(sample, "area")]]$labels,
                                 ref.cell.type.v2 = combined.results[[paste(sample, "cell.type.v2")]]$labels,
                                 ref.custom.clusterv2 = combined.results[[paste(sample, "custom.clusterv2")]]$labels)
    # get ordered and unique label names from reference data
    colors.celltype <- unique(meta$celltype)[order(unique(meta$celltype))]
    colors.age <- unique(meta$age)[order(unique(meta$age))]
    colors.individual <- unique(meta$individual)[order(unique(meta$individual))]
    colors.structure <- unique(meta$structure)[order(unique(meta$structure))]
    colors.area <- unique(meta$area)[order(unique(meta$area))]
    colors.cell.type.v2 <- unique(meta$cell.type.v2)[order(unique(meta$cell.type.v2))]
    colors.custom.clusterv2 <- unique(meta$custom.clusterv2)[order(unique(meta$custom.clusterv2))]
    # set colors to each unique label from reference data
    names(colors.celltype) <- my.color.palettes(type = 'mixed', n = length(colors.celltype))
    names(colors.age) <- my.color.palettes(type = 'mixed', n = length(colors.age))
    names(colors.individual) <- my.color.palettes(type = 'mixed', n = length(colors.individual))
    names(colors.structure) <- my.color.palettes(type = 'mixed', n = length(colors.structure))
    names(colors.area) <- my.color.palettes(type = 'mixed', n = length(colors.area))
    names(colors.cell.type.v2) <- my.color.palettes(type = 'mixed', n = length(colors.cell.type.v2))
    names(colors.custom.clusterv2) <- my.color.palettes(type = 'mixed', n = length(colors.custom.clusterv2))
    # swap names and values of named vector to proper formatting
    colors.celltype <- setNames(names(colors.celltype), colors.celltype)
    colors.age <- setNames(names(colors.age), colors.age)
    colors.individual <- setNames(names(colors.individual), colors.individual)
    colors.structure <- setNames(names(colors.structure), colors.structure)
    colors.area <- setNames(names(colors.area), colors.area)
    colors.cell.type.v2 <- setNames(names(colors.cell.type.v2), colors.cell.type.v2)
    colors.custom.clusterv2 <- setNames(names(colors.custom.clusterv2), colors.custom.clusterv2)
    # set annotation colors for all transferred labels
    annotation_colors <- list(ref.celltype = colors.celltype[unique(annotation_col$ref.celltype)],
                              ref.age = colors.age[unique(annotation_col$ref.age)],
                              ref.individual = colors.individual[unique(annotation_col$ref.individual)],
                              ref.structure = colors.structure[unique(annotation_col$ref.structure)],
                              ref.area = colors.area[unique(annotation_col$ref.area)],
                              ref.cell.type.v2 = colors.cell.type.v2[unique(annotation_col$ref.cell.type.v2)],
                              ref.custom.clusterv2 = colors.custom.clusterv2[unique(annotation_col$ref.custom.clusterv2)])
    # order annotation columns for visualization
    annotation_col <- annotation_col[, c("ref.custom.clusterv2", "ref.cell.type.v2", "ref.celltype",
                                         "ref.structure", "ref.area", "ref.individual", "ref.age")]
    # set rownames for identification of rows during plotting
    rownames(combined.results[[paste(sample, anno)]]$scores) <- levels((data.list[[sample]]$seurat_clusters))
    # plot pretty heatmap
    p <- pheatmap::pheatmap(t(combined.results[[paste(sample, anno)]]$scores),
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
