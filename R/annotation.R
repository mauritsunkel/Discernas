#' Run Kriegstein reference data chunking for annotation.
#'
#' Need to chunk data before training and annoation as to reduce memory needed
#' per interation. If no chunking needed set n_chunkcs to 1.
#'
#' @param n_chunks default: 25, integer, amount of data chunks, less memory needed per chunk
#' @param kriegstein_data_dir String with Kriegstein folder path, containing meta.tsv and exprMatrix.tsv.gz.
#' @param kriegstein_chunks_output_dir String with Kriegstein output directory for data chunks.
#' @param kriegstein_custom_annotation .csv with Kriegstein clusterv2 reference data
#'
#' @export
#'
#' @examplesIf FALSE
#' kriegstein_custom_annotation <- system.file(
#' "extdata",
#' "kriegstein_custom_annotation.txt", package = 'EMC.SKlab.scRNAseq')
#' chunk_kriegstein_data(
#'   n_chunks = 25,
#'   kriegstein_data_dir = "path/to/Kriegstein_data/",
#'   kriegstein_chunks_output_dir = paste0(kriegstein_data_dir, "RData/chunks/"),
#'   kriegstein_custom_annotation = kriegstein_custom_annotation
#' )
#'
#' @note
#' Added Kriegstein mapping based on visual inspection of their UCSC
#' CellBrowser UMAP, see custom_annotation.txt.
#'
#' Can't train first and then intersect to common genes space to classify after.
#'
#' Run EMC.SKlab.scRNAseq::chunk_kriegstein_data() to have chunked reference data before annotating!
#' No need to run if chunked data already exists!
chunk_kriegstein_data <- function(n_chunks, kriegstein_data_dir, kriegstein_chunks_output_dir, kriegstein_custom_annotation) {
  # get start time to measure run time
  start_run_time <- Sys.time()

  # set working directory
  if (!file.exists(file.path(kriegstein_data_dir, "meta.tsv")) || !file.exists(file.path(kriegstein_data_dir, "exprMatrix.tsv.gz"))) {
    stop("Kriegstein data folder should contain meta.tsv and exprMatrix.tsv.gz")
  }
  dir.create(kriegstein_chunks_output_dir, recursive = T)

  # get custom meta data from file, create from original if doesn't exist yet
  if (file.exists(file.path(kriegstein_data_dir, "custom.meta.tsv"))) {
    meta <- utils::read.table(file.path(kriegstein_data_dir, "custom.meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)
  } else {
    meta <- utils::read.table(file.path(kriegstein_data_dir, "meta.tsv"), header=T, sep="\t", as.is=T, row.names=1)
    # create custom merged cluster annotations for clusterv2 from reference data
    anno_df <- read.csv2(kriegstein_custom_annotation)
    meta$custom.clusterv2 <- plyr::mapvalues(meta$clusterv2, anno_df$from, anno_df$to)

    # save custom annotation for use in processing and visualization
    utils::write.table(meta, file = file.path(kriegstein_data_dir, "custom.meta.tsv"), sep="\t", row.names=T)
  }

  # get genes, if non-existent create from reference data
  if (file.exists(file.path(kriegstein_data_dir, "kriegstein_genes.csv")) & file.exists(file.path(kriegstein_data_dir, "n_cols.txt"))) {
    genes = utils::read.table(file.path(kriegstein_data_dir, "kriegstein_genes.csv"), sep = "\t", col.names = 'gene')
    genes = genes$gene[2:length(genes$gene)]
    n_cols <- read.csv2(file.path(kriegstein_data_dir, "n_cols.txt"), header = F)
    n_cols <- n_cols$V1
  } else {
    # read only columns for efficiency
    initial <- data.table::fread(file.path(kriegstein_data_dir, "exprMatrix.tsv.gz"), nrows=0, colClasses=c("gene"="character"))
    n_cols <- dim(initial)[2]
    # read first row only for efficiency
    initial <- data.table::fread(file.path(kriegstein_data_dir, "exprMatrix.tsv.gz"), select=c(1), colClasses=c("gene"="character"))
    genes = initial[,1][[1]]
    genes = gsub(".+[|]", "", genes)
    utils::write.csv2(genes, file = file.path(kriegstein_data_dir, "kriegstein_genes.csv"), row.names = FALSE)
    write(n_cols, file.path(kriegstein_data_dir, "n_cols.txt"))
  }

  # set constants for data chunking
  n_cells <- n_cols - 1
  cell_ids <- 2:n_cols
  n_rows <- length(genes)
  n_sample <- ceiling(n_cells/n_chunks)
  # set seed for randomness in SingleR aggregation
  set.seed(0)

  # print initial run time and memory usage
  message('pre-chunk run time: ', Sys.time()-start_run_time)
  message('point 1 mem ', utils::memory.size(), ' - ', utils::memory.size(max=TRUE))
  # get run time before iterating
  start_run_time <- Sys.time()
  # initialize sampled_cells (to store sampled cell_ids)
  sampled_cells <- vector()

  # iterate to chunk reference data
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

    cell_data <- data.table::fread(file.path(kriegstein_data_dir, "exprMatrix.tsv.gz"), select = c(unique(sample_cells)), colClasses=c("gene"="character"))
    rownames(cell_data) <- genes
    # convert table to SO and add metadata matching by cell name
    cell_data <- SeuratObject::CreateSeuratObject(counts = cell_data, meta.data = meta[names(cell_data),])
    cell_data <- Seurat::as.SingleCellExperiment(cell_data)
    # SingleR() expects REFERENCE datasets to be normalized and log-transformed
    cell_data <- scuttle::logNormCounts(cell_data)

    save(cell_data, file = file.path(kriegstein_chunks_output_dir, paste0("chunk.", i, ".RData")))
    # remove cell data before next iteration to save memory
    rm(cell_data)

    message("Iteration ", i ," runtime was ", Sys.time()-start_iter_time, " minutes")
  }
  message("Total run time: ", Sys.time() - start_run_time)
}

#' Perform annotation with chunked Kriegstein reference data
#'
#' @param sample_names Character vector with sample names.
#' @param sample_files Character vector with sample .qs files to be annotated.
#' @param output_dir String with output directory for results.
#' @param kriegstein_data_dir String with Kriegstein folder path, containing (custom.)meta.tsv and exprMatrix.tsv.gz.
#' @param kriegstein_chunks_input_dir String with Kriegstein input directory for data chunks.
#' @param kriegstein_annotated_output_dir String with Kriegstein output directory for annotated data chunks.
#' @param annotations default: c("age", "structure", "custom.clusterv2"), see (custom.)meta.tsv for selectable features.
#' @param annotations_to_plot default: c("custom.clusterv2"), annotations from meta features used for individual heatmaps.
#' @param ref_aggr_strategy default: "max", choose one of "max" or "mean"
#' @param run_only_visualization default: FALSE, if TRUE, skip/assume SingleR already ran, need to adjust visualisations only
#'
#' @export
#'
#' @examplesIf FALSE
#' kriegstein_data_dir <-  "path/to/Kriegstein_data/"
#'
#' annotate_visualize_with_kriegstein_data(
#'   sample_names = c("sample_a", "sample_b"),
#'   sample_files = c("path/to/sample_a.qs", "path/to/sample_b.qs"),
#'   kriegstein_data_dir = kriegstein_data_dir,
#'   kriegstein_chunks_input_dir = paste0(kriegstein_data_dir, "/Kriegstein_chunks/"),
#'   kriegstein_annotated_output_dir = paste0(kriegstein_data_dir, "Kriegstein_annotated_RData/"),
#'   annotations = c("age", "structure", "custom.clusterv2")
#' )
#'
#' @note
#' Run EMC.SKlab.scRNAseq::chunk_kriegstein_data() to have chunked reference data before annotating!
#' No need to run if chunked data already exists!
annotate_visualize_with_kriegstein_data <- function(
    sample_names, sample_files, output_dir,
    kriegstein_data_dir,
    kriegstein_chunks_input_dir,
    kriegstein_annotated_output_dir,
    annotations = c("age", "structure", "custom.clusterv2"),
    annotations_to_plot = c("custom.clusterv2"),
    ref_aggr_strategy = "max",
    run_only_visualization = FALSE) {

  dir.create(kriegstein_annotated_output_dir)

  if (!run_only_visualization) {
    genes <- getGenes(kriegstein_data_dir)

    # iterate files and perform SingleR for annotation with Pearson correlation
    for (j in 1:length(sample_files)) {
      # initialize iteration timer
      start_iter_time <- Sys.time()

      sample_name <- sample_names[j]
      message("start iteration of ", sample_name)

      sample_output_dir <- file.path(output_dir, sample_name, 'annotation_kriegstein')
      dir.create(sample_output_dir, recursive = T)

      sample_data <- qs::qread(file = sample_files[j])
      # set RNA assay to have genes in same feature space as the reference data
      SeuratObject::DefaultAssay(sample_data) <- "SCT"

      # get overlapping genes between data and reference
      sample_genes <- rownames(sample_data)
      overlapping_genes <- intersect(sample_genes, genes)
      openxlsx::write.xlsx(x = overlapping_genes, file = file.path(sample_output_dir, "overlapping_genes.xlsx"))

      message("create Seurat sample -> SCE")
      sample_data <- Seurat::as.SingleCellExperiment(sample_data, assay = 'SCT')
      sample_data <- sample_data[overlapping_genes,]
      # perform transformation to have genes in the same feature space as the reference data
      sample_data <- scuttle::logNormCounts(sample_data)

      for (i in 1:length(list.files(path = kriegstein_chunks_input_dir))) {
        message("start iteration of ", sample_name, " at cell RData chunk:  ", kriegstein_chunks_input_dir, "/iter.", i, ".RData")

        load(file.path(kriegstein_chunks_input_dir, paste0("chunk.", i, ".RData"))) # cell_data (R object name)
        cell_data <- cell_data[overlapping_genes,]

        # set filename base for RData saving
        filename_base <- file.path(kriegstein_annotated_output_dir, paste0(sample_name, ".iter.", i))
        # run SingleR for each annotation and save RData
        for (annotation in annotations) {
          message("sample=", sample_name, " iter=", i, " annotation=", annotation)
          result <- SingleR::SingleR(
            test = sample_data,
            ref = cell_data,
            labels = SummarizedExperiment::colData(cell_data)[, annotation],
            clusters = SummarizedExperiment::colData(sample_data)[, "seurat_clusters"],
            de.method = 'wilcox',
            aggr.ref = FALSE)
          save(result, file = paste0(filename_base, ".", annotation, ".RData"))
          # remove result before next iteration to save memory
          rm(result)
        }
      }
      message("Iteration", i ,"runtime was", Sys.time()-start_iter_time, "minutes - memory in use:", as.data.frame(gc())$'(Mb)'[2])
      message('point 1 mem ', utils::memory.size(), ' ', utils::memory.size(max=TRUE))
    }
  }

  # visualize UMAPs and heatmap after annotation
  visualize_kriegstein_annotated_data(
    sample_names = sample_names, sample_files = sample_files, output_dir = output_dir,
    kriegstein_data_dir, kriegstein_annotated_input_dir = kriegstein_annotated_output_dir,
    annotations = annotations,
    annotations_to_plot = annotations_to_plot,
    ref_aggr_strategy = ref_aggr_strategy)
}

#' Perform visualizations of annotated data from Kriegstein reference data.
#'
#' @param sample_names Character vector with sample names.
#' @param sample_files Character vector with sample .qs files.
#' @param output_dir String with output directory for results.
#' @param kriegstein_data_dir String with Kriegstein folder path, containing (custom.)meta.tsv and exprMatrix.tsv.gz.
#' @param kriegstein_annotated_input_dir String with Kriegstein input directory containing annotated data chunks.
#' @param annotations default: c("age", "structure", "custom.clusterv2"), annotations from meta features used for heatmap metadata.
#' @param annotations_to_plot default: c("custom.clusterv2"), annotations from meta features used for individual heatmaps.
#' @param ref_aggr_strategy default: "max", choose one of "max" or "mean".
#' "max" as SingleR intended (combineCommonResults), max.scores/max.labels across references.
#' "mean" custom for averaging scores (and labels) across references.
#'
#' @export
#'
#' @examplesIf FALSE
#' kriegstein_data_dir <-  "path/to/Kriegstein_data/"
#'
#' visualize_kriegstein_annotated_data(
#'   sample_names = c("A", "B"),
#'   sample_files = c(file.path("path", "to", "sampleA.qs"), file.path("path", "to", "sampleB.qs")),
#'   output_dir = file.path("path", "to", "results"),
#'   kriegstein_data_dir = kriegstein_data_dir,
#'   kriegstein_annotated_input_dir = file.path("path", "to", "kriegstein_annotated_Rdata"),
#'   annotations = c("age", "structure", "custom.clusterv2"),
#'   annotations_to_plot = c("custom.clusterv2"),
#'   ref_aggr_strategy = "max"
#' )
visualize_kriegstein_annotated_data <- function(
    sample_names, sample_files, output_dir,
    kriegstein_data_dir, kriegstein_annotated_input_dir,
    annotations = c("age", "structure", "custom.clusterv2"),
    annotations_to_plot = c("custom.clusterv2"),
    ref_aggr_strategy = "max") {

  library(Seurat) # added because of error
  # Error: package or namespace load failed for ‘Seurat’ in .doLoadActions(where, attach):
  #   error in load action .__A__.1 for package RcppAnnoy: loadModule(module = "AnnoyAngular", what = TRUE, env = ns, loadNow = TRUE): Unable to load module "AnnoyAngular": attempt to apply non-function
  # Error in .requirePackage(package) :
  #   unable to find required package ‘Seurat’

  names(sample_files) <- sample_names

  # get Kriegstein custom feature metadata with custom clusterv2 celltype mapping
  meta <- getMeta(kriegstein_data_dir)

  # read and get .qs data
  data.list <- lapply(X = sample_files, FUN = function(x) {
    return(qs::qread(x))
  })
  names(data.list) <- sample_names

  # initialize list to store SingleR results
  results.list <- list()
  # iterate sample_names and annotations: grab all corresponding files, then load and get corresponding data
  for (sample in sample_names) {
    for (anno in annotations) {
      # list all files in kriegstein_annotated_input_dir based on sample and annotation
      files <- list.files(path = kriegstein_annotated_input_dir, pattern = paste0("^",  stringr::str_replace(sample, " \\+ ", " .*"), ".iter.*", anno), full.names=T)

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
  combined.results <- lapply(seq_along(results.list), function(i) {

    # get name from names list
    # name <- n[[i]]
    name <- names(results.list)[i]
    sample <- strsplit(name, " ")[[1]][1]

    # get scores from combined results
    # combined <- SingleR::combineCommonResults(x)
    combined <- SingleR::combineCommonResults(results.list[[i]])
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
    if (grepl('postSelect', sample_files[[sample]])) {
      filename <- file.path(output_dir, sample, 'postSelect', 'annotation_kriegstein', paste0("Pearson.correlation.max_", name ,".csv"))
    } else {
      filename <- file.path(output_dir, sample, 'annotation_kriegstein', paste0("Pearson.correlation.max_", name ,".csv"))
    }
    utils::write.csv2(combined$max.scores, file = filename)

    ## get mean score and label of all references
    # for each unique column name select all columns
    combined$mean.scores <- sapply(unique(names(df)), function(names) {
      rowMeans(df[names(df) == names]) # mean per cluster (row) per set of named columns (scores)
    })
    combined$mean.labels <- colnames(combined$mean.scores)[max.col(combined$mean.scores)]

    # write scores for figure reference
    if (grepl('postSelect', sample_files[[sample]])) {
      filename <- file.path(output_dir, sample, 'postSelect', 'annotation_kriegstein', paste0("Pearson.correlation.mean_", name, ".csv"))
    } else {
      filename <- file.path(output_dir, sample, 'annotation_kriegstein', paste0("Pearson.correlation.mean_", name, ".csv"))
    }
    utils::write.csv2(combined$mean.scores, file = filename)

    return(combined)
  })
  names(combined.results) <- names(results.list)

  # save Kriegstein cluster labels into Seurat object --> .qs
  for (sample in sample_names) {
    for (anno in annotations) {
      # Kriegstein data and labels from sample-annotation correlation overlap
      misc_data <- combined.results[[paste(sample, anno)]]
      # add as miscellaneous data to Seurat object
      SeuratObject::Misc(object = data.list[[sample]], slot = paste0("Kriegstein.SingleR.", anno)) <- misc_data
    }

    SeuratObject::Misc(object = data.list[[sample]], slot = "Kriegstein.SingleR.ref_aggr_strategy") <- ref_aggr_strategy
    # copy seurat clusters metadata to aggregate with Kriegstein labels
    data.list[[sample]]$kriegstein.seurat.custom.clusters.mean <- data.list[[sample]]$seurat_clusters
    # overwrite copied metadata to an aggregation of Kriegstein.Seurat custom cluster labels
    levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean) <- paste0(combined.results[[paste(sample, "custom.clusterv2")]]$mean.labels, ".", levels(data.list[[sample]]$kriegstein.seurat.custom.clusters.mean))

    # plot composition
    composition_df <- as.data.frame(table(data.list[[sample]]@meta.data[,c("orig.ident", "kriegstein.seurat.custom.clusters.mean")]))
    composition_df <- composition_df %>%
      dplyr::arrange(orig.ident, kriegstein.seurat.custom.clusters.mean) %>%
      dplyr::group_by(orig.ident) %>%
      dplyr::mutate(csum = cumsum(Freq))
    openxlsx::write.xlsx(composition_df, file = file.path(output_dir, sample, 'annotation_kriegstein', 'samples_composition.xlsx'))
    composition_df$kriegstein.seurat.custom.clusters.mean <- factor(composition_df$kriegstein.seurat.custom.clusters.mean, levels = rev(levels(composition_df$kriegstein.seurat.custom.clusters.mean)))
    # remove labels with less than 10 cells
    composition_df$Freq[composition_df$Freq < 10] <- NA
    p <- ggplot2::ggplot(composition_df, ggplot2::aes(x = orig.ident, y = Freq, fill = kriegstein.seurat.custom.clusters.mean)) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      # ggrepel::geom_text_repel(ggplot2::aes(x = orig.ident, y = csum, label=csum), color="white", size=3, max.overlaps = Inf) +
      ggplot2::geom_text(ggplot2::aes(x = orig.ident, y = csum, label=Freq), nudge_y = -75, color="white", size=3) + #
      # ggplot2::geom_text() +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "Kriegstein.Seurat cluster")) +
      ggplot2::labs(x = "", y = "# cells") +
      ggplot2::theme_minimal()
    p
    # p + ggplot2::scale_fill_manual()
    ggplot2::ggsave(plot = p, file = file.path(output_dir, sample, 'annotation_kriegstein', 'samples_composition.png'), width = 30, height = 20, units = "cm")

    # overwrite .qs file with new misc annotation
    qs::qsave(data.list[[sample]], sample_files[[sample]], preset = 'custom', algorithm = "zstd_stream", compress_level = 4, shuffle_control = 15, nthreads = 1)
  }



  # custom visualizations per sample-reference comparison for each annotation
  for (sample in sample_names) {
    # first plot umaps of Kriegstein annotated data
    p <- Seurat::DimPlot(data.list[[sample]], reduction = "umap", group.by = "kriegstein.seurat.custom.clusters.mean", pt.size = .75)
    p + ggplot2::labs(title = "Kriegstein annotated clusters")
    ggplot2::ggsave(file = file.path(output_dir, sample, "annotation_kriegstein", paste0("UMAP_", sample, ".png")), width = 30, height = 20, units = "cm")
    p <- Seurat::DimPlot(data.list[[sample]], reduction = "umap", group.by = "kriegstein.seurat.custom.clusters.mean", label = T) + Seurat::NoLegend()
    p + ggplot2::labs(title = "Kriegstein annotated clusters")
    ggplot2::ggsave(file = file.path(output_dir, sample, "annotation_kriegstein", paste0("UMAP_noLegend_", sample, ".png")), width = 30, height = 20, units = "cm")

    # then prepare data for plotting Kriegstein annotation heatmap
    annotation_col <- data.frame(row.names = levels(data.list[[sample]]$seurat_clusters))
    annotation_colors <- list()

    for (anno in annotations) {
      # set annotation column for transferred labels from reference data
      annotation_col[ , ncol(annotation_col) + 1] <- data.list[[sample]]@misc[[paste0("Kriegstein.SingleR.", anno)]][[paste0(ref_aggr_strategy, ".labels")]]
      colnames(annotation_col)[ncol(annotation_col)] <- paste0("ref.", anno)

      # get ordered and unique label names from reference data
      labels <- unique(meta[which(colnames(meta) == anno)][,1])
      annotation_colors[[length(annotation_colors) + 1]] <- labels[order(labels)]
      names(annotation_colors)[length(annotation_colors)] <- paste0("ref.", anno)

      # set colors to each unique label from reference data
      names(annotation_colors[[paste0("ref.", anno)]]) <- generate_color_palette(type = 'mixed', n = length(annotation_colors[[paste0("ref.", anno)]]))

      # swap names and values of named vector to proper formatting
      annotation_colors[[paste0("ref.", anno)]] <- stats::setNames(names(annotation_colors[[paste0("ref.", anno)]]), annotation_colors[[paste0("ref.", anno)]])

      # set annotation colors for all transferred labels
      annotation_colors[[paste0("ref.", anno)]] <- annotation_colors[[paste0("ref.", anno)]][unique(annotation_col[which(colnames(annotation_col) == paste0("ref.", anno))][,1])]
    }

    # order annotation columns for visualization
    names(annotation_col) <- c("fetal.brain.age", "fetal.brain.structure", "fetal.brain.celltype")
    annotation_col <- annotation_col[, c("fetal.brain.celltype", "fetal.brain.structure", "fetal.brain.age")]

    for (anno in annotations_to_plot) {
      if (grepl('postSelect', sample_files[[sample]])) {
        filename <- file.path(output_dir, sample, 'postSelect', 'annotation_kriegstein', paste0("Heatmap_", sample, "_", anno, ".png"))
      } else {
        filename <- file.path(output_dir, sample, 'annotation_kriegstein', paste0("Heatmap_", sample, "_", anno, ".png"))
      }
      message("plotting: ", filename)

      # set rownames for identification of rows during plotting
      rownames(combined.results[[paste(sample, anno)]][[paste0(ref_aggr_strategy, ".scores")]]) <- levels((data.list[[sample]]$seurat_clusters))

      # plot pretty heatmap
      p <- pheatmap::pheatmap(t(combined.results[[paste(sample, anno)]][[paste0(ref_aggr_strategy, ".scores")]]),
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
}

#' Get Kriegstein genes
#'
#' @param kriegstein_data_dir folder path to Kriegstein data
#'
#' @return Kriegstein genes character vector.
getGenes <- function(kriegstein_data_dir) {
  genesFile <- file.path(kriegstein_data_dir, "kriegstein_genes.csv")
  if (!file.exists(genesFile)) {
    stop(genesFile, " does not exist, generate this file with EMC.SKlab.scRNAseq::chunk_kriegstein_data().")
  } else {
    genes <- utils::read.table(genesFile, sep = "\t", col.names = 'gene')
    return(genes$gene[2:length(genes$gene)])
  }
}

#' Get Kriegstein meta with custom annotation mapping
#'
#' @param kriegstein_data_dir folder path to Kriegstein data
#'
#' @return Kriegstein meta features with custom clusterv2 mapping.
getMeta <- function(kriegstein_data_dir) {
  metaFile <- file.path(kriegstein_data_dir, "custom.meta.tsv")
  if (!file.exists(metaFile)) {
    stop(metaFile, " does not exist, generate this file with EMC.SKlab.scRNAseq::chunk_kriegstein_data().")
  } else {
    return(utils::read.table(metaFile, header = TRUE, sep = "\t", as.is = TRUE, row.names = 1))
  }
}
