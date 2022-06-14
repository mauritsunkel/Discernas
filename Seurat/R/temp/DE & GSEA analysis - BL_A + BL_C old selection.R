### start initialization
library(Seurat)
library(future)
library(patchwork)
library(ggplot2)
library(chron)
library(tidyr)
library(dplyr)


## TODO build in SCT based DEG: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
# SCT DEG: we make use of 'corrected counts' that are stored in the data slot of the the SCT assay


## use SCT assay for visualization
# DefaultAssay(immune.combined.sct) <- "SCT"
# visualization_functions_here()



# SCTransform combines NormalizeData (better normalization), ScaleData and FindVariableFeatures
## v2 upgrades: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#perform-integration-using-pearson-residuals-1
# normally SCT used for dadta transformations that are then used for: PCA, clustering, UMAP
# before, RNA assay was used for DE, as SCT Pearson residuals could not resolve batch effects in comparison
## now: this is the optimized way of doing DE analysis
# can use top VariableFeatureGenes for marker identifictation: https://www.biostars.org/p/406388/
# RNA assay: counts@rawcounts / TPMs, data@normalized data, scale.data@scaled data
# SCT assay: counts@corrected counts, data@log1p(counts), scale.data@pearson residuals


# TODO commented out FGSEA analysis for now, focus on DEGA

# TODO put in it's own file and source from here or pipeline
## TODO tryout sourcing in test file by setting working directory and seeing how sourcing behaves with a mock function
# TODO add rm cleanup environment to end of function
# TODO check maximum Windows file length and output pathway (partial) name instead of order for enrichment plot
# FGSEA package: vignette: http://127.0.0.1:31440/library/fgsea/doc/fgsea-tutorial.html
FGSEA_analysis <- function(markers, working_directory, marker_type, cluster) {
  library(biomaRt)
  library(fgsea)
  library(data.table)
  library(ggplot2)

  dir.create(paste0(work_dir, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/"), recursive = T)

  ## fix infinite values later by applying -log10 function
  markers$p_val[markers$p_val == 0] <- min(markers$p_val[markers$p_val != 0])

  ## calculate metric by FoldChangeSign and -LogPvalue
  markers$fcsign <- sign(markers$avg_log2FC)
  markers$logPval <- -log10(markers$p_val)

  ## create ranked vector
  fgsea_ranks <- markers$logPval/markers$fcsign

  ## get Entrez IDs by HGNC symbol to match gene names and provide a translation map
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
  mapping <- getBM(
    attributes = c('entrezgene_id', 'hgnc_symbol'),
    filters = 'hgnc_symbol',
    values = rownames(markers),
    mart = hsmart
  )
  names(fgsea_ranks) <- match(rownames(markers), mapping$hgnc_symbol)

  ## get Reactome pathways by Entrez IDs
  pathways <- reactomePathways(names(fgsea_ranks))

  fgsea_results <- fgsea(pathways = pathways,
                         stats    = fgsea_ranks,
                         eps      = 0.0,
                         minSize  = 15,
                         maxSize  = 500)

  topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

  png(filename=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/overview_table.png"), width = 1600)
  plotGseaTable(pathways[topPathways], fgsea_ranks, fgsea_results,
                     gseaParam=0.5)
  dev.off()

  ## can try to collapse pathways if there are many seemlingly alike in the plot above
  collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < 0.01],
                                        pathways, fgsea_ranks)
  mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
    order(-NES), pathway]

  ## check if mainPathways is empty (likely collapsedPathways is empty too)
  if (length(mainPathways) > 0) {
    png(filename=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/collapsed_table.png"), width = 1600)
    p <- plotGseaTable(pathways[mainPathways], fgsea_ranks, fgsea_results,
                       gseaParam = 0.5)
    dev.off()
  }

  fwrite(fgsea_results, file=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/overview.xls"), sep="\t", sep2=c("", ",", ""))

  for (i in seq_along(topPathways)) {
    # png(filename=paste0(working_directory, "GSEA/cluster_", cluster, "/enriched_", i, ".png"), width = 1600)
    p <- plotEnrichment(pathways[[topPathways[i]]],
                        fgsea_ranks) + labs(title=topPathways[[i]])
    ggsave(file = paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/enriched_", i, ".png"), width = 30, height = 20, units = "cm")
  }
}



## OVERRIDE SEURAT DE FUNCTIONS WITH ADDITIONAL FUNCTIONALITY
# fixInNamespace(FindMarkers.Assay, pos = "package:Seurat") # use this to copy source code
# trace(Seurat:::FindMarkers.Assay, edit = T)
## but then need: untrace(Seurat:::FindMarkers.Assay) # to undo edits by tracing, will also be undone on reloading R
### all these functions call edit() under the hood
FoldChange.default.adjusted <- function (object, cells.1, cells.2, mean.fxn, fc.name, mean.fxn.adj = mean.fxn.adj, features = NULL, ...)
{
  features <- features %||% rownames(x = object)
  thresh.min <- 0
  pct.1 <- round(x = rowSums(x = object[features, cells.1,
                                        drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  pct.2 <- round(x = rowSums(x = object[features, cells.2,
                                        drop = FALSE] > thresh.min)/length(x = cells.2), digits = 3)
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  fc <- (data.1 - data.2)

  ### MY INJECTED CUSTOM CODE
  n_nonzero.1 <- rowSums(x = object[features, cells.1, drop = FALSE] > 0)
  n_nonzero.2 <- rowSums(x = object[features, cells.2, drop = FALSE] > 0)
  object[object==0] <- NA

  # Seurat code
  mean.fxn <- mean.fxn %||% switch(
    EXPR = slot,
    'data' = function(x) {
      return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
    },
    'scale.data' = rowMeans,
    function(x) {
      return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
    }
  )
  ## Custom code for temporary custom mean.fxn.adj to calculate proper rowMeans for non-zero expressing cells
  # added: na.rm = TRUE, such that 0's not taken into account for row means
  # pseudocount.use = 1, base = 2, hardcoded because of namespace issues
  mean.fxn.adj <- function (x)
  {
    return(log(x = mean(x = expm1(x = x), na.rm = T) + 1,
               base = 2))
  }
  mean.fxn.adj <- mean.fxn.adj %||% switch(
    EXPR = slot,
    'data' = function(x) {
      return(log(x = mean(x = expm1(x = x), na.rm = T) + 1, base = 2))
    },
    'scale.data' = rowMeans,
    function(x) {
      return(log(x = mean(x = x, na.rm = T) + 1, base = 2))
    }
  )
  nonzero_data.1 <- apply(object[features, cells.1, drop = FALSE], 1, mean.fxn.adj)
  nonzero_data.2 <- apply(object[features, cells.2, drop = FALSE], 1, mean.fxn.adj)
  nonzero_fc <- (nonzero_data.1 - nonzero_data.2)
  nonzero_fc[is.na(nonzero_fc)] <- 0

  # ADJUSTED by adding: nonzero_fc, n_nonzero.1, n_nonzero.2
  fc.results <- as.data.frame(x = cbind(fc, pct.1, pct.2, data.1, data.2, nonzero_fc, n_nonzero.1, n_nonzero.2,
                                        nonzero_data.1, nonzero_data.2))
  # ADJUSTED by adding: "nz_log_fc", "n_nz.1", "n_nz.2"
  colnames(fc.results) <- c(fc.name, "pct.1", "pct.2", "meanExpr.1", "meanExpr.2",
                            paste0("nz_", fc.name), "nz_n.1", "nz_n.2", "nz_meanExpr.1", "nz_meanExpr.2")

  return(fc.results)
}
# WORKS, but, if it fails after all continue searching
## https://stackoverflow.com/questions/8204008/redirect-intercept-function-calls-within-a-package-function
# namespace of customFunction is R_globalenv, where it is defined, bur should be Seurat as that is ns of my targeted function
environment(FoldChange.default.adjusted) <- asNamespace("Seurat")
# then with assignInNameSpace I can basically inject my code their copied function and then substitute it back in their environment
assignInNamespace("FoldChange.default", FoldChange.default.adjusted, ns = "Seurat")

WilcoxDETest.adjusted <- function (data.use, cells.1, cells.2, verbose = TRUE, ...)
{
  # save data.use for non-zero calculation in case it is filtered too much beforehand
  data.use.orig <- data.use


  # use only data needed for group comparison
  data.use <- data.use[, c(cells.1, cells.2), drop = FALSE]
  # create sequential index of group 1 cells
  j <- seq_len(length.out = length(x = cells.1))
  # use ProgressBarSApply or FutureSApply (sequential/parallel processing)
  my.sapply <- ifelse(test = verbose && future::nbrOfWorkers() ==
                        1, yes = pbsapply, no = future_sapply)
  # check if not overflowing (data is NaN)
  overflow.check <- ifelse(test = is.na(x = suppressWarnings(length(x = data.use[1,
  ]) * length(x = data.use[1, ]))), yes = FALSE, no = TRUE)
  # check if limma package available in R session
  limma.check <- PackageCheck("limma", error = FALSE)
  # if data not overflowing and limma package available in R session
  if (limma.check[1] && overflow.check) {
    # calculate p-value with defined sapply function
    p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
      return(min(2 * min(limma::rankSumTestWithCorrelation(index = j,
                                                           statistics = data.use[x, ])), 1))
    })
  }
  else {
    if (getOption("Seurat.limma.wilcox.msg", TRUE) && overflow.check) {
      message("For a more efficient implementation of the Wilcoxon Rank Sum Test,",
              "\n(default method for FindMarkers) please install the limma package",
              "\n--------------------------------------------",
              "\ninstall.packages('BiocManager')", "\nBiocManager::install('limma')",
              "\n--------------------------------------------",
              "\nAfter installation of limma, Seurat will automatically use the more ",
              "\nefficient implementation (no further action necessary).",
              "\nThis message will be shown once per session")
      options(Seurat.limma.wilcox.msg = FALSE)
    }
    group.info <- data.frame(row.names = c(cells.1, cells.2))
    group.info[cells.1, "group"] <- "Group1"
    group.info[cells.2, "group"] <- "Group2"
    group.info[, "group"] <- factor(x = group.info[, "group"])
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
      return(wilcox.test(data.use[x, ] ~ group.info[,
                                                    "group"], ...)$p.value)
    })
  }

  ### CUSTOM CODE for calculating p-value for non-zero expression cells
  # use only non-zero expression data needed for group comparison
  data.use <- data.use.orig[, c(cells.1, cells.2), drop = FALSE]

  if (limma.check[1] && overflow.check) {
    # calculate p-value with defined sapply function
    nz_p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
      return(min(2 * min(limma::rankSumTestWithCorrelation(index = seq_len(sum(names(data.use[x,][data.use[x,] > 0]) %in% cells.1)),
                                                           statistics = data.use[x,][data.use[x,] > 0]
      )), 1))
    })
  }
  else {
    data.use <- data.use[, rownames(x = group.info), drop = FALSE]
    nz_p_val <- my.sapply(X = 1:nrow(x = data.use), FUN = function(x) {
      return(wilcox.test(data.use[x, ] ~ group.info[,
                                                    "group"], ...)$p.value)
    })
  }

  return(data.frame(p_val, nz_p_val, row.names = rownames(x = data.use)))
}
# change namespace of adjusted function to target function
environment(WilcoxDETest.adjusted) <- asNamespace("Seurat")
# now overrice target function within that namespace with my custom function
assignInNamespace("WilcoxDETest", WilcoxDETest.adjusted, ns = "Seurat")

FindMarkers.default.adjusted <- function(
    object,
    slot = "data",
    counts = numeric(),
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    logfc.threshold = 0.25,
    test.use = 'wilcox',
    min.pct = 0.1,
    min.diff.pct = -Inf,
    verbose = TRUE,
    only.pos = FALSE,
    max.cells.per.ident = Inf,
    random.seed = 1,
    latent.vars = NULL,
    min.cells.feature = 3,
    min.cells.group = 3,
    pseudocount.use = 1,
    fc.results = NULL,
    densify = FALSE,
    ...
) {
  ValidateCellGroups(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  features <- features %||% rownames(x = object)
  # reset parameters so no feature filtering is performed
  if (test.use %in% DEmethods_noprefilter()) {
    features <- rownames(x = object)
    min.diff.pct <- -Inf
    logfc.threshold <- 0
  }
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  # feature selection (based on logFC)
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1] #first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }
  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  de.results <- PerformDE(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    verbose = verbose,
    min.cells.feature = min.cells.feature,
    latent.vars = latent.vars,
    densify = densify,
    ...
  )
  de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% DEmethods_nocorrect()) {
    de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -de.results[, 1]), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
    # MY CUSTOM CODE: also perform Bonferroni correction for non-zero expression DE data
    de.results$nz_p_val_adj = p.adjust(
      p = de.results$nz_p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
    # sort table column names
    de.results <- de.results[,c(1, 13, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 14)]
  }
  return(de.results)
}
# change namespace of adjusted function to target function
environment(FindMarkers.default.adjusted) <- asNamespace("Seurat")
# now overrice target function within that namespace with my custom function
assignInNamespace("FindMarkers.default", FindMarkers.default.adjusted, ns = "Seurat")
### END DE INITIALIZATION ###



### USER PARAMETERS
# read an integrated saved RDS file
sample_name <- "BL_A + BL_C"
integrated <- readRDS("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_corrected_13-06/integrated - old selection/BL_A + BL_C/BL_A + BL_C.rds")

# set amount of cells used for 'downsampling' clusters during FindMarkers function (max amount of cells per cluster)
nCellsDownsampling <- Inf

# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"

# load future library and set plan to run certain functions with multiprocessing
plan("multisession", workers = 1) # DEVNOTE: n_workers > 1 for parallelization (for me, 5 is max, 4 is safe)

# play system sounds, call function to alarm user that run is done!
beep <- function(n = 5){
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
}

start_time <- format(Sys.time(), "%F %H-%M-%S")
dir.create(paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/DE_analysis/"), recursive = T)
work_dir <- paste0(work_dir, 'results/', start_time, '/integrated/', sample_name, "/DE_analysis/")
setwd(work_dir)
dir.create(paste0(work_dir, "../GSE_analysis/"))
dir.create(paste0(work_dir, "sample_markers/"))
dir.create(paste0(work_dir, "markers/"))
dir.create(paste0(work_dir, "conserved_markers/"))
dir.create(paste0(work_dir, "condition_markers/"))
### END USER PARAMETERS

### create marker dfs to count N cells used in comparisons
## sample markers
sample_markers_columns <- c(paste0("n_cells_", names(table(integrated$orig.ident))[1]),
                            paste0("n_cells_", names(table(integrated$orig.ident))[2]))
sample_markers_df <- data.frame(matrix(nrow = 0, ncol = length(sample_markers_columns)))
colnames(sample_markers_df) <- sample_markers_columns
sample_markers_df[nrow(sample_markers_df) + 1,] = c(table(integrated$orig.ident)[1],
                                                    table(integrated$orig.ident)[2])

# set idents to compare cells at sample level instead of cluster level
Idents(integrated) <- integrated$orig.ident
# get sample vs sample markers (now: monoculture vs coculture (for neurons and astrocytes))
## note: p_val_adj = Adjusted p-value, based on Bonferroni correction using all genes (including non-zero expression) in the dataset
### adjusted both defaults: logfc.threshold = 0.25, min.pct = 0.1    to 0
sample_markers <- FindMarkers(integrated, assay = "SCT", ident.1 = names(table(integrated$orig.ident))[1], only.pos = FALSE, verbose = T,
                              logfc.threshold = 0, min.pct = 0)

# set pct variable based on BL_C orig.identity index
if (names(table(integrated$orig.ident))[1] == "BL_C") {
  pct <- "pct.1"
} else {
  pct <- "pct.2"
}

# filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferronu correction
sample_markers <- sample_markers[!(sample_markers$p_val_adj > 0.05 & sample_markers$nz_p_val_adj > 0.05),]
# filter on pct.ref and order by avg_log2FC
sample_markers_pval_adj <- sample_markers %>% arrange(desc(avg_log2FC)) # DEPRECATED: filter(pct > 0.1)

# sample_markers_pval_adj <- sample_markers %>% filter(p_val_adj <= 0.05) %>% filter(pct > 0.1) %>% arrange(desc(avg_log2FC))
write.csv2(sample_markers_pval_adj, file = paste0("sample_markers/pct1=", names(table(integrated$orig.ident))[1], "-pct2=", names(table(integrated$orig.ident))[2], " - (nz-)p-val st 0.05.csv"), row.names = TRUE)
# FGSEA_analysis(markers = sample_markers, working_directory = work_dir, marker_type = 'sample_markers', cluster = pct)

# set idents back to cluster level
Idents(integrated) <- integrated$seurat_clusters

# Note: Custom Bonferroni correction based on amount of genes TESTED would be way less stringent then using n_genes in dataset, as we filter features for downstream processing before DE
## Note: alternatives to this stringent Bonferroni are not supported by Seurat but can be custom made on the output p-values of FindMarkers functions
### Note: answer to Steven specific discussion: https://github.com/satijalab/seurat/issues/4112
# dim(integrated[which(rowSums(integrated)!=0),]) # 22056, number of nonzero expression genes
# dim(integrated) # 23173, number of total genes, also if they have 0 expression, used for DE Bonferroni correction



## markers
markers_columns <- c('cluster_ID', 'n_cells_cluster', 'n_all_other_cells')
markers_df <- data.frame(matrix(nrow = 0, ncol = length(markers_columns)))
colnames(markers_df) <- markers_columns
## conserved markers
conserved_markers_columns <- c('cluster_ID',
                               paste0('n_cells_cluster_identity_', names(table(integrated$orig.ident))[1]),
                               paste0('n_cells_cluster_identity_', names(table(integrated$orig.ident))[2]),
                               paste0('n_all_other_cells_identity_', names(table(integrated$orig.ident))[1]),
                               paste0('n_all_other_cells_identity_', names(table(integrated$orig.ident))[2]))
conserved_markers_df <- data.frame(matrix(nrow = 0, ncol = length(conserved_markers_columns)))
colnames(conserved_markers_df) <- conserved_markers_columns
## conition markers
condition_markers_columns <- c('cluster_id',
                               paste0('n_cells_', names(table(integrated$orig.ident))[1]),
                               paste0('n_cells_', names(table(integrated$orig.ident))[2]))
condition_markers_df <- data.frame(matrix(nrow = 0, ncol = length(condition_markers_columns)))
colnames(condition_markers_df) <- condition_markers_columns

## get cluster IDs to loop over
cluster_ids <- levels(integrated$seurat_clusters)


for (i in cluster_ids) {
  print(paste('Cluster ID:', i))

  ## add amount of cells used for markers comparison to df
  markers_df[nrow(markers_df) + 1,] = c(i,
                                        table(integrated$seurat_clusters)[i],
                                        sum(table(integrated$seurat_clusters)[-as.integer(i)]))
  ## create markers for integrated data for each cluster vs all other clusters
  markers <- FindMarkers(integrated, assay = "SCT", ident.1 = i, max.cells.per.ident = nCellsDownsampling, only.pos = FALSE, verbose = T)
  # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
  markers <- markers[!(markers$p_val_adj > 0.05 & markers$nz_p_val_adj > 0.05),]
  write.csv2(markers, file = paste0("markers/all_cluster", i, "_m.csv"))
  # FGSEA_analysis(markers = markers, working_directory = work_dir, marker_type = 'markers', cluster = i)




  ## add amount of cells used for conserved_markers comparison to df
  df <- data.frame('orig.ident' = integrated$orig.ident, 'seurat_clusters' = integrated$seurat_clusters)
  # condition 1 & match cluster
  cm_val1 <- nrow(df %>% filter(orig.ident == names(table(integrated$orig.ident))[1] & seurat_clusters == i))
  # condition 2 & match cluster
  cm_val2 <- nrow(df %>% filter(orig.ident == names(table(integrated$orig.ident))[2] & seurat_clusters == i))
  # condition 1 & no match cluster
  cm_val3 <- nrow(df %>% filter(orig.ident == names(table(integrated$orig.ident))[1] & seurat_clusters != i))
  # condition 2 & no match cluster
  cm_val4 <- nrow(df %>% filter(orig.ident == names(table(integrated$orig.ident))[2] & seurat_clusters != i))
  conserved_markers_df[nrow(conserved_markers_df) + 1,] = c(i, cm_val1, cm_val2, cm_val3, cm_val4)
  ## create markers conserved between groups (conditions) for integrated data for each cluster vs all other clusters
  conserved_markers <- FindConservedMarkers(integrated, assay = "SCT", ident.1 = i, only.pos = FALSE,
                                            max.cells.per.ident = nCellsDownsampling,
                                            grouping.var = "orig.ident", verbose = T)
  # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
  conserved_markers <- conserved_markers[!(conserved_markers$p_val_adj > 0.05 & conserved_markers$nz_p_val_adj > 0.05),]
  # TODO check if filters are still correct now that order of column names is different etc
  head(conserved_markers)
  # pos_conserved_markers <- conserved_markers %>% filter((.[[2]] > 0) & (.[[7]] > 0))
  # neg_conserved_markers <- conserved_markers %>% filter((.[[2]] < 0) | (.[[7]] < 0))
  write.csv2(conserved_markers, file = paste0("conserved_markers/all_cluster", i, "_cm.csv"))
  # FGSEA_analysis(markers = conserved_markers, working_directory = work_dir, marker_type = 'conserved_markers', cluster = i)








  # ## create condition markers for integrated data within each cluster between each condition
  # ## DEV NOTE: this is not pairwise if more than 2 conditions are integrated at the same time
  subset <- subset(integrated, seurat_clusters == i)
  # ##  change cluster identity to original identity to find markers between conditions
  Idents(subset) <- subset$orig.ident
  ## check if subset contains cells for at least 2 condtions/samples for comparison
  ### DEVNOTE: also need to check if a condition has ENOUGH cells for comparison?
  if (length(names(table(subset$orig.ident))) == 1) {
    print(paste0('In cluster ', i, ' only cells for condition/sample ',
                 names(table(subset$orig.ident)), ' were found, cannot create condition markers for this cluster.'))
    next
  }
  ## add amount of cells used for condition_markers comparison to df
  condition_markers_df[nrow(condition_markers_df) + 1,] = c(i, table(subset$orig.ident)[1], table(subset$orig.ident)[2])
  ## create condition_markers for subset data for within each cluster to compare conditions
  condition_markers <- FindMarkers(subset, assay = "SCT", recorrect_umi = FALSE, ident.1 = "BL_C", verbose = T, only.pos = FALSE)
  # filters rows (genes) if they are >0.05 for both p_val and non-zero p_val with Bonferroni correction
  condition_markers <- condition_markers[!(condition_markers$p_val_adj > 0.05 & condition_markers$nz_p_val_adj > 0.05),]
  # pos_condition_markers <- condition_markers %>% filter(avg_log2FC > 0)
  # neg_condition_markers <- condition_markers %>% filter(avg_log2FC < 0)
  write.csv2(condition_markers, file = paste0("condition_markers/all_cluster", i, ".csv"))
  print(paste('Cluster ID:', i, ' before condition_markers FGSEA call'))
  # FGSEA_analysis(markers = condition_markers, working_directory = work_dir, marker_type = 'condition_markers', cluster = i)

  # DEVNOTE if want to assign each table to its own variable, use assign() and get()
  # assign(paste0("cluster", i, "_markers"), markers)
  ## get(paste0("cluster", i, "_markers"))
}
## write n cells for comparison to CSV files
write.csv2(sample_markers_df, file = "sample_markers/n_cells_for_comparison_m.csv", row.names = FALSE)
write.csv2(markers_df, file = "markers/n_cells_for_comparison_m.csv", row.names = FALSE)
write.csv2(conserved_markers_df, file = "conserved_markers/n_cells_for_comparison_cm.csv", row.names = FALSE)
write.csv2(condition_markers_df, file = "condition_markers/n_cells_for_comparison.csv", row.names = FALSE)

beep()



# topGO vignette: http://127.0.0.1:25748/library/topGO/doc/topGO.pdf


# # DEVNOTE: load markers and prep for topGoData class object
# go <- markers$p_val
# hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
# mapping <- getBM(
#   attributes = c('entrezgene_id', 'hgnc_symbol'),
#   filters = 'hgnc_symbol',
#   values = rownames(markers),
#   mart = hsmart
# )
# names(go) <- paste0(match(rownames(markers), mapping$hgnc_symbol), '_at')
#
# # DEVNOTE used this to find Affymetrix glossary for _at _f _g _i _r _s
# table(sapply(strsplit(names(geneList), '_'), "[[", 2))
#
# library(topGO)
# library(ALL) # Acute Lymphoblastic Leukemia
# data(ALL) # Acute Lymphoblastic Leukemia
# data(geneList)
#
# affyLib <- paste(annotation(ALL), "db", sep = ".")
# library(package = affyLib, character.only = TRUE)
# # sum(topDiffGenes(geneList)) # same as sum(geneList < 0.01)
#
# # create topGoData class object
# ## DEVNOTE I don't see how I can do this for our data, as they use a certain database for their specific
# ## affymetrix probes that they used in the ALL experiment, I can add _at behind our integer Entrez IDs which
# ## seems to work but I cannot be sure that ID 1000 in Entrez is same gene in their affylib (hgu95av2.db) database
# ### code in comments  below is from vignette chapter 4.3 Custom Annotation, but this doesn't seem applicable for us
# #### geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
# #### str(head(geneID2GO))
# sampleGOdata <- new("topGOdata",
#                     description = "Simple session", ontology = "BP",
#                     allGenes = go, geneSel = topDiffGenes,
#                     nodeSize = 10,
#                     annot = annFUN.db, affyLib = affyLib)
# sampleGOdata
#
# # TODO which tests would be proper for our data?
# resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
# resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
# resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")
#
# allRes <- GenTable(sampleGOdata,
#                    classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim,
#                    orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)
# allRes
#
# colMap <- function(x) {
#   .col <- rep(rev(heat.colors(length(unique(x)), alpha = 1)), time = table(x))
#   return(.col[match(1:length(x), order(x))])
# }
#
# pValue.classic <- score(resultKS)
# pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
# gstat <- termStat(sampleGOdata, names(pValue.classic))
# gSize <- gstat$Annotated / max(gstat$Annotated) * 4
# gCol <- colMap(gstat$Significant)
# plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
#        pch = 19, cex = gSize, col = gCol)
#
# sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
# cbind(termStat(sampleGOdata, sel.go),
#         elim = pValue.elim[sel.go],
#         classic = pValue.classic[sel.go])
#
# showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')


### TESTING BAS CUSTOM GENE LISTS FOR DE
# bas <- read.table("C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/data/Bas Genes of interest scRNA seq DE.txt",
#                   sep = '\t')
# bas_markers <- c(bas$V1)
# sample_markers[rownames(sample_markers) %in% bas_markers,]
#
# bas[bas_markers %in% rownames(sample_markers),] # which genes are in DE at all
# bas[!bas_markers %in% rownames(sample_markers),] # which genes are not in DE at all
# bas[bas_markers %in% rownames(integrated@assays$RNA@meta.features),] # which genes are in all features
# bas[!bas_markers %in% rownames(integrated@assays$RNA@meta.features),] # which genes are not in all features - GJB6

# SYN synonyms: SYN1 & SYN3
# PSD95 synonym: DLG4 (SAP90)
# ---
# GJB6 synonyms: CX30, EDH, HED, DFNA3
###
