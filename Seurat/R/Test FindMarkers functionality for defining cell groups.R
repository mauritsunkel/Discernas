
library(Seurat)
library(future) # https://satijalab.org/seurat/articles/future_vignette.html#futurized-functions-in-seurat

### USER PARAMETERS
# read an integrated saved RDS file
sample_name <- "BL_A + BL_C"
integrated <- readRDS(paste0("C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/Exploration results/SCTransform + Leiden - Cellcycle/integrated/", sample_name, "/integrated.rds"))

# work dir should contain forward slashes (/) on Windows
work_dir <- "C:/Users/mauri/Desktop/M/Work & Education/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/"
### END USER PARAMETERS

work_dir <- paste0(work_dir, 'results/')
dir.create(work_dir)
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- paste0(work_dir, start_time, '/')
dir.create(work_dir)
work_dir <- paste0(work_dir, 'integrated/')
dir.create(work_dir)
work_dir <- paste0(work_dir, sample_name, "/")
dir.create(work_dir)
work_dir <- paste0(work_dir, 'DE_analysis/')
dir.create(work_dir)
setwd(work_dir)
### end initialization ###



# fixInNamespace(FindMarkers.Assay, pos = "package:Seurat") # use this to copy source code
# trace(Seurat:::FindMarkers.Assay, edit = T)
## but then need: untrace(Seurat:::FindMarkers.Assay) # to undo edits by tracing, will also be undone on reloading R
### all these functions call edit() under the hood
FoldChange.default.adjusted <- function (object, cells.1, cells.2, mean.fxn, fc.name, mean.fxn.adj = mean.fxn.adj, features = NULL, ...)
{
  features <- features %||% rownames(x = object)
  features <- c('HES6') # TODO REMOVE AFTER TESTING
  thresh.min <- 0
  # print(paste("HES6-cells.1:", head(object["HES6", cells.1, drop = FALSE])))
  # print(paste("HES6-cells.2:", head(object["HES6", cells.2, drop = FALSE])))
  pct.1 <- round(x = rowSums(x = object[features, cells.1,
                                        drop = FALSE] > thresh.min)/length(x = cells.1), digits = 3)
  pct.2 <- round(x = rowSums(x = object[features, cells.2,
                                        drop = FALSE] > thresh.min)/length(x = cells.2), digits = 3)
  print(paste("HES6-pct.1:", pct.1))
  print(paste("HES6-pct.2:", pct.2))
  data.1 <- mean.fxn(object[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(object[features, cells.2, drop = FALSE])
  print(paste("HES6-meanExpression.1:", data.1))
  print(paste("HES6-meanExpression.2:", data.2))
  fc <- (data.1 - data.2)
  print(paste("HES6-fc:", fc))

  ### MY INJECTED CUSTOM CODE
  n_nonzero.1 <- rowSums(x = object[features, cells.1, drop = FALSE] > 0)
  n_nonzero.2 <- rowSums(x = object[features, cells.2, drop = FALSE] > 0)
  print(paste("HES6-n_nz.1:", n_nonzero.1))
  print(paste("HES6-n_nz.2:", n_nonzero.2))
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
  print(paste("HES6-nz_meanExpression.1:", nonzero_data.1))
  print(paste("HES6-nz_meanExpression.2:", nonzero_data.2))
  nonzero_fc <- (nonzero_data.1 - nonzero_data.2)
  print(paste("HES6-nz_fc:", nonzero_fc))
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
  data.use <- data.use['HES6', c(cells.1, cells.2), drop = FALSE] # TODO remove HES6 after testing
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
  data.use <- data.use.orig['HES6', c(cells.1, cells.2), drop = FALSE] # TODO remove HES6 after testing
  print(paste("HES6-p_val_in-cells.1:", table(data.use['HES6',cells.1] > 0)))
  print(paste("HES6-p_val_in-cells.2:", table(data.use['HES6',cells.2] > 0)))
  print(paste("HES6-p_val-statistic:", table(data.use['HES6', ] > 0)))
  print(paste("HES6-nz_p_val_in-cells.1:", table(data.use['HES6',cells.1][data.use['HES6',cells.1] > 0] > 0)))
  print(paste("HES6-nz_p_val_in-cells.2:", table(data.use['HES6',cells.2][data.use['HES6',cells.2] > 0] > 0)))
  print(paste("HES6-nz_p_val-statistic:", table(data.use['HES6',][data.use['HES6',] != 0] > 0)))


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
  print(paste("HES6-p_val:", p_val))
  print(paste("HES6-nz_p_val:", nz_p_val))

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

# play system sounds, call function to alarm user than running is done!
beep <- function(n = 5){
  for(i in seq(n)){
    system("rundll32 user32.dll, MessageBeep -1")
    Sys.sleep(.5)
  }
}

# set idents to compare cells at sample level instead of cluster level
Idents(integrated) <- integrated$orig.ident
# load future library and set plan to run certain functions with multiprocessing
plan("multisession", workers = 1) # n_workers > 1 for parallelization (voor mij, 5 is max, 4 is safe)
# get sample vs sample markers (now: monoculture vs coculture (for neurons and astrocytes))
## note: p_val_adj = Adjusted p-value, based on Bonferroni correction using all genes (including non-zero expression) in the dataset
### adjusted both defaults: logfc.threshold = 0.25, min.pct = 0.1    to 0
sample_markers <- FindMarkers(integrated, ident.1 = names(table(integrated$orig.ident))[1], only.pos = FALSE, verbose = T,
                              logfc.threshold = 0, min.pct = 0)



## where to find Seurat source code
# https://github.com/satijalab/seurat/blob/master/R/differential_expression.R
# https://github.com/satijalab/seurat/blob/master/R/utilities.R
