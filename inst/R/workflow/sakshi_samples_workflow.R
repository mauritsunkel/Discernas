# devtools::install_github("kushnerlab/scRNAseqR")
library(EMC.SKlab.scRNAseq)



#####  USER INITIALIZATION -----
### USER CONFIG ###
# general
run_name <- "sakshi_pipeV8"
project_dir <- "C:/SynologyDrive/Projects/scRNAseqR"
samples_dir <- file.path(project_dir, "data/samples/")

# sample analysis
sample_names <- c("NS", "M", "NC", "NSM")
if (any(duplicated(sample_names))) stop("Make sure no sample_names are duplicated")
# sample integration
sample_integrations <- list(
  c("NSM", "NS", "NC", "M")
)

integration_method <- "RPCA"

features_of_interest <- list(
  "astrocyte" = c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                           "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE"),
  "astrocyte_maturity" = c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                           "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE"),
  "neuron" = c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1"),
  "neuron_maturity" = c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                        "SYN1", "SNAP25", "SYT1", "APOE"),
  "schema_psych" = c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                              "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11"),
  "sloan_2017" = c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                            "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                            "HOPX", "FAM107A", "AGT"),
  "interneuron" = c("SST", "PVALB", "GAD1"),
  "microglia" = c("IBA1", "TMEM119", "P2RY12", "CXCR1", "ITGAM", "PTPRC", "SALL1", "TREM2", "SPI1", "CSF1", "CSF1R",
                  "AIF1", "C1QA", "C1QB", "CX3CR1", "TGFB1", "CX3CL1", "CSF2"),
  "microglia_absence" = c("CD163", "CCL2", "MRC1", "IL34"),
  "proliferating" = c("MKI67", "SOX2", "HOPX", "NES", "POU5F1"),
  "supplement2" = c("TGFB1", "CSF1", "IL34", "LEFTY2")
)
pseudotime_root_markers <- list(
  "Microglia" = c("AIF1", "CSF1R", "SPI1"),
  "Astrocyte" = c("VIM", "S100B", "SOX9"),
  "Neuron"    = c("MAP2", "DCX", "NEUROG2"),
  "Dividing"  = c("MKI67"),
  "other"     = c("FOXJ1") # based on Siletti et al ..., gene used for both choroid plexus and ependymal cells, seemingly early brain cells
)

kriegstein_data_dir <- file.path(project_dir, "data/Kriegstein")
kriegstein_chunks_output_dir <- file.path(kriegstein_data_dir, "RData", "chunks_25")

file.copy(EMC.SKlab.scRNAseq::thisFilePath(), results_dir, overwrite = TRUE)
### END USER CONFIG ###



# create results/run_start_time/ directories from project directories
if (run_name == "") run_name <- format(Sys.time(), "%F_%H-%M-%S")
results_dir <- file.path(project_dir, "results", run_name)

integrated_sample_names <- unlist(sapply(sample_integrations, simplify = F, function(integrated_sample_name) {
  paste(integrated_sample_name, collapse = "-")
}))
integrated_sample_files <- unname(unlist(sapply(integrated_sample_names, simplify = F, function(integrated_sample_name) {
  file.path(project_dir, "results", run_name, integrated_sample_name, paste0(integrated_sample_name, ".qs"))
})))
#### END USER INITIALIZATION #####



## SAMPLE ANALYSIS ----
# for (sample_name in sample_names) {
#   message("RUNNING sample_analysis sample: ", sample_name)
#   sample_analysis(
#     samples_dir = samples_dir,
#     sample_name = sample_name,
#     output_dir = results_dir,
#     features_of_interest = features_of_interest,
#     run_cell_cycle_regression = FALSE
#   )
# }
#

## SAMPLES INTEGRATION ----
# message("RUNNING samples_integration")
# samples_integration(
#   sample_files = c(
#     file.path(results_dir, sample_integrations[[1]][1], paste0(sample_integrations[[1]][1], ".qs")),
#     file.path(results_dir, sample_integrations[[1]][2], paste0(sample_integrations[[1]][2], ".qs")),
#     file.path(results_dir, sample_integrations[[1]][3], paste0(sample_integrations[[1]][3], ".qs")),
#     file.path(results_dir, sample_integrations[[1]][4], paste0(sample_integrations[[1]][4], ".qs"))
#   ),
#   sample_names = sample_integrations[[1]],
#   output_dir = results_dir,
#   features_of_interest = features_of_interest,
#   integration_method = integration_method
# )
#

## ANNOTATE AND VISUALIZE (KRIEGSTEIN) ----
# message("RUNNING annotate_with_kriegstein_data")
# annotate_visualize_with_kriegstein_data(
#   sample_names = integrated_sample_names,
#   sample_files = integrated_sample_files,
#   output_dir = results_dir,
#   kriegstein_data_dir = kriegstein_data_dir,
#   kriegstein_chunks_input_dir = kriegstein_chunks_output_dir,
#   kriegstein_annotated_output_dir = file.path(kriegstein_data_dir, "RData", run_name),
#   run_only_visualization = FALSE # DEVNOTE: check if TRUE, only when testing
# )
#
## PSEUDOTIME ----
# message("RUNNING pseudotime")
# pseudotime(
#   input_files = integrated_sample_files,
#   input_names = integrated_sample_names,
#   output_dir = results_dir,
#   pseudotime_root_markers = pseudotime_root_markers
# )


#### MANUAL selections ----
## Subset selection microglia ----
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "harmony",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "microglia"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  exclude_samples = c('NS', 'NC'),
  selection_markers = c("AIF1", "CSF1R", "SPI1"), percent_expressed = 30, reference_annotations = NULL)
## Subset selection astrocytes ----
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "harmony",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  exclude_samples = c('M'),
  selection_markers = c("VIM", "S100B", "SOX9"), percent_expressed = 30, reference_annotations = NULL)
## Subset selection neurons ----
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "harmony",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  exclude_samples = c('M'),
  selection_markers = c("MAP2", "DCX", "NEUROG2"), percent_expressed = 30, reference_annotations = NULL)

## selections PSEUDOTIME ----
message("RUNNING pseudotime selections")
pseudotime(
  input_files = file.path(results_dir, integrated_sample_names[[1]], "subset", "microglia", basename(integrated_sample_files[[1]])),
  input_names = integrated_sample_names[[1]],
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "microglia"),
  pseudotime_root_markers = pseudotime_root_markers,
  single_partition = TRUE
)
pseudotime(
  input_files = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes", basename(integrated_sample_files[[1]])),
  input_names = integrated_sample_names[[1]],
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes"),
  pseudotime_root_markers = pseudotime_root_markers,
  single_partition = TRUE
)
pseudotime(
  input_files = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons", basename(integrated_sample_files[[1]])),
  input_names = integrated_sample_names[[1]],
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons"),
  pseudotime_root_markers = pseudotime_root_markers,
  single_partition = TRUE
)







# TODO check if /sample_name/ in output_dir at every R script file: if (!grepl(paste0("/", sample_name, "/"), output_dir))
## if (!grepl(paste0("/", sample_name, "/"), output_dir)) output_dir <- file.path(output_dir, sample_name)




# sample_celltype_options <- unique(paste(integrated$orig.ident, integrated@meta.data[, metadata_celltype_annotation], sep = "_"))
#
#
# ### MANUAL USER DEA ###
# ## to provide sample(s)-sample(s) or sample(s)-celltype(s) DEA comparisons in sample_celltype_DEA list, follow and adjust example below
# ## first get seurat_object to perform DEA on
# integrated <- qs::qread(list(integrated_sample_names, integrated_sample_files)[[2]][1])
# ## then check metadata column names to see which can be used for DE comparisons, look specifically for annotation column names: MapMyCells & Kriegstein
# colnames(integrated@meta.data)
# ## MapMyCells default: mapmycells_supercluster
# ## Kriegstein default: kriegstein.seurat.custom.clusters.mean
# sample_options <- unique(integrated$orig.ident)
# sample_celltype_mmc_options <- unique(paste(integrated$orig.ident, integrated@meta.data[, "mapmycells_supercluster"], sep = "_"))
# sample_celltype_kriegstein.seurat_options <- unique(paste(integrated$orig.ident, integrated@meta.data[, "kriegstein.seurat.custom.clusters.mean"], sep = "_"))
# ## see above options for MapMyCells annotated single-cells or Kriegstein.Seurat annotated clusters for DE comparisons
# ## setup sample_celltype_DEA list, with named list by metadata column and then (named) list(s) of comparisons, follow and adjust the sample_celltype format below
# sample_celltype_DEA <- list(
#   orig.ident = list(
#     ## comparing whole samples
#     name = list(ref = "NC", vs = "NS")
#   ),
#   mapmycells_supercluster = list(
#     ## comparing microglia/astrocyte/neuron pairs
#     name = list(ref = "NSM_Astrocyte", vs = "NS_Astrocyte"),
#     name = list(ref = "NSM_Neuron", vs = "NS_Neuron"),
#     name = list(ref = "NSM_Astrocyte", vs = "NC_Astrocyte"),
#     name = list(ref = "NSM_Neuron", vs = "NC_Neuron"),
#     name = list(ref = "NC_Astrocyte", vs = "NS_Astrocyte"),
#     name = list(ref = "NC_Neuron", vs = "NS_Neuron"),
#     name = list(ref = "NSM_Microglia", vs = "M_Microglia")
#   ),
#   kriegstein.seurat.custom.clusters.mean = list(
#     ## comparing microglia pairs
#     name = list(ref = "M_Microglia.4", vs = "M_Microglia.5"),
#     name = list(ref = "M_Microglia.4", vs = "M_Microglia.6"),
#     name = list(ref = "M_Microglia.4", vs = "M_Microglia.9"),
#     name = list(ref = "M_Microglia.5", vs = "M_Microglia.6"),
#     name = list(ref = "M_Microglia.5", vs = "M_Microglia.9"),
#     name = list(ref = "M_Microglia.6", vs = "M_Microglia.9"),
#     name = list(ref = "NSM_Microglia.4", vs = "NSM_Microglia.5"),
#     name = list(ref = "NSM_Microglia.4", vs = "NSM_Microglia.6"),
#     name = list(ref = "NSM_Microglia.4", vs = "NSM_Microglia.9"),
#     name = list(ref = "NSM_Microglia.5", vs = "NSM_Microglia.6"),
#     name = list(ref = "NSM_Microglia.5", vs = "NSM_Microglia.9"),
#     name = list(ref = "NSM_Microglia.6", vs = "NSM_Microglia.9"),
#     name = list(ref = "NSM_Microglia.4", vs = "M_Microglia.4"),
#     name = list(ref = "NSM_Microglia.5", vs = "M_Microglia.5"),
#     name = list(ref = "NSM_Microglia.6", vs = "M_Microglia.6"),
#     name = list(ref = "NSM_Microglia.9", vs = "M_Microglia.9"),
#     ## comparing microglia multiples
#     NSM_Microglia.4.5.6_vs_M_Microglia.4.5.6 = list(ref = c("NSM_Microglia.4", "NSM_Microglia.5", "NSM_Microglia.6"), vs = c("M_Microglia.4", "M_Microglia.5", "M_Microglia.6")),
#     NSM_Microglia.4.5_vs_M_Microglia.4.5 = list(ref = c("NSM_Microglia.4", "NSM_Microglia.5"), vs = c("M_Microglia.4", "M_Microglia.5")),
#
#     ## for "progenitor" cells, differently labeled by MapMyCells and Kriegstein, manually inspect
#     M_Microglia.9_vs_NSM_NS_NC_Microglia.9 = list(ref = "M_Microglia.9", vs = c("NSM_Microglia.9", "NS_Microglia.9", "NC_Microglia.9")),
#     M_Dividing.8_vs_NSM_NS_NC_Dividing.8 = list(ref = "M_Dividing.8", vs = c("NSM_Dividing.8", "NS_Dividing.8", "NC_Dividing.8")),
#     NSM_Progenitors.0.7.8.9.14.17_vs_NS_Progenitors.0.7.8.9.14.17 = list(ref = c("NSM_Dividing.0", "NSM_IPC.7", "NSM_Dividing.8", "NSM_Microglia.9", "NSM_Endo 1.14", "NSM_Dividing.17"), vs = c("NS_Dividing.0", "NS_IPC.7", "NS_Dividing.8", "NS_Microglia.9", "NS_Endo 1.14", "NS_Dividing.17")),
#     NSM_Progenitors.0.7.8.9.14.17_vs_NC_Progenitors.0.7.8.9.14.17 = list(ref = c("NSM_Dividing.0", "NSM_IPC.7", "NSM_Dividing.8", "NSM_Microglia.9", "NSM_Endo 1.14", "NSM_Dividing.17"), vs = c("NC_Dividing.0", "NC_IPC.7", "NC_Dividing.8", "NC_Microglia.9", "NC_Endo 1.14", "NC_Dividing.17")),
#     NS_Progenitors.0.7.8.9.14.17_vs_NC_Progenitors.0.7.8.9.14.17 = list(ref = c("NS_Dividing.0", "NS_IPC.7", "NS_Dividing.8", "NS_Microglia.9", "NS_Endo 1.14", "NS_Dividing.17"), vs = c("NC_Dividing.0", "NC_IPC.7", "NC_Dividing.8", "NC_Microglia.9", "NC_Endo 1.14", "NC_Dividing.17")),
#     ## only positive markers
#     M_Microglia.4_vs_M_Microglia.5.6.9 = list(ref = c("M_Microglia.4"), vs = c("M_Microglia.5", "M_Microglia.6", "M_Microglia.9")),
#     M_Microglia.5_vs_M_Microglia.4.6.9 = list(ref = c("M_Microglia.5"), vs = c("M_Microglia.4", "M_Microglia.6", "M_Microglia.9")),
#     M_Microglia.6_vs_M_Microglia.4.5.9 = list(ref = c("M_Microglia.6"), vs = c("M_Microglia.4", "M_Microglia.5", "M_Microglia.9")),
#     M_Microglia.9_vs_M_Microglia.4.5.6 = list(ref = c("M_Microglia.9"), vs = c("M_Microglia.4", "M_Microglia.5", "M_Microglia.6")),
#     NSM_Microglia.4_vs_NSM_Microglia.5.6 = list(ref = c("NSM_Microglia.4"), vs = c("NSM_Microglia.5", "NSM_Microglia.6")),
#     NSM_Microglia.5_vs_NSM_Microglia.4.6 = list(ref = c("NSM_Microglia.5"), vs = c("NSM_Microglia.4", "NSM_Microglia.6")),
#     NSM_Microglia.6_vs_NSM_Microglia.4.5 = list(ref = c("NSM_Microglia.6"), vs = c("NSM_Microglia.4", "NSM_Microglia.5"))
#   )
# )
# ### END MANUAL USER DEA
#
#
#
# # differential expression analysis
# sample_info <- list(integrated_sample_names, integrated_sample_files)
# for (i in seq_along(integrated_sample_names)) {
#   differential_expression_analysis(
#     sample_name = sample_info[[1]][i],
#     qs_file = sample_info[[2]][i],
#     output_dir = dirname(sample_info[[2]][i]),
#     sample_celltype_DEA = sample_celltype_DEA,
#     features_of_interest = features_of_interest
#   )
# }
