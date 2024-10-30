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
### END USER CONFIG ###



# create results/run_start_time/ directories from project directories
if (run_name == "") run_name <- format(Sys.time(), "%F_%H-%M-%S")
results_dir <- file.path(project_dir, "results", run_name)

integrated_sample_names <- unlist(sapply(sample_integrations, simplify = F, function(integrated_sample_name) {
  paste(integrated_sample_name, collapse = "-")
}))
integrated_sample_files <- unname(unlist(sapply(integrated_sample_names, simplify = F, function(integrated_sample_name) {
  file.path(results_dir, integrated_sample_name, paste0(integrated_sample_name, ".qs"))
})))
#### END USER INITIALIZATION #####



## START ANALYSIS ----

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
# ## Subset selection microglia ----
# selection_reintegration(
#   so_filename = integrated_sample_files[[1]],
#   integration_method = "harmony",
#   output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "microglia"),
#   sample_name = integrated_sample_names[[1]],
#   features_of_interest = features_of_interest,
#   exclude_samples = c('NS', 'NC'),
#   selection_markers = c("AIF1", "CSF1R", "SPI1"), percent_expressed = 30, reference_annotations = NULL)
# ## Subset selection astrocytes ----
# selection_reintegration(
#   so_filename = integrated_sample_files[[1]],
#   integration_method = "harmony",
#   output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes"),
#   sample_name = integrated_sample_names[[1]],
#   features_of_interest = features_of_interest,
#   exclude_samples = c('M'),
#   selection_markers = c("VIM", "S100B", "SOX9"), percent_expressed = 30, reference_annotations = NULL)
# ## Subset selection neurons ----
# selection_reintegration(
#   so_filename = integrated_sample_files[[1]],
#   integration_method = "harmony",
#   output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons"),
#   sample_name = integrated_sample_names[[1]],
#   features_of_interest = features_of_interest,
#   exclude_samples = c('M'),
#   selection_markers = c("MAP2", "DCX", "NEUROG2"), percent_expressed = 30, reference_annotations = NULL)

## Selections KRIEGSTEIN ANNOTATION ----
# message("RUNNING subset annotate_with_kriegstein_data")
# subset_names <- c("microglia", "astrocytes", "neurons")
# for (subset_name in subset_names) {
#   annotate_visualize_with_kriegstein_data(
#     sample_names = integrated_sample_names[[1]],
#     sample_files = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name, basename(integrated_sample_files[[1]])),
#     output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name),
#     kriegstein_data_dir = kriegstein_data_dir,
#     kriegstein_chunks_input_dir = kriegstein_chunks_output_dir,
#     kriegstein_annotated_output_dir = file.path(kriegstein_data_dir, "RData", run_name, integrated_sample_names[[1]], "subset", subset_name),
#     run_only_visualization = FALSE # DEVNOTE: check if TRUE, only when testing
#   )
# }


## selections PSEUDOTIME ----
# message("RUNNING subset pseudotime selections")
# subset_names <- c("microglia", "astrocytes", "neurons")
# for (subset_name in subset_names) {
#   pseudotime(
#     input_files = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name, basename(integrated_sample_files[[1]])),
#     input_names = integrated_sample_names[[1]],
#     output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name),
#     pseudotime_root_markers = pseudotime_root_markers,
#     single_partition = TRUE
#   )
# }










#### MANUAL DEA ----
## to provide sample(s)-sample(s) or sample(s)-celltype(s) DEA comparisons in sample_celltype_DEA list, follow and adjust example below
## first get seurat_object to perform DEA on
## then check metadata column names to see which can be used for DE comparisons, look specifically for annotation column names: MapMyCells & Kriegstein
# library(Seurat)
# DEA_qs <- file.path(results_dir, integrated_sample_names[[1]], "subset", "microglia", basename(integrated_sample_files[[1]]))
# integrated <- qs::qread(DEA_qs)
#
# colnames(integrated@meta.data) # MapMyCells default: mapmycells_supercluster - Kriegstein default: kriegstein.seurat.custom.clusters.mean
# sample_options <- sort(table(factor(integrated$orig.ident)), decreasing = T)
# sample_celltype_mmc_options <- sort(table(paste(integrated$orig.ident, factor(integrated@meta.data[, "mapmycells_supercluster"]), sep = "_")), decreasing = T)
# sample_celltype_kriegstein.seurat_options <- sort(table(paste(integrated$orig.ident, factor(integrated@meta.data[, "kriegstein.seurat.custom.clusters.mean"]), sep = "_")), decreasing = T)
## sample_cluster_options <- sort(table(paste(integrated$orig.ident, integrated$seurat_clusters, sep = "_")), decreasing = T)
## see above options for MapMyCells annotated single-cells or Kriegstein.Seurat annotated clusters for DE comparisons
## look at the (subset) UMAPs for sample original identity, MapMyCells and Kriegstein to see which cells/clusters to compare
## setup sample_celltype_DEA list, with custom named list by metadata and list(s) of comparisons, recreate and add to the sample_celltype format below
###
### DEA SETUP
###
sample_celltype_DEA <- list(
  'microglia' = list(
    orig.ident = list(
      ## comparing whole samples
      name = list(ref = "M", vs = "NSM")
    ),
    mapmycells_supercluster = list(
      ## comparing microglia pairs
      name = list(ref = "M_Microglia", vs = "NSM_Microglia")
    ),
    kriegstein.seurat.custom.clusters.mean = list(
      ## comparing microglia pairs between M and NSM
      name = list(ref = "M_Microglia.0", vs = "NSM_Microglia.0"),
      name = list(ref = "M_Microglia.1", vs = "NSM_Microglia.1"),
      name = list(ref = "M_Microglia.2", vs = "NSM_Microglia.2"),
      name = list(ref = "M_Microglia.3", vs = "NSM_Microglia.3"),
      name = list(ref = "M_Microglia.4", vs = "NSM_Microglia.4"),
      name = list(ref = "M_Microglia.5", vs = "NSM_Microglia.5"),
      name = list(ref = "M_Microglia.6", vs = "NSM_Microglia.6"),
      name = list(ref = "M_Microglia.7", vs = "NSM_Microglia.7"),
      name = list(ref = "M_Microglia.8", vs = "NSM_Microglia.8"),
      ## comparing vs all other cells, with focus on positive markers
      "M_NSM_Microglia.0_vs_rest" = list(ref = c("M_Microglia.0", "NSM_Microglia.0"), vs = "rest"),
      "M_NSM_Microglia.1_vs_rest" = list(ref = c("M_Microglia.1", "NSM_Microglia.1"), vs = "rest"),
      "M_NSM_Microglia.2_vs_rest" = list(ref = c("M_Microglia.2", "NSM_Microglia.2"), vs = "rest"),
      "M_NSM_Microglia.3_vs_rest" = list(ref = c("M_Microglia.3", "NSM_Microglia.3"), vs = "rest"),
      "M_NSM_Microglia.4_vs_rest" = list(ref = c("M_Microglia.4", "NSM_Microglia.4"), vs = "rest"),
      "M_NSM_Microglia.5_vs_rest" = list(ref = c("M_Microglia.5", "NSM_Microglia.5"), vs = "rest"),
      "M_NSM_Microglia.6_vs_rest" = list(ref = c("M_Microglia.6", "NSM_Microglia.6"), vs = "rest"),
      "M_NSM_Microglia.7_vs_rest" = list(ref = c("M_Microglia.7", "NSM_Microglia.7"), vs = "rest"),
      "M_NSM_Microglia.8_vs_rest" = list(ref = c("M_Microglia.8", "NSM_Microglia.8"), vs = "rest")
    )
 ),
 'astrocytes' = list(
   orig.ident = list(
     ## comparing whole samples
     name = list(ref = "NC", vs = "NS"),
     name = list(ref = "NC", vs = "NSM"),
     name = list(ref = "NS", vs = "NSM"),
     "NS_NC_vs_NSM" = list(ref = c("NS", "NC"), vs = "NSM")
   ),
   mapmycells_supercluster = list(
     ## comparing astrocyte pairs
     name = list(ref = "NC_Astrocyte", vs = "NS_Astrocyte"),
     name = list(ref = "NC_Astrocyte", vs = "NSM_Astrocyte"),
     name = list(ref = "NS_Astrocyte", vs = "NSM_Astrocyte"),
     "NC_NS_Astrocyte_vs_NSM_Astrocyte" = list(ref = c("NC_Astrocyte", "NS_Astrocyte"), vs = "NSM_Astrocyte"),
     "Ependymal_vs_rest" = list(ref = c("NC_Ependymal", "NS_Ependymal", "NSM_Ependymal"), vs = "rest"),
     "Neuron_vs_rest" = list(ref = c("NC_Neuron", "NS_Neuron", "NSM_Neuron"), vs = "rest")
   ),
   kriegstein.seurat.custom.clusters.mean = list(
     ## comparing astrocyte pairs
     'NC_NS_Astrocyte 2.0_vs_NSM_Astrocyte 2.0' = list(ref = c("NC_Astrocyte 2.0", "NS_Astrocyte 2.0"), vs = "NSM_Astrocyte 2.0"),
     'NC_NS_Radial glia 2.1_vs_NSM_Radial glia 2.1' = list(ref = c("NC_Radial glia 2.1", "NS_Radial glia 2.1"), vs = "NSM_Radial glia 2.1"),
     'NC_NS_Astrocyte 2.2_vs_NSM_Astrocyte 2.2' = list(ref = c("NC_Astrocyte 2.2", "NS_Astrocyte 2.2"), vs = "NSM_Astrocyte 2.2"),
     'NC_NS_Astrocyte 2.3_vs_NSM_Astrocyte 2.3' = list(ref = c("NC_Astrocyte 2.3", "NS_Astrocyte 2.3"), vs = "NSM_Astrocyte 2.3"),
     'NC_NS_Radial glia 2.4_vs_NSM_Radial glia 2.4' = list(ref = c("NC_Radial glia 2.4", "NS_Radial glia 2.4"), vs = "NSM_Radial glia 2.4"),
     'NC_NS_Radial glia 2.5_vs_NSM_Radial glia 2.5' = list(ref = c("NC_Radial glia 2.5", "NS_Radial glia 2.5"), vs = "NSM_Radial glia 2.5"),
     'NC_NS_Radial glia 2.6_vs_NSM_Radial glia 2.6' = list(ref = c("NC_Radial glia 2.6", "NS_Radial glia 2.6"), vs = "NSM_Radial glia 2.6"),
     'NC_NS_Radial glia 2.7_vs_NSM_Radial glia 2.7' = list(ref = c("NC_Radial glia 2.7", "NS_Radial glia 2.7"), vs = "NSM_Radial glia 2.7"),
     ## comparing astrocytes & microglia
     'NC_NS_Astrocyte 2_vs_NSM_Astrocyte 2' = list(ref = c("NC_Astrocyte 2.0", "NS_Astrocyte 2.0", "NC_Astrocyte 2.2", "NS_Astrocyte 2.2", "NC_Astrocyte 2.3", "NS_Astrocyte 2.3"), vs = c("NSM_Astrocyte 2.0", "NSM_Astrocyte 2.2", "NSM_Astrocyte 2.3")),
     'NC_NS_Radial glia 2_vs_NSM_Radial glia 2' = list(ref = c("NC_Radial glia 2.1", "NS_Radial glia 2.1", "NC_Radial glia 2.4", "NS_Radial glia 2.4", "NC_Radial glia 2.5", "NS_Radial glia 2.5", "NC_Radial glia 2.6", "NS_Radial glia 2.6", "NC_Radial glia 2.7", "NS_Radial glia 2.7"), vs = c("NSM_Radial glia 2.1", "NSM_Radial glia 2.4", "NSM_Radial glia 2.5", "NSM_Radial glia 2.6", "NSM_Radial glia 2.7")),
     'Astrocyte 2_vs_Radial glia 2' = list(ref = c("NC_Astrocyte 2.0", "NS_Astrocyte 2.0", "NC_Astrocyte 2.2", "NS_Astrocyte 2.2", "NC_Astrocyte 2.3", "NS_Astrocyte 2.3", "NSM_Astrocyte 2.0", "NSM_Astrocyte 2.2", "NSM_Astrocyte 2.3"), vs = c("NC_Radial glia 2.1", "NS_Radial glia 2.1", "NC_Radial glia 2.4", "NS_Radial glia 2.4", "NC_Radial glia 2.5", "NS_Radial glia 2.5", "NC_Radial glia 2.6", "NS_Radial glia 2.6", "NC_Radial glia 2.7", "NS_Radial glia 2.7", "NSM_Radial glia 2.1", "NSM_Radial glia 2.4", "NSM_Radial glia 2.5", "NSM_Radial glia 2.6", "NSM_Radial glia 2.7")),
     ## comparing vs rest
     'NC_NS_NSM_Astrocyte 2.0_vs_rest' = list(ref = c("NC_Astrocyte 2.0", "NS_Astrocyte 2.0", "NSM_Astrocyte 2.0"), vs = "rest"),
     'NC_NS_NSM_Astrocyte 2.2_vs_rest' = list(ref = c("NC_Astrocyte 2.2", "NS_Astrocyte 2.2", "NSM_Astrocyte 2.2"), vs = "rest"),
     'NC_NS_NSM_Astrocyte 2.3_vs_rest' = list(ref = c("NC_Astrocyte 2.3", "NS_Astrocyte 2.3", "NSM_Astrocyte 2.3"), vs = "rest"),
     'NC_NS_NSM_Radial glia 2.1_vs_rest' = list(ref = c("NC_Radial glia 2.1", "NS_Radial glia 2.1", "NSM_Radial glia 2.1"), vs = "rest"),
     'NC_NS_NSM_Radial glia 2.4_vs_rest' = list(ref = c("NC_Radial glia 2.4", "NS_Radial glia 2.4", "NSM_Radial glia 2.4"), vs = "rest"),
     'NC_NS_NSM_Radial glia 2.5_vs_rest' = list(ref = c("NC_Radial glia 2.5", "NS_Radial glia 2.5", "NSM_Radial glia 2.5"), vs = "rest"),
     'NC_NS_NSM_Radial glia 2.6_vs_rest' = list(ref = c("NC_Radial glia 2.6", "NS_Radial glia 2.6", "NSM_Radial glia 2.6"), vs = "rest"),
     'NC_NS_NSM_Radial glia 2.7_vs_rest' = list(ref = c("NC_Radial glia 2.7", "NS_Radial glia 2.7", "NSM_Radial glia 2.7"), vs = "rest")
   )
 ),
 'neurons' = list(
   orig.ident = list(
     ## comparing whole samples
     name = list(ref = "NC", vs = "NS"),
     name = list(ref = "NC", vs = "NSM"),
     name = list(ref = "NS", vs = "NSM"),
     "NS_NC_vs_NSM" = list(ref = c("NS", "NC"), vs = "NSM")
   ),
   mapmycells_supercluster = list(
     ## comparing neuron pairs
     name = list(ref = "NC_Neuron", vs = "NS_Neuron"),
     name = list(ref = "NC_Neuron", vs = "NSM_Neuron"),
     name = list(ref = "NS_Neuron", vs = "NSM_Neuron"),
     "NC_NS_Neuron_vs_NSM_Neuron" = list(ref = c("NC_Neuron", "NS_Neuron"), vs = "NSM_Neuron")
   ),
   kriegstein.seurat.custom.clusters.mean = list(
     ## comparing neuron pairs (>100 or more cells per group)
     'NC_NS_Excitatory neuron 1.0_vs_NSM_Excitatory neuron 1.0' = list(ref = c("NC_Excitatory neuron 1.0", "NS_Excitatory neuron 1.0"), vs = c("NSM_Excitatory neuron 1.0")),
     'NC_NS_Excitatory neuron 4.1_vs_NSM_Excitatory neuron 4.1' = list(ref = c("NC_Excitatory neuron 4.1", "NS_Excitatory neuron 4.1"), vs = c("NSM_Excitatory neuron 4.1")),
     'NC_NS_Excitatory neuron 4.2_vs_NSM_Excitatory neuron 4.2' = list(ref = c("NC_Excitatory neuron 4.2", "NS_Excitatory neuron 4.2"), vs = c("NSM_Excitatory neuron 4.2")),
     'NC_NS_Excitatory neuron 1.3_vs_NSM_Excitatory neuron 1.3' = list(ref = c("NC_Excitatory neuron 1.3", "NS_Excitatory neuron 1.3"), vs = c("NSM_Excitatory neuron 1.3")),
     'NC_NS_Excitatory neuron 1.4_vs_NSM_Excitatory neuron 1.4' = list(ref = c("NC_Excitatory neuron 1.4", "NS_Excitatory neuron 1.4"), vs = c("NSM_Excitatory neuron 1.4")),
     'NC_NS_Excitatory neuron 4.5_vs_NSM_Excitatory neuron 4.5' = list(ref = c("NC_Excitatory neuron 4.5", "NS_Excitatory neuron 4.5"), vs = c("NSM_Excitatory neuron 4.5")),
     'NC_NS_Excitatory neuron 1.6_vs_NSM_Excitatory neuron 1.6' = list(ref = c("NC_Excitatory neuron 1.6", "NS_Excitatory neuron 1.6"), vs = c("NSM_Excitatory neuron 1.6")),
     'NC_NS_Excitatory neuron 1.7_vs_NSM_Excitatory neuron 1.7' = list(ref = c("NC_Excitatory neuron 1.7", "NS_Excitatory neuron 1.7"), vs = c("NSM_Excitatory neuron 1.7")),
     'NC_NS_Excitatory neuron 1.8_vs_NSM_Excitatory neuron 1.8' = list(ref = c("NC_Excitatory neuron 1.8", "NS_Excitatory neuron 1.8"), vs = c("NSM_Excitatory neuron 1.8")),
     'NC_NS_Excitatory neuron 4.9_vs_NSM_Excitatory neuron 4.9' = list(ref = c("NC_Excitatory neuron 4.9", "NS_Excitatory neuron 4.9"), vs = c("NSM_Excitatory neuron 4.9")),
     'NC_NS_Excitatory neuron 1.10_vs_NSM_Excitatory neuron 1.10' = list(ref = c("NC_Excitatory neuron 1.10", "NS_Excitatory neuron 1.10"), vs = c("NSM_Excitatory neuron 1.10")),
     'NC_NS_Excitatory neuron 1.12_vs_NSM_Excitatory neuron 1.12' = list(ref = c("NC_Excitatory neuron 1.12", "NS_Excitatory neuron 1.12"), vs = c("NSM_Excitatory neuron 1.12")),
     ## comparing En1 vs En4
     'Excitatory neuron 1_vs_Excitatory neuron 4' = list(ref = c('Excitatory neuron 1'), vs = c('Excitatory neuron 4'), regexp = TRUE),
     ## comparing vs rest
     'Excitatory neuron 1.0_vs_rest' = list(ref = c('Excitatory neuron 1.0'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.3_vs_rest' = list(ref = c('Excitatory neuron 1.3'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.4_vs_rest' = list(ref = c('Excitatory neuron 1.4'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.6_vs_rest' = list(ref = c('Excitatory neuron 1.6'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.7_vs_rest' = list(ref = c('Excitatory neuron 1.7'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.8_vs_rest' = list(ref = c('Excitatory neuron 1.8'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.10_vs_rest' = list(ref = c('Excitatory neuron 1.10'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 1.12_vs_rest' = list(ref = c('Excitatory neuron 1.12'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 4.1_vs_rest' = list(ref = c('Excitatory neuron 4.1'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 4.2_vs_rest' = list(ref = c('Excitatory neuron 4.2'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 4.5_vs_rest' = list(ref = c('Excitatory neuron 4.5'), vs = 'rest', regexp = TRUE),
     'Excitatory neuron 4.9_vs_rest' = list(ref = c('Excitatory neuron 4.9'), vs = 'rest', regexp = TRUE)
   )
 )
)
### END MANUAL USER DEA


## DEA EXECUTION ----
message("RUNNING subset DEA")
subset_names <- c("microglia", "astrocytes", "neurons")
for (subset_name in subset_names) {
  differential_expression_analysis(
    sample_name = integrated_sample_names[[1]],
    qs_file = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name, basename(integrated_sample_files[[1]])),
    output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", subset_name),
    sample_celltype_DEA = sample_celltype_DEA[[subset_name]],
    features_of_interest = NULL, # features_of_interest, # if NULL, don't plot features etc
    pct.both = 0.01,
    pct.either = 0.05,
    DE_test = 'wilcox'
  )
}

## end with adding this file as config file.R to results folder
file.copy(EMC.SKlab.scRNAseq::thisFilePath(), results_dir, overwrite = TRUE)
