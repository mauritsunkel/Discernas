# devtools::install_github("mauritsunkel/EMC-SKlab-scRNAseq")
library(EMC.SKlab.scRNAseq)

#####  USER INITIALIZATION #####
### USER CONFIG ###
# general
run_name <- "Bas_pipe_V8" # may be empty string: ""
project_dir <- "C:/SynologyDrive/Projects/scRNAseqR"
samples_dir <- file.path(project_dir, "data/samples/")

# sample analysis
sample_names <- c("A", "C", "N")
if (any(duplicated(sample_names))) stop("Make sure no sample_names are duplicated")
# sample integration
sample_integrations <- list(
  c("A", "C", "N")
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
  "microglia_absence" = c("CD163", "CCL2", "MRC1"),
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
  file.path(project_dir, "results", run_name, integrated_sample_name, paste0(integrated_sample_name, ".rds"))
})))
##### END USER INITIALIZATION #####



## INDIVIDUAL SAMPLE ANALYSIS ----
message("\n RUNNING sample analysis \n")
for (sample_name in sample_names) {
  sample_analysis(
    samples_dir = samples_dir,
    sample_name = sample_name,
    output_dir = results_dir,
    features_of_interest = features_of_interest,
    run_cell_cycle_regression = FALSE
  )
}

## SAMPLES INTEGRATION ----
message("\n RUNNING samples integration \n")
samples_integration(
  sample_files = c(
    file.path(results_dir, sample_integrations[[1]][1], paste0(sample_integrations[[1]][1], ".rds")),
    file.path(results_dir, sample_integrations[[1]][2], paste0(sample_integrations[[1]][2], ".rds")),
    file.path(results_dir, sample_integrations[[1]][3], paste0(sample_integrations[[1]][3], ".rds"))
  ),
  sample_names = sample_integrations[[1]],
  output_dir = results_dir,
  features_of_interest = features_of_interest,
  integration_method = integration_method
)



## ANNOTATE AND VISUALIZE (KRIEGSTEIN) ----
### MANUALLY PERFORM MAPMYCELLS SILETTI ANNOTATION ###
message("RUNNING annotate_visualize_with_kriegstein_data")
annotate_visualize_with_kriegstein_data(
  sample_names = integrated_sample_names,
  sample_files = integrated_sample_files,
  output_dir = results_dir,
  kriegstein_data_dir = kriegstein_data_dir,
  kriegstein_chunks_input_dir = kriegstein_chunks_output_dir,
  kriegstein_annotated_output_dir = file.path(kriegstein_data_dir, "RData", run_name),
  run_only_visualization = FALSE # DEVNOTE: check if TRUE, only when testing
)

## Subset selection neurons ----
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "RPCA",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons_RPCA"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("MAP2", "DCX", "NEUROG2"), percent_expressed = 30, reference_annotations = NULL)
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "CCA",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons_CCA"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("MAP2", "DCX", "NEUROG2"), percent_expressed = 30, reference_annotations = NULL)
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "harmony",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "neurons_harmony"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("MAP2", "DCX", "NEUROG2"), percent_expressed = 30, reference_annotations = NULL)
## Subset selection astrocytes ----
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "RPCA",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes_RPCA"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("VIM", "S100B", "SOX9"), percent_expressed = 30, reference_annotations = NULL)
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "CCA",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes_CCA"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("VIM", "S100B", "SOX9"), percent_expressed = 30, reference_annotations = NULL)
selection_reintegration(
  so_filename = integrated_sample_files[[1]],
  integration_method = "harmony",
  output_dir = file.path(results_dir, integrated_sample_names[[1]], "subset", "astrocytes_harmony"),
  sample_name = integrated_sample_names[[1]],
  features_of_interest = features_of_interest,
  selection_markers = c("VIM", "S100B", "SOX9"), percent_expressed = 30, reference_annotations = NULL)



# TODO test percent expressed in histograms
# TODO check if /sample_name/ in output_dir at every R script file: if (!grepl(paste0("/", sample_name, "/"), output_dir))

## PSEUDOTIME ----
message("RUNNING pseudotime")
pseudotime(
  input_files = integrated_sample_files,
  input_names = integrated_sample_names,
  output_dir = results_dir,
  pseudotime_root_markers = pseudotime_root_markers,
  single_partition = FALSE
)

















# annotate integrated samples
## chunk Kriegstein reference data
# chunk_kriegstein_data(
#   n_chunks = 25,
#   kriegstein_data_dir = kriegstein_data_dir,
#   kriegstein_chunks_output_dir = kriegstein_chunks_output_dir,
#   kriegstein_custom_annotation = system.file("extdata", "kriegstein_custom_annotation.txt", package = 'EMC.SKlab.scRNAseq')
# )







# ## annotate and visualize with kriegstein reference data
# annotate_visualize_with_kriegstein_data(
#   sample_names = integrated_sample_names,
#   sample_files = integrated_sample_files,
#   output_dir = results_dir,
#   kriegstein_data_dir = kriegstein_data_dir,
#   kriegstein_chunks_input_dir = kriegstein_chunks_output_dir,
#   kriegstein_annotated_output_dir = file.path(kriegstein_data_dir, "RData", run_name)
# )
#


## pseudotime
# message("RUNNING pseudotime")
# pseudotime(
#   input_files = integrated_sample_files[1],
#   input_names = integrated_sample_names[1],
#   output_dir = results_dir,
#   pseudotime_root_markers = pseudotime_root_markers
# )
# pseudotime(
#   input_files = integrated_sample_files[2],
#   input_names = integrated_sample_names[2],
#   output_dir = results_dir,
#   genes_of_interest = pseudotime_root_markers
# )











# # differential expression analysis
# sample_info <- list(integrated_sample_names, integrated_sample_files)
# for (i in seq_along(integrated_sample_names)) {
#   differential_expression_analysis(
#     sample_name = sample_info[[1]][i],
#     rds_file = sample_info[[2]][i],
#     output_dir = dirname(sample_info[[2]][i])
#   )
# }
#

#
# # FGSEA: functional gene set enrichment- and over-representation analysis
# for (path in integrated_sample_files) {
#   run_fgsea(
#     output_dir = dirname(path),
#     dea_result_file = file.path(dirname(path), 'DEA', 'sample_markers', 'dea.csv'),
#     cellRanger_ensembl_features = system.file("extdata", "ensembl_genes.tsv", package = 'EMC.SKlab.scRNAseq')
#   )
# }

