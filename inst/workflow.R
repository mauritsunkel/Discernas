# devtools::install_github("mauritsunkel/EMC-SKlab-scRNAseq")

library(EMC.SKlab.scRNAseq)

#####  USER INITIALIZATION #####
### USER CONFIG ###
# general
run_name <- "package_workflow_test_8" # may be empty string: ""
project_dir <- "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq"
samples_dir <- "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/data/samples/bas cultures"

features_of_interest <- list(
  "astrocyte_interest" = c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                          "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE"),
  "astrocyte_maturity" = c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                          "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE"),
  "neuron_interest" = c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1"),
  "neuron_maturity" = c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                       "SYN1", "SNAP25", "SYT1", "APOE"),
  "schema_psych_interest" = c("SETD1A", "CUL1", "XPO7", "TRIO", "CACNA1G", "SP4",
                             "GRIA3", "GRIN2A", "HERC1", "RB1CC1", "HCN4", "AKAP11"),
  "sloan_2017_interest" = c("AQP4", "ALDH1L1", "RANBP3L", "IGFBP7", "TOP2A", "TMSB15A", "NNAT", "HIST1H3B",
                           "STMN2", "SYT1", "SNAP25", "SOX9", "CLU", "SLC1A3", "UBE2C", "NUSAP1", "PTPRZ1",
                           "HOPX", "FAM107A", "AGT"),
  "interneuron_interest" = c("SST", "PVALB", "GAD1")
)

# sample analysis
sample_names <- c("BL_A", "BL_C", "BL_N")

# sample integration
run_selection <- TRUE
sample_integrations <- list(
  c("BL_A", "BL_C"),
  c("BL_N", "BL_C")
)
marker_selection_panels <- list(
  c("VIM", "S100B", "SOX9"), # SOX9 <-> FABP7
  c("MAP2", "DCX", "NEUROG2") # RBFOX3 <-> DCX
)
integration_reference_sample <- "BL_C"

# annotation
## chunk Kriegstein reference
kriegstein_data_dir <- "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/data/Kriegstein"
kriegstein_chunks_output_dir <- file.path(kriegstein_data_dir, "RData", "chunks_25")
### END USER CONFIG ###

# create results/run_start_time/ directories from project directories
if (run_name == "") {
  run_name <- format(Sys.time(), "%F_%H-%M-%S")
}
results_dir <- file.path(project_dir, "results", run_name)

integrated_sample_names <- unlist(sapply(sample_integrations, simplify = F, function(integrated_sample_name) {
  paste(integrated_sample_name, collapse = "-")
}))
integrated_sample_files <- unname(unlist(sapply(integrated_sample_names, simplify = F, function(integrated_sample_name) {
  if (run_selection) {
    file.path(project_dir, "results", run_name, "integrated", integrated_sample_name, "postSelect", paste0(integrated_sample_name, ".rds"))
  } else {
    file.path(project_dir, "results", run_name, "integrated", integrated_sample_name, paste0(integrated_sample_name, ".rds"))
  }
})))
##### END USER INITIALIZATION #####





# initial sample analysis
for (sample_name in sample_names) {
  sample_analysis(
    samples_dir = samples_dir,
    sample_name = sample_name,
    output_dir = results_dir,
    features_of_interest = features_of_interest,
    run_cell_cycle_regression = FALSE
  )
}

# integrate samples
## BL_A + BL_C
samples_integration(
  sample_files = c(
    file.path(results_dir, sample_integrations[[1]][1], paste0(sample_integrations[[1]][1], ".rds")),
    file.path(results_dir, sample_integrations[[1]][2], paste0(sample_integrations[[1]][2], ".rds"))
  ),
  sample_names = sample_integrations[[1]],
  output_dir = results_dir,
  selection_panel = marker_selection_panels[[1]],
  features_of_interest = features_of_interest
)
## BL_C + BL_N
samples_integration(
  sample_files = c(
    file.path(results_dir, sample_integrations[[2]][1], paste0(sample_integrations[[2]][1], ".rds")),
    file.path(results_dir, sample_integrations[[2]][2], paste0(sample_integrations[[2]][2], ".rds"))
  ),
  sample_names = sample_integrations[[2]],
  output_dir = results_dir,
  selection_panel = marker_selection_panels[[2]],
  features_of_interest = features_of_interest
)

# annotate integrated samples
## chunk Kriegstein reference data
chunk_kriegstein_data(
  n_chunks = 25,
  kriegstein_data_dir = kriegstein_data_dir,
  kriegstein_chunks_output_dir = kriegstein_chunks_output_dir,
  kriegstein_custom_annotation = system.file("extdata", "kriegstein_custom_annotation.txt", package = 'EMC.SKlab.scRNAseq')
)

## annotate with kriegstein reference data
annotate_with_kriegstein_data(
  sample_names = integrated_sample_names,
  sample_files = integrated_sample_files,
  kriegstein_data_dir = kriegstein_data_dir,
  kriegstein_chunks_input_dir = kriegstein_chunks_output_dir,
  kriegstein_annotated_output_dir = file.path(kriegstein_data_dir, "RData", "workflow_annotated_Rdata")
)

## visualize Kriegstein reference annotated data
visualize_kriegstein_annotated_data(
  sample_names = integrated_sample_names,
  sample_files = integrated_sample_files,
  kriegstein_data_dir = kriegstein_data_dir,
  kriegstein_annotated_input_dir = file.path(kriegstein_data_dir, "RData", "workflow_annotated_Rdata"),
  output_dir = results_dir
)

# differential expression analysis
sample_info <- list(integrated_sample_names, integrated_sample_files)
for (i in seq_along(sample_info)) {
  differential_expression_analysis(
    sample_name = sample_info[[1]][i],
    rds_file = sample_info[[2]][i],
    output_dir = dirname(sample_info[[2]][i])
  )
}

DE2volcanoPlot(
  DE_csv = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_A-BL_C/postSelect/DEA/sample_markers/dea.csv",
  output_filename = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_A-BL_C/postSelect/DEA/sample_markers/dea",
  log2FC_threshold = 1
)
DE2volcanoPlot(
  DE_csv = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_N-BL_C/postSelect/DEA/sample_markers/dea.csv",
  output_filename = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_N-BL_C/postSelect/DEA/sample_markers/dea",
  log2FC_threshold = 1
)

# FGSEA: functional gene set enrichment- and over-representation analysis
for (path in integrated_sample_files) {
  run_fgsea(
    output_dir = dirname(path),
    dea_result_file = file.path(dirname(path), 'DEA', 'sample_markers', 'dea.csv'),
    cellRanger_ensembl_features = system.file("extdata", "ensembl_genes.tsv", package = 'EMC.SKlab.scRNAseq')
  )
}

# pseudotime
## astrocytes
pseudotime(
  input_files = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_A-BL_C/postSelect/BL_A-BL_C.rds",
  input_names = "BL_A-BL_C",
  output_dir = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_A-BL_C/postSelect/",
  genes_of_interest = c(marker_selection_panels[[1]])
)
# neurons
pseudotime(
  input_files = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_N-BL_C/postSelect/BL_N-BL_C.rds",
  input_names = "BL_N-BL_C",
  output_dir = "C:/Users/mauri/Desktop/EMC-SKlab-scRNAseq/results/package_workflow_test_8/integrated/BL_N-BL_C/postSelect/",
  genes_of_interest = c(marker_selection_panels[[1]])
)
