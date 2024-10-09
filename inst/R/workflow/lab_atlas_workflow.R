# devtools::install_github("kushnerlab/scRNAseqR")
library(EMC.SKlab.scRNAseq)



#####  USER INITIALIZATION #####
### USER CONFIG ###
# general
run_name <- "atlas_V1"
project_dir <- "C:/SynologyDrive/Projects/scRNAseqR"
samples_dir <- file.path(project_dir, "data/samples/")

# sample analysis
sample_names <- c("NS", "M", "NC", "NSM", "A", "C", "N", "t90", "t149", "t275")
if (any(duplicated(sample_names))) stop("Make sure no sample_names are duplicated")
# sample integration
sample_integrations <- list(
  c("NS", "M", "NC", "NSM", "A", "C", "N", "t90", "t149", "t275")
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
                  "AIF1", "C1QA", "C1QB", "CX3CR1"),
  "microglia_absence" = c("CD163", "CCL2", "MRC1"),
  "proliferating" = c("MKI67", "SOX2", "HOPX", "NES", "POU5F1"),
  "supplement2" = c("TGFB1", "CSF1", "IL34", "LEFTY2")
)
pseudotime_root_markers <- list(
  "Microglia" = c("AIF1"),
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

integrated_sample_names <- "atlas"
# integrated_sample_names <- unlist(sapply(sample_integrations, simplify = F, function(integrated_sample_name) {
#   paste(integrated_sample_name, collapse = "-")
# }))
integrated_sample_files <- unname(unlist(sapply(integrated_sample_names, simplify = F, function(integrated_sample_name) {
  file.path(project_dir, "results", run_name, integrated_sample_name, paste0(integrated_sample_name, ".rds"))
})))
#### END USER INITIALIZATION #####



# initial sample analysis
for (sample_name in sample_names) {
  message("RUNNING sample_analysis sample: ", sample_name)
  sample_analysis(
    samples_dir = samples_dir,
    sample_name = sample_name,
    output_dir = results_dir,
    features_of_interest = features_of_interest,
    run_cell_cycle_regression = FALSE
  )
}

# integrate samples
message("RUNNING samples_integration")
samples_integration(
  sample_files = c(
    file.path(results_dir, sample_integrations[[1]][1], paste0(sample_integrations[[1]][1], ".rds")),
    file.path(results_dir, sample_integrations[[1]][2], paste0(sample_integrations[[1]][2], ".rds")),
    file.path(results_dir, sample_integrations[[1]][3], paste0(sample_integrations[[1]][3], ".rds")),
    file.path(results_dir, sample_integrations[[1]][4], paste0(sample_integrations[[1]][4], ".rds")),
    file.path(results_dir, sample_integrations[[1]][5], paste0(sample_integrations[[1]][5], ".rds")),
    file.path(results_dir, sample_integrations[[1]][6], paste0(sample_integrations[[1]][6], ".rds")),
    file.path(results_dir, sample_integrations[[1]][7], paste0(sample_integrations[[1]][7], ".rds")),
    file.path(results_dir, sample_integrations[[1]][8], paste0(sample_integrations[[1]][8], ".rds")),
    file.path(results_dir, sample_integrations[[1]][9], paste0(sample_integrations[[1]][9], ".rds")),
    file.path(results_dir, sample_integrations[[1]][10], paste0(sample_integrations[[1]][10], ".rds"))
  ),
  sample_names = sample_integrations[[1]],
  output_dir = results_dir,
  features_of_interest = features_of_interest,
  integration_method = integration_method,
  integrated_sample_name = "atlas"
)
