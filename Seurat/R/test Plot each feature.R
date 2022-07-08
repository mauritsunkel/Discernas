library(Seurat)
library(ggplot2)

### USER PARAMETERS
sample_names <- c('BL_N + BL_C')

rds.files <- c("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_N + BL_C/after_selection_new/BL_N + BL_C.rds")

# work dir should contain forward slashes (/) on Windows
start_time <- format(Sys.time(), "%F %H-%M-%S")
work_dir <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/"




plot_DEF <- function(data, features, name) {
  dir.create(paste0("DE_analysis/", name, "/"))
  dir.create(paste0("DE_analysis/", name, "/Feature/"))

  for (i in seq_along(features)) {
    tryCatch({
      p <- Seurat::FeaturePlot(data, features = features[i])
      ggplot2::ggsave(file=paste0("DE_analysis/", name ,"/Feature/", features[i], ".png"), width = 30, height = 20, units = "cm")
    },
      error=function(e) {
        message(features[i], ' plot is skipped, as it was not found with FetchData')
      })
  }
}



for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  message(sample_name)
  data <- readRDS(rds.files[i])

  dir.create(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"), recursive = T)
  setwd(paste0(work_dir, 'results/', start_time, "/", sample_name, "/"))
  dir.create('DE_analysis/')



  astrocyte_interest <- c("GFAP", "VIM", "S100B", "SOX9", "CD44", "AQP4", "ALDH1L1",
                          "HIST1H4C", "FABP7", "SLC1A2", "SLC1A3", "GJA1", "APOE")
  neuron_interest <- c("TUBB3", "MAP2", "CAMK2A", "GAD2", "NEUROG2", "SYN1", "RBFOX3", "GJA1")
  astrocyte_maturity <- c("CD44", "FABP7", "VIM", "SOX9", "TOP2A", "S100B",
                          "GJA", "SLC1A3", "IGFBP7", "ALDH1L1", "APOE")
  neuron_maturity <- c("NEUROG2", "DCX", "MAP2", "RBFOX3",
                       "SYN1", "SNAP25", "SYT1", "APOE")
  plot_DEF(data = data, features = astrocyte_interest, name = "astrocyte")
  plot_DEF(data = data, features = neuron_interest, name = "neuron")
  plot_DEF(data = data, features = astrocyte_maturity, name = "astrocyte_maturity")
  plot_DEF(data = data, features = neuron_maturity, name = "neuron_maturity")
}
