

setwd("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/after_selection/Plots/astrocyte")

data <- readRDS("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe +SCT +Leiden -Cellcycle +SingleR +Autoselection +05-05-2022/integrated/BL_A + BL_C/after_selection/BL_A + BL_C.rds")

Seurat::FeaturePlot(data, features = "APOE")
ggplot2::ggsave(file=paste0("APOE.png"), width = 30, height = 20, units = "cm")

Seurat::FeaturePlot(data, features = "APOE", split.by = "orig.ident", cols = c("grey", "red"))
ggplot2::ggsave(file=paste0("APOE-split.png"), width = 30, height = 20, units = "cm")
