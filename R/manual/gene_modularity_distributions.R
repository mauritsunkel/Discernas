data <- readRDS("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/BL_A + BL_C.rds")

gene <- "VIM"

blc <- data$orig.ident %in% "BL_C"
bla <- data$orig.ident %in% "BL_A"

plot(density(as.numeric(data@assays$RNA[gene,blc])), col = "green")
lines(density(as.numeric(data@assays$RNA[gene,bla])), col = "red")

# TODO would be interesting to automatically compare distributions in order to define modularity type per gene
## and then globally define how modular the dataset is
