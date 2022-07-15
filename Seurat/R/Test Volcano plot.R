library(ggplot2)
library(ggrepel)

### USER PARAMETERS
## A+C
# file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv"
## N+C
file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_N + BL_C/after_selection_old/DE_analysis/sample_markers/pct1=BL_C-pct2=BL_N - (nz-)p-val st 0.05.csv"

log2FC.threshold <- 1
### END USER PARAMETERS


# read data
data <- read.csv2(file)

# add up/down/no differential expression flags
data$diffexpressed <- "NO"
data$diffexpressed[data$avg_log2FC > log2FC.threshold & data$p_val_adj < 0.05] <- "UP"
data$diffexpressed[data$avg_log2FC < -log2FC.threshold & data$p_val_adj < 0.05] <- "DOWN"

# color for diffexpressed
diffcolors <- c("blue", "red", "black")
names(diffcolors) <- c("DOWN", "UP", "NO")

# label up and down regulated genes
data$delabel <- NA
data$delabel[data$diffexpressed != "NO"] <- data$X[data$diffexpressed != "NO"]

# plot data
ggplot(data=data, aes(x=avg_log2FC, y=-log10(p_val_adj), col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = diffcolors) +
  geom_vline(xintercept=c(-log2FC.threshold, log2FC.threshold), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  xlab("Average log2 Fold Change") +
  ylab("-log10 P-val (Bonferroni adjusted)")
