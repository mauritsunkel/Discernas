library(ggplot2)
library(ggrepel)



### USER PARAMETERS
## A+C
file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv"
## N+C
# file <- "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_N + BL_C/after_selection_old/DE_analysis/sample_markers/pct1=BL_C-pct2=BL_N - (nz-)p-val st 0.05.csv"

# log2FC threshold for labelling differentially expressed genes
log2FC.threshold <- 1
# amount of top genes to plot
top_n <- 5
### END USER PARAMETERS



# read data
data <- read.csv2(file)

# set up/down/no differential expression conditions
data$diffexpressed <- "NO"
data$diffexpressed[data$avg_log2FC > log2FC.threshold & data$p_val_adj < 0.05] <- "DOWN"
data$diffexpressed[data$avg_log2FC < -log2FC.threshold & data$p_val_adj < 0.05] <- "UP"

# color for diff expressed
diffcolors <- c("blue", "red", "black")
names(diffcolors) <- c("DOWN", "UP", "NO")

# define up and down thresholds for which top n genes to label in plot
up <- data$avg_log2FC[data$p_val_adj == 0 & data$avg_log2FC > log2FC.threshold][top_n]
down <- data$avg_log2FC[data$p_val_adj == 0 & data$avg_log2FC < log2FC.threshold]
down <- down[length(down)-top_n+1]
# label up and down regulated genes
data$delabel <- NA
data$delabel[data$diffexpressed != "NO" & data$p_val_adj != 0] <- data$X[data$diffexpressed != "NO" & data$p_val_adj != 0] #  & data$p_val_adj != 0
data$delabel[data$p_val_adj == 0 & data$avg_log2FC >= up]  <- data$X[data$p_val_adj == 0 & data$avg_log2FC >= up]
data$delabel[data$p_val_adj == 0 & data$avg_log2FC <= down]  <- data$X[data$p_val_adj == 0 & data$avg_log2FC <= down]

# plot data
ggplot(data=data, aes(x=-avg_log2FC, y=-log10(p_val_adj), col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_classic() +
  geom_text_repel(
    box.padding = .25,
    max.overlaps = 15,
    label.size = .5
  ) +
  scale_color_manual(values = diffcolors) +
  geom_vline(xintercept=c(-log2FC.threshold, log2FC.threshold), col="red") +
  geom_hline(yintercept=-log10(0.01), col="red") +
  xlab("Average log2 Fold Change") +
  ylab("-log10 P-val (Bonferroni adjusted)")

# save plot
ggplot2::ggsave(file=paste0("A+C DEG volcano plot v2.png"), width = 20, height = 20, units = "cm")
