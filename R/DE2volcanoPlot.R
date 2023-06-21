#' Visualize DEG in Volcano plot.
#'
#' @param DE_csv Input DEA as .csv format.
#' @param output_filename Plot output filename.
#' @param log2FC_threshold log2 fold-change threshold for defining DEG.
#' @param top_n default: 5. Amount of top up/down regulated genes to plot.
#'
#' @export
#'
#' @examplesIf FALSE
#' DE_csv <- file.path("DEA_results", "DEG.csv")
#' output_filename <- "A+C_DEG_volcano-plot_log2FC-threshold=1.png"
#' log2FC_threshold <- 1
#' top_n <- 5
#' DE2volcanoPlot(DE_csv, output_filename, log2FC_threshold, top_n = top_n)
DE2volcanoPlot <- function(DE_csv, output_filename, log2FC_threshold, top_n = 5) {
  # read data
  data <- utils::read.csv2(DE_csv)

  # set up/down/no differential expression conditions
  data$diffexpressed <- "NO"
  data$diffexpressed[data$avg_log2FC > log2FC_threshold & data$p_val_adj < 0.05] <- "DOWN"
  data$diffexpressed[data$avg_log2FC < -log2FC_threshold & data$p_val_adj < 0.05] <- "UP"

  # color for diff expressed
  diffcolors <- c("blue", "red", "black")
  names(diffcolors) <- c("DOWN", "UP", "NO")

  # define up and down thresholds for which top n genes to label in plot
  up <- data$avg_log2FC[data$p_val_adj == 0 & data$avg_log2FC > log2FC_threshold]
  if (length(up) < top_n) {
    up <- up[length(up)]
  } else {
    up[top_n]
  }
  down <- data$avg_log2FC[data$p_val_adj == 0 & data$avg_log2FC < log2FC_threshold]
  if(length(down) < top_n) {
    down <- down[1]
  } else {
    down <- down[length(down)-top_n+1]
  }
  # label up and down regulated genes
  data$delabel <- NA
  data$delabel[data$diffexpressed != "NO" & data$p_val_adj != 0] <- data$X[data$diffexpressed != "NO" & data$p_val_adj != 0]
  data$delabel[data$p_val_adj == 0 & data$avg_log2FC >= up]  <- data$X[data$p_val_adj == 0 & data$avg_log2FC >= up]
  data$delabel[data$p_val_adj == 0 & data$avg_log2FC <= down]  <- data$X[data$p_val_adj == 0 & data$avg_log2FC <= down]

  # set 0 to highest value to allow labels to be plotted above plot y-margin
  data$p_val_adj[data$p_val_adj == 0] <- min(data$p_val_adj[data$p_val_adj != 0], na.rm=TRUE)

  # plot data
  ggplot2::ggplot(data=data, ggplot2::aes(x=-avg_log2FC, y=-log10(p_val_adj), col = diffexpressed, label = delabel)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggrepel::geom_text_repel(
      box.padding = .25,
      max.overlaps = 15,
      label.size = .3) +
    ggplot2::scale_color_manual(values = diffcolors) +
    ggplot2::geom_vline(xintercept=c(-log2FC_threshold, log2FC_threshold), col="red") +
    ggplot2::geom_hline(yintercept=-log10(0.01), col="red") +
    ggplot2::xlab("Average log2 Fold Change") +
    ggplot2::ylab("-log10 P-val (Bonferroni adjusted)")

  # save plot
  ggplot2::ggsave(file=paste0(output_filename, ".png"), width = 20, height = 20, units = "cm")
}
