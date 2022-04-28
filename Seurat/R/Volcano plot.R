
neg <- read.csv2(file = "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-11-02 13-06-19/integrated/BL_A + BL_C/DE_analysis/sample_markers/negative_pct1=BL_A-pct2=BL_C_p-val-adj=bonferroni.csv")
pos <- read.csv2(file = "C:/Users/mauri/Desktop/M/Erasmus MC PhD/Projects/Single Cell RNA Sequencing/Seurat/results/2021-11-02 13-06-19/integrated/BL_A + BL_C/DE_analysis/sample_markers/positive_pct1=BL_A-pct2=BL_C_p-val-adj=bonferroni.csv")

data <- rbind(neg, pos)

plot(x = data$avg_log2FC, y = -log(data$p_val), xlab = "Average log-2 Fold Change", ylab = "p-val (Bonferroni adjusted)", ylim = c(min(-log(data$p_val)), round(max(-log(data$p_val_adj)[is.finite(-log(data$p_val_adj))]), -1)))
abline(v = -.25, col = "red")
abline(v = .25, col = "red")
abline(h = -log(0.01), col = "red")
