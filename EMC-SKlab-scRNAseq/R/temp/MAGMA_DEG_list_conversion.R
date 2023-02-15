

# df <- read.csv2("C:/Users/mauri/Desktop/astrocytes.csv")
df <- read.csv2("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/method=DESeq2-pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv")

## if .csv was saved in Excel, some gene names e.g. SEP2 become dates
# get index of these genes
excel_genes_ind <- grepl(x = df$X..., pattern = "^[A-Z][a-z]{2}/[0-9]{2}$")
# get these genes
excel_genes <- df$X...[excel_genes_ind]
# convert back to original names from date names
df$X...[excel_genes_ind] <- sapply(X = excel_genes, USE.NAMES = FALSE, FUN = function(X) {
  ss <- unlist(strsplit(X, split = '/'))
  paste0(toupper(ss[1]), as.integer(ss[2]))
})


df$avg_log2FC <- df$avg_log2FC*-1
df$meanExpr.astrocytes.monoculture <- rev(df$meanExpr.astrocytes.monoculture)
df$meanExpr.astrocytes.coculture <- rev(df$meanExpr.astrocytes.coculture)
rownames(df) <- df$X...


write.csv2(x=df, file = "C:/Users/mauri/Desktop/astrocytes2.csv")






