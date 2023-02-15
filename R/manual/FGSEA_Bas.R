# get from manually running run_fgsea function ...
deg_ids
deg_ids_positive
deg_ids_negative

setwd("C:/Users/mauri/Desktop")
p_adjust_cutoff <- 0.05
pval_cutoff <- 0.05
qval_cutoff <- 0.05

res_enrich_go <- clusterProfiler::enrichGO(
  gene = names(deg_ids),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = 'BP',
  pvalueCutoff = pval_cutoff,
  pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
  universe = go_universe,
  qvalueCutoff = qval_cutoff,
  minGSSize = 10, # default: 10
  maxGSSize = 500, # default: 500
  readable = FALSE, # default: FALSE
  pool = FALSE # default: FALSE
)
res_enrich_go_pos <- clusterProfiler::enrichGO(
  gene = deg_ids_positive,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = pval_cutoff,
  pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
  universe = go_universe,
  qvalueCutoff = qval_cutoff,
  minGSSize = min_gene_set_size_ora, # default: 10
  maxGSSize = max_gene_set_size_ora, # default: 500
  readable = FALSE, # default: FALSE
  pool = FALSE # default: FALSE
)
res_enrich_go_neg <- clusterProfiler::enrichGO(
  gene = deg_ids_negative,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = 'BP',
  pvalueCutoff = pval_cutoff,
  pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
  universe = go_universe,
  qvalueCutoff = qval_cutoff,
  minGSSize = min_gene_set_size_ora, # default: 10
  maxGSSize = max_gene_set_size_ora, # default: 500
  readable = FALSE, # default: FALSE
  pool = FALSE # default: FALSE
)

# filter on p.adjust, pvalue, qvalue (devnote: thresholds do not get applied in enrichGO results)
res_enrich_go@result <- res_enrich_go@result %>% dplyr::filter(p.adjust < p_adjust_cutoff)
res_enrich_go_pos@result <- res_enrich_go_pos@result %>% dplyr::filter(p.adjust < p_adjust_cutoff)
res_enrich_go_neg@result <- res_enrich_go_neg@result %>% dplyr::filter(p.adjust < p_adjust_cutoff)
# set readable for plots
res_enrich_go <- setReadable(res_enrich_go, 'org.Hs.eg.db', 'ENTREZID')
res_enrich_go_pos <- setReadable(res_enrich_go_pos, 'org.Hs.eg.db', 'ENTREZID')
res_enrich_go_neg <- setReadable(res_enrich_go_neg, 'org.Hs.eg.db', 'ENTREZID')
# aggregate selected terms if needed
term_select <- c(
  "antigen processing and presentation of peptide",
  "antigen processing and presentation of peptide antigen via MHC class II",
  "antigen processing and presentation of exogenous peptide antigen via MHC class II",
  "antigen processing and presentation of exogenous peptide antigen",
  "antigen processing and presentation of exogenous antigen",
  "antigen processing and presentation of peptide antigen",
  "antigen processing and presentation of endogenous peptide antigen",
  "antigen processing and presentation of endogenous antigen",
  "antigen processing and presentation")
collapse_terms <- function(term_select, ORA_result, new_description) {
  # rename for clarity
  df <- ORA_result@result

  # get index of selected terms
  index <- df$Description %in% term_select

  # summarize terms selected to single row:
  ## max(p.adjust, pvalue, qvalue), first (ID, GEneRatio, BgRatio, Count), re-overlap for geneID
  new_row <- data.frame(
    "ID" = df[index,][1,]$ID,
    "GeneRatio" = df[index,][1,]$GeneRatio,
    "BgRatio" = df[index,][1,]$BgRatio,
    "Count" = df[index,][1,]$Count,
    "pvalue" = min(df[index,]$pvalue),
    "qvalue" = min(df[index,]$qvalue),
    "p.adjust" = min(df[index,]$p.adjust),
    "geneID" = paste0(unique(unlist(strsplit(df[index,]$geneID, "/"))), collapse = "/"),
    "Description" = new_description
  )

  # replace geneSet
  res_enrich_go_pos@geneSets[new_row[,"ID"]] <- new_row[,"geneID"]

  # get rowname
  rowname <- df[index,][1,]$ID
  # remove rows by index
  df <- df[-which(index),]
  # add new row
  df <- rbind(new_row, df)
  # set rowname of new row
  rownames(df)[nrow(df)] <- rowname
  # sort
  df %>% arrange(desc(Count), pvalue, p.adjust, qvalue)

  # overwrite
  ORA_result@result <- df
  # remove index terms from semantic similarity
  termsim_index <- which(rownames(ORA_result@termsim) %in% term_select)
  ORA_result@termsim <- ORA_result@termsim[-termsim_index,]
  # recalculate termsim
  ORA_result <- pairwise_termsim(ORA_result)

  return(ORA_result)
}
res_enrich_go_pos <- collapse_terms(
  term_select = term_select,
  ORA_result = res_enrich_go_pos,
  new_description = "antigen processing and presentation")
# calculate pairwise_termsim
res_enrich_go <- pairwise_termsim(res_enrich_go)
res_enrich_go_pos <- pairwise_termsim(res_enrich_go_pos)
res_enrich_go_neg <- pairwise_termsim(res_enrich_go_neg)
# plot and manually save
enrichplot::treeplot(
  res_enrich_go,
  showCategory = 40,
  label_format = 14)
enrichplot::treeplot(
  res_enrich_go_pos,
  showCategory = 20,
  label_format = 14)
enrichplot::treeplot(
  res_enrich_go_neg,
  showCategory = 20,
  label_format = 14)
# save the data to excel
write.csv2(res_enrich_go@result, "ORA_astrocyte_mix_results_sort=p.adjust.csv")
write.csv2(res_enrich_go_pos@result, "ORA_astrocyte_pos_results_sort=p.adjust.csv")
write.csv2(res_enrich_go_neg@result, "ORA_astrocyte_neg_results_sort=p.adjust.csv")
























# res_r_pt <- res_enrich_go
# res_r_pt_pos <- res_enrich_go_pos
# res_r_pt_neg <- res_enrich_go_neg
# nsig <- sum(res_r_pt@result$p.adjust < 0.05)
# genesets <- res_r_pt@geneSets[match(res_r_pt@result$ID[1:nsig], names(res_r_pt@geneSets))]
# geneset_length <- sapply(genesets, function(x) {
#   length(x)
# })
# upcount <- sapply(genesets, function(x) {
#   sum(deg_ids_positive %in% unname(unlist(x)))
# })
# upratio <- round((upcount / res_r_pt@result$Count[1:nsig]) * 100, 1)
# downcount <- sapply(genesets, function(x) {
#   sum(deg_ids_negative %in% unname(unlist(x)))
# })
# downratio <- round((downcount / res_r_pt@result$Count[1:nsig]) * 100, 1)
# df <- data.frame(
#   "GO term ID" = res_r_pt@result$ID[1:nsig],
#   "GO term name" = res_r_pt@result$Description[1:nsig],
#   "GO term nGenesTotal" = geneset_length,
#   "GO term nGenesMatch" = res_r_pt@result$Count[1:nsig],
#   "GO term geneSet" = res_r_pt@result$geneID[1:nsig],
#   "p-val-adj" = res_r_pt@result$p.adjust[1:nsig],
#   "geneRatio" = res_r_pt@result$GeneRatio[1:nsig],
#   "upcount" = upcount,
#   "upRatio" = upratio,
#   "downcount" = downcount,
#   "downRatio" = downratio
# )
# df
# write.csv2(df, file = "ORA_astrocyte_mix_upAndDownRatio_GO-BP_info.csv")

# library(VennDiagram)
# threshold <- 0.05
# n_pos <- sum(res_r_pt_pos@result$p.adjust < threshold)
# n_neg <- sum(res_r_pt_neg@result$p.adjust < threshold)
# n_mix <- sum(res_r_pt@result$p.adjust < threshold)
# terms_pos <- res_r_pt_pos@result$Description[1:n_pos]
# terms_neg <- res_r_pt_neg@result$Description[1:n_neg]
# terms_mix <- res_r_pt@result$Description[1:n_mix]
# grid.newpage()
# venn.plot <- draw.triple.venn(
#   area1 = n_pos,
#   area2 = n_neg,
#   area3 = n_mix,
#   n12 = length(intersect(terms_pos, terms_neg)),
#   n23 = length(intersect(terms_neg, terms_mix)),
#   n13 = length(intersect(terms_pos, terms_mix)),
#   n123 = length(Reduce(intersect, list(terms_pos, terms_neg, terms_mix))),
#   category = c("Upregulated", "Downregulated", "Mixed"),
#   fill = c("blue", "red", "green"),
#   lty = "blank",
#   cex = 2,
#   cat.cex = 2,
#   cat.col = c("blue", "red", "green")
# )
# write.csv2(setdiff(terms_pos, c(terms_neg, terms_mix)), "terms_only_in_pos.csv")
# write.csv2(setdiff(terms_neg, c(terms_pos, terms_mix)), "terms_only_in_neg.csv")
# write.csv2(setdiff(terms_mix, c(terms_neg, terms_pos)), "terms_only_in_mix.csv")
# write.csv2(Reduce(intersect, list(terms_pos, terms_neg, terms_mix)), "terms_in_all.csv")




# ## CHECKING GENES IN UP/DOWN BUT NOT MIXED ANALYSIS   &   CHECKING 9 OVERLAPPING TERMS IN EACH ANALYSES
# # check why some terms are in up/down but not mixed (based on nGenes, and nOverlapGenes)
# t <- setdiff(intersect(terms_pos, terms_neg), terms_mix)
# # check geneset length of 9 terms that are in each analysis
# t <- Reduce(intersect, list(terms_pos, terms_neg, terms_mix))
# t2 <- res_enrich_go_pos@result$Description
# indices <- match(t,t2)
# ids <- res_enrich_go_pos@result$ID[indices]
# df <- sapply(res_enrich_go_pos@geneSets[ids], function(x) {
#   data.frame(
#     neg = sum(deg_ids_negative %in% x),
#     pos = sum(deg_ids_positive %in% x),
#     negRatio = round((sum(deg_ids_negative %in% x)/length(x))*100, digits = 2),
#     posRatio = round((sum(deg_ids_positive %in% x)/length(x))*100, digits = 2),
#     totRatio = round(((sum(deg_ids_negative %in% x) + sum(deg_ids_positive %in% x)) / length(x))*100, digits = 2),
#     nGenes = length(x),
#     negGenes = paste(deg_ids_negative[deg_ids_negative %in% x], collapse="/"),
#     posGenes = paste(deg_ids_positive[deg_ids_negative %in% x], collapse="/")
#   )
# })
# colnames(df) <- t
# df
#
# write.csv2(t(df), "upAndDown_notMixed.csv")
# write.csv2(t(df), "9_allOverlapping.csv")
