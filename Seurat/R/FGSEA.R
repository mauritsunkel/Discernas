library(fgsea)
library(org.Hs.eg.db) # Organism.HomoSapiens.EntrezGene.DataBase
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)



set.seed(42)
lfc_threshold <- 1

# load DEA result
dea_result <- read.csv2("C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/method=DESeq2-pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv")

## if .csv was saved in Excel, some gene names e.g. SEP2 become dates
# get index of these genes
excel_genes_ind <- grepl(x = dea_result$X, pattern = "^[A-Z][a-z]{2}/[0-9]{2}$")
# get these genes
excel_genes <- dea_result$X[excel_genes_ind]
# convert back to original names from date names
dea_result$X[excel_genes_ind] <- sapply(X = excel_genes, USE.NAMES = FALSE, FUN = function(X) {
  ss <- unlist(strsplit(X, split = '/'))
  paste0(toupper(ss[1]), as.integer(ss[2]))
})


# TODO test swapping signs
# dea_result$avg_log2FC <- dea_result$avg_log2FC * -1
# filter by log2 fold change threshold
dea_result <- dea_result[dea_result$avg_log2FC > lfc_threshold,]
## DEVNOTE: some tens of genes like SEPT10 became a date in excel, hard to read back in properly
## DEVNOTE: don't remove '.' and everything after because of matching ENSEMBL IDs with ENSEMBL symbols at origin
## DEVNOTE: Normally remove for proper Ensembl gene version matching to Entrez
## dea_result$X <- sub("\\.[^.]*$", "", dea_result$X)
# load CellRanger features - Ensembl genes and IDs
ensembl_genes_ids <- read.delim(file = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Ensembl genes.tsv",
                                col.names = c("ensembl_gene_id", "ensembl_gene_name", "gene_type"))
dea_result$ensembl_id <- ensembl_genes_ids$ensembl_gene_id[match(dea_result$X, ensembl_genes_ids$ensembl_gene_name)]
entrez_ids <- mapIds(org.Hs.eg.db, keys=dea_result$ensembl_id, keytype="ENSEMBL", column="ENTREZID", multiVals="first")
# set length 0 values to NA, otherwise lose them while unlisting
entrez_ids[lengths(entrez_ids) == 0] <- NA
dea_result$entrez_id <- unname(unlist(entrez_ids))
# get DEG names and lfc
deg_ids_lfc <- dea_result$avg_log2FC
names(deg_ids_lfc) <- dea_result$entrez_id
# remove NA
deg_ids_lfc <- deg_ids_lfc[!is.na(names(deg_ids_lfc))]
# sort from high to low
deg_ids_lfc <- sort(deg_ids_lfc, decreasing = TRUE)




gene_ontology_types <- c("biological_process" = "BP")
# gene_ontology_types <- c("biological_process" = "BP",
#                          "molecular_function" = "MF",
                         # "cellular_component" = "CC")
for (ind in 1:length(gene_ontology_types)) {
  # GOEA: provide a set of genes, compare against predefined gene set
  ## look for over representtion, but not by rank!
  gse <- clusterProfiler::gseGO(
    geneList = deg_ids_lfc,
    ont = gene_ontology_types[ind],
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    minGSSize = 3, # default: 10
    maxGSSize = 800, # default: 500
    pvalueCutoff = 1, # default: 0.05
    eps = 0, # default: 1e-10
    verbose = TRUE,
    pAdjustMethod = "none", # default: "BH" (Benjamini-Hochberg)
    seed = TRUE,
    by = "fgsea",
    scoreType = "pos"
  )

  enrichplot::dotplot(
    gse,
    title=paste0("Gene Ontology: ", names(gene_ontology_types[ind])),
    # split=".sign", # uncomment to show suppressed GO terms
    showCategory=50) +
    # facet_grid(.~.sign) +
    labs(subtitle = paste0("Log fold-change threshold: >", lfc_threshold))

  # TODO uncomment
  # ggsave(file = paste0("enrichmentPlot_geneOntology_", names(gene_ontology_types[ind]), ".png"), width = 30, height = 20, units = "cm")
}

# TODO create issue on clusterProfiler, mention how certain genelist with only pos and length 55 works
## other with only pos and length only works for MF, not BP CC or ALL
## with the error: Error in FUN(X[[i]], ...) : GSEA statistic is not defined when all genes are selected
### on tracing back through package code the error is generated once length of 'stats' (geneList) vs 'selectedStats' (geneSets/(filtered) pathways) are compared
#### is this a bug?
#### Anyhow, how to overcome this?
# based on the fsea::preparePathwaysAndStats() function I assume the filtered pathways become the same length as the geneList

# should I post this here or to fgsea @ctlab Github?

# clusterProfiler::gseGO --> clusterProfiler::GSEA_internal --> GSEA_internal <- DOSE:::GSEA_internal
## https://github.com/YuLab-SMU/clusterProfiler/blob/f164715c2eb696542ba31141c964cfd407201c2b/R/gseAnalyzer.R
# .GSEA <- GSEA_fgsea
## https://github.com/YuLab-SMU/DOSE/blob/d94b63f50c3ac8c238800e4052576aa7f6af63ef/R/gsea.R
# DOSE::fgsea
## https://github.com/YuLab-SMU/DOSE/blob/d94b63f50c3ac8c238800e4052576aa7f6af63ef/R/gsea.R
# fgsea::fgseaMultilevel
## https://github.com/ctlab/fgsea/blob/c4fb22b2deae6d8808b90439ea3d7297dfc7b349/R/fgsea.R
# fgsea::preparePathwaysAndStats & fgsea::calcGseaStat
## https://github.com/ctlab/fgsea/blob/master/R/fgsea.R
### https://github.com/ctlab/fgsea/blob/c4fb22b2deae6d8808b90439ea3d7297dfc7b349/R/fgsea.R preparePathwaysAndStats

# TODO can switch x and size by passing 'count' and 'geneRatio' to one and the other


# internl function used by clusterProfiler to select data from gse object
df <- fortify(gse, showCategory = 10, split=".sign")
df <- fortify(gse, showCategory = 10, split=NULL)
df


gse <- df[df$.sign == "activated"]







enrichGO() # for subset analysis, and use full set for gseGO: https://github.com/YuLab-SMU/clusterProfiler/issues/509


gseKEGG()
gseMKEGG()


# GSEA works a bit differently.
# For GSEA, you will provide ALL the genes in your analysis without a cutoff.
# You will rank them and weigh them by some statistic (fold change, or log(p) for example).
# The software will walk through from top to bottom, investigating a pathway.
# If it runs into a gene in the pathway it'll boost your enrichment score,
# if it doesn't it'll subtract some points.
# The number of points it gives depends on the rank/weight,
# such that the higher the rank and fold change, the bigger the boost in enrichment score.
# At the end you're left with a single score that is high if the pathway genes are clustered at the top/bottom,
# and low if they are randomly distributed throughout your list.





# use Entrez IDs with ClusterProfiler etc to perform GSEA/GOEA/pathway analyses with FGSEA








# FGSEA_analysis <- function(markers, working_directory, marker_type, cluster) {
#   library(biomaRt)
#   library(fgsea)
#   library(data.table)
#   library(ggplot2)
#
#   dir.create(paste0(work_dir, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/"), recursive = T)
#
#   ## fix infinite values later by applying -log10 function
#   markers$p_val[markers$p_val == 0] <- min(markers$p_val[markers$p_val != 0])
#
#   ## calculate metric by FoldChangeSign and -LogPvalue
#   markers$fcsign <- sign(markers$avg_log2FC)
#   markers$logPval <- -log10(markers$p_val)
#
#   ## create ranked vector
#   fgsea_ranks <- markers$logPval/markers$fcsign
#
#   ## get Entrez IDs by HGNC symbol to match gene names and provide a translation map
#   hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
#   mapping <- getBM(
#     attributes = c('entrezgene_id', 'hgnc_symbol'),
#     filters = 'hgnc_symbol',
#     values = rownames(markers),
#     mart = hsmart
#   )
#   names(fgsea_ranks) <- match(rownames(markers), mapping$hgnc_symbol)
#
#   ## get Reactome pathways by Entrez IDs
#   pathways <- reactomePathways(names(fgsea_ranks))
#
#   fgsea_results <- fgsea(pathways = pathways,
#                          stats    = fgsea_ranks,
#                          eps      = 0.0,
#                          minSize  = 15,
#                          maxSize  = 500)
#
#   topPathwaysUp <- fgsea_results[ES > 0][head(order(pval), n=10), pathway]
#   topPathwaysDown <- fgsea_results[ES < 0][head(order(pval), n=10), pathway]
#   topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#
#   png(filename=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/overview_table.png"), width = 1600)
#   plotGseaTable(pathways[topPathways], fgsea_ranks, fgsea_results,
#                 gseaParam=0.5)
#   dev.off()
#
#   ## can try to collapse pathways if there are many seemlingly alike in the plot above
#   collapsedPathways <- collapsePathways(fgsea_results[order(pval)][padj < 0.01],
#                                         pathways, fgsea_ranks)
#   mainPathways <- fgsea_results[pathway %in% collapsedPathways$mainPathways][
#     order(-NES), pathway]
#
#   ## check if mainPathways is empty (likely collapsedPathways is empty too)
#   if (length(mainPathways) > 0) {
#     png(filename=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/collapsed_table.png"), width = 1600)
#     p <- plotGseaTable(pathways[mainPathways], fgsea_ranks, fgsea_results,
#                        gseaParam = 0.5)
#     dev.off()
#   }
#
#   fwrite(fgsea_results, file=paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/overview.xls"), sep="\t", sep2=c("", ",", ""))
#
#   for (i in seq_along(topPathways)) {
#     # png(filename=paste0(working_directory, "GSEA/cluster_", cluster, "/enriched_", i, ".png"), width = 1600)
#     p <- plotEnrichment(pathways[[topPathways[i]]],
#                         fgsea_ranks) + labs(title=topPathways[[i]])
#     ggsave(file = paste0(working_directory, "../GSE_analysis/", marker_type, "/cluster ", cluster, "/enriched_", i, ".png"), width = 30, height = 20, units = "cm")
#   }
# }
