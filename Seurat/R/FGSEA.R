library(fgsea)
library(org.Hs.eg.db) # Organism.HomoSapiens.EntrezGene.DataBase
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)


### USER PARAMETERS
set.seed(42)
lfc_threshold <- 1
### END USER PARAMETERS ###




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

## DEVNOTE: don't remove '.' and everything after because of matching ENSEMBL IDs with ENSEMBL symbols at origin
## DEVNOTE: Normally remove for proper Ensembl gene version matching to Entrez
## dea_result$X <- sub("\\.[^.]*$", "", dea_result$X)
# load CellRanger features - Ensembl genes and IDs
ensembl_genes_ids <- read.delim(file = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Ensembl genes.tsv",
                                col.names = c("ensembl_gene_id", "ensembl_gene_name", "gene_type"))
dea_result$ensembl_id <- ensembl_genes_ids$ensembl_gene_id[match(dea_result$X, ensembl_genes_ids$ensembl_gene_name)]
# map ENSEMBL to ENTREZIDs
entrez_ids <- bitr(geneID = dea_result$ensembl_id, OrgDb = org.Hs.eg.db, fromType = "ENSEMBL", toType = "ENTREZID")
## DEVNOTE: clusterProfiler::bitr mapped more non-NA than annotationDbi::mapIds
# entrez_ids <- mapIds(org.Hs.eg.db, keys=dea_result$ensembl_id, keytype="ENSEMBL", column="ENTREZID", multiVals="first")
# set length 0 values to NA, otherwise lose them while unlisting
# entrez_ids[lengths(entrez_ids) == 0] <- NA

# get DEG names and lfc
deg_ids_lfc <- dea_result$avg_log2FC
names(deg_ids_lfc) <- entrez_ids$ENTREZID
# remove NA
deg_ids_lfc <- deg_ids_lfc[!is.na(names(deg_ids_lfc))]
# sort decreasingly
deg_ids_lfc <- sort(deg_ids_lfc, decreasing = TRUE)




# TODO check which genes are not mapped?
## check of those which are in head/tail 100 positions, because they would be very weighty in GSEA
## and check specifically which ones are in abs(log2fc) > threshold positions

## TODO perform gseXXX functions first




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
    pAdjustMethod = "none", # TODO try (always used by GuidoHooiberg fromclusterProfiler) the default: "BH" (Benjamini-Hochberg)
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



## TODO then subset with dplyr::filter(abs(log2fc) > threshold)
# TODO test swapping signs
# dea_result$avg_log2FC <- dea_result$avg_log2FC * -1
# filter by log2 fold change threshold
dea_result <- dea_result[dea_result$avg_log2FC > lfc_threshold,]

## TODO then run enrichXXX functions






# clusterProfiler vignette: https://yulab-smu.top/biomedical-knowledge-mining-book/useful-utilities.html
# general intro ORA vs GSEA: https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/
# Over-representation (enrichment) analysis (ORA): check DEG vs known biological relationship (via gene set)
## will find genes where the difference is large and will fail where the difference is small, but evidenced in coordinated way in a set of related genes
## limitations: 1) arbitrary user-defined inclusion criteria, 2) equal importance to each gene, 3) assume gene independance
# Gene set enrichment analysis (GSEA): unordered collection of related genes (sets)
## Pathway: type of gene set if functional relations are ignored
## gene set enrichment analysis (GSEA): aggregates the per gene statistics across genes within a gene set
### overcoming ORA limitation, by detecting situations where all genes in a predefined set change in a small but coordinated way
## steps
### local gene-level statistic: adj-p-val / log fold change
### global gene-set statistics: weighted running sum by ranked-list gene prevalence, Enrichment Score (ES) = max score
#### Normalized ES (NES): accounting for gene set size
### significance testing & multiple testing correction
#### Type 1 error: family of hypotheses, chance of false positive (significant p-val)
##### False Discovery Rate (FDR): fraction of rejected null hypothesis (Type 1 error) that are true (false negatives)
###### https://www.pathwaycommons.org/guide/primers/statistics/multiple_testing/
###### https://www.nature.com/articles/nbt1209-1135
###### https://jtd.amegroups.com/article/view/13609/11598

## OVER REPRESENTATION ANALYSIS (one-sided version of Fishers exact test, input DEG subsets)
# enrichGO, enrichKEGG, enrichDO, enrichWP, enrichPathway, and the universal one enrichr
## GO=Gene ontology, DOSE=Disease Ontology, WP=WikiPathways
enrichGO()
enrichKEGG()
enrichWP()
enrichPathway() # ReactomePA package
enrichDO() # Disease Ontology
enrichDGN() # Disease Gene Network
enricher() # custom annotation, like MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)

## GENE SET ENRICHMENTS ANALYSIS (full DEG)
# gseGO, gseKEGG, gseWP, and the universal one GSEA
gseGO()
gseKEGG()
gseWP()
gsePathway() # ReactomePA package
gseDO() # Disease Ontology
gseDGN() # Disease Gene Network
GSEA() # custom annotation, like MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)


## LEADING EDGE ANALYSIS
# https://guangchuangyu.github.io/2016/07/leading-edge-analysis/



## GENERAL
goplot(gse)
setReadable() # get symbols with IDs
bitr() # convert gene IDs
compareCluster() # try when have multiple gene lists from different experiments/conditions

## VISUALIZATION
enrichplot() # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#enrichplot
browseKEGG()
pathview()
viewPathway() # reactome pathway analysis







# use Entrez IDs with ClusterProfiler etc to perform GSEA/GOEA/pathway analyses with FGSEA





## internal function used by clusterProfiler to select data from gse object
# df <- fortify(gse, showCategory = 10, split=".sign")
# df <- fortify(gse, showCategory = 10, split=NULL)
# df
# gse <- df[df$.sign == "activated"]



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
