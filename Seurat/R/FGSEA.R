# clusterProfiler vignette: https://yulab-smu.top/biomedical-knowledge-mining-book
# general intro ORA vs GSEA: https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/
# use Entrez IDs with clusterProfiler etc to perform GSEA/ORA analyses with FGSEA
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

# TODO remove these when function incorporated in package
library(fgsea)
library(clusterProfiler)
library(DOSE)
library(ReactomePA)
library(org.Hs.eg.db) # Organism.HomoSapiens.EntrezGene.DataBase
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(ggupset)

# TODO remove these when function incorporated in package
### USER PARAMETERS
set.seed(42)
lfc_threshold <- 1
# set TRUE if up/down regulated GSEA visualization results need to be swapped
swap_GSEA_groups <- TRUE

use_internal_universe <- TRUE
p_adjust_method = "BH"
pval_cutoff <- 1 # default: 0.05
qval_cutoff <- 1 # default: 0.2
# Gene Set Size (GSS) is the number of genes match a term and are present in the gene universe (background)
min_gene_set_size_gsea <- 10 # default: 10
max_gene_set_size_gsea <- 500 # default: 500
min_gene_set_size_ora <- 3 # default: 10
max_gene_set_size_ora <- NA # default: 500
plot_n_category <- 30
gsea_plot_folder = "gene_set_enrichment_analysis"
ora_plot_folder = "over_representation_analysis"
### END USER PARAMETERS ###



# TODO testing inputs
## "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/method=DESeq2-pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv"
## "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/2022-11-01 14-02-50/integrated/BL_A + BL_C/DE_analysis/sample_markers/method=DESeq2-pct1=BL_C-pct2=BL_A - (nz-)p-val st 0.05.csv








#' Run GSEA and ORA.
#'
#' Run gene set enrichment and over representation analyses with clusterProfiler using fgsea.
#'
#' @param output_dir output directory for plots, string
#' @param dea_result_file Seurat differential expression analysis result object in .csv
#' @param cellRanger_ensembl_features CellRanger .csv file with mapping of Ensembl IDs to Ensembl gene symbols
#' @param lfc_threshold ORA differentially expressed genes log fold-change threshold, default: 1
#' @param seed set pseudorandomness, default: 42
#' @param swap_GSEA_groups invert log fold changes to visually swap upregulation/overrepresentation and downregulation/underrepresentation for GSEA, this is independent of positive/negative with ORA, default: FALSE
#' @param use_internal_universe TRUE: use DEA measured genes, FALSE: use clusterProfiler::ORA_db() specific genes, default: TRUE
#' @param p_adjust_method p-value multiple testing correction method, default: "BH" (Benjamini-Hochberg)
#' @param pval_cutoff GSEA/ORA term taken into account by p-val statistic
#' @param qval_cutoffGSEA/ORA term taken into account by q-val statistic
#' @param min_gene_set_size_gsea minimum gene set size for GSEA, default: 10
#' @param max_gene_set_size_gsea maximum gene set size for GSEA, default: 500
#' @param min_gene_set_size_ora minimum gene set size for ORA, default: 3
#' @param max_gene_set_size_ora maximum gene set size for ORA, default: NA
#' @param plot_n_category plot top n most significant terms, default: 30
#'
#' @return Nothing
#'
#' @import fgsea
#' @import clusterProfiler
#' @import DOSE
#' @import ReactomePA
#' @import org.Hs.eg.db
#' @import enrichplot
#' @import ggplot2
#' @import ggnewscale
#' @import ggupset
#'
#' @export
run_fgsea <- function(
    output_dir, dea_result_file, cellRanger_ensembl_features,
    lfc_threshold = 1, seed = 42,
    swap_GSEA_groups = FALSE, use_internal_universe = TRUE,
    p_adjust_method = "BH", pval_cutoff = 1, qval_cutoff = 1,
    min_gene_set_size_gsea = 10, max_gene_set_size_gsea = 500,
    min_gene_set_size_ora = 3, max_gene_set_size_ora = NA,
    plot_n_category = 30,
    gsea_plot_folder = "gene_set_enrichment_analysis",
    ora_plot_folder = "over_representation_analysis") {

  ## set and create output directories
  dir.create(output_dir, recursive = TRUE)
  setwd(output_dir)

  fgsea_dbs <- c(
    "gene_ontology", "disease_ontology", "disease_gene_network",
    "wiki_pathways", "reactome_pathways", "KEGG")
  for (db in fgsea_dbs) {
    dir.create(paste(gsea_plot_folder, db, sep = "/"), recursive = TRUE)
    dir.create(paste(ora_plot_folder, "upregulated", db, sep = "/"), recursive = TRUE)
    dir.create(paste(ora_plot_folder, "downregulated", db, sep = "/"), recursive = TRUE)
  }

  ## prep GSEA/ORA input
  # load DEA result
  dea_result <- read.csv2(dea_result_file)

  # DEVNOTE: non-unique/duplicate, as reconverting dates have multiple gene name options (e.g. MAR2/MARC2/MARCH2 --> 02/03/year)
  dea_result$X <- date2gene(gene_names = dea_result$X)

  ## DEVNOTE: don't remove '.' and everything after because of matching ENSEMBL IDs with ENSEMBL symbols at origin
  ## DEVNOTE: Normally remove for proper Ensembl gene version matching to Entrez
  ## dea_result$X <- sub("\\.[^.]*$", "", dea_result$X)
  # load CellRanger features - Ensembl genes and IDs
  ensembl_genes_ids <- read.delim(file = cellRanger_ensembl_features,
                                  col.names = c("ensembl_gene_id", "ensembl_gene_name", "gene_type"))
  # map enseml gene names to ensembl gene IDs from DEA result gene names
  dea_result$ensembl_id <- plyr::mapvalues(
    x = dea_result$X,
    from = ensembl_genes_ids$ensembl_gene_name,
    to = ensembl_genes_ids$ensembl_gene_id,
    warn_missing = FALSE)
  # map ENSEMBL IDs to ENTREZ IDs
  entrez_ids <- bitr(geneID = dea_result$ensembl_id, OrgDb = org.Hs.eg.db, fromType = "ENSEMBL", toType = "ENTREZID")

  fgsea_examine_unmapped_genes(
    gene_names = dea_result$X,
    dea_ensembl_ids = dea_result$ensembl_id,
    bitr_ensembl_ids = entrez_ids$ENSEMBL)

  ## prep GSEA input
  # get log fold-change (lfc)
  dea_ids_lfc <- dea_result$avg_log2FC
  # get Entrez names by matching Ensembl IDs
  names(dea_ids_lfc) <- plyr::mapvalues(x = dea_result$ensembl_id, from = entrez_ids$ENSEMBL, to = entrez_ids$ENTREZID)
  # remove NA
  dea_ids_lfc <- dea_ids_lfc[!is.na(names(dea_ids_lfc))]
  # swap signs to swap up/down regulation for GSEA visualization
  if (swap_GSEA_groups) {
    dea_ids_lfc <- dea_ids_lfc * -1
  }
  # sort decreasingly
  dea_ids_lfc <- sort(dea_ids_lfc, decreasing = TRUE)


  # TODO put this in a function
  ## prep ORA input
  # subset DEA result to DEG result based on log fold-change threshold
  deg_ids_positive <- names(dea_ids_lfc[dea_ids_lfc >= lfc_threshold])
  deg_ids_negative <- names(dea_ids_lfc[dea_ids_lfc <= -lfc_threshold])
  # save ORA deg subset names
  deg_names_positive <- bitr(geneID = deg_ids_positive, OrgDb = org.Hs.eg.db, fromType = "ENTREZID", toType = "ENSEMBL")
  deg_names_negative <- bitr(geneID = deg_ids_negative, OrgDb = org.Hs.eg.db, fromType = "ENTREZID", toType = "ENSEMBL")

  deg_names_positive <- plyr::mapvalues(
    x = deg_names_positive$ENSEMBL,
    from = ensembl_genes_ids$ensembl_gene_id,
    to = ensembl_genes_ids$ensembl_gene_name,
    warn_missing = FALSE)
  deg_names_positive <- deg_names_positive[!grepl("^ENSG", deg_names_positive)]
  deg_names_negative <- plyr::mapvalues(
    x = deg_names_negative$ENSEMBL,
    from = ensembl_genes_ids$ensembl_gene_id,
    to = ensembl_genes_ids$ensembl_gene_name,
    warn_missing = FALSE)
  deg_names_negative <- deg_names_negative[!grepl("^ENSG", deg_names_negative)]
  write.csv2(deg_names_positive, file = paste(ora_plot_folder, "upregulated", "deg_subset_pos.csv", sep = "/"))
  write.csv2(deg_names_negative, file = paste(ora_plot_folder, "downregulated", "deg_subset_neg.csv", sep = "/"))

  # define background gene universe, either internal (all measured/measurable genes) or external (package function built-in)
  if (use_internal_universe) {
    universe <- names(dea_ids_lfc)
    go_universe = unique(sort(as.data.frame(org.Hs.egGO)$gene_id))

    fgsea_compare_DEG_GO_universe(
      deg_universe = universe,
      go_universe = go_universe)
  } else {
    universe <- NULL
  }




  # TODO build in msigdb
  # msigdb.hs.h <- get_msigdb_term2gene(collection = 'h')
  # msigdb.hs.c2.cp <- get_msigdb_term2gene(collection = 'c2', subcollection = 'CP')





  # Gene Set Enrichment Analyses (GSEA. full DEA list)
  gsea_results <- run_fgsea_gsea(
    dea_ids_lfc, min_gene_set_size_gsea, max_gene_set_size_gsea,
    pval_cutoff, p_adjust_method
  )
  # Over Representation Analyses (ORA, subset DEG list)
  ora_results_positive <- run_fgsea_ora(
    deg_ids_positive, posneg = "upregulated",
    min_gene_set_size_ora, max_gene_set_size_ora,
    pval_cutoff, p_adjust_method, qval_cutoff, universe
  )
  ora_results_negative <- run_fgsea_ora(
    deg_ids_negative, posneg = "downregulated",
    min_gene_set_size_ora, max_gene_set_size_ora,
    pval_cutoff, p_adjust_method, qval_cutoff, universe
  )

  # plot results
  all_results <- c(gsea_results, ora_results_positive, ora_results_negative)
  for (res_name in names(all_results)) {
    plot_fgsea_result(
      res = all_results[res_name][[1]],
      res_name = res_name,
      plot_n_category,
      dea_ids_lfc,
      gsea_plot_folder,
      ora_plot_folder
    )
  }
}

# TODO optimize GSEA/ORA parameters
run_fgsea(
  output_dir = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/fgsea_GOuniverse_treeplot_KEGG_pvalCutoff=1/",
  dea_result_file = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/2022-11-01 14-02-50/integrated/BL_A + BL_C/DE_analysis/sample_markers/method=DESeq2-pct1=BL_C-pct2=BL_A - (nz-)p-val st 0.05.csv",
  cellRanger_ensembl_features = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Ensembl genes.tsv",
  pval_cutoff = 1
)

### TODO SET FUNCTIONS TO UTILS ###
#' gene_names: get gene names in character vector
#' return: gene names in character vector
date2gene <- function(gene_names) {
  ## if .csv was saved in Excel, some gene names e.g. SEP2 become dates
  # get index of these genes
  excel_genes_ind <- grepl(x = gene_names, pattern = "^[A-Z][a-z]{2}/[0-9]{2}$")
  # get these genes
  excel_genes <- gene_names[excel_genes_ind]
  # convert back to original names from date names
  gene_names[excel_genes_ind] <- sapply(X = excel_genes, USE.NAMES = FALSE, FUN = function(X) {
    ss <- unlist(strsplit(X, split = '/'))
    paste0(toupper(ss[1]), as.integer(ss[2]))
  })
  return(unlist(gene_names))
}
### END SET FUNCTIONS TO UTILS


fgsea_examine_unmapped_genes <- function(gene_names, dea_ensembl_ids, bitr_ensembl_ids) {
  dir.create("gene_mapping_bitr")
  d <- dea_ensembl_ids
  b <- bitr_ensembl_ids

  ratio_mapped <- round(length(b)/(length(d)-sum(is.na(d)))*100, digits = 2)

  genes_not_mapped <- gene_names[d %in% d[!d %in% b]]
  write.csv2(x = genes_not_mapped, file = "gene_mapping_bitr/genes_not_mapped.csv")

  rank_genes_not_mapped <- which(d %in% d[!d %in% b])
  png("gene_mapping_bitr/distribution_rank_genes_not_mapped.png")
  hist(rank_genes_not_mapped, breaks = 200, main = paste0("bitr rank unmapped genes - %mapped: ", ratio_mapped))
  dev.off()
}

fgsea_compare_DEG_GO_universe <- function(deg_universe, go_universe) {
  overlapping_gene_ids <- go_universe[go_universe %in% deg_universe]
  gene_ids_not_in_go <- go_universe[!go_universe %in% deg_universe]

  overlapping_gene_ids <- bitr(geneID = overlapping_gene_ids, OrgDb = org.Hs.eg.db, fromType = "ENTREZID", toType = "ENSEMBL")
  gene_ids_not_in_go <- bitr(geneID = gene_ids_not_in_go, OrgDb = org.Hs.eg.db, fromType = "ENTREZID", toType = "ENSEMBL")

  overlapping_gene_names <- plyr::mapvalues(
    x = overlapping_gene_ids$ENSEMBL,
    from = ensembl_genes_ids$ensembl_gene_id,
    to = ensembl_genes_ids$ensembl_gene_name,
    warn_missing = FALSE)
  gene_names_not_in_go <- plyr::mapvalues(
    x = gene_ids_not_in_go$ENSEMBL,
    from = ensembl_genes_ids$ensembl_gene_id,
    to = ensembl_genes_ids$ensembl_gene_name,
    warn_missing = FALSE)
  write.csv2(overlapping_gene_names, file = "DEG_GO_universe_overlapping_genes.csv")
  write.csv2(gene_names_not_in_go, file = "DEG_GO_universe_nonOverlapping_genes.csv")
}

run_fgsea_gsea <- function(
    dea_ids_lfc, min_gene_set_size_gsea, max_gene_set_size_gsea,
    pval_cutoff, p_adjust_method) {

  # initialize results
  gsea_results <- list()

  gene_ontology_types <- c("GO-BP" = "BP",
                           "GO-MF" = "MF",
                           "GO-CC" = "CC")
  for (i in 1:length(gene_ontology_types)) {
    res_gse_go <- clusterProfiler::gseGO(
      geneList = dea_ids_lfc,
      ont = gene_ontology_types[i],
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      exponent = 1, # default: 1
      minGSSize = min_gene_set_size_gsea, # default: 10
      maxGSSize = max_gene_set_size_gsea, # default: 500
      pvalueCutoff = pval_cutoff, # default: 0.05
      eps = 0, # default: 1e-10
      verbose = TRUE,
      pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
      seed = TRUE,
      by = "fgsea"
    )
    gsea_results[[paste0("GSEA-", names(gene_ontology_types)[i])]] <- res_gse_go
  }

  # Kyoto Encyclopedia of Genes and Genomes (KEGG)
  res_gse_kegg <- clusterProfiler::gseKEGG(
    geneList = dea_ids_lfc,
    organism = "hsa",
    keyType = "kegg", # "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
    exponent = 1, # default: 1
    minGSSize = min_gene_set_size_gsea, # default: 10
    maxGSSize = max_gene_set_size_gsea, # default: 500
    pvalueCutoff = pval_cutoff, # default: 0.05
    pAdjustMethod = p_adjust_method,
    eps = 0, # default: 1e-10
    verbose = TRUE,
    use_internal_data = FALSE, # default: FALSE
    seed = TRUE,
    by = "fgsea"
  )
  gsea_results[["GSEA-KEGG"]] <- res_gse_kegg

  # WikiPathways
  res_gse_wp <- clusterProfiler::gseWP(
    geneList = dea_ids_lfc,
    organism = "Homo sapiens",
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size_gsea,
    maxGSSize = max_gene_set_size_gsea
  )
  gsea_results[["GSEA-WP"]] <- res_gse_wp

  # Reactome Pathway Analysis
  res_gse_rp <- ReactomePA::gsePathway(
    geneList = dea_ids_lfc,
    organism = "human",
    exponent = 1, # default: 1
    minGSSize = min_gene_set_size_gsea, # default: 10
    maxGSSize = max_gene_set_size_gsea, # default: 500
    pvalueCutoff = pval_cutoff, # default: 0.05
    eps = 0, # default: 1e-10
    verbose = TRUE,
    pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
    seed = TRUE,
    by = "fgsea"
  )
  gsea_results[["GSEA-RP"]] <- res_gse_rp

  # DiseaseOntology (DO)
  res_gse_do <- gseDO(
    geneList = dea_ids_lfc,
    exponent = 1, # default: 1
    minGSSize = min_gene_set_size_gsea, # default: 10
    maxGSSize = max_gene_set_size_gsea, # default: 500
    pvalueCutoff = pval_cutoff, # default: 0.05
    verbose = TRUE,
    pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
    seed = TRUE,
    by = "fgsea"
  )
  gsea_results[["GSEA-DO"]] <- res_gse_do

  # Disease Gene Network (DGN)
  res_gse_dgn <- gseDGN(
    geneList = dea_ids_lfc,
    exponent = 1, # default: 1
    minGSSize = min_gene_set_size_gsea, # default: 10
    maxGSSize = max_gene_set_size_gsea, # default: 500
    pvalueCutoff = pval_cutoff, # default: 0.05
    verbose = TRUE,
    pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
    seed = TRUE,
    by = "fgsea"
  )
  gsea_results[["GSEA-DGN"]] <- res_gse_dgn

  return(gsea_results)
}

run_fgsea_ora <- function(
    deg_names, posneg, min_gene_set_size_ora, max_gene_set_size_ora,
    pval_cutoff, p_adjust_method, qval_cutoff, universe) {

  # initialize results
  ora_results <- list()

  go_universe = unique(sort(as.data.frame(org.Hs.egGO)$gene_id))
  gene_ontology_types <- c("GO-BP" = "BP",
                           "GO-MF" = "MF",
                           "GO-CC" = "CC")
  for (i in 1:length(gene_ontology_types)) {
    # DEG universe
    res_enrich_go <- clusterProfiler::enrichGO(
      gene = deg_names,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = gene_ontology_types[i],
      pvalueCutoff = pval_cutoff,
      pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
      universe = universe,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gene_set_size_ora, # default: 10
      maxGSSize = max_gene_set_size_ora, # default: 500
      readable = FALSE, # default: FALSE
      pool = FALSE # default: FALSE
    )
    ora_results[[paste0("ORA-", names(gene_ontology_types)[i], "-DEGuniverse-", posneg)]] = res_enrich_go

    # GO gene universe
    res_enrich_go <- clusterProfiler::enrichGO(
      gene = deg_names,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = gene_ontology_types[i],
      pvalueCutoff = pval_cutoff,
      pAdjustMethod = p_adjust_method, # default: "BH" (Benjamini-Hochberg)
      universe = go_universe,
      qvalueCutoff = qval_cutoff,
      minGSSize = min_gene_set_size_ora, # default: 10
      maxGSSize = max_gene_set_size_ora, # default: 500
      readable = FALSE, # default: FALSE
      pool = FALSE # default: FALSE
    )
    ora_results[[paste0("ORA-", names(gene_ontology_types)[i], "-GOuniverse-", posneg)]] = res_enrich_go
  }

  # Kyoto Encyclopedia of Genes and Genomes (KEGG)
  res_enrich_kegg <- enrichKEGG(
    gene = deg_names,
    organism = "hsa",
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    universe = universe,
    minGSSize = min_gene_set_size_ora, # default: 10
    maxGSSize = max_gene_set_size_ora, # default: 500
    qvalueCutoff = qval_cutoff, # default: 0.2
    use_internal_data = FALSE
  )
  ora_results[[paste0("ORA-KEGG-", posneg)]] = res_enrich_kegg

  # WikiPathways
  res_enrich_wp <- enrichWP(
    gene = deg_names,
    organism = "Homo sapiens",
    universe = universe,
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size_ora,
    maxGSSize = max_gene_set_size_ora,
    qvalueCutoff = qval_cutoff)
  ora_results[[paste0("ORA-WP-", posneg)]] = res_enrich_wp

  # Reactome Pathway
  res_enrich_rp <- enrichPathway(
    gene = deg_names,
    organism = "human",
    universe = universe,
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size_ora,
    maxGSSize = max_gene_set_size_ora,
    qvalueCutoff = qval_cutoff)
  ora_results[[paste0("ORA-RP-", posneg)]] = res_enrich_rp

  # Disease Ontology
  res_enrich_do <- enrichDO(
    gene = deg_names,
    universe = universe,
    ont = "DO",
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size_ora,
    maxGSSize = max_gene_set_size_ora,
    qvalueCutoff = qval_cutoff,
    readable = FALSE)
  ora_results[[paste0("ORA-DO-", posneg)]] = res_enrich_do

  # Disease Gene Network
  res_enrich_dgn <- enrichDGN(
    gene = deg_names,
    universe = universe,
    pvalueCutoff = pval_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gene_set_size_ora,
    maxGSSize = max_gene_set_size_ora,
    qvalueCutoff = qval_cutoff,
    readable = FALSE)
  ora_results[[paste0("ORA-DGN-", posneg)]] = res_enrich_dgn

  return(ora_results)
}

plot_fgsea_result <- function(
    res, res_name, plot_n_category, dea_ids_lfc,
    gsea_plot_folder, ora_plot_folder) {
  message("plotting: ", res_name)
  # get db folder name by res name
  fgsea_dbs <- c(
    "gene_ontology" = "GO",
    "disease_ontology" = "DO",
    "disease_gene_network" = "DGN",
    "wiki_pathways" = "WP",
    "reactome_pathways" = "RP",
    "KEGG" = "KEGG")
  db_name <- strsplit(res_name, "-")[[1]][2]
  db_folder <- names(fgsea_dbs)[which(fgsea_dbs %in% db_name)]

  # prep results
  res_r <- setReadable(res, 'org.Hs.eg.db', 'ENTREZID')
  res_r_pt <- pairwise_termsim(res_r)
  res_df <- as.data.frame(res)

  if (class(res) == "enrichResult") {
    # set plot_folder name
    ora_posneg_folder <- tail(strsplit(res_name, "-")[[1]], n = 1)
    plot_folder <- paste(ora_plot_folder, ora_posneg_folder, db_folder, sep = "/")
    if (grepl("^GO:", as.data.frame(res)$ID[1])) {
      plot_folder <- paste(plot_folder, res@ontology, sep = "/")

      if (grepl("-GOuniverse-", res_name)) {
        plot_folder <- paste(plot_folder, "GOuniverse", sep = "/")
      } else {
        plot_folder <- paste(plot_folder, "DEGuniverse", sep = "/")
      }
    }
    dir.create(plot_folder, recursive = TRUE)

    # check heatmap of terms with associated genes and associated log fold change
    p <- enrichplot::heatplot(res_r, foldChange = dea_ids_lfc, showCategory = plot_n_category) + ggplot2::theme_bw()
    p + ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggplot2::ggsave(paste(plot_folder, "heatmap_genes+lfc.png", sep = "/"), width = 30, height = 20, units = "cm")

    # GOplot specifically for Gene Ontology ORA
    if (grepl("^GO:", as.data.frame(res)$ID[1])) {
      enrichplot::goplot(res) + ggplot2::theme_bw()
      ggplot2::ggsave(paste(plot_folder, "goplot.png", sep = "/"), width = 30, height = 20, units = "cm")
    }
  }

  if (class(res) == "gseaResult") {
    # set plot_folder name
    plot_folder <- paste(gsea_plot_folder, db_folder, sep = "/")
    if (grepl("^GO:", as.data.frame(res)$ID[1])) {
      plot_folder <- paste(plot_folder, res@setType, sep = "/")
    }
    dir.create(plot_folder, recursive = TRUE)

    enrichplot::ridgeplot(res) + ggplot2::theme_bw()
    ggplot2::ggsave(paste(plot_folder, "ridgeplot.png", sep = "/"), width = 30, height = 20, units = "cm")
    p <- enrichplot::dotplot(
      res,
      showCategory = plot_n_category/2,
      split=".sign",
      label_format = Inf)
    p + ggplot2::theme_bw()
    p + facet_grid(.~.sign)
    ggplot2::ggsave(paste(plot_folder, "dotplot_split.png", sep = "/"), width = 30, height = 20, units = "cm")

    gsea_shown <- 0
    while(gsea_shown < plot_n_category) {
      enrichplot::gseaplot2(res, geneSetID = res_df$ID[(gsea_shown+1):(gsea_shown+10)], pvalue_table = TRUE) + ggplot2::theme_bw()
      ggplot2::ggsave(paste(plot_folder, paste0("gseaplot2_categories_", gsea_shown+1,"-", gsea_shown+10, ".png"), sep = "/"), width = 30, height = 20, units = "cm")
      gsea_shown <- gsea_shown + 10
    }
  }

  ## plots for GSEA and ORA results
  # term vs gene ratio & p-adjust & gene count
  enrichplot::dotplot(
    res,
    showCategory = plot_n_category,
    label_format = Inf) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "dotplot_top.png", sep = "/"), width = 30, height = 20, units = "cm")
  # check network of terms with associated genes and amount of genes
  enrichplot::cnetplot(res_r, showCategory = plot_n_category) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "cnetplot_terms+genes.png", sep = "/"), width = 30, height = 20, units = "cm")
  # check relations between amount of genes on terms and term association
  enrichplot::upsetplot(res, n = 10)
  ggplot2::ggsave(paste(plot_folder, "upsetplot_terms+ngenes.png", sep = "/"))
  # check network of associated terms with gene amounts but not gene associations
  enrichplot::emapplot(res_r_pt) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "emapplot_terms+ngenes.png", sep = "/"), width = 30, height = 20, units = "cm")
  # terms grouped and semantically summarized
  enrichplot::treeplot(res_r_pt, showCategory = plot_n_category)
  ggplot2::ggsave(paste(plot_folder, "treeplot.png", sep = "/"), width = 30, height = 20, units = "cm")
}

#' Get TERM2GENE dataframe for custom MSigDB collection/subcollection terms.
#'
#' The molecular signatures database (MSigDB) contains curated collections and
#' subcollection of gene sets. These can be used for comparison in GSEA/ORA
#' analyses.
#'
#' @param organism Either Homo Sapiems 'hs' or Mus Musculus 'mm'
#' @param collection The collection(s) must be one or more from
#' msigdb::listCollections()
#' @param subcollection The subcollection(s) must be one or more from
#' msigdb::listSubCollections()
#' @param id Either gene names 'SYM' or Entrez IDs 'EZID'
#' @param version 'latest' or one of msigdb::getMsigdbVersions()
#'
#' @return Dataframe with 'term' and 'gene' columns.
#'
#' @import msigdb
#'
#' @export
get_msigdb_term2gene <- function(
    organism = 'hs',
    collection = 'h',
    subcollection = c(),
    id = 'EZID',
    version = 'latest') {

  if (version == "latest") {
    version <- msigdb::getMsigdbVersions()[1]
  }

  msigdb <- msigdb::getMsigdb(org = organism, id = id, version = version)
  msigdb.sub <- msigdb::subsetCollection(msigdb, collection, subcollection)

  # list: name = term, values = EntrezIDs
  msigdb.sub = GSEABase::geneIds(msigdb.sub)

  # create term to gene dataframe
  term2gene <- DataFrame()
  for (name in names(msigdb.sub)) {
    term2gene <- rbind(term2gene, as.data.frame(cbind(name, msigdb.sub[[name]])))
  }
  colnames(term2gene) <- c("term", "gene")

  return(term2gene)
}
