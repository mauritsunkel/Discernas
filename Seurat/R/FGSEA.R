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
set.seed(seed)
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



# TODO build in that ORA runs with separate + and - lfc analysis, GSEA with complete list, comparing biological samples

# TODO optimize GSEA/ORA parameters
run_fgsea(
  output_dir = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/fgsea_test/",
  dea_result_file = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/results/Pipe_SCTv2_23-06 cluster_level_selection/integrated/BL_A + BL_C/after_selection_old/DE_analysis/sample_markers/method=DESeq2-pct1=BL_A-pct2=BL_C - (nz-)p-val st 0.05.csv",
  cellRanger_ensembl_features = "C:/Users/mauri/Desktop/Single Cell RNA Sequencing/Seurat/data/Ensembl genes.tsv"
)





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
    "wiki_pathways", "reactome_pathways")
  for (db in fgsea_dbs) {
    dir.create(paste(gsea_plot_folder, db, sep = "/"), recursive = TRUE)
    dir.create(paste(ora_plot_folder, "positive", db, sep = "/"), recursive = TRUE)
    dir.create(paste(ora_plot_folder, "negative", db, sep = "/"), recursive = TRUE)
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

  ## prep ORA input
  # subset DEA result to DEG result based on log fold-change threshold
  deg_names_positive <- names(dea_ids_lfc[dea_ids_lfc >= lfc_threshold])
  deg_names_negative <- names(dea_ids_lfc[dea_ids_lfc <= -lfc_threshold])

  # define background gene universe, either internal (all measured/measurable genes) or external (package function built-in)
  if (use_internal_universe) {
    universe <- names(dea_ids_lfc)
  } else {
    universe <- NULL
  }

  # Gene Set Enrichment Analyses (GSEA. full DEA list)
  gsea_results <- run_fgsea_gsea(
    dea_ids_lfc, min_gene_set_size_gsea, max_gene_set_size_gsea,
    pval_cutoff, p_adjust_method)

  # Over Representation Analyses (ORA, subset DEG list)
  ora_results_positive <- run_fgsea_ora(
    deg_names_positive, posneg = "positive",
    min_gene_set_size_ora, max_gene_set_size_ora,
    pval_cutoff, p_adjust_method, qval_cutoff, universe
  )
  ora_results_negative <- run_fgsea_ora(
    deg_names_negative, posneg = "negative",
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
      ora_plot_folder)
  }
}

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

  ## BUGFIX 1
  # library(R.utils)
  # R.utils::setOption("clusterProfiler.download.method","auto")
  ## BUGFIX 2: update R --> update Bioconductor --> update clusterProfiler
  ## BUGFIX 3: restart R (lol)
  # res_gsekegg <- clusterProfiler::gseKEGG(
  #   geneList = dea_ids_lfc,
  #   organism = "hsa",
  #   keyType = "kegg", # "kegg", 'ncbi-geneid', 'ncib-proteinid' and 'uniprot'
  #   exponent = 1, # default: 1
  #   minGSSize = min_gene_set_size_gsea, # default: 10
  #   maxGSSize = max_gene_set_size_gsea, # default: 500
  #   pvalueCutoff = pval_cutoff, # default: 0.05
  #   pAdjustMethod = p_adjust_method,
  #   eps = 0, # default: 1e-10
  #   verbose = TRUE,
  #   use_internal_data = FALSE, # default: FALSE
  #   seed = TRUE,
  #   by = "fgsea"
  # )
  # TODO check clusterProfiler::download_KEGG("hsa") # to skip gse part and just try download hsa pathways from KEGG
  # TODO check if bitr_kegg() needed for keyType
  ## TODO process error messages
  # Reading KEGG annotation online:
  #
  #   fail to download KEGG data...
  # Error in download.KEGG.Path(species) :
  #   'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...
  # In addition: Warning message:
  #   In utils::download.file(url, quiet = quiet, method = "libcurl",  :
  #                             URL 'https://rest.kegg.jp/link/hsa/pathway': status was 'Failure when receiving data from the peer'

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

  # DEPRECATED: define GO specific gene universe
  # go_gene_list = unique(sort(as.data.frame(org.Hs.egGO)$gene_id))
  gene_ontology_types <- c("GO-BP" = "BP",
                           "GO-MF" = "MF",
                           "GO-CC" = "CC")
  for (i in 1:length(gene_ontology_types)) {
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
    ora_results[[paste0("ORA-", names(gene_ontology_types)[i], "-", posneg)]] = res_enrich_go
  }

  ## Kyto Encyclopedia of Genes and Genomes
  # res_enrich_kegg <- enrichKEGG(
  #   gene = deg_names,
  #   organism = "hsa",
  #   pvalueCutoff = pval_cutoff,
  #   pAdjustMethod = p_adjust_method,
  #   universe = universe,
  #   minGSSize = min_gene_set_size_ora, # default: 10
  #   maxGSSize = max_gene_set_size_ora, # default: 500
  #   qvalueCutoff = qval_cutoff, # default: 0.2
  #   use_internal_data = FALSE
  # )
  # TODO try KEGG.REST (as above)/Pathway Commons/reactome.db/fixing error messages
  # Reading KEGG annotation online:
  #
  #   fail to download KEGG data...
  # Error in download.KEGG.Path(species) :
  #   'species' should be one of organisms listed in 'http://www.genome.jp/kegg/catalog/org_list.html'...
  # In addition: Warning message:
  #   In utils::download.file(url, quiet = quiet, method = "libcurl",  :
  #                             URL 'https://rest.kegg.jp/link/hsa/pathway': status was 'Failure when receiving data from the peer'

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
    "gene_ontology" = "GO", "disease_ontology" = "DO",
    "disease_gene_network" = "DGN",
    "wiki_pathways" = "WP", "reactome_pathways" = "RP")
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

    gsea_shown <- 0
    while(gsea_shown < plot_n_category) {
      enrichplot::gseaplot2(res, geneSetID = res_df$ID[(gsea_shown+1):(gsea_shown+10)], pvalue_table = TRUE) + ggplot2::theme_bw()
      ggplot2::ggsave(paste(plot_folder, paste0("gseaplot2_categories_", gsea_shown+1,"-", gsea_shown+10, ".png"), sep = "/"), width = 30, height = 20, units = "cm")
      gsea_shown <- gsea_shown + 10
    }
  }

  ## plots for GSEA and ORA results
  # term vs gene ratio & p-adjust & gene count
  enrichplot::dotplot(res, showCategory = plot_n_category) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "dotplot.png", sep = "/"), width = 30, height = 20, units = "cm")
  # check network of terms with associated genes and amount of genes
  enrichplot::cnetplot(res_r, showCategory = plot_n_category) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "cnetplot_terms+genes.png", sep = "/"), width = 30, height = 20, units = "cm")
  # check relations between amount of genes on terms and term association
  enrichplot::upsetplot(res, n = 10)
  ggplot2::ggsave(paste(plot_folder, "upsetplot_terms+ngenes.png", sep = "/"))
  # check network of associated terms with gene amounts but not gene associations
  enrichplot::emapplot(res_r_pt) + ggplot2::theme_bw()
  ggplot2::ggsave(paste(plot_folder, "emapplot_terms+ngenes.png", sep = "/"), width = 30, height = 20, units = "cm")
}

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
  return(gene_names)
}
### END SET FUNCTIONS TO UTILS












# TODO upgrade R version and then BioCOnductor version then packages to test this function
enrichplot::treeplot(res_enrich_do_r2)



# TODO check how to get genes associated with a certain term from results object



# TODO try msigdb with custom GSEA/ORA?
enricher() # custom annotation, like MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)
GSEA() # custom annotation, like MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/human/collection_details.jsp#H)

## user input: visualize specific pathways
# browseKEGG() # KEGG pathways
# pathview::pathview() # general
# viewPathway() # reactome pathway analysis
