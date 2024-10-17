
library(Libra)

data <- qs::qread("C:/SynologyDrive/Projects/scRNAseqR/results/sakshi_4/NSM-NS-NC-M/NSM-NS-NC-M.qs")
# SeuratObject::DefaultAssay(data) <- "RNA"
SeuratObject::DefaultAssay(data) <- "RNA"

comp_data <- data[, data$orig.ident %in% c("NS", "NC")]
comp_data <- comp_data[, comp_data$mapmycells_supercluster == "Neuron"]

meta_df <- data.frame(
  cell_type = comp_data$mapmycells_supercluster,
  replicate = comp_data$orig.ident, # TODO we have no replicates right?
  label = comp_data$orig.ident
)

# meta_df$cell_type

# mtrix <- Libra::to_pseudobulk(comp_data@assays$RNA$counts, meta_df)


DE <- Libra::run_de(
  input = comp_data@assays$RNA$counts, # comp_data@assays$RNA$counts,
  meta = meta_df,
  de_family = "pseudobulk",
  de_method = "edgeR",
  de_type = "LRT",
  n_threads = 2,
  min_reps = 1)
DE
min(DE$p_val_adj)
hist(DE$p_val)


data("hagai_toy")
hagai_toy <- Seurat::UpdateSeuratObject(hagai_toy)
DE <- Libra::run_de(hagai_toy)
factor(hagai_toy$label)
factor(hagai_toy$cell_type)
factor(hagai_toy$replicate)
DE

colnames(hagai_toy@assays$RNA)




data
data$orig.ident
Seurat::Idents(data)

data

res3 <- Seurat::FindMarkers(data, ident.1 = 2, ident.2 = 3, assay = "RNA", test.use = "DESeq2")
identical(res$p_val, res2$p_val)


mean(res$p_val)
mean(res2$p_val)
mean(res$p_val_adj)
mean(res2$p_val_adj)


# TODO get old code to work?
## markers, conserved markers, condition markers, sample markers
## Seurat defaults: wilcox --> presto
## Seurat defaults: DESeq2
# TODO via Muscat: check pseudobulk with edge R sum.counts (see Muscat paper)
# TODO read paper series to motivate choice for or against pseudobulk/scDE (MAST/DESeq2)
# TODO can still try aggregation pre-pseudobulk, but then there's no DE test afterwards







## Junttilla et al 2022 (summarizes both Zimmerman et al 2021/2022 and Squair et al 2021!)
### pseudobulk (best: ROTS with sum aggregation) > pseudobulk mean > mixed models (best: MAST_RE) > naive models > latent variable models
#### AUROC suggests that p-values are generally in order for all methods, however large amount of false positives for naive/latent variable models
### https://github.com/elolab/multisubjectDSanalysis
### https://academic.oup.com/bib/article/23/5/bbac286/6649780?login=true



## Zimmerman 2021
### propose generalized linear mixed models with a random effect for individual (GLMM-RE)
### pseudobulk: aggregate of cell-type specific expression within individual by sum (favorite) or mean, controls for zero-inflation and within-sample correlation, but are conservative, lose information, and lose power
#### pseudobulking should be balanced, i.e. comparable number of cells summed/averaged
#### implemented in DESeq2?
#### after custom pseudobulk, should to DESeq2 or other DEA?
### MAST (two-part Hurdle): type 1 error remains inflated if random effect for individual is not accounted for
#### MAST with fixed-effect term for individual may be considered next to MAST-RE, especially for low sample size (see ref 33)
### ComBat: batch effect correction prior to DE within cell type, to improve type 1 error
### with low amount of samples (individuals), pseudobulking is a lot better for power during statistics

## Murphy & Skene 2022: a response on Zimmerman 2021
### MCC (Matthew's Correlation Coefficient): to take into account type-1 and type-2 error simultaneously instead of in isolation
### pseudobulk mean outperforms everything, maybe a safe choice since this was consistently performant in Zimmerman 2021 as well
### Tweedie: GLMM als consistent and performant in both!
### weird that MAST-RE drops eavily in this paper
### pseudobulk sum should outperfrom pseudobulk mean: IMPORTANT to have a normalisation step when using sum > mean, sum can account for intra-sample variance which is lost with averaging, however this should be tested especially in imbalanced datasets
#### mean is the safe option, for now, because unbalanced datasets in single-cell are frequent and the performance loss is minimal comparatively
#### mean is also a better option for fewer samples and fewer cells per sample

## Zimmerman 2022: a response on response of Murphy & Skene 2022
### they have significantly better statistical understanding, they argue well and correctly and counter a lot of the issued points
### mix interpretations from both papers to be safe


## Squair et al 2021: Confronting false discoveries in single-cell differential expression
### edgeR-LRT for pseudobulk?
### Negative Binomial Generalized Linear Mixed Model with offset-LRT: highest performance but significant runtine and computational resource
### "the central principle underlying valid DE analysis is the ability of statistical methods to account for the intrinsic variability of biological replicates"
### they propose pseudobulk methods (not any specific pseudobulk method)
#### see their R package
#### runtime and memory can definitely play a role in choosing an applicable DE method

## Muscat paper
### edgeR.sum.counts proposed as best pseudobulk method in Muscat


# DEVNOTE: pseudobulking or correcting for intra-sample variation is needed for correct single-cell approaches and pseudobulking, however, you NEED sample replicates for this
## pseudobulk: mean is the safe choice, choose sum for optimal, if properly setup balanced comparison - edgeR-LRT or edgeR.sum.counts - generally speed and memory efficient
## scDE: Tweedie-GLMM safe, MAST-RE optimal in theory stats but lower in ML stats - generally speed and memory intensive
### or NB-GLMM with offset-LRT from R Libra package (Squair et al 2021) -> R Libra implements edgeR-LRT and others!
### MAST-RE custom/direct implmentation
#### https://github.com/satijalab/seurat/issues/7868
#### https://github.com/satijalab/seurat/issues/6056
#### https://github.com/satijalab/seurat/issues/3712
## pseudobulk Marijn Groningen: glmQLF custom implement per celltype
# Tuurlijk, dit is een loop die DE doet op elk celtype in een seurat object en even checkt of beide patientengroepen wel in de subset zitten.
#
# subjectContrastTable = unique(tcells@meta.data[,c("DonorID", "PatientGroup")])
# rownames(subjectContrastTable) = subjectContrastTable[,1]
#
# # Cell types and levels to use
# rm(resultsTable)
# for(celltype in levels(tcells)){
#   celltypeSubset = subset(tcells, idents = celltype)
#   genes.percent.expression <- rowMeans(celltypeSubset@assays$RNA@counts >0 )*100
#   pseudobulk = as.data.frame(AggregateExpression(celltypeSubset, assays = "RNA", group.by = "DonorID"))
#
#   subjects = str_replace(colnames(pseudobulk), "RNA\\.", "")
#
#   contrast = subjectContrastTable[subjects, 2]
#
#   subjects = data.frame(subjects, contrast)
#   rownames(subjects) = subjects$subjects
#
#   if(length(table(contrast)) == 2){
#     deObject = DGEList(pseudobulk, samples = subjects, group = contrast)
#
#     keep = filterByExpr(deObject)
#     # table(keep)
#     deObject = deObject[keep,]
#     deObject = calcNormFactors(deObject)
#
#     design = model.matrix(~contrast)
#     deObject = estimateDisp(deObject, design)
#
#     # use QLFT
#     fit = glmQLFit(deObject, design)
#     qlf = glmQLFTest(fit, coef=2)
#
#     res = qlf$table
#     res$celltype = celltype
#     res$gene = rownames(res)
#
#     res = res[order(res$PValue),]
#     res$p.adjust <- p.adjust(res$PValue, method = 'bonferroni')
#
#     if(exists("resultsTable")){
#       resultsTable = rbind(resultsTable, res)
#     }else{
#       resultsTable = res
#     }
#   }
# }
