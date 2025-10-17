# Introduction
Single-cell RNA sequencing workflow used at the Erasmus Medical Center, department of Psychiatry in the lab of Steven Kushner and Femke de Vrij. 
This repository is linked with github.com/kushnerlab/scRNAseqR.

Preprint (BioRxiv) paper: https://www.biorxiv.org/content/10.1101/2022.08.25.505166v1
Project data shared with UCSC Cell Browser: https://cells.ucsc.edu/?ds=ipsc-astrocyte-neuron

Using custom fork of Seurat v4.3.0 at github.com/kushnerlab/Seurat.

# Installation
install.packages(c("devtools", "BiocManager"))

BiocManager::install(c("SingleCellExperiment", "SingleR", "scuttle", "enrichplot", "msigdb", "org.Hs.eg.db", "EnhancedVolcano", "HDF5Array", "scDblFinder", "glmGamPoi"))

devtools::install_github('cole-trapnell-lab/monocle3')

devtools::install_github('immunogenomics/presto')

devtools::install_github("kushnerlab/scRNAseqR")

remotes::install_github('satijalab/seurat-wrappers')

Kriegstein reference data used for annotation: https://cells.ucsc.edu/?ds=dev-brain-regions+wholebrain
- 'Info & Download' -> 'Data Download' tab -> 'exprMatrix.tsv.gz' & 'meta.tsv'
-- put in project/data/Kriegstein/ folder

# Usage
See workflow.py for reproducing our results, this simultaneously serves as to how to work with the package.

# Notes
Louvain (instead of Leiden) algorithm in package for now, Leiden algorithm user for data clustering in paper.
