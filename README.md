# Introduction
Single-cell RNA sequencing workflow used at the Erasmus Medical Center, department of Psychiatry at the Steven Kushner lab. 

Preprint (BioRxiv) paper: https://www.biorxiv.org/content/10.1101/2022.08.25.505166v1
Project data shared with UCSC Cell Browser: https://cells.ucsc.edu/?ds=ipsc-astrocyte-neuron

Using custom fork of Seurat v4.3.0 at github.com/kushnerlab/Seurat.

# Installation
install.packages(c("devtools", "BiocManager"))
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install(c("SingleCellExperiment", "SingleR", "scuttle", "enrichplot", "msigdb", "org.Hs.eg.db", "EnhancedVolcano", "HDF5Array"))
devtools::install_github('cole-trapnell-lab/monocle3')
devtools::install_github("kushnerlab/scRNAseqR")

Kriegstein reference data used for annotation: https://cells.ucsc.edu/?ds=dev-brain-regions+wholebrain
- 'Info & Download' -> 'Data Download' tab -> 'exprMatrix.tsv.gz' & 'meta.tsv'
-- put in project/data/Kriegstein/ folder

# Usage
See workflow.py for reproducing our results, this simultaneously serves as to how to work with the package.

# Notes
Louvain (instead of Leiden) algorithm in package for now, Leiden algorithm user for data clustering in paper.
