library(EnsDb.Hsapiens.v86)
library(Signac)


EnsDb.Hsapiens.v86

# extract gene annotations from EnsDb
annotations <- signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg38
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(seur) <- annotations

