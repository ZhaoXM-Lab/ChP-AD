#!/usr/bin/env Rscript
# Prepare Signac and ArchR gene annotations from the matched GRCh38 GTF.

library(ArchR)
library(rtracklayer)
library(GenomicFeatures)
library(org.Hs.eg.db)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir,"data/reference")
gtf <- import(file.path(out,"genes.gtf.gz"))
GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
genome(gtf) <- "hg38"
saveRDS(gtf,file.path(out,"hg38_gene_annotations.rds"))
txdb <- makeTxDbFromGRanges(gtf)
gene <- createGeneAnnotation(TxDb=txdb,OrgDb=org.Hs.eg.db,annoStyle="ENSEMBL",singleStrand=TRUE)
saveRDS(gene,file.path(out,"archr_gene_annotations.rds"))
