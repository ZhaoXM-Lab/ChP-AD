#!/usr/bin/env Rscript
# Infer high-confidence peak-gene associations within 250 kb.

library(ArchR)
library(GenomicRanges)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/epi")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
project <- loadArchRProject(file.path(proj_dir, "data/ATAC/output/ccres/ArchR_peaks"))
project <- addPeak2GeneLinks(project, useMatrix="GeneIntegrationMatrix", reducedDims="IterativeLSI", maxDist=250000)
links <- getPeak2GeneLinks(project, corCutOff=0.5, FDRCutOff=1e-4, resolution=1, returnLoops=FALSE)
genes <- metadata(links)$geneSet
peaks <- metadata(links)$peakSet
tab <- as.data.frame(links)
tab$gene <- genes$name[tab$idxRNA]
tab$peak <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))[tab$idxATAC]
marker <- read.csv(file.path(proj_dir, "data/ATAC/output/ccres/specific_ccres.csv"))
tab <- inner_join(tab, unique(marker[, c("peak", "cell_type")]), by="peak")
write.csv(distinct(tab, cell_type, peak, gene, .keep_all=TRUE), file.path(out, "archr_links.csv"), row.names=FALSE)
saveRDS(getGenes(project), file.path(out, "genes.rds"))
saveArchRProject(project, outputDirectory=file.path(out, "ArchR_links"), load=FALSE)
