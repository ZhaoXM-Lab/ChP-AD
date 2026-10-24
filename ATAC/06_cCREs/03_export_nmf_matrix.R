#!/usr/bin/env Rscript
# Aggregate accessibility by cell type and diagnosis for NMF.

library(ArchR)
library(Matrix)
library(SummarizedExperiment)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/ccres/nmf")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
project <- loadArchRProject(file.path(proj_dir, "data/ATAC/output/ccres/ArchR_peaks"))
se <- getMatrixFromProject(project, useMatrix="PeakMatrix")
m <- assay(se)
meta <- as.data.frame(getCellColData(project))[colnames(m), ]
keep <- !meta$cell_type %in% c("Neuron", "OPC", "CD8+T")
m <- m[, keep, drop=FALSE]
meta <- meta[keep, , drop=FALSE]
group <- factor(paste(meta$cell_type, meta$disease, sep="|"))
design <- sparse.model.matrix(~ group - 1)
avg <- m %*% design %*% Diagonal(x=1/as.numeric(table(group)))
colnames(avg) <- levels(group)
peaks <- getPeakSet(project)
rownames(avg) <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
write.table(as.matrix(avg), file.path(out, "accessibility.tsv"), sep="\t", quote=FALSE, col.names=NA)
