#!/usr/bin/env Rscript
# Call reproducible 501-bp peaks and add peak and RNA integration matrices.

library(ArchR)
library(Seurat)
library(GenomicRanges)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/ccres")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
addArchRGenome("hg38")
addArchRThreads(threads=as.integer(Sys.getenv("THREADS", "8")))
project <- loadArchRProject(file.path(proj_dir, "data/ATAC/output/qc/ArchR_singlets"))
atac <- readRDS(file.path(proj_dir, "data/ATAC/output/annotation/atac_annotated.rds"))
meta <- atac[[]][getCellNames(project), , drop=FALSE]
stopifnot(!anyNA(meta$cell_subtype))
project$cell_subtype <- meta$cell_subtype
project$cell_type <- meta$cell_type
project <- addIterativeLSI(project, useMatrix="TileMatrix", name="IterativeLSI")
project <- addGroupCoverages(project, groupBy="cell_subtype")
project <- addReproduciblePeakSet(project, groupBy="cell_subtype", peakMethod="Macs2", pathToMacs2=Sys.getenv("MACS2", "macs2"), extendSummits=250, reproducibility="2", promoterRegion=c(2000,100), maxPeaks=300000, cutOff=0.05)
project <- addPeakMatrix(project)
rna <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
groups <- lapply(intersect(unique(project$cell_subtype), unique(rna$cell_subtype)), function(ct) list(ATAC=getCellNames(project)[project$cell_subtype == ct], RNA=colnames(rna)[rna$cell_subtype == ct]))
names(groups) <- intersect(unique(project$cell_subtype), unique(rna$cell_subtype))
project <- addGeneIntegrationMatrix(project, useMatrix="GeneScoreMatrix", matrixName="GeneIntegrationMatrix", reducedDims="IterativeLSI", seRNA=rna, groupRNA="cell_subtype", groupList=groups, addToArrow=TRUE)
saveRDS(getPeakSet(project), file.path(out, "union_peaks.rds"))
saveRDS(getMatrixFromProject(project, useMatrix="PeakMatrix"), file.path(out, "peak_matrix.rds"))
saveArchRProject(project, outputDirectory=file.path(out, "ArchR_peaks"), load=FALSE)
