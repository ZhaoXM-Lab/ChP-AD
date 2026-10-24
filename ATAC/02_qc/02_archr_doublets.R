#!/usr/bin/env Rscript
# Retain the intersection of Signac QC cells and ArchR singlets.

library(ArchR)
library(Seurat)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/qc")
addArchRGenome("hg38")
addArchRThreads(threads=as.integer(Sys.getenv("THREADS", "8")))
samples <- read.delim(file.path(proj_dir, "ATAC/00_reference/samples.tsv"), colClasses="character")
samples <- samples[samples$include == "1", ]
files <- file.path(proj_dir, samples$outs_dir, "fragments.tsv.gz")
names(files) <- samples$sample
arrow_dir <- file.path(out, "arrows")
dir.create(arrow_dir, recursive=TRUE, showWarnings=FALSE)
setwd(arrow_dir)
arrows <- createArrowFiles(inputFiles=files, sampleNames=names(files), minTSS=2, minFrags=1000, addTileMat=TRUE, addGeneScoreMat=TRUE, geneAnnotation=readRDS(file.path(proj_dir,"data/reference/archr_gene_annotations.rds")))
addDoubletScores(input=arrows, k=10, knnMethod="UMAP", LSIMethod=1)
project <- ArchRProject(ArrowFiles=arrows, outputDirectory=file.path(out, "ArchR"), copyArrows=TRUE, geneAnnotation=readRDS(file.path(proj_dir,"data/reference/archr_gene_annotations.rds")))
project <- filterDoublets(project)
objects <- readRDS(file.path(out, "signac_samples.rds"))
keep <- intersect(unlist(lapply(objects, colnames), use.names=FALSE), getCellNames(project))
stopifnot(length(keep) > 0)
project <- project[keep, ]
project$disease <- samples$disease[match(project$Sample, samples$sample)]
objects <- lapply(objects, function(x) subset(x, cells=intersect(colnames(x), keep)))
saveRDS(objects, file.path(out, "singlet_samples.rds"))
saveArchRProject(project, outputDirectory=file.path(out, "ArchR_singlets"), load=FALSE)
