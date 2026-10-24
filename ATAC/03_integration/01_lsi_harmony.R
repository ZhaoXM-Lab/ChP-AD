#!/usr/bin/env Rscript
# Recount a shared peak space and integrate ATAC LSI coordinates.

library(Signac)
library(Seurat)
library(GenomicRanges)
library(harmony)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/integration")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
objects <- readRDS(file.path(proj_dir, "data/ATAC/output/qc/singlet_samples.rds"))
peaks <- reduce(do.call(c, unname(lapply(objects, granges))))
peaks <- peaks[width(peaks) > 20 & width(peaks) < 10000]
objects <- lapply(objects, function(x) {
    m <- FeatureMatrix(fragments=Fragments(x), features=peaks, cells=colnames(x))
    x[["peaks"]] <- CreateChromatinAssay(m, fragments=Fragments(x), annotation=Annotation(x))
    x
})
x <- merge(objects[[1]], y=objects[-1])
x <- RunTFIDF(x)
x <- FindTopFeatures(x, min.cutoff="q5")
x <- RunSVD(x, n=50)
x <- RunHarmony(x, group.by.vars="sample", reduction.use="lsi", dims.use=2:30, project.dim=FALSE)
saveRDS(x, file.path(out, "atac_integrated.rds"))
