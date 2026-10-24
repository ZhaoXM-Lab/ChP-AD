#!/usr/bin/env Rscript
# Visualize the batch-corrected ATAC embedding.

library(Seurat)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
input <- file.path(proj_dir,"data/ATAC/output/integration/atac_integrated.rds")
out <- file.path(proj_dir,"data/ATAC/output/umap")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(input)
x <- RunUMAP(x,reduction="harmony",dims=1:28,seed.use=1)
saveRDS(x,file.path(out,"atac_umap.rds"))
