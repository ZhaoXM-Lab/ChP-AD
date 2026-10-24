#!/usr/bin/env Rscript
# Expand FDR 0.1 cell-specific cCREs by 5 kb and harmonize the reference build.

library(GenomicRanges)
library(rtracklayer)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
build <- Sys.getenv("LDSC_BUILD")
stopifnot(build %in% c("hg19", "hg38"))
out <- file.path(proj_dir, "data/ATAC/output/ldsc/beds")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
markers <- read.csv(file.path(proj_dir, "data/ATAC/output/ccres/ldsc_ccres_fdr01.csv"))
for (ct in setdiff(unique(markers$cell_type), c("Neuron","OPC","CD8+T"))) {
    gr <- GRanges(unique(markers$peak[markers$cell_type == ct]))
    start(gr) <- pmax(1, start(gr)-5000)
    end(gr) <- end(gr)+5000
    if (build == "hg19") {
        lifted <- liftOver(gr, import.chain(file.path(proj_dir,"data/reference/hg38ToHg19.over.chain")))
        gr <- unlist(lifted[elementNROWS(lifted) == 1], use.names=FALSE)
    }
    gr <- reduce(gr[as.character(seqnames(gr)) %in% paste0("chr",1:22)])
    export.bed(gr, file.path(out, paste0(ct,".bed")))
}
