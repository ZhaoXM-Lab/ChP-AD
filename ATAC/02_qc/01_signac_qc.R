#!/usr/bin/env Rscript
# Build Signac objects and filter ATAC nuclei.

library(Signac)
library(Seurat)
library(GenomicRanges)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/qc")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
samples <- read.delim(file.path(proj_dir, "ATAC/00_reference/samples.tsv"), colClasses="character")
samples <- samples[samples$include == "1", ]
stopifnot(nrow(samples) > 0)
annotations <- readRDS(file.path(proj_dir, "data/reference/hg38_gene_annotations.rds"))
objects <- list()
for (i in seq_len(nrow(samples))) {
    sample <- samples$sample[i]
    input <- file.path(proj_dir, samples$outs_dir[i])
    counts <- Read10X_h5(file.path(input, "filtered_peak_bc_matrix.h5"))
    meta <- read.csv(file.path(input, "singlecell.csv"), row.names=1)
    meta <- meta[colnames(counts), , drop=FALSE]
    assay <- CreateChromatinAssay(counts, sep=c(":", "-"), genome="hg38", fragments=file.path(input, "fragments.tsv.gz"), annotation=annotations)
    x <- CreateSeuratObject(assay, assay="peaks", meta.data=meta)
    x$sample <- sample
    x$disease <- samples$disease[i]
    x$barcode <- colnames(x)
    standard <- as.character(seqnames(granges(x))) %in% GenomeInfoDb::standardChromosomes(granges(x))
    x <- x[standard, ]
    x <- NucleosomeSignal(x)
    x <- TSSEnrichment(x, fast=FALSE)
    x$pct_reads_in_peaks <- x$peak_region_fragments / x$passed_filters * 100
    x$blacklist_ratio <- x$blacklist_region_fragments / x$peak_region_fragments
    x <- subset(x, subset=nCount_peaks > 1000 & nCount_peaks < 30000 & TSS.enrichment >= 2 & nucleosome_signal < 4 & blacklist_ratio < 0.05 & pct_reads_in_peaks > 15)
    x <- RenameCells(x, new.names=paste0(sample, "#", x$barcode))
    objects[[sample]] <- x
}
saveRDS(objects, file.path(out, "signac_samples.rds"))
