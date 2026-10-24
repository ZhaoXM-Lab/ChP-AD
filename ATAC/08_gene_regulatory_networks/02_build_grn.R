#!/usr/bin/env Rscript
# Connect matched motifs through cell-specific cCREs and EPIs to target genes.

library(ArchR)
library(Seurat)
library(SummarizedExperiment)
library(Matrix)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/grn")
project <- loadArchRProject(file.path(out, "ArchR_motifs"))
matches <- getMatches(project, name="VierstraMotif")
mat <- assays(matches)$matches
peaks <- getPeakSet(project)
rownames(mat) <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
motif_names <- colnames(mat)
tf_names <- sub("_.*$", "", motif_names)
epi <- read.csv(file.path(proj_dir, "data/ATAC/output/epi/archr_links.csv"))
markers <- read.csv(file.path(out, "specific_tfs.csv"))
rna <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
expr <- GetAssayData(rna, assay="RNA", layer="data")
deg <- read.csv(file.path(proj_dir, "data/RNA/output/differential_expression/consensus_degs.csv"))
all <- list()
for (ct in setdiff(unique(epi$cell_type), c("Neuron","OPC","CD8+T"))) {
    edges <- list()
    ct_epi <- epi[epi$cell_type == ct, ]
    motif_idx <- which(tf_names %in% markers$TF[markers$cell_type == ct])
    for (j in motif_idx) {
        matched <- rownames(mat)[which(mat[, j] > 0)]
        z <- ct_epi[ct_epi$peak %in% matched, c("peak", "gene")]
        if (nrow(z)) { z$TF <- tf_names[j]; edges[[length(edges)+1]] <- z }
    }
    z <- unique(bind_rows(edges))
    if (!nrow(z)) next
    cells <- colnames(rna)[(if (ct %in% c("CPEC-1","CPEC-2")) rna$cell_subtype == ct else rna$cell_type == ct) & rna$disease == "AD"]
    z <- z[z$TF %in% rownames(expr) & z$gene %in% rownames(expr), ]
    pair <- unique(z[,c("TF","gene")])
    pair$correlation <- vapply(seq_len(nrow(pair)), function(i) suppressWarnings(cor(as.numeric(expr[pair$TF[i],cells]), as.numeric(expr[pair$gene[i],cells]))), numeric(1))
    pair <- pair[is.finite(pair$correlation), ]
    pair$weight <- abs(pair$correlation)
    pair$cell_type <- ct
    d <- deg[deg$cell_type == ct, ]
    pair$direction <- d$direction[match(pair$gene, d$gene)]
    pair$direction[is.na(pair$direction)] <- "NS"
    write.csv(z, file.path(out, paste0(ct,"_tf_peak_gene.csv")), row.names=FALSE)
    write.csv(pair, file.path(out, paste0(ct,"_grn.csv")), row.names=FALSE)
    mcl <- pair[pair$weight >= 0.1, c("TF","gene","weight")]
    mcl$TF <- paste0("TF_", mcl$TF)
    mcl$gene <- paste0("G_", mcl$gene)
    write.table(mcl, file.path(out, paste0(ct,".abc")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    all[[ct]] <- pair
}
write.csv(bind_rows(all), file.path(out, "all_grns.csv"), row.names=FALSE)
