#!/usr/bin/env Rscript
# Create the full-gene Seurat reference from aligned sparse inputs.

library(Seurat)
library(Matrix)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
species <- Sys.getenv("SPECIES", "human")
input <- file.path(proj_dir, "data/RNA/output/seurat_input", species)
output <- file.path(proj_dir, "data/RNA/output/annotation", paste0(species, "_final.rds"))
counts <- readMM(file.path(input, "matrix.mtx"))
rownames(counts) <- readLines(file.path(input, "genes.tsv"))
meta <- read.csv(file.path(input, "cells.csv"), row.names=1, check.names=FALSE)
colnames(counts) <- rownames(meta)
obj <- CreateSeuratObject(counts=as(counts, "CsparseMatrix"), meta.data=meta)
obj <- NormalizeData(obj, scale.factor=10000, verbose=FALSE)
VariableFeatures(obj) <- readLines(file.path(proj_dir, "data/RNA/output/integration", paste0(species, "_hvg.txt")))
for (key in c("X_pca_harmony", "X_umap")) {
    emb <- as.matrix(read.csv(file.path(input, paste0(key, ".csv")), row.names=1, check.names=FALSE))
    emb <- emb[colnames(obj), , drop=FALSE]
    nm <- if (key == "X_umap") "umap" else "harmony"
    colnames(emb) <- paste0(nm, "_", seq_len(ncol(emb)))
    obj[[nm]] <- CreateDimReducObject(embeddings=emb, key=paste0(nm, "_"), assay="RNA")
}
Idents(obj) <- "cell_subtype"
saveRDS(obj, output)
