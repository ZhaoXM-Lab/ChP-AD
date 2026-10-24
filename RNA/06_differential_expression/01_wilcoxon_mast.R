#!/usr/bin/env Rscript
# Run unadjusted Wilcoxon and MAST tests in each major lineage and CPEC subtype.

library(Seurat)
library(MAST)
library(edgeR)
library(data.table)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/RNA/output/differential_expression")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
obj <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
sex <- readLines(file.path(proj_dir, "RNA/00_reference/human_sex_genes.txt"))
types <- c("CPEC", "Ependymal", "Astrocyte", "Oligodendrocyte", "Macrophage", "Endothelial", "Fibroblast", "Pericyte", "V-SMC", "M-SMC", "CPEC-1", "CPEC-2")
results <- list()
for (ct in types) {
    cells <- colnames(obj)[if (ct %in% c("CPEC-1", "CPEC-2")) obj$cell_subtype == ct else obj$cell_type == ct]
    x <- subset(obj, cells=cells)
    counts <- GetAssayData(x, assay="RNA", layer="counts")
    meta <- x[[]][colnames(counts), , drop=FALSE]
    stopifnot(all(c("AD", "WT") %in% meta$disease), identical(rownames(meta), colnames(counts)))
    detectable <- sapply(c("AD", "WT"), function(d) Matrix::rowMeans(counts[, meta$disease == d, drop=FALSE] > 0))
    genes <- rownames(counts)[apply(detectable, 1, max) >= 0.05 & !rownames(counts) %in% sex]
    x <- x[genes, ]
    x <- NormalizeData(x, normalization.method="LogNormalize", scale.factor=10000, verbose=FALSE)
    Idents(x) <- "disease"
    wil <- FindMarkers(x, ident.1="AD", ident.2="WT", test.use="wilcox", logfc.threshold=0, min.pct=0)
    wil$gene <- rownames(wil)
    wil$padj <- p.adjust(wil$p_val, "bonferroni")
    counts <- GetAssayData(x, assay="RNA", layer="counts")
    cpms <- edgeR::cpm(edgeR::calcNormFactors(edgeR::DGEList(counts=counts)))
    meta <- x[[]][colnames(counts), , drop=FALSE]
    meta$wellKey <- rownames(meta)
    meta$disease <- relevel(factor(meta$disease), "WT")
    feature <- data.frame(primerid=rownames(counts), row.names=rownames(counts))
    sca <- MAST::FromMatrix(log2(cpms + 1), meta, feature)
    fit <- MAST::zlm(~ disease, sca)
    dt <- summary(fit, doLRT="diseaseAD")$datatable
    mast <- merge(dt[contrast == "diseaseAD" & component == "H", .(primerid, pvalue=`Pr(>Chisq)`)],
                  dt[contrast == "diseaseAD" & component == "logFC", .(primerid, log2FC=coef)], by="primerid")
    mast$padj <- p.adjust(mast$pvalue, "bonferroni")
    write.csv(wil, file.path(out, paste0(ct, "_wilcoxon.csv")), row.names=FALSE)
    write.csv(mast, file.path(out, paste0(ct, "_mast.csv")), row.names=FALSE)
    w <- wil[wil$padj < 0.05 & abs(wil$avg_log2FC) > 0.3, ]
    m <- mast[mast$padj < 0.05 & abs(mast$log2FC) > 0.3, ]
    consensus <- merge(w[, c("gene", "avg_log2FC", "padj")], m[, c("primerid", "log2FC", "padj")], by.x="gene", by.y="primerid", suffixes=c("_wilcoxon", "_mast"))
    consensus <- consensus[sign(consensus$avg_log2FC) == sign(consensus$log2FC), ]
    consensus$cell_type <- ct
    consensus$direction <- ifelse(consensus$avg_log2FC > 0, "Up", "Down")
    results[[ct]] <- consensus
}
all <- do.call(rbind, results)
all$specific <- FALSE
for (family in list(setdiff(types, c("CPEC-1", "CPEC-2")), c("CPEC-1", "CPEC-2"))) {
    idx <- which(all$cell_type %in% family)
    key <- paste(all$gene[idx], all$direction[idx], sep="|")
    all$specific[idx] <- !duplicated(key) & !duplicated(key, fromLast=TRUE)
}
write.csv(all, file.path(out, "consensus_degs.csv"), row.names=FALSE)
