#!/usr/bin/env Rscript
# Validate comparator links against matched eQTL and pcHi-C references.

library(GenomicRanges)
library(PRROC)
library(dplyr)
label_epi_links <- function(x, gold, assay, genes) {
    if (assay == "eQTL") {
        gr <- makeGRangesFromDataFrame(gold, seqnames.field="Chromosome", start.field="Start", end.field="End")
        gr <- resize(gr, width=width(gr)+4000, fix="center")
        start(gr) <- pmax(start(gr), 1)
        hits <- findOverlaps(GRanges(x$peak), gr)
        ids <- setNames(sub("\\..*$", "", genes$gene_id), genes$symbol)
        valid <- unname(ids[x$gene[queryHits(hits)]]) == sub("\\..*$", "", gold$gold_gene_id[subjectHits(hits)])
        positives <- unique(queryHits(hits)[which(valid)])
    } else if (assay == "pcHiC") {
        enh <- makeGRangesFromDataFrame(gold, seqnames.field="chr2_hg38", start.field="start2_hg38", end.field="end2_hg38")
        prom <- makeGRangesFromDataFrame(gold, seqnames.field="chr1_hg38", start.field="start1_hg38", end.field="end1_hg38")
        hits <- findOverlaps(GRanges(x$peak), enh)
        gene_idx <- match(x$gene[queryHits(hits)], genes$symbol)
        stopifnot(!anyNA(gene_idx))
        pairs <- findOverlaps(genes[gene_idx], prom[subjectHits(hits)])
        same <- queryHits(pairs) == subjectHits(pairs)
        positives <- unique(queryHits(hits)[queryHits(pairs)[same]])
    } else stop("Assay must be eQTL or pcHiC.")
    as.integer(seq_len(nrow(x)) %in% positives)
}
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/epi")
refs <- read.delim(file.path(proj_dir, "ATAC/00_reference/epi_validation.tsv"))
stopifnot(nrow(refs) > 0)
genes <- readRDS(file.path(out, "genes.rds"))
gene_ids <- setNames(sub("\\..*$", "", genes$gene_id), genes$symbol)
archr <- read.csv(file.path(out, "archr_links.csv"))
archr$method <- "ArchR"
archr$score <- archr$Correlation
archr$iteration <- 0L
baseline <- read.csv(file.path(out, "baseline_links.csv"))
pred <- bind_rows(archr[, c("cell_type","gene","peak","method","score","iteration")], baseline[, c("cell_type","gene","peak","method","score","iteration")])
result <- list()
for (i in seq_len(nrow(refs))) {
    gold <- read.csv(file.path(proj_dir, refs$file[i]))
    x <- pred[pred$cell_type == refs$cell_type[i], ]
    x$label <- label_epi_links(x, gold, refs$assay[i], genes)
    for (method in unique(x$method)) for (iteration in unique(x$iteration[x$method == method])) {
        z <- x[x$method == method & x$iteration == iteration & is.finite(x$score), ]
        auc <- if (length(unique(z$label)) == 2) pr.curve(scores.class0=z$score[z$label == 1], scores.class1=z$score[z$label == 0])$auc.integral else NA_real_
        result[[length(result)+1]] <- data.frame(cell_type=refs$cell_type[i], assay=refs$assay[i], method=method, iteration=iteration, AUPRC=auc, n=nrow(z), positives=sum(z$label))
    }
}
write.csv(bind_rows(result), file.path(out, "validation_auprc.csv"), row.names=FALSE)
