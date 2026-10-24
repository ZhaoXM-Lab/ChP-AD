#!/usr/bin/env Rscript
# Build distance, ranked-activity and 30 random EPI comparators.

library(ArchR)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/epi")
project <- loadArchRProject(file.path(out, "ArchR_links"))
peak_se <- getMatrixFromProject(project, useMatrix="PeakMatrix")
gene_se <- getMatrixFromProject(project, useMatrix="GeneScoreMatrix")
peak_mat <- assay(peak_se)
gene_mat <- assay(gene_se)
peaks <- getPeakSet(project)
peak_names <- paste0(seqnames(peaks), ":", start(peaks), "-", end(peaks))
rownames(peak_mat) <- peak_names
rownames(gene_mat) <- rowData(gene_se)$name
genes <- getGenes(project)
genes <- genes[!duplicated(genes$symbol), ]
markers <- read.csv(file.path(proj_dir, "data/ATAC/output/ccres/specific_ccres.csv"))
all <- list()
set.seed(0)
for (ct in unique(markers$cell_type)) {
    cells <- getCellNames(project)[if (ct %in% c("CPEC-1","CPEC-2")) project$cell_subtype == ct else project$cell_type == ct]
    candidate <- unique(markers$peak[markers$cell_type == ct])
    candidate <- match(candidate, peak_names)
    gm <- getMarkerFeatures(project, useMatrix="GeneScoreMatrix", groupBy=if (ct %in% c("CPEC-1","CPEC-2")) "cell_subtype" else "cell_type", useGroups=ct, bias=c("TSSEnrichment","log10(nFrags)"), testMethod="wilcoxon")
    gtab <- getMarkers(gm, cutOff="Pval <= 0.01 & Log2FC >= 0.5")[[ct]]
    for (gene in intersect(gtab$name, rownames(gene_mat))) {
        g <- genes[genes$symbol == gene]
        if (!length(g)) next
        tss <- ifelse(as.character(strand(g))[1] == "-", end(g)[1], start(g)[1])
        dist <- abs(start(resize(peaks[candidate], 1, "center")) - tss)
        eligible <- which(as.character(seqnames(peaks[candidate])) == as.character(seqnames(g))[1] & dist <= 250000)
        if (!length(eligible)) next
        idx <- candidate[eligible]
        distance <- dist[eligible]
        add <- function(j, method, score, iteration=0L) data.frame(cell_type=ct, gene=gene, peak=peak_names[idx[j]], method=method, score=score, iteration=iteration)
        j <- which.min(distance)
        all[[length(all)+1]] <- add(j, "distance", if (distance[j] == 0) 1 else NA_real_)
        y <- rank(as.numeric(gene_mat[gene, cells]))
        correlations <- vapply(idx, function(i) suppressWarnings(cor(rank(as.numeric(peak_mat[i, cells])), y)), numeric(1))
        statistic <- correlations * sqrt((length(cells)-2)/pmax(1-correlations^2, .Machine$double.eps))
        p <- 2*pt(abs(statistic), df=length(cells)-2, lower.tail=FALSE)
        valid <- which(is.finite(correlations) & correlations > 0 & p < 0.05)
        if (length(valid)) {
            j <- valid[which.max(correlations[valid])]
            all[[length(all)+1]] <- add(j, "activity", correlations[j])
        }
        for (iteration in 1:30) {
            j <- sample.int(length(idx), 1)
            all[[length(all)+1]] <- add(j, "random", 1, iteration)
        }
    }
}
tab <- bind_rows(all)
gr <- GRanges(tab$peak)
gene_idx <- match(tab$gene, genes$symbol)
tss <- ifelse(as.character(strand(genes[gene_idx])) == "-", end(genes[gene_idx]), start(genes[gene_idx]))
tab$distance <- abs(start(resize(gr, 1, "center")) - tss)
for (ct in unique(tab$cell_type)) for (method in c("distance","random")) for (iteration in unique(tab$iteration[tab$method == method])) {
    idx <- which(tab$cell_type == ct & tab$method == method & tab$iteration == iteration)
    positive <- idx[tab$distance[idx] > 0]
    raw <- log10(1000/tab$distance[positive])
    tab$score[positive] <- if (length(unique(raw)) > 1) (raw-min(raw))/(max(raw)-min(raw)) else 1
    tab$score[idx[tab$distance[idx] == 0]] <- 1
}
write.csv(tab, file.path(out, "baseline_links.csv"), row.names=FALSE)
