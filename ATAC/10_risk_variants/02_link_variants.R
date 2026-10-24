#!/usr/bin/env Rscript
# Map PIP-prioritized hg38 variants to cell-specific cCRE-gene links.

library(GenomicRanges)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir,"data/ATAC/output/finemap")
loci <- read.delim(file.path(proj_dir,"ATAC/00_reference/loci.tsv"))
epi <- read.csv(file.path(proj_dir,"data/ATAC/output/epi/archr_links.csv"))
all <- list()
for (i in seq_len(nrow(loci))) {
    gene <- loci$gene[i]; ct <- loci$cell_type[i]; lead <- loci$lead_snp[i]
    prefix <- file.path(base,paste0(gene,"_",ct))
    variants <- read.csv(paste0(prefix,"_variants.csv"))
    r2 <- as.matrix(read.delim(paste0(prefix,"_r2.tsv"),row.names=1,check.names=FALSE))
    stopifnot(lead %in% colnames(r2))
    variants$r2 <- r2[variants$SNP_ID,lead]
    links <- epi[epi$cell_type == ct & epi$gene == gene, ]
    hits <- findOverlaps(GRanges(variants$chr,IRanges(variants$position_hg38,variants$position_hg38)),GRanges(links$peak))
    tab <- cbind(variants[queryHits(hits), ],links[subjectHits(hits),c("peak","Correlation","FDR"),drop=FALSE])
    tab$gene <- gene; tab$cell_type <- ct; tab$lead_snp <- lead
    tab$prioritized <- tab$PIP_ATAC > 0.05
    all[[length(all)+1]] <- tab
    write.csv(variants,paste0(prefix,"_variants_annotated.csv"),row.names=FALSE)
}
write.csv(bind_rows(all),file.path(base,"variant_peak_gene_links.csv"),row.names=FALSE)
