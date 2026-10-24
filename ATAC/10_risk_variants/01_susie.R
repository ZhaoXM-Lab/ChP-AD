#!/usr/bin/env Rscript
# Fine-map AD loci using aligned GWAS and cell-type accessibility priors.

library(data.table)
library(susieR)
library(genio)
library(GenomicRanges)
library(rtracklayer)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
build <- Sys.getenv("FINEMAP_BUILD")
stopifnot(build %in% c("hg19","hg38"))
out <- file.path(proj_dir, "data/ATAC/output/finemap")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
reference <- file.path(proj_dir, "data/reference/ldsc/1000G_EUR_Phase3_plink")
gwas <- fread(file.path(proj_dir,"data/reference/gwas/AD.finemap.tsv"), data.table=FALSE)
stopifnot(all(c("SNP","A1","A2","Z","N") %in% names(gwas)), !anyDuplicated(gwas$SNP))
loci <- read.delim(file.path(proj_dir, "ATAC/00_reference/loci.tsv"))
peaks <- read.csv(file.path(proj_dir,"data/ATAC/output/ccres/specific_ccres.csv"))
for (i in seq_len(nrow(loci))) {
    gene <- loci$gene[i]; ct <- loci$cell_type[i]; lead <- loci$lead_snp[i]
    found <- list()
    for (chr in 1:22) {
        bim <- fread(file.path(reference,paste0("1000G.EUR.QC.",chr,".bim")), header=FALSE, data.table=FALSE)
        hit <- which(bim[[2]] == lead)
        if (length(hit)) found <- list(chr=chr, position=bim[[4]][hit])
    }
    if (length(found) == 0) stop(paste("Lead SNP missing from reference:",lead))
    plink <- read_plink(file.path(reference,paste0("1000G.EUR.QC.",found$chr)))
    bim <- plink$bim
    snps <- bim$id[abs(bim$pos-found$position) <= 500000]
    snps <- intersect(snps,gwas$SNP)
    summary <- gwas[match(snps,gwas$SNP), ]
    variant <- bim[match(snps,bim$id), ]
    direct <- summary$A1 == variant$alt & summary$A2 == variant$ref
    reverse <- summary$A1 == variant$ref & summary$A2 == variant$alt
    summary$Z[reverse] <- -summary$Z[reverse]
    keep <- (direct | reverse) & is.finite(summary$Z) & is.finite(summary$N) & summary$N > 0
    summary <- summary[keep, ]; variant <- variant[keep, ]
    X <- t(plink$X[match(variant$id,bim$id), , drop=FALSE])
    for (j in seq_len(ncol(X))) X[is.na(X[,j]),j] <- mean(X[,j],na.rm=TRUE)
    variable <- apply(X,2,function(x) all(is.finite(x)) && sd(x) > 0)
    X <- X[,variable,drop=FALSE]; summary <- summary[variable, ]; variant <- variant[variable, ]
    stopifnot(ncol(X) >= 10, identical(as.character(summary$SNP),as.character(variant$id)))
    gr <- GRanges(paste0("chr",found$chr), IRanges(variant$pos,variant$pos))
    if (build == "hg19") {
        lifted <- liftOver(gr,import.chain(file.path(proj_dir,"data/reference/hg19ToHg38.over.chain")))
        ok <- elementNROWS(lifted) == 1
        gr <- unlist(lifted[ok],use.names=FALSE)
        X <- X[,ok,drop=FALSE]; summary <- summary[ok, ]; variant <- variant[ok, ]
    }
    marker <- GRanges(peaks$peak[peaks$cell_type == ct])
    stopifnot(length(marker) > 0)
    prior <- ifelse(countOverlaps(gr,marker) > 0,1,0.1)
    prior <- prior/sum(prior)
    ld <- cor(X)
    n <- mean(summary$N)
    lambda <- estimate_s_rss(summary$Z,ld,n=n)
    reg <- (1-lambda)*ld+lambda*diag(ncol(ld))
    fits <- list(original=susie_rss(z=summary$Z,R=reg,n=n,L=5), atac=susie_rss(z=summary$Z,R=reg,n=n,L=5,prior_weights=prior))
    result <- data.frame(SNP_ID=summary$SNP, chr=as.character(seqnames(gr)), position_hg38=start(gr), reference_position=variant$pos,
                         PIP_original=fits$original$pip, PIP_ATAC=fits$atac$pip, prior=prior)
    prefix <- file.path(out,paste0(gene,"_",ct))
    write.csv(result,paste0(prefix,"_variants.csv"),row.names=FALSE)
    dimnames(ld) <- list(summary$SNP,summary$SNP)
    write.table(ld^2,paste0(prefix,"_r2.tsv"),sep="\t",quote=FALSE,col.names=NA)
    saveRDS(list(fits=fits,lambda=lambda,sample_size=n,build=build),paste0(prefix,"_fits.rds"))
}
