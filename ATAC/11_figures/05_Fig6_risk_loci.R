#!/usr/bin/env Rscript
# Render chromatin, peak-gene links and fine-mapped LD-colored variants.

library(ArchR)
library(GenomicRanges)
library(ggplot2)
library(gridExtra)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir,"data/ATAC/output/finemap")
out <- file.path(proj_dir,"figures/ATAC")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
project <- loadArchRProject(file.path(proj_dir,"data/ATAC/output/grn/ArchR_motifs"))
project <- addGroupCoverages(project,groupBy="cell_subtype")
loci <- read.delim(file.path(proj_dir,"ATAC/00_reference/loci.tsv"))
links <- read.csv(file.path(base,"variant_peak_gene_links.csv"))
for (i in seq_len(nrow(loci))) {
    gene <- loci$gene[i]; ct <- loci$cell_type[i]
    z <- read.csv(file.path(base,paste0(gene,"_",ct,"_variants_annotated.csv")))
    lead <- z[z$SNP_ID == loci$lead_snp[i], ]
    stopifnot(nrow(lead) == 1)
    region <- GRanges(lead$chr,IRanges(max(1,lead$position_hg38-500000),lead$position_hg38+500000))
    loops <- getPeak2GeneLinks(project,corCutOff=0.5,FDRCutOff=1e-4,returnLoops=TRUE)
    track <- plotBrowserTrack(project,groupBy="cell_subtype",useGroups=ct,region=region,loops=loops,plotSummary=c("bulkTrack","loopTrack","geneTrack"))[[1]]
    p <- ggplot(z,aes(position_hg38,PIP_ATAC,color=r2))+geom_point(size=1)+scale_color_gradient(low="grey80",high="#B02425",limits=c(0,1))+theme_classic()+labs(x="hg38 position",y="PIP",color="LD r2")
    p <- p+geom_point(data=lead,shape=23,fill="red",size=2)
    grob <- if (inherits(track,"ggplot")) ggplotGrob(track) else track
    combined <- arrangeGrob(grobs=list(grob,ggplotGrob(p)),ncol=1,heights=c(3,1))
    ggsave(file.path(out,paste0("Fig6_",gene,"_",ct,".pdf")),combined,width=8,height=6)
}
