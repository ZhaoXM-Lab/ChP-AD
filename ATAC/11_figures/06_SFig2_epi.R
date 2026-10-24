#!/usr/bin/env Rscript
# Plot EPI distributions, linked-feature comparisons and validation AUPRC.
library(ArchR)
library(Seurat)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir,"data/ATAC/output/epi")
out <- file.path(proj_dir,"figures/ATAC")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
links <- read.csv(file.path(base,"archr_links.csv"))
project <- loadArchRProject(file.path(base,"ArchR_links"))
genes <- getGenes(project)
g <- genes[match(links$gene,genes$symbol)]
tss <- ifelse(as.character(strand(g)) == "-",end(g),start(g))
peak <- GRanges(links$peak)
links$distance <- abs(start(resize(peak,1,"center"))-tss)
degree <- links %>% group_by(cell_type,gene) %>% summarise(n_enhancers=n_distinct(peak),.groups="drop")
p <- ggplot(degree,aes(n_enhancers))+geom_histogram(binwidth=1)+geom_vline(xintercept=median(degree$n_enhancers),linetype=2)+theme_classic()
ggsave(file.path(out,"SFig2f_enhancers_per_gene.pdf"),p,width=5,height=3)
p <- ggplot(links,aes(distance,fill=cell_type))+geom_histogram(bins=50)+theme_classic()
ggsave(file.path(out,"SFig2h_EPI_distance.pdf"),p,width=6,height=3)
rna <- readRDS(file.path(proj_dir,"data/RNA/output/annotation/human_final.rds"))
rna$group <- rna$cell_type
avg <- AverageExpression(rna,assays="RNA",group.by="group",layer="data")$RNA
degree$expression <- mapply(function(g,ct) if (g %in% rownames(avg) && ct %in% colnames(avg)) avg[g,ct] else NA_real_,degree$gene,degree$cell_type)
p <- ggplot(degree,aes(factor(pmin(n_enhancers,5)),log1p(expression),fill=cell_type))+geom_boxplot(outlier.shape=NA)+theme_classic()+labs(x="Linked enhancers (5 = 5 or more)")
ggsave(file.path(out,"SFig2g_j_link_degree_expression.pdf"),p,width=7,height=4)
pk <- getGroupSE(project,useMatrix="PeakMatrix",groupBy="cell_type",divideN=TRUE)
pm <- assay(pk); gr <- getPeakSet(project)
rownames(pm) <- paste0(seqnames(gr),":",start(gr),"-",end(gr))
comparison <- list()
for (ct in intersect(colnames(pm),colnames(avg))) {
    for (modality in c("RNA","ATAC")) {
        values <- if (modality == "RNA") avg[,ct] else pm[,ct]
        linked <- if (modality == "RNA") links$gene[links$cell_type == ct] else links$peak[links$cell_type == ct]
        comparison[[length(comparison)+1]] <- data.frame(cell_type=ct,modality=modality,value=as.numeric(values),linked=names(values) %in% linked)
    }
}
table <- bind_rows(comparison)
p <- ggplot(table,aes(factor(linked),log1p(value),fill=factor(linked)))+geom_violin(scale="width")+facet_grid(modality~cell_type,scales="free_y")+theme_classic()+labs(x="Linked feature",y="log1p activity/expression")
ggsave(file.path(out,"SFig2i_linked_features.pdf"),p,width=12,height=5)
baseline <- read.csv(file.path(base,"baseline_links.csv"))
counts <- bind_rows(links %>% transmute(cell_type,peak,gene,method="ArchR",iteration=0),baseline[,c("cell_type","peak","gene","method","iteration")]) %>% group_by(cell_type,method,iteration) %>% summarise(n=n(),.groups="drop")
p <- ggplot(counts,aes(cell_type,n,fill=method))+geom_col(position="dodge")+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(file.path(out,"SFig2k_EPI_counts.pdf"),p,width=8,height=4)
validation <- read.csv(file.path(base,"validation_auprc.csv"))
p <- ggplot(validation,aes(cell_type,AUPRC,fill=method))+geom_boxplot()+facet_wrap(~assay)+theme_classic()
ggsave(file.path(out,"SFig2l_external_validation.pdf"),p,width=8,height=4)
