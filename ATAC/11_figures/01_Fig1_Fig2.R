#!/usr/bin/env Rscript
# Plot co-embedding, ATAC markers and cross-modal feature concordance.

library(Seurat)
library(ArchR)
library(SummarizedExperiment)
library(ComplexHeatmap)
library(ggplot2)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir,"figures/ATAC")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
coembed <- readRDS(file.path(proj_dir,"data/ATAC/output/annotation/coembedding.rds"))
atac <- readRDS(file.path(proj_dir,"data/ATAC/output/annotation/atac_annotated.rds"))
rna <- readRDS(file.path(proj_dir,"data/RNA/output/annotation/human_final.rds"))
for (modality in c("RNA","ATAC")) {
    x <- subset(coembed,cells=colnames(coembed)[coembed$modality == modality])
    p <- DimPlot(x,group.by="cell_type",label=TRUE,repel=TRUE)+NoLegend()
    ggsave(file.path(out,paste0("SFig2_coembedding_",modality,".pdf")),p,width=6,height=5)
}
ggsave(file.path(out,"Fig1c_ATAC.pdf"),DimPlot(atac,group.by="cell_type",label=TRUE,repel=TRUE),width=7,height=5)
ggsave(file.path(out,"SFig2a_joint_embedding.pdf"),DimPlot(coembed,group.by="modality"),width=6,height=5)
cpec <- subset(coembed,cells=colnames(coembed)[coembed$modality == "ATAC" & coembed$cell_type == "CPEC"])
ggsave(file.path(out,"Fig2b_ATAC_CPEC.pdf"),DimPlot(cpec,group.by="cell_subtype"),width=5,height=4)
project <- loadArchRProject(file.path(proj_dir,"data/ATAC/output/grn/ArchR_motifs"))
gs <- getMatrixFromProject(project,useMatrix="GeneScoreMatrix")
gsm <- assay(gs); rownames(gsm) <- rowData(gs)$name
common <- intersect(colnames(atac),colnames(gsm))
atac <- subset(atac,cells=common)
atac[["GeneScore"]] <- CreateAssayObject(counts=gsm[,common,drop=FALSE])
atac <- NormalizeData(atac,assay="GeneScore")
markers <- c("OTX2","LMX1A","GFAP","AQP4","CFAP299","UGT8","CYBB","CIITA","FLT1","CDH5","SLC4A4","RGS5","ABCA9","KCNK2")
markers <- intersect(markers,rownames(gsm))
ggsave(file.path(out,"Fig1e_gene_activity.pdf"),DotPlot(atac,features=markers,group.by="cell_type",assay="GeneScore",col.min=-2,col.max=2)+RotatedAxis(),width=7,height=4)
ggsave(file.path(out,"Fig2d_CPEC_gene_activity.pdf"),DotPlot(subset(atac,cells=colnames(atac)[atac$cell_type == "CPEC"]),features=intersect(c("DNAH5","SLC24A4","NDUFA1","NDUFAB1"),rownames(gsm)),group.by="cell_subtype",assay="GeneScore")+RotatedAxis(),width=5,height=3)
gs_avg <- assay(getGroupSE(project,useMatrix="GeneScoreMatrix",groupBy="cell_type",divideN=TRUE))
rownames(gs_avg) <- rowData(getGroupSE(project,useMatrix="GeneScoreMatrix",groupBy="cell_type",divideN=TRUE))$name
rna_avg <- AverageExpression(rna,assays="RNA",group.by="cell_type",layer="data")$RNA
Idents(rna) <- "cell_type"
rna_markers <- FindAllMarkers(rna,only.pos=TRUE,min.pct=0.25,logfc.threshold=0.25)
write.csv(rna_markers,file.path(out,"RNA_cell_type_markers.csv"),row.names=FALSE)
gs_marker <- getMarkerFeatures(project,useMatrix="GeneScoreMatrix",groupBy="cell_type",bias=c("TSSEnrichment","log10(nFrags)"),testMethod="wilcoxon")
gm <- getMarkers(gs_marker,cutOff="FDR <= 0.05 & Log2FC >= 0.5")
shared <- unique(unlist(lapply(names(gm),function(ct) intersect(gm[[ct]]$name,rna_markers$gene[rna_markers$cluster == ct]))))
genes <- Reduce(intersect,list(shared,rownames(gs_avg),rownames(rna_avg)))
cts <- intersect(colnames(gs_avg),colnames(rna_avg))
scale_rows <- function(m) { z <- t(scale(t(as.matrix(m)))); z[!is.finite(z)] <- 0; z }
if (length(genes)) {
    pdf(file.path(out,"Fig1f_shared_gene_concordance.pdf"),width=10,height=7)
    draw(Heatmap(scale_rows(rna_avg[genes,cts]),name="RNA",show_row_names=FALSE,cluster_rows=FALSE,use_raster=TRUE)+Heatmap(scale_rows(gs_avg[genes,cts]),name="GeneScore",show_row_names=FALSE,cluster_rows=FALSE,use_raster=TRUE))
    dev.off()
}
peak <- getGroupSE(project,useMatrix="PeakMatrix",groupBy="cell_type",divideN=TRUE)
pm <- assay(peak)
gr <- getPeakSet(project)
rownames(pm) <- paste0(seqnames(gr),":",start(gr),"-",end(gr))
epi <- read.csv(file.path(proj_dir,"data/ATAC/output/epi/archr_links.csv"))
epi <- unique(epi[,c("peak","gene")])
epi <- epi[epi$peak %in% rownames(pm) & epi$gene %in% rownames(rna_avg), ]
if (nrow(epi)) {
    pdf(file.path(out,"Fig1f_linked_cCRE_gene_concordance.pdf"),width=10,height=7)
    draw(Heatmap(scale_rows(pm[epi$peak,cts]),name="cCRE",show_row_names=FALSE,cluster_rows=FALSE,use_raster=TRUE)+Heatmap(scale_rows(rna_avg[epi$gene,cts]),name="Linked_RNA",show_row_names=FALSE,cluster_rows=FALSE,use_raster=TRUE))
    dev.off()
}
pdf(file.path(out,"SFig2d_peak_annotation.pdf"),width=5,height=4)
barplot(table(gr$peakType),las=2,ylab="cCRE count")
dev.off()
marker_ccres <- read.csv(file.path(proj_dir,"data/ATAC/output/ccres/specific_ccres.csv"))
marker_ccres <- marker_ccres[order(marker_ccres$cell_type,-marker_ccres$Log2FC), ]
marker_ccres <- marker_ccres[!duplicated(marker_ccres$peak), ]
marker_peaks <- intersect(marker_ccres$peak,rownames(pm))
pdf(file.path(out,"SFig2e_specific_cCRE_heatmap.pdf"),width=8,height=6)
draw(Heatmap(scale_rows(pm[marker_peaks,,drop=FALSE]),name="Accessibility",show_row_names=FALSE,cluster_rows=FALSE,use_raster=TRUE))
dev.off()
composition <- rbind(RNA=prop.table(table(factor(rna$cell_type,levels=cts))),ATAC=prop.table(table(factor(atac$cell_type,levels=cts))))
pdf(file.path(out,"Fig1f_modality_composition.pdf"),width=6,height=4)
barplot(t(composition),beside=FALSE,legend.text=colnames(composition),args.legend=list(cex=0.5),ylab="Cell fraction")
dev.off()
