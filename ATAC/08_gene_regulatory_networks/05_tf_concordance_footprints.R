#!/usr/bin/env Rscript
# Compare motif, gene-score and RNA profiles and plot OTX2/CEBPA footprints.

library(ArchR)
library(Seurat)
library(SummarizedExperiment)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/grn")
project <- loadArchRProject(file.path(out, "ArchR_motifs"))
rna <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
mot <- getGroupSE(project, useMatrix="VierstraMotifMatrix", groupBy="cell_type", divideN=TRUE)
mot <- mot[rowData(mot)$seqnames == "z", ]
mm <- assay(mot)
rownames(mm) <- rowData(mot)$name
gs <- getGroupSE(project, useMatrix="GeneScoreMatrix", groupBy="cell_type", divideN=TRUE)
gsm <- assay(gs)
rownames(gsm) <- rowData(gs)$name
rm <- AverageExpression(rna, group.by="cell_type", assays="RNA", layer="data")$RNA
cts <- Reduce(intersect, list(colnames(mm), colnames(gsm), colnames(rm)))
cts <- setdiff(cts, c("Neuron","OPC","CD8+T"))
specific <- read.csv(file.path(out, "specific_tfs.csv"))
res <- list()
for (i in seq_len(nrow(mm))) {
    tf <- sub("_.*$", "", rownames(mm)[i])
    if (!tf %in% intersect(rownames(gsm), rownames(rm))) next
    a <- as.numeric(mm[i,cts]); b <- as.numeric(gsm[tf,cts]); c <- as.numeric(rm[tf,cts])
    if (min(sd(a),sd(b),sd(c)) == 0) next
    ab <- cor.test(a,b); bc <- cor.test(b,c)
    res[[length(res)+1]] <- data.frame(motif=rownames(mm)[i], TF=tf, motif_GS=ab$estimate, p_motif_GS=ab$p.value, GS_RNA=bc$estimate, p_GS_RNA=bc$p.value,
                                      motif_ct=cts[which.max(a)], GS_ct=cts[which.max(b)], RNA_ct=cts[which.max(c)])
}
tab <- bind_rows(res)
tab$FDR_motif_GS <- p.adjust(tab$p_motif_GS,"BH")
tab$FDR_GS_RNA <- p.adjust(tab$p_GS_RNA,"BH")
tab <- tab %>% arrange(desc(motif_GS)) %>% distinct(TF, .keep_all=TRUE)
write.csv(tab, file.path(out,"tf_concordance_all.csv"), row.names=FALSE)
candidate <- tab[tab$motif_ct == tab$GS_ct & tab$GS_ct == tab$RNA_ct & tab$FDR_motif_GS < 0.01 & tab$FDR_GS_RNA < 0.01, ]
write.csv(candidate, file.path(out,"tf_concordance_candidates.csv"), row.names=FALSE)
write.csv(mm, file.path(out,"motif_group_matrix.csv"))
sub <- getGroupSE(project, useMatrix="VierstraMotifMatrix", groupBy="cell_subtype", divideN=TRUE)
sub <- sub[rowData(sub)$seqnames == "z", ]
submat <- assay(sub)
rownames(submat) <- rowData(sub)$name
write.csv(submat, file.path(out,"motif_subtype_group_matrix.csv"))
write.csv(gsm, file.path(out,"gene_score_group_matrix.csv"))
write.csv(as.matrix(rm), file.path(out,"rna_group_matrix.csv"))
z <- getMatrixFromProject(project, useMatrix="VierstraMotifMatrix", useSeqnames="z")
zmat <- as.matrix(assay(z))
rownames(zmat) <- rowData(z)$name
zmat <- zmat[rownames(zmat) %in% unique(specific$name), , drop=FALSE]
nonconstant <- apply(zmat,1,sd) > 0
zmat <- zmat[nonconstant, , drop=FALSE]
hc <- hclust(as.dist(1-cor(t(zmat))), method="ward.D2")
write.csv(data.frame(motif=names(cutree(hc,k=10)), module=cutree(hc,k=10)), file.path(out,"motif_modules.csv"), row.names=FALSE)
saveRDS(hc, file.path(out,"motif_hclust.rds"))
project <- addGroupCoverages(project, groupBy="cell_type")
positions <- getPositions(project, name="VierstraMotif")
names_use <- unique(unlist(lapply(c("OTX2", "CEBPA"), function(tf) grep(paste0("^",tf,"_"), names(positions), value=TRUE))))
footprints <- getFootprints(project, positions=positions[names_use], groupBy="cell_type")
plotFootprints(footprints, ArchRProj=project, normMethod="Subtract", plotName="OTX2_CEBPA", smoothWindow=5, addDOC=FALSE)
saveRDS(footprints, file.path(out,"OTX2_CEBPA_footprints.rds"))
