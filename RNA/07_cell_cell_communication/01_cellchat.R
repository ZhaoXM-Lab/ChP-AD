#!/usr/bin/env Rscript
# Infer matched AD and control networks with population-size weighting.

library(CellChat)
library(Seurat)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/RNA/output/cellchat")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
obj <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
cells <- colnames(obj)[!obj$cell_type %in% c("Neuron", "OPC", "CD8+T")]
obj <- subset(obj, cells=cells)
objects <- list()
for (condition in c("WT", "AD")) {
    x <- subset(obj, cells=colnames(obj)[obj$disease == condition])
    x <- NormalizeData(x, normalization.method="LogNormalize", scale.factor=10000)
    meta <- x[[]]
    meta$samples <- factor(meta$sample)
    meta$cell_subtype <- droplevels(factor(meta$cell_subtype))
    cc <- createCellChat(object=GetAssayData(x, layer="data", assay="RNA"), meta=meta, group.by="cell_subtype")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    set.seed(0)
    cc <- computeCommunProb(cc, type="triMean", population.size=TRUE)
    cc <- filterCommunication(cc, min.cells=10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name="netP")
    objects[[condition]] <- cc
    saveRDS(cc, file.path(out, paste0(condition, ".rds")))
}
stopifnot(identical(levels(objects$WT@idents), levels(objects$AD@idents)))
saveRDS(mergeCellChat(objects, add.names=names(objects)), file.path(out, "merged.rds"))
