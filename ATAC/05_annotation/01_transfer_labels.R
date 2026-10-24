#!/usr/bin/env Rscript
# Transfer RNA labels through CCA anchors and create a joint embedding.

library(Signac)
library(Seurat)
library(harmony)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/annotation")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
rna <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
atac <- readRDS(file.path(proj_dir, "data/ATAC/output/umap/atac_umap.rds"))
DefaultAssay(rna) <- "RNA"
activities <- GeneActivity(atac, features=VariableFeatures(rna))
atac[["Activity"]] <- CreateAssayObject(counts=activities)
DefaultAssay(atac) <- "Activity"
atac <- NormalizeData(atac)
atac <- FindVariableFeatures(atac, nfeatures=2000)
genes <- intersect(VariableFeatures(rna), VariableFeatures(atac))
stopifnot(length(genes) >= 30)
atac <- ScaleData(atac, features=genes)
set.seed(12355)
anchors <- FindTransferAnchors(reference=rna, query=atac, reference.assay="RNA", query.assay="Activity", features=genes, reduction="cca", normalization.method="LogNormalize")
labels <- TransferData(anchorset=anchors, refdata=rna$cell_subtype, weight.reduction=atac[["harmony"]], dims=1:20, k.weight=5)
atac <- AddMetaData(atac, labels)
atac$cell_subtype <- atac$predicted.id
atac$cell_type <- ifelse(atac$cell_subtype %in% c("CPEC-1", "CPEC-2"), "CPEC", atac$cell_subtype)
imputed <- TransferData(anchorset=anchors, refdata=GetAssayData(rna, assay="RNA", layer="data")[genes, ], weight.reduction=atac[["harmony"]], dims=2:28, k.weight=5)
atac[["RNA"]] <- imputed
rna$modality <- "RNA"
atac$modality <- "ATAC"
coembed <- merge(rna, y=atac)
DefaultAssay(coembed) <- "RNA"
coembed <- JoinLayers(coembed)
coembed <- ScaleData(coembed, features=genes, do.scale=FALSE)
coembed <- RunPCA(coembed, features=genes)
coembed <- RunHarmony(coembed, group.by.vars="sample", reduction.use="pca", project.dim=FALSE)
coembed <- RunUMAP(coembed, reduction="harmony", dims=1:15, seed.use=1)
write.csv(aggregate(atac$prediction.score.max, list(cell_type=atac$cell_subtype), mean), file.path(out, "label_prediction_scores.csv"), row.names=FALSE)
saveRDS(atac, file.path(out, "atac_annotated.rds"))
saveRDS(coembed, file.path(out, "coembedding.rds"))
