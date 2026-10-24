#!/usr/bin/env Rscript
# Infer cell-type trajectories and export aligned Monocle3 coordinates.

library(Seurat)
library(monocle3)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
ct <- Sys.getenv("CELL_TYPE", "Macrophage")
out <- file.path(proj_dir, "data/ATAC/output/celloracle", ct)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
rna <- readRDS(file.path(proj_dir, "data/RNA/output/annotation/human_final.rds"))
cells <- colnames(rna)[if (ct %in% c("CPEC-1","CPEC-2")) rna$cell_subtype == ct else rna$cell_type == ct]
x <- subset(rna, cells=cells)
x <- FindVariableFeatures(x, nfeatures=2000)
counts <- GetAssayData(x, layer="counts")
cds <- new_cell_data_set(counts, cell_metadata=x[[]], gene_metadata=data.frame(gene_short_name=rownames(counts), row.names=rownames(counts)))
set.seed(0)
cds <- preprocess_cds(cds, num_dim=25, use_genes=VariableFeatures(x), norm_method="log")
cds <- reduce_dimension(cds, preprocess_method="PCA", reduction_method="UMAP", umap.n_neighbors=15, umap.min_dist=0.5, max_components=2)
cds <- cluster_cells(cds, resolution=5e-5)
cds <- learn_graph(cds, use_partition=FALSE, close_loop=FALSE)
cds <- order_cells(cds, root_cells=colnames(cds)[cds$disease == "WT"])
emb <- SingleCellExperiment::reducedDims(cds)$UMAP
write.csv(data.frame(cell=rownames(emb), UMAP1=emb[,1], UMAP2=emb[,2], pseudotime=pseudotime(cds)[rownames(emb)]), file.path(out,"trajectory.csv"), row.names=FALSE)
writeLines(VariableFeatures(x), file.path(out,"hvg.txt"))
saveRDS(cds, file.path(out,"monocle3.rds"))
for (color in c("pseudotime","disease")) {
    p <- plot_cells(cds,color_cells_by=color,label_cell_groups=FALSE,label_leaves=FALSE,label_branch_points=FALSE,label_roots=FALSE,label_principal_points=FALSE)
    if (color == "disease") p <- p+ggplot2::scale_color_manual(values=c(AD="#cc0000",WT="#00167a"))
    ggplot2::ggsave(file.path(out,paste0("monocle3_",color,".pdf")),p,width=5,height=5)
}
