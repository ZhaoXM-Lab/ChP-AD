#!/usr/bin/env Rscript
# Plot motif specificity, TF families, concordance and prioritized GRN modules.

library(ggplot2)
library(ComplexHeatmap)
library(igraph)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir,"data/ATAC/output/grn")
out <- file.path(proj_dir,"figures/ATAC")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tf <- read.csv(file.path(base,"specific_tfs.csv"))
mot <- as.matrix(read.csv(file.path(base,"motif_group_matrix.csv"),row.names=1,check.names=FALSE))
gs <- as.matrix(read.csv(file.path(base,"gene_score_group_matrix.csv"),row.names=1,check.names=FALSE))
rna <- as.matrix(read.csv(file.path(base,"rna_group_matrix.csv"),row.names=1,check.names=FALSE))
scale_rows <- function(m) { z <- t(scale(t(m))); z[!is.finite(z)] <- 0; z }
top <- tf %>% group_by(cell_type) %>% arrange(FDR) %>% slice_head(n=10)
motifs <- intersect(unique(top$name),rownames(mot))
pdf(file.path(out,"Fig5a_motif_specificity.pdf"),width=8,height=7)
draw(Heatmap(scale_rows(mot[motifs,,drop=FALSE]),name="Motif")); dev.off()
sets <- split(tf$TF,tf$cell_type)
pdf(file.path(out,"Fig5b_TF_upset.pdf"),width=8,height=5)
draw(UpSet(make_comb_mat(sets),set_order=names(sets))); dev.off()
candidate <- read.csv(file.path(base,"tf_concordance_candidates.csv"))
cts <- Reduce(intersect,list(colnames(mot),colnames(gs),colnames(rna)))
if (nrow(candidate)) {
    pdf(file.path(out,"Fig5c_TF_concordance.pdf"),width=12,height=6)
    draw(Heatmap(scale_rows(mot[candidate$motif,cts,drop=FALSE]),name="Motif")+Heatmap(scale_rows(gs[candidate$TF,cts,drop=FALSE]),name="GeneScore")+Heatmap(scale_rows(rna[candidate$TF,cts,drop=FALSE]),name="RNA"))
    dev.off()
}
p <- ggplot(tf,aes(cell_type,fill=family))+geom_bar()+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(file.path(out,"SFig8c_TF_families.pdf"),p,width=8,height=4)
hc <- readRDS(file.path(base,"motif_hclust.rds"))
pdf(file.path(out,"SFig8a_TF_clustering.pdf"),width=8,height=5)
plot(hc,labels=FALSE); dev.off()
stats <- read.csv(file.path(base,"module_statistics.csv"))
nodes <- read.csv(file.path(base,"module_nodes.csv"))
edges <- read.csv(file.path(base,"all_grns.csv"))
for (ct in unique(stats$cell_type)) {
    selected <- stats[stats$cell_type == ct & stats$n_targets >= 3 & stats$n_TFs > 0, ]
    selected <- head(selected[order(-selected$DEG_ratio), ],if (ct == "Macrophage") 5 else 3)
    for (module in selected$module) {
        n <- nodes[nodes$cell_type == ct & nodes$module == module, ]
        e <- edges[edges$cell_type == ct & edges$TF %in% n$name[n$type == "TF"] & edges$gene %in% n$name[n$type == "Gene"], ]
        e$from <- paste0("TF_",e$TF); e$to <- paste0("G_",e$gene)
        g <- graph_from_data_frame(e[,c("from","to","weight")],directed=TRUE)
        is_tf <- grepl("^TF_",V(g)$name)
        gene <- sub("^(TF_|G_)","",V(g)$name)
        direction <- e$direction[match(gene,e$gene)]
        colors <- ifelse(direction == "Up","#B02425",ifelse(direction == "Down","#5355A2","#bbbbbb"))
        colors[is.na(colors)] <- "#bbbbbb"
        pdf(file.path(out,paste0("GRN_",ct,"_module_",module,".pdf")),width=6,height=6)
        set.seed(0)
        plot(g,vertex.label=gene,vertex.size=ifelse(is_tf,12,5),vertex.color=colors,edge.width=E(g)$weight*4,edge.arrow.size=0.2)
        dev.off()
    }
}
parameter <- read.csv(file.path(base,"mcl_parameter_statistics.csv"))
p <- ggplot(parameter,aes(inflation,largest_size,color=cell_type))+geom_line()+theme_classic()
ggsave(file.path(out,"SFig9a_MCL_parameters.pdf"),p,width=6,height=4)
for (file in list.files(file.path(proj_dir,"data/ATAC/output/celloracle"),pattern="TF_perturb_score_by_inhibition.csv",recursive=TRUE,full.names=TRUE)) {
    z <- read.csv(file)
    p <- ggplot(z,aes(sum_neg_PS,mean_PS,label=TF))+geom_point()+theme_classic()
    ggsave(file.path(out,paste0("SFig10_",basename(dirname(file)),"_PS.pdf")),p,width=5,height=4)
}
modules <- read.csv(file.path(base,"motif_modules.csv"))
modules <- modules[modules$motif %in% rownames(mot), ]
activity <- do.call(rbind,lapply(seq_len(nrow(modules)),function(i) data.frame(motif=modules$motif[i],module=modules$module[i],cell_type=colnames(mot),value=as.numeric(mot[modules$motif[i],]))))
p <- ggplot(activity,aes(cell_type,value,fill=factor(module)))+geom_boxplot(outlier.shape=NA)+facet_wrap(~module)+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(file.path(out,"SFig8b_TF_module_activity.pdf"),p,width=12,height=8)
ranked <- tf %>% group_by(cell_type) %>% arrange(desc(MeanDiff)) %>% slice_head(n=15) %>% ungroup()
p <- ggplot(ranked,aes(reorder(TF,MeanDiff),MeanDiff,fill=cell_type))+geom_col()+coord_flip()+facet_wrap(~cell_type,scales="free")+theme_classic()
ggsave(file.path(out,"SFig8e_f_TF_ranking.pdf"),p,width=12,height=8)
burden <- edges[edges$direction != "NS", ] %>% distinct(cell_type,TF,gene) %>% count(cell_type,TF,name="linked_DEGs") %>% group_by(cell_type) %>% slice_max(linked_DEGs,n=3,with_ties=FALSE) %>% ungroup()
burden <- left_join(burden,unique(tf[,c("TF","family")]),by="TF")
p <- ggplot(burden,aes(reorder(TF,linked_DEGs),linked_DEGs,fill=family))+geom_col()+coord_flip()+facet_wrap(~cell_type,scales="free")+theme_classic()
ggsave(file.path(out,"SFig8d_TF_DEG_burden.pdf"),p,width=10,height=7)
for (metric in c("cluster_count","largest_TFs","largest_genes")) {
    p <- ggplot(parameter,aes(inflation,.data[[metric]],color=cell_type))+geom_line()+theme_classic()
    ggsave(file.path(out,paste0("SFig9a_",metric,".pdf")),p,width=6,height=4)
}
summary <- stats %>% filter(n_targets >= 3,n_TFs > 0) %>% mutate(type=ifelse(n_TFs == 1,"Single TF","Multiple TFs"))
p <- ggplot(summary[summary$cell_type == "Macrophage", ],aes(reorder(factor(module),DEG_ratio),DEG_ratio,fill=type))+geom_col()+coord_flip()+theme_classic()
ggsave(file.path(out,"SFig9b_macrophage_module_ranking.pdf"),p,width=7,height=6)
top_modules <- summary %>% group_by(cell_type) %>% slice_max(DEG_ratio,n=3,with_ties=FALSE) %>% ungroup()
p <- ggplot(top_modules,aes(cell_type,factor(module),size=DEG_ratio,color=n_targets))+geom_point()+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(file.path(out,"SFig9c_top_modules.pdf"),p,width=7,height=4)
sorl1 <- edges[edges$cell_type == "Macrophage" & edges$gene == "SORL1", ]
if (nrow(sorl1)) {
    g <- graph_from_data_frame(sorl1[,c("TF","gene","weight")],directed=TRUE)
    pdf(file.path(out,"SFig9d_SORL1_network.pdf"),width=6,height=6)
    set.seed(0)
    plot(g,vertex.size=ifelse(V(g)$name == "SORL1",8,14),vertex.color=ifelse(V(g)$name == "SORL1","#5355A2","#bbbbbb"),edge.width=E(g)$weight*4,edge.arrow.size=0.2)
    dev.off()
}
ps_files <- list.files(file.path(proj_dir,"data/ATAC/output/celloracle"),pattern="TF_perturb_score_by_inhibition.csv",recursive=TRUE,full.names=TRUE)
if (length(ps_files)) {
    ps <- bind_rows(lapply(ps_files,function(file) { x <- read.csv(file); x$cell_type <- basename(dirname(file)); x }))
    counts <- ps %>% count(cell_type,name="n_TFs")
    write.csv(counts,file.path(out,"Fig5g_retained_TF_counts.csv"),row.names=FALSE)
    p <- ggplot(counts,aes(cell_type,n_TFs))+geom_col()+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
    ggsave(file.path(out,"Fig5g_retained_TFs.pdf"),p,width=7,height=4)
    top <- ps %>% group_by(cell_type) %>% slice_min(sum_neg_PS,n=3,with_ties=FALSE) %>% ungroup()
    p <- ggplot(top,aes(reorder(TF,-sum_neg_PS),-sum_neg_PS,fill=cell_type))+geom_col()+coord_flip()+facet_wrap(~cell_type,scales="free")+theme_classic()+labs(y="Negative perturbation score magnitude",x="TF")
    ggsave(file.path(out,"SFig10a_top3_TFs.pdf"),p,width=10,height=6)
}
