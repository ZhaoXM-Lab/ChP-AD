#!/usr/bin/env Rscript
# Plot CellChat differences, APP-SORL1 signaling and outgoing patterns.

library(CellChat)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(circlize)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir,"data/RNA/output/cellchat")
out <- file.path(proj_dir,"figures/RNA")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
merged <- readRDS(file.path(base,"merged.rds"))
pdf(file.path(out,"Fig4a_CCC_difference.pdf"),width=6,height=6)
netVisual_diffInteraction(merged,weight.scale=TRUE,measure="weight")
dev.off()
delta <- as.matrix(read.csv(file.path(base,"weight_difference.csv"),row.names=1,check.names=FALSE))
summary <- data.frame(cell_type=rownames(delta),delta=rowSums(delta))
p <- ggplot(summary,aes(reorder(cell_type,delta),delta,fill=delta > 0))+geom_col()+coord_flip()+theme_classic()+labs(x="",y="Differential outgoing strength")
ggsave(file.path(out,"Fig4b_outgoing.pdf"),p,width=5,height=4)
pdf(file.path(out,"Fig4c_CPEC1_chord.pdf"),width=6,height=6)
chordDiagram(data.frame(source="CPEC-1",target=colnames(delta),value=abs(delta["CPEC-1",])),col=ifelse(delta["CPEC-1",] > 0,"#B02425","#5355A2"))
circos.clear(); dev.off()
top <- read.csv(file.path(base,"cpec1_top20.csv"))
p <- ggplot(top,aes(target,interaction_name,size=abs(delta),color=delta))+geom_point()+scale_color_gradient2(low="#5355A2",mid="white",high="#B02425")+theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1))
ggsave(file.path(out,"Fig4d_CPEC1_top20_LR.pdf"),p,width=6,height=5)
for (d in c("WT","AD")) {
    cc <- readRDS(file.path(base,paste0(d,".rds")))
    pdf(file.path(out,paste0("SFig6_",d,"_network.pdf")),width=6,height=6)
    netVisual_circle(cc@net$weight,vertex.weight=as.numeric(table(cc@idents)),weight.scale=TRUE)
    dev.off()
    pdf(file.path(out,paste0("Fig4e_APP_SORL1_",d,".pdf")),width=6,height=6)
    netVisual_chord_gene(cc,pairLR.use=data.frame(interaction_name="APP_SORL1"))
    dev.off()
    cells <- read.csv(file.path(base,paste0(d,"_cell_patterns.csv")))
    pathways <- read.csv(file.path(base,paste0(d,"_pathway_patterns.csv")))
    for (nm in c("cells","pathways")) {
        x <- if (nm == "cells") cells else pathways
        x$label <- if (nm == "cells") x$CellGroup else x$Signaling
        p <- ggplot(x,aes(axis1=label,axis2=Pattern,y=Contribution))+geom_alluvium(aes(fill=factor(Pattern)),alpha=0.6)+geom_stratum()+geom_text(stat="stratum",aes(label=after_stat(stratum)),size=2)+theme_void()
        ggsave(file.path(out,paste0("Fig4f_",d,"_",nm,".pdf")),p,width=6,height=5)
    }
}
ggsave(file.path(out,"SFig6c_interaction_count.pdf"),compareInteractions(merged,measure="count"),width=4,height=3)
ggsave(file.path(out,"SFig6d_interaction_weight.pdf"),compareInteractions(merged,measure="weight"),width=4,height=3)
