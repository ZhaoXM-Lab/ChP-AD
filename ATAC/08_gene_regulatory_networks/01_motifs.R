#!/usr/bin/env Rscript
# Annotate Vierstra v2.1 motifs and identify cell-specific motif deviations.

library(ArchR)
library(SummarizedExperiment)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/grn")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
project <- loadArchRProject(file.path(proj_dir, "data/ATAC/output/epi/ArchR_links"))
motifs <- readRDS(file.path(proj_dir, "data/reference/Vierstra2.1_Motif_individual.rds"))
project <- addMotifAnnotations(project, motifPWMs=motifs, annoName="VierstraMotif", species="Homo sapiens", cutOff=5e-5)
project <- addBgdPeaks(project)
project <- addDeviationsMatrix(project, peakAnnotation="VierstraMotif")
se <- getMatrixFromProject(project, useMatrix="VierstraMotifMatrix")
saveRDS(se, file.path(out, "motif_deviations.rds"))
all <- list()
for (group in c("cell_type", "cell_subtype")) {
    m <- getMarkerFeatures(project, useMatrix="VierstraMotifMatrix", useSeqnames="z", groupBy=group, bias=c("TSSEnrichment","log10(nFrags)"), testMethod="wilcoxon")
    saveRDS(m, file.path(out, paste0(group,"_motif_markers.rds")))
    markers <- getMarkers(m, cutOff="FDR <= 0.05 & MeanDiff >= 0.3")
    for (ct in names(markers)) {
        if (group == "cell_subtype" && !ct %in% c("CPEC-1", "CPEC-2")) next
        z <- as.data.frame(markers[[ct]])
        if (!nrow(z)) next
        z$cell_type <- ct
        z$TF <- sub("_.*$", "", z$name)
        all[[length(all)+1]] <- z
    }
}
tab <- bind_rows(all)
family <- read.delim(file.path(proj_dir, "data/reference/tf_families.tsv"))
tab <- left_join(tab, family, by="TF")
write.csv(tab, file.path(out, "specific_tfs.csv"), row.names=FALSE)
saveArchRProject(project, outputDirectory=file.path(out, "ArchR_motifs"), load=FALSE)
