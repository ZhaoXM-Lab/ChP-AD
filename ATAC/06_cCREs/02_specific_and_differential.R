#!/usr/bin/env Rscript
# Identify cell-type markers and AD-associated cCREs with binomial tests.

library(ArchR)
library(SummarizedExperiment)
library(dplyr)
library(rtracklayer)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
out <- file.path(proj_dir, "data/ATAC/output/ccres")
project <- loadArchRProject(file.path(out, "ArchR_peaks"))
flatten <- function(x) bind_rows(lapply(names(x), function(ct) {
    z <- as.data.frame(x[[ct]])
    if (!nrow(z)) return(NULL)
    z$cell_type <- ct
    z$peak <- paste0(z$seqnames, ":", z$start, "-", z$end)
    z
}))
markers <- list()
diffs <- list()
ldsc <- list()
for (group in c("cell_type", "cell_subtype")) {
    m <- getMarkerFeatures(project, useMatrix="PeakMatrix", groupBy=group, bias=c("TSSEnrichment", "log10(nFrags)"), testMethod="binomial", binarize=TRUE)
    saveRDS(m, file.path(out, paste0(group, "_markers.rds")))
    tab <- flatten(getMarkers(m, cutOff="FDR <= 0.05 & Log2FC >= 0.5"))
    ld <- flatten(getMarkers(m, cutOff="FDR <= 0.1 & Log2FC >= 0.5"))
    if (group == "cell_subtype") {
        tab <- tab[tab$cell_type %in% c("CPEC-1", "CPEC-2"), ]
        ld <- ld[ld$cell_type %in% c("CPEC-1", "CPEC-2"), ]
    }
    markers[[group]] <- tab
    ldsc[[group]] <- ld
    types <- if (group == "cell_type") setdiff(unique(project$cell_type), c("Neuron","OPC","CD8+T")) else c("CPEC-1","CPEC-2")
    for (ct in types) {
        cells <- getCellNames(project)[getCellColData(project)[[group]] == ct]
        x <- project[cells, ]
        m <- getMarkerFeatures(x, useMatrix="PeakMatrix", groupBy="disease", useGroups="AD", bgdGroups="WT", bias=c("TSSEnrichment", "log10(nFrags)"), testMethod="binomial", binarize=TRUE)
        saveRDS(m, file.path(out, paste0(ct, "_differential.rds")))
        tab <- flatten(getMarkers(m, cutOff="FDR <= 0.05 & abs(Log2FC) >= 0.5"))
        if (nrow(tab)) {
            tab$cell_type <- ct
            tab$direction <- ifelse(tab$Log2FC > 0, "AD", "WT")
            diffs[[ct]] <- tab
        }
    }
}
marker <- bind_rows(markers) %>% distinct(cell_type, peak, .keep_all=TRUE)
diff <- bind_rows(diffs)
diff$cell_type_specific <- paste(diff$cell_type, diff$peak) %in% paste(marker$cell_type, marker$peak)
write.csv(marker, file.path(out, "specific_ccres.csv"), row.names=FALSE)
write.csv(diff, file.path(out, "differential_ccres.csv"), row.names=FALSE)
write.csv(bind_rows(ldsc) %>% distinct(cell_type, peak, .keep_all=TRUE), file.path(out, "ldsc_ccres_fdr01.csv"), row.names=FALSE)
