#!/usr/bin/env Rscript
# Compute control-reference epigenomic information and module accessibility per nucleus.

library(ArchR)
library(GenomicRanges)
library(dplyr)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir, "data/ATAC/output/ccres")
project <- loadArchRProject(file.path(base, "ArchR_peaks"))
diff <- read.csv(file.path(base, "differential_ccres.csv"))
scores <- list()
for (ct in unique(diff$cell_type)) {
    panel <- unique(diff$peak[diff$cell_type == ct & diff$direction == "WT" & diff$FDR < 0.05])
    if (!length(panel)) next
    field <- "ControlReference"
    x <- addFeatureCounts(project, features=GRanges(panel), name=field)
    meta <- as.data.frame(getCellColData(x))
    keep <- if (ct %in% c("CPEC-1", "CPEC-2")) meta$cell_subtype == ct else meta$cell_type == ct
    meta <- meta[keep, ]
    scores[[ct]] <- data.frame(cell=rownames(meta), cell_type=ct, disease=meta$disease, score=meta[[paste0(field,"Counts")]]/meta$nFrags)
}
write.csv(bind_rows(scores), file.path(base, "epigenomic_information.csv"), row.names=FALSE)
modules <- read.csv(file.path(base, "nmf/module_peaks.csv"))
values <- list()
for (m in unique(modules$module)) {
    x <- addFeatureCounts(project, features=GRanges(modules$peak[modules$module == m]), name="Module")
    meta <- as.data.frame(getCellColData(x))
    values[[m]] <- data.frame(cell=rownames(meta), cell_type=meta$cell_type, disease=meta$disease, module=m, accessibility=meta$ModuleCounts/meta$nFrags)
}
values <- bind_rows(values)
write.csv(values, file.path(base, "nmf/module_cell_accessibility.csv"), row.names=FALSE)
tests <- values %>% group_by(module, cell_type) %>% summarise(p=if (all(c("AD","WT") %in% disease)) wilcox.test(accessibility[disease == "AD"], accessibility[disease == "WT"], alternative="two.sided")$p.value else NA_real_, .groups="drop")
tests$p_adj <- p.adjust(tests$p, "BH")
write.csv(tests, file.path(base, "nmf/module_disease_tests.csv"), row.names=FALSE)
