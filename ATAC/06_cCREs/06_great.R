#!/usr/bin/env Rscript
# Run GREAT on marker, differential and NMF peak sets against the union background.

library(rGREAT)
library(GenomicRanges)
proj_dir <- Sys.getenv("PROJ_DIR")
stopifnot(nzchar(proj_dir))
base <- file.path(proj_dir, "data/ATAC/output/ccres")
out <- file.path(base, "great")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
bg <- readRDS(file.path(base, "union_peaks.rds"))
run_set <- function(peaks, name) {
    if (!length(peaks)) return(invisible(NULL))
    result <- great(GRanges(unique(peaks)), gene_sets="GO:BP", tss_source="TxDb.Hsapiens.UCSC.hg38.knownGene", min_gene_set_size=5, background=bg)
    write.csv(getEnrichmentTable(result, min_region_hits=2), file.path(out, paste0(name, ".csv")), row.names=FALSE)
}
marker <- read.csv(file.path(base, "specific_ccres.csv"))
diff <- read.csv(file.path(base, "differential_ccres.csv"))
for (ct in unique(marker$cell_type)) run_set(marker$peak[marker$cell_type == ct], paste0(ct, "_markers"))
for (ct in unique(diff$cell_type)) for (d in c("AD", "WT")) run_set(diff$peak[diff$cell_type == ct & diff$direction == d & diff$cell_type_specific], paste0(ct, "_", d))
module_file <- file.path(base, "nmf/module_peaks.csv")
if (file.exists(module_file)) {
    modules <- read.csv(module_file)
    for (m in unique(modules$module)) run_set(modules$peak[modules$module == m], paste0("NMF_", m))
}
