#!/usr/bin/env bash
set -euo pipefail
: "${PROJ_DIR:?Set PROJ_DIR}"
: "${LDSC_DIR:?Set LDSC_DIR to the LDSC source directory}"
python_bin="${LDSC_PYTHON:-python}"
base="$PROJ_DIR/data/ATAC/output/ldsc"
ref="$PROJ_DIR/data/reference/ldsc"
mkdir -p "$base/ldscores"
shopt -s nullglob
beds=("$base/beds/"*.bed)
((${#beds[@]})) || { echo "No annotation BED files found" >&2; exit 1; }
for bed in "${beds[@]}"; do
    cell=$(basename "$bed" .bed)
    for chr in $(seq 1 22); do
        "$python_bin" "$LDSC_DIR/ldsc.py" --l2 \
            --bfile "$ref/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr" \
            --annot "$base/annotations/$cell.$chr.annot.gz" \
            --ld-wind-cm 1 --print-snps "$ref/hm3.snplist" \
            --out "$base/ldscores/$cell.$chr"
    done
done
