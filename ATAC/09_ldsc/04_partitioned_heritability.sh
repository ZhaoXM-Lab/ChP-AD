#!/usr/bin/env bash
set -euo pipefail
: "${PROJ_DIR:?Set PROJ_DIR}"
: "${LDSC_DIR:?Set LDSC_DIR}"
python_bin="${LDSC_PYTHON:-python}"
base="$PROJ_DIR/data/ATAC/output/ldsc"
ref="$PROJ_DIR/data/reference/ldsc"
shopt -s nullglob
beds=("$base/beds/"*.bed)
((${#beds[@]})) || { echo "No annotation BED files found" >&2; exit 1; }
while IFS=$'\t' read -r trait sumstats; do
    [[ "$trait" == trait || -z "$trait" ]] && continue
    for bed in "${beds[@]}"; do
        cell=$(basename "$bed" .bed)
        mkdir -p "$base/results/$cell"
        "$python_bin" "$LDSC_DIR/ldsc.py" --h2 "$PROJ_DIR/$sumstats" \
            --ref-ld-chr "$base/ldscores/$cell.,$ref/baselineLD_v2.1/baselineLD." \
            --w-ld-chr "$ref/weights/weights.hm3_noMHC." \
            --frqfile-chr "$ref/frequencies/1000G.EUR.QC." \
            --overlap-annot --print-coefficients --out "$base/results/$cell/$trait"
    done
done < "$PROJ_DIR/ATAC/00_reference/gwas.tsv"
