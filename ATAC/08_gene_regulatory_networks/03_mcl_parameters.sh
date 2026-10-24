#!/usr/bin/env bash
set -euo pipefail
: "${PROJ_DIR:?Set PROJ_DIR}"
base="$PROJ_DIR/data/ATAC/output/grn"
shopt -s nullglob
files=("$base/"*.abc)
((${#files[@]})) || { echo "No MCL inputs found" >&2; exit 1; }
for input in "${files[@]}"; do
    cell=$(basename "$input" .abc)
    output="$base/mcl/$cell"
    mkdir -p "$output"
    for inflation in $(seq 1.2 0.2 7.0); do
        mcl "$input" --abc -I "$inflation" -o "$output/Clusters_I${inflation}.txt"
    done
done
