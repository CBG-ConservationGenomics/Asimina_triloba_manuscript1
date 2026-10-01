#!/usr/bin/env bash
set -euo pipefail

CAFE_BIN="$HOME/miniconda3/envs/cafe5/bin/cafe5"

INPUT="cafe_gene_families_filtered.tsv"
TREE="SpeciesTree_ultrametric.txt"
OUT="cafe_zeroRoot"

echo "[1/2] Running CAFE5..."

"$CAFE_BIN" \
    -i "$INPUT" \
    -t "$TREE" \
    -o "$OUT" \
    -k 2 \
    -z
echo
echo "[2/2] Done!"
echo "Output directory: $OUT"
