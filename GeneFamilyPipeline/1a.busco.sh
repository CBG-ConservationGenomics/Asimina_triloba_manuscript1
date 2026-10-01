#!/bin/bash
set -e

# Set lineage and mode
LINEAGE="embryophyta_odb10"
MODE="protein"
THREADS=30

# Input + output folders
SEQ_DIR="protein_sequences"
OUT_PARENT="busco_results"

mkdir -p "$OUT_PARENT"

# Loop over each FASTA file in protein_sequences
for fasta in "$SEQ_DIR"/*.fa "$SEQ_DIR"/*.fasta; do
    [ -e "$fasta" ] || continue

    # Get file name only (strip path)
    fname=$(basename "$fasta")

    # Strip extension
    BASENAME="${fname%%.*}"

    echo "Running BUSCO on $fasta (basename: $BASENAME)..."

    busco \
        -i "$fasta" \
        -l "$LINEAGE" \
        -m "$MODE" \
        -c "$THREADS" \
        -o "$BASENAME" \
        --out_path "$OUT_PARENT" \
        -f
done

echo "All BUSCO analyses complete."

