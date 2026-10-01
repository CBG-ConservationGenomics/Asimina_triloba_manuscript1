# to run outside of script: conda activate orthofinder_env
#!/bin/bash
set -e

INPUT_DIR="protein_sequences"
OUTPUT_DIR="orthofinder"
THREADS=30

[ -d "$INPUT_DIR" ] || { echo "Missing folder: $INPUT_DIR"; exit 1; }

orthofinder \
  -f "$INPUT_DIR" \
  -o "$OUTPUT_DIR" \
  -t "$THREADS" \
  -a "$THREADS" \
  -S diamond

echo "=== OrthoFinder finished! Check: $OUTPUT_DIR ==="
