
#!/usr/bin/env bash
set -e

# ==========================================
# InterProScan input preparation
#
# TEST MODE:
#   set SUBSET=50 to prepare only 50 random orthogroups
#
# FULL MODE:
#   change SUBSET=50 to SUBSET=ALL
#   to prepare every orthogroup FASTA
# ==========================================

SUBSET=ALL
OUT_DIR="interproscan_input"

mkdir -p "$OUT_DIR"

# Find newest OrthoFinder results
OG_DIR=$(ls -td orthofinder/Results_*/Orthogroup_Sequences | head -n 1)
RESULTS_DIR=$(dirname "$OG_DIR")
ORTHOGROUPS_TSV="$RESULTS_DIR/Orthogroups/Orthogroups.tsv"

echo "Preparing InterProScan input"
echo "Source FASTAs:     $OG_DIR"
echo "Orthogroups table: $ORTHOGROUPS_TSV"
echo "Output directory:  $OUT_DIR"
echo

# Remove old contents so reruns are clean
echo "Clearing old contents from $OUT_DIR"
rm -f "$OUT_DIR"/*.fa "$OUT_DIR"/Orthogroups.tsv "$OUT_DIR"/subset.list

# Get orthogroup IDs from Orthogroups.tsv
tail -n +2 "$ORTHOGROUPS_TSV" | cut -f1 > og.list

# Select subset or all orthogroups
if [ "$SUBSET" != "ALL" ]; then
    echo "Selecting random subset of $SUBSET orthogroups..."
    shuf og.list | head -n "$SUBSET" > subset.list
else
    echo "Using all orthogroups..."
    cp og.list subset.list
fi

echo
echo "Cleaning selected FASTAs..."

while read -r og; do
    in_file="$OG_DIR/${og}.fa"
    out_file="$OUT_DIR/${og}.fa"

    if [ -f "$in_file" ]; then
        echo "Cleaning $(basename "$in_file")"
        awk '
            /^>/ { print; next }
            { gsub(/\*/, "", $0); print }
        ' "$in_file" > "$out_file"
    else
        echo "Warning: missing FASTA for $og"
    fi
done < subset.list

echo
echo "Copying Orthogroups.tsv and subset list..."
cp "$ORTHOGROUPS_TSV" "$OUT_DIR/"
cp subset.list "$OUT_DIR/"

rm -f og.list subset.list

echo
echo "=== Done! ==="
if [ "$SUBSET" != "ALL" ]; then
    echo "Prepared $SUBSET cleaned orthogroup FASTAs in: $OUT_DIR"
    echo "Saved selected orthogroups in: $OUT_DIR/subset.list"
else
    echo "Prepared all cleaned orthogroup FASTAs in: $OUT_DIR"
fi
