#!/bin/bash
set -euo pipefail

# -------- paths --------
BASE_DIR=$(pwd)                                  # where you run the script from
INTERPRO_OUT="interproscan_output"         # relative path
GO_COL=14                                        # adjust ONLY if needed

OUT_MAP="all_go_unique_clean.tsv"

echo "Running TopGO prep from: $BASE_DIR"
echo "Using InterPro output from: $INTERPRO_OUT"

cd "$INTERPRO_OUT"

# ---- output files ----
OUT_EXT="all_go_extended.tsv"
OUT_GO="all_go.tsv"
OUT_GO_UNIQ="all_go_unique.tsv"

: > "$OUT_EXT"
: > "$OUT_GO"

# ---- 1. merge cleaned interpro outputs ----
for f in *.tsv; do
  [ -f "$f" ] || continue
  og="${f%%_cleaned.tsv}"
  awk -v og="$og" 'BEGIN{FS=OFS="\t"}
    $0 !~ /^#/ && tolower($0) !~ /transpos/ { print $0, og }' "$f" >> "$OUT_EXT"
done

# ---- 2. extract OG + GO column ----
awk -v c="$GO_COL" 'BEGIN{FS=OFS="\t"}
  $c != "-" && index($c,"GO:") { print $NF, $c }' "$OUT_EXT" > "$OUT_GO"

# ---- 3. deduplicate ----
LC_ALL=C sort -u "$OUT_GO" > "$OUT_GO_UNIQ"

# ---- 4. clean GO formatting for topGO ----
awk 'BEGIN{FS=OFS="\t"}
{
  gsub(/\r/,"",$2)
  gsub(/\(InterPro\)|\(PANTHER\)/,"",$2)
  gsub(/\|/,",",$2)
  gsub(/[[:space:]]+/,"",$2)
  print $1, $2
}' "$OUT_GO_UNIQ" > "$OUT_MAP"

# ---- 5. copy mapping back to main directory ----
cp "$OUT_MAP" "$BASE_DIR/"

echo "TopGO mapping written to:"
echo "  $BASE_DIR/$OUT_MAP"
