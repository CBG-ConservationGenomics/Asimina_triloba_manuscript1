#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 3b.filter_interpro_output.sh
#
# PURPOSE
# -------
# This script performs quality control on InterProScan output files and prepares
# annotation tables for downstream analyses.
#
# Specifically, it:
#
#   1. Scans every TSV file in interproscan_output/.
#   2. Removes ONLY empty TSV files (no InterPro annotations).
#   3. Keeps ALL annotated orthogroups, including those with transposon-
#      related annotations.
#   4. Creates:
#        - include_orthogroups.txt
#            List of orthogroups retained for downstream analyses.
#
#        - excluded_orthogroups_with_reason.tsv
#            Lists only orthogroups excluded because they contained no
#            annotations (empty TSV).
#
#        - transposon_associated_orthogroups.tsv
#            Lists retained orthogroups whose annotations contain the word
#            "transpos". This file is informative only and allows transposon-
#            associated orthogroups to be filtered later if desired.
#
#        - all_go_unique.tsv
#            Table containing orthogroup IDs and their associated GO terms for
#            every retained orthogroup.
#
# NOTE
# ----
# Transposon-associated orthogroups are intentionally retained so that they can
# be evaluated separately during downstream analyses (e.g. GO enrichment,
# expanded gene family interpretation, candidate gene exploration).
###############################################################################

# =========================
# SETTINGS
# =========================
BASE_DIR="$(pwd)"
INPUT_DIR="$BASE_DIR/interproscan_output"
MOVE_BAD=true   # change to true if you want bad TSVs moved into rejected_tsv/

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

REJECT_DIR="$INPUT_DIR/rejected_tsv"
if [ "$MOVE_BAD" = true ]; then
    mkdir -p "$REJECT_DIR"
fi

INCLUDE_LIST="include_orthogroups.txt"
EXCLUDE_REPORT="excluded_orthogroups_with_reason.tsv"
TRANSPOS_REPORT="transposon_associated_orthogroups.tsv"
OUT="all_go_unique.tsv"
GO_COL=14

: > "$INCLUDE_LIST"
: > "$EXCLUDE_REPORT"
: > "$TRANSPOS_REPORT"
: > "$OUT"

echo -e "Orthogroup\tReason\tFile" > "$EXCLUDE_REPORT"
echo -e "Orthogroup\tFile" > "$TRANSPOS_REPORT"

#############################################
### PART 1 — IDENTIFY GOOD/BAD TSV FILES  ###
#############################################

total=0
excluded=0
kept=0
transpos=0

shopt -s nullglob
for tsv in "$INPUT_DIR"/*.tsv; do
    base=$(basename "$tsv")

    # skip output/report files if rerunning
    case "$base" in
        all_go_unique.tsv|excluded_orthogroups_with_reason.tsv|transposon_associated_orthogroups.tsv)
            continue
            ;;
    esac

    og="${base%.tsv}"
    ((total+=1))

    # empty file?
    if [ ! -s "$tsv" ]; then
        echo -e "${og}\tempty\t${tsv}" >> "$EXCLUDE_REPORT"
        ((excluded+=1))

        if [ "$MOVE_BAD" = true ]; then
            mv "$tsv" "$REJECT_DIR/"
        fi

        continue
    fi

    # keep all annotated orthogroups
    echo "$og" >> "$INCLUDE_LIST"
    ((kept+=1))

    # record transposon-associated orthogroups for later reference
    if grep -Iqi 'transpos' "$tsv"; then
        echo -e "${og}\t${tsv}" >> "$TRANSPOS_REPORT"
        ((transpos+=1))
    fi
done

sort -u "$INCLUDE_LIST" -o "$INCLUDE_LIST"

# sort exclusion report but keep header on top
{
    head -n 1 "$EXCLUDE_REPORT"
    tail -n +2 "$EXCLUDE_REPORT" | sort -u
} > "${EXCLUDE_REPORT}.tmp" && mv "${EXCLUDE_REPORT}.tmp" "$EXCLUDE_REPORT"

# sort transpos report but keep header
{
    head -n 1 "$TRANSPOS_REPORT"
    tail -n +2 "$TRANSPOS_REPORT" | sort -u
} > "${TRANSPOS_REPORT}.tmp" && mv "${TRANSPOS_REPORT}.tmp" "$TRANSPOS_REPORT"

echo "Total .tsv files checked: $total"
echo "Excluded (empty): $excluded"
echo "Kept: $kept"
echo "Transposon-associated (retained): $transpos"

#############################################
### PART 2 — BUILD GO TABLE               ###
#############################################

echo "Building GO table using $INCLUDE_LIST ..."

while read -r og; do
    [ -n "$og" ] || continue
    tsv="$INPUT_DIR/${og}.tsv"
    [ -f "$tsv" ] || continue

    awk -v og="$og" -v c="$GO_COL" '
        BEGIN {FS=OFS="\t"}
        $c ~ /GO:/ { print og, $c }
    ' "$tsv" >> "$OUT"
done < "$INCLUDE_LIST"

sort -u "$OUT" -o "$OUT"

echo "Created: $OUT"
echo "Created keep list: $INCLUDE_LIST"
echo "Created exclusion report: $EXCLUDE_REPORT"
echo "Created transposon report: $TRANSPOS_REPORT"

if [ "$MOVE_BAD" = true ]; then
    echo "Empty TSVs were moved to: $REJECT_DIR"
fi

echo "All done!"
