#!/bin/bash

# ============================================================
# Summarize gene-family expansions and contractions for each
# terminal species branch from CAFE5.
#
# Run this script from the directory containing cafe_zeroRoot/
#
# Inputs:
#   cafe_zeroRoot/Gamma_change.tab
#   cafe_zeroRoot/Gamma_branch_probabilities.tab
#
# Output:
#   cafe_zeroRoot/species_expansion_contraction_summary.tsv
#
# Definitions:
#
#   Total expanded family:
#       Change > 0 in Gamma_change.tab
#
#   Total contracted family:
#       Change < 0 in Gamma_change.tab
#
#   Rapidly changing expanded family:
#       Change > 0 AND CAFE5 branch probability < 0.05
#
#   Rapidly changing contracted family:
#       Change < 0 AND CAFE5 branch probability < 0.05
#
# IMPORTANT:
#   Values in Gamma_branch_probabilities.tab are CAFE5 branch
#   probabilities and should NOT be described as conventional
#   p-values.
#
#   A branch probability < 0.05 corresponds to branches marked
#   with "*" in Gamma_asr.tre. This was independently verified
#   for the Asimina_triloba branch in this analysis:
#
#       291 starred families total
#       = 15 expansions + 276 contractions
#
# Only terminal species branches are summarized.
# Internal CAFE5 nodes are excluded.
# ============================================================

CHANGE="cafe_zeroRoot/Gamma_change.tab"
PROB="cafe_zeroRoot/Gamma_branch_probabilities.tab"
OUTPUT="cafe_zeroRoot/species_expansion_contraction_summary.tsv"

# CAFE5 branch-probability threshold
THRESHOLD=0.05


# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

for file in "$CHANGE" "$PROB"; do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Cannot find $file"
        exit 1
    fi
done


# ------------------------------------------------------------
# Count total and rapidly changing expansions/contractions
# ------------------------------------------------------------

awk -v threshold="$THRESHOLD" '

BEGIN {
    FS = OFS = "\t"
}


# ============================================================
# FIRST FILE: Gamma_change.tab
# ============================================================

ARGIND == 1 {

    # Read header and store terminal species names
    if (FNR == 1) {

        # Columns 2-7 correspond to the six terminal species
        for (i = 2; i <= 7; i++) {

            species[i] = $i

            # Remove CAFE5 node labels such as <1>, <2>, etc.
            sub(/<[0-9]+>$/, "", species[i])
        }

        next
    }

    family = $1

    for (i = 2; i <= 7; i++) {

        c = $i

        # Store change for comparison with branch probabilities
        change[family,i] = c

        # Count all inferred changes
        if (c > 0)
            total_expanded[i]++

        else if (c < 0)
            total_contracted[i]++
    }

    next
}


# ============================================================
# SECOND FILE: Gamma_branch_probabilities.tab
# ============================================================

ARGIND == 2 {

    if (FNR == 1)
        next

    family = $1

    for (i = 2; i <= 7; i++) {

        probability = $i
        c = change[family,i]

        # Ignore unavailable branch probabilities
        if (probability == "N/A" || probability == "")
            continue

        # CAFE5 rapidly changing branch
        if (probability < threshold) {

            if (c > 0)
                rapid_expanded[i]++

            else if (c < 0)
                rapid_contracted[i]++
        }
    }
}


# ============================================================
# OUTPUT
# ============================================================

END {

    print "Species", \
          "Total_Expanded_Families", \
          "Total_Contracted_Families", \
          "Rapidly_Expanded_Families", \
          "Rapidly_Contracted_Families"

    for (i = 2; i <= 7; i++) {

        print species[i], \
              total_expanded[i] + 0, \
              total_contracted[i] + 0, \
              rapid_expanded[i] + 0, \
              rapid_contracted[i] + 0
    }
}

' "$CHANGE" "$PROB" > "$OUTPUT"


# ------------------------------------------------------------
# Display results
# ------------------------------------------------------------

echo
echo "CAFE5 gene-family expansion/contraction summary"
echo
echo "Rapidly changing families:"
echo "CAFE5 branch probability < $THRESHOLD"
echo
column -t -s $'\t' "$OUTPUT"
echo
echo "Results written to:"
echo "$OUTPUT"


#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Summarize CAFE5 gene-family expansions and contractions
# across all terminal and internal branches of the phylogeny.
#
# Input:
#   cafe_zeroRoot/Gamma_change.tab
#
# Output:
#   phylogeny_changevalues.tsv
#
# Definitions:
#   Expanded   = gene family change > 0
#   Contracted = gene family change < 0
#   Unchanged  = gene family change = 0
#
# Includes:
#   - Six terminal species branches
#   - Five internal branches/nodes
# ============================================================


# ------------------------------------------------------------
# 1. Input and output files
# ------------------------------------------------------------

CHANGE="cafe_zeroRoot/Gamma_change.tab"
OUTPUT="phylogeny_changevalues.tsv"


# ------------------------------------------------------------
# 2. Check that input exists
# ------------------------------------------------------------

if [[ ! -f "$CHANGE" ]]; then
    echo "ERROR: Cannot find $CHANGE"
    exit 1
fi


# ------------------------------------------------------------
# 3. Count expansions and contractions
# ------------------------------------------------------------

awk '
BEGIN {
    FS = OFS = "\t"
}

# Save branch/node names from header
NR == 1 {
    for (i = 2; i <= NF; i++) {
        branch[i] = $i
    }
    next
}

# Count changes for every branch
{
    for (i = 2; i <= NF; i++) {

        if ($i > 0) {
            expanded[i]++
        }

        else if ($i < 0) {
            contracted[i]++
        }

        else if ($i == 0) {
            unchanged[i]++
        }
    }
}

# Write summary
END {

    print "Branch", "Expanded", "Contracted", "Unchanged"

    for (i = 2; i <= NF; i++) {

        print branch[i],
              expanded[i] + 0,
              contracted[i] + 0,
              unchanged[i] + 0
    }
}
' "$CHANGE" > "$OUTPUT"


# ------------------------------------------------------------
# 4. Report completion
# ------------------------------------------------------------

echo "Done!"
echo "Output written to: $OUTPUT"
echo
column -t -s $'\t' "$OUTPUT"
~
