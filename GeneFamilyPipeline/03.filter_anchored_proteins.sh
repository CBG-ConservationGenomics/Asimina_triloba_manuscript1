#!/bin/bash

# ============================================================
# Filter protein FASTAs to chromosome/pseudomolecule-anchored
# proteins only.
#
# IMPORTANT:
# Original files in protein_sequences/ are NEVER modified.
#
# Outputs:
#   protein_sequences_anchored/
#   preprocessing_results/protein_filter_summary.tsv
#   preprocessing_results/protein_filter_by_chromosome.tsv
# ============================================================


# ----------------------------
# Set up directories
# ----------------------------

INPUT_DIR="unfiltered_protein_sequences"
OUTPUT_DIR="protein_sequences"
RESULTS_DIR="preprocessing_results"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$RESULTS_DIR"

SUMMARY="$RESULTS_DIR/protein_filter_summary.tsv"
CHROM_COUNTS="$RESULTS_DIR/protein_filter_by_chromosome.tsv"


# ----------------------------
# Create summary files
# ----------------------------

echo -e "Species\tTotal_proteins\tAnchored_proteins\tRemoved_proteins\tPercent_retained" \
    > "$SUMMARY"

echo -e "Species\tChromosome\tProteins" \
    > "$CHROM_COUNTS"


# ============================================================
# Annona cherimola
#
# Protein IDs:
#   Anche102Chr... = chromosome anchored
#   Anche102Scf... = scaffold/unanchored
#
# Keep: Chr
# Remove: Scf
# ============================================================

echo "=== Annona cherimola ==="

input="$INPUT_DIR/Annona_cherimola.fa"
output="$OUTPUT_DIR/Annona_cherimola.fa"

# Filter FASTA
awk '
    /^>/ {
        keep = ($1 ~ /^>Anche102Chr/)
    }
    keep
' "$input" > "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Annona_cherimola\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each chromosome
grep "^>" "$output" \
    | sed -E 's/^>Anche102Chr([0-9]+).*/Chr\1/' \
    | sort -V \
    | uniq -c \
    | awk '{print "Annona_cherimola\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Annona montana
#
# Protein headers contain OriSeqID:
#   OriSeqID=LG01-LG07 = chromosome/linkage-group anchored
#   OriSeqID=Contig... = unanchored contigs
#
# Keep: LG01-LG07
# Remove: Contig
# ============================================================

echo "=== Annona montana ==="

input="$INPUT_DIR/Annona_montana.fa"
output="$OUTPUT_DIR/Annona_montana.fa"

# Filter FASTA
awk '
    /^>/ {
        keep = ($0 ~ /OriSeqID=LG[0-9]+([[:space:]]|$)/)
    }
    keep
' "$input" > "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Annona_montana\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each linkage group
grep "^>" "$output" \
    | grep -o 'OriSeqID=LG[0-9]*' \
    | cut -d= -f2 \
    | sort -V \
    | uniq -c \
    | awk '{print "Annona_montana\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Asimina triloba
#
# Protein headers contain seq_id:
#   seq_id=Atri1.0C1-C8 = chromosome/pseudomolecule anchored
#   seq_id=Atri1.0S...   = unanchored scaffold
#
# Keep: C1-C8
# Remove: S
# ============================================================

echo "=== Asimina triloba ==="

input="$INPUT_DIR/Asimina_triloba.fa"
output="$OUTPUT_DIR/Asimina_triloba.fa"

# Filter FASTA
awk '
    /^>/ {
        keep = ($0 ~ /seq_id=Atri1\.0C[1-8]([[:space:]]|$)/)
    }
    keep
' "$input" > "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Asimina_triloba\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each chromosome
grep "^>" "$output" \
    | grep -o 'seq_id=Atri1\.0C[0-9]*' \
    | cut -d= -f2 \
    | sort -V \
    | uniq -c \
    | awk '{print "Asimina_triloba\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Lindera megaphylla
#
# Protein headers contain OriSeqID:
#   OriSeqID=chr01-chr12 = nuclear chromosomes
#   OriSeqID=ctg...      = unanchored contigs
#   OriSeqID=Mt          = mitochondrial
#   OriSeqID=Pt          = plastid
#
# Keep ONLY: chr01-chr12
# ============================================================

echo "=== Lindera megaphylla ==="

input="$INPUT_DIR/Lindera_megaphylla.fa"
output="$OUTPUT_DIR/Lindera_megaphylla.fa"

# Filter FASTA
awk '
    /^>/ {
        keep = ($0 ~ /OriSeqID=chr(0[1-9]|1[0-2])([[:space:]]|$)/)
    }
    keep
' "$input" > "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Lindera_megaphylla\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each chromosome
grep "^>" "$output" \
    | grep -o 'OriSeqID=chr[0-9]*' \
    | cut -d= -f2 \
    | sort -V \
    | uniq -c \
    | awk '{print "Lindera_megaphylla\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Magnolia kwangsiensis
#
# All proteins are chromosome-anchored (Chr1-Chr19)
# Retained: 100%
# ============================================================

echo "=== Magnolia kwangsiensis ==="

input="$INPUT_DIR/Magnolia_kwangsiensis.fa"
output="$OUTPUT_DIR/Magnolia_kwangsiensis.fa"

# All proteins are chromosome-anchored, so copy original
cp "$input" "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Magnolia_kwangsiensis\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each chromosome
grep "^>" "$output" \
    | grep -o 'OriSeqID=Chr[0-9]*' \
    | cut -d= -f2 \
    | sort -V \
    | uniq -c \
    | awk '{print "Magnolia_kwangsiensis\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Persea americana
#
# All proteins correspond to chromosomes Pa01-Pa12
# Retained: 100%
# ============================================================

echo "=== Persea americana ==="

input="$INPUT_DIR/Persea_americana.fa"
output="$OUTPUT_DIR/Persea_americana.fa"

# All proteins are chromosome-anchored, so copy original
cp "$input" "$output"

# Count proteins
total=$(grep -c "^>" "$input")
anchored=$(grep -c "^>" "$output")
removed=$((total - anchored))

percent=$(awk -v a="$anchored" -v t="$total" \
    'BEGIN {printf "%.2f", (a/t)*100}')

# Write overall summary
echo -e "Persea_americana\t$total\t$anchored\t$removed\t$percent" \
    >> "$SUMMARY"

# Count proteins on each chromosome
grep "^>" "$output" \
    | sed 's/^>//' \
    | sed -E 's/(Pa[0-9]+)g.*/\1/' \
    | sort -V \
    | uniq -c \
    | awk '{print "Persea_americana\t"$2"\t"$1}' \
    >> "$CHROM_COUNTS"

echo "Original proteins: $total"
echo "Anchored proteins: $anchored"
echo "Removed proteins:  $removed"
echo "Percent retained:  $percent%"
echo


# ============================================================
# Finished
# ============================================================

echo "========================================"
echo "Filtering complete."
echo "========================================"
echo

echo "Anchored FASTAs:"
echo "  $OUTPUT_DIR/"
echo

echo "Results:"
echo "  $SUMMARY"
echo "  $CHROM_COUNTS"
echo

echo "Protein filtering summary:"
column -t -s $'\t' "$SUMMARY"
echo
