#!/bin/bash

# ==============================================================================
# Purpose:
#   Calculates basic genome assembly statistics for each species to summarize
#   assembly size and contiguity before downstream gene family analysis.
#
# Inputs:
#   - genome_sequences/*.fa   Genome assembly FASTA files
#   - FastaSeqStats           Perl script used to calculate assembly statistics
#
# Outputs:
#   - preprocessing_results/genome_stats.tsv
#     Tab-separated summary containing sequence count, total assembly size (Gb),
#     longest sequence (Mb), N50 (Mb), and L50 for each genome assembly.
# ==============================================================================

# Create preprocessing results directory if it does not already exist
mkdir -p preprocessing_results

# Output file
OUTPUT="preprocessing_results/genome_stats.tsv"

# Create table header
echo -e "Species\tSequences\tTotal_Gb\tLongest_Mb\tN50_Mb\tL50" > "$OUTPUT"

# Calculate genome assembly statistics for each species
for fasta in genome_sequences/*.fa
do
    species=$(basename "$fasta" .fa)

    perl FastaSeqStats -i "$fasta" > temp_stats.txt

    seqs=$(awk '/Sequence Count:/ {print $3}' temp_stats.txt)
    total=$(awk '/Total Length:/ {print $3}' temp_stats.txt)
    longest=$(awk '/Longest sequence:/ {print $3}' temp_stats.txt)
    n50=$(awk '/N50 length:/ {print $3}' temp_stats.txt)
    l50=$(awk '/N50 index:/ {print $3}' temp_stats.txt)

    # Convert bp to Gb/Mb and round to 2 decimal places
    total_gb=$(awk -v x="$total" 'BEGIN {printf "%.2f", x/1000000000}')
    longest_mb=$(awk -v x="$longest" 'BEGIN {printf "%.2f", x/1000000}')
    n50_mb=$(awk -v x="$n50" 'BEGIN {printf "%.2f", x/1000000}')

    # Add species to results table
    echo -e "$species\t$seqs\t$total_gb\t$longest_mb\t$n50_mb\t$l50" >> "$OUTPUT"
done

# Remove temporary file
rm -f temp_stats.txt

echo
echo "Done! Results saved to:"
echo "$OUTPUT"
echo

# Display results nicely
column -t -s $'\t' "$OUTPUT"
