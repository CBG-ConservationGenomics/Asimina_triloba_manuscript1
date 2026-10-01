#!/usr/bin/env bash

set -u
set -o pipefail

INPUT_DIR="interproscan_input"
OUTPUT_DIR="interproscan_output"
IPS_PATH="/opt/interproscan/interproscan-5.75-106.0/interproscan.sh"

# Number of InterProScan analyses running simultaneously.
# Start with 10. You can increase it later if the server remains stable.
JOBS=40

# CPUs assigned to each individual InterProScan analysis.
CPU=1

RUN_ID=$(date +"%Y-%m-%d_%H-%M-%S")
JOBLOG="interproscan_parallel_joblog_${RUN_ID}.tsv"
MISSING_LIST="interproscan_missing_${RUN_ID}.list"
FAILED_LIST="interproscan_failed_${RUN_ID}.txt"

# Confirm required directories and programs exist.
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory does not exist: $INPUT_DIR" >&2
    exit 1
fi

if [[ ! -x "$IPS_PATH" ]]; then
    echo "ERROR: InterProScan executable was not found or is not executable:" >&2
    echo "$IPS_PATH" >&2
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "ERROR: GNU Parallel is not available." >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

TOTAL_INPUTS=$(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.fa" |
    wc -l
)

EXISTING_OUTPUTS=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.tsv" |
    wc -l
)

# Build a null-delimited list containing only FASTA files
# that do not already have a corresponding TSV output.
: > "$MISSING_LIST"

while IFS= read -r -d '' file; do
    base=$(basename "$file" .fa)
    out="$OUTPUT_DIR/${base}.tsv"

    # An empty TSV can be a valid result, so existence is used,
    # rather than requiring the TSV to be nonempty.
    if [[ ! -e "$out" ]]; then
        printf '%s\0' "$file" >> "$MISSING_LIST"
    fi
done < <(
    find "$INPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.fa" \
        -print0
)

MISSING_COUNT=$(tr -cd '\0' < "$MISSING_LIST" | wc -c)

echo "=================================================="
echo "InterProScan resume run"
echo "=================================================="
echo "Start time:             $(date)"
echo "Input directory:        $INPUT_DIR"
echo "Output directory:       $OUTPUT_DIR"
echo "Total FASTA files:      $TOTAL_INPUTS"
echo "Existing TSV outputs:   $EXISTING_OUTPUTS"
echo "Files to process:       $MISSING_COUNT"
echo "Parallel jobs:          $JOBS"
echo "CPUs per job:           $CPU"
echo "Job log:                $JOBLOG"
echo "=================================================="
echo

if [[ "$MISSING_COUNT" -eq 0 ]]; then
    echo "Nothing to run. Every FASTA already has a TSV output."
    rm -f "$MISSING_LIST"
    exit 0
fi

export OUTPUT_DIR
export IPS_PATH
export CPU

run_interproscan() {
    file="$1"
    base=$(basename "$file" .fa)
    out="$OUTPUT_DIR/${base}.tsv"

    # Include the GNU Parallel sequence and shell PID to prevent
    # simultaneous jobs from using the same temporary filename.
    temp_out="${out}.tmp.${PARALLEL_SEQ}.$$.tsv"

    if [[ -e "$out" ]]; then
        echo "[SKIP] $base — output appeared before job started"
        return 0
    fi

    rm -f "$temp_out"

    echo "[START] $base"

    "$IPS_PATH" \
        -i "$file" \
        -f tsv \
        -appl Pfam,PANTHER \
        --iprlookup \
        --goterms \
        -cpu "$CPU" \
        -o "$temp_out"

    status=$?

    if [[ "$status" -eq 0 && -e "$temp_out" ]]; then
        mv "$temp_out" "$out"
        echo "[DONE] $base"
        return 0
    fi

    echo "[FAILED] $base — InterProScan exit status: $status" >&2
    rm -f "$temp_out"
    return "$status"
}

export -f run_interproscan

parallel \
    -0 \
    -j "$JOBS" \
    --joblog "$JOBLOG" \
    --line-buffer \
    run_interproscan {} \
    < "$MISSING_LIST"

PARALLEL_STATUS=$?

FINAL_OUTPUTS=$(
    find "$OUTPUT_DIR" \
        -maxdepth 1 \
        -type f \
        -name "*.tsv" |
    wc -l
)

REMAINING=$((TOTAL_INPUTS - FINAL_OUTPUTS))

# Column 7 of a standard GNU Parallel job log is Exitval.
awk -F '\t' '
    NR > 1 && $7 != 0 {
        print $9
    }
' "$JOBLOG" > "$FAILED_LIST"

FAILED_COUNT=$(wc -l < "$FAILED_LIST")

echo
echo "=================================================="
echo "InterProScan resume run ended"
echo "=================================================="
echo "End time:               $(date)"
echo "Parallel exit status:   $PARALLEL_STATUS"
echo "Total input FASTAs:     $TOTAL_INPUTS"
echo "Current TSV outputs:    $FINAL_OUTPUTS"
echo "Remaining without TSV:  $REMAINING"
echo "Failed jobs this run:   $FAILED_COUNT"
echo "Job log:                $JOBLOG"
echo "Failed-job list:        $FAILED_LIST"
echo "=================================================="

rm -f "$MISSING_LIST"

if [[ "$PARALLEL_STATUS" -ne 0 ]]; then
    exit "$PARALLEL_STATUS"
fi

exit 0
