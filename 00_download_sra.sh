#!/usr/bin/env bash
# =============================================================================
# 00_download_sra.sh
#
# Downloads raw sequencing data from NCBI SRA and converts it to FASTQ.
#
# FASTQ files are kept uncompressed for simplicity and speed early in the pipeline.
#
# All samples are single-end. `fasterq-dump` usually outputs `_1.fastq`,
# so we check that naming pattern (plus a few common alternatives).
#
# Expected outputs in raw/:
#   SRR7898026_1.fastq
#   SRR7897501_1.fastq
#   SRR1039508_1.fastq
#
# Optional compression after pipeline:
#   pigz -p 8 raw/*.fastq
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

RAW_DIR="${PROJECT_DIR}/raw"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$RAW_DIR" "$LOG_DIR"

# Show files created for a given accession
# Helps confirm naming and output after conversion
show_new_fastqs() {
    local acc="$1"
    echo "[$(date '+%F %T')] Files for ${acc}:"
    find "$RAW_DIR" -maxdepth 2 -name "*${acc}*" \
        \( -name "*.fastq" -o -name "*.fastq.gz" \) 2>/dev/null | sort | \
        while read -r f; do echo "    $(du -sh "$f" | cut -f1)  $f"; done
}

# Check if FASTQ output already exists and is non-empty
files_ready() {
    local acc="$1"
    [[ -f "${RAW_DIR}/${acc}_1.fastq"    && -s "${RAW_DIR}/${acc}_1.fastq"    ]] ||
    [[ -f "${RAW_DIR}/${acc}_1.fastq.gz" && -s "${RAW_DIR}/${acc}_1.fastq.gz" ]] ||
    [[ -f "${RAW_DIR}/${acc}.fastq"      && -s "${RAW_DIR}/${acc}.fastq"      ]] ||
    [[ -f "${RAW_DIR}/${acc}.fastq.gz"   && -s "${RAW_DIR}/${acc}.fastq.gz"   ]]
}

# Download + convert one SRA accession
download_accession() {
    local acc="$1"

    echo ""
    echo "[$(date '+%F %T')] ---- $acc ----"

    # Skip if already done
    if files_ready "$acc"; then
        echo "[$(date '+%F %T')] Already exists — skipping"
        return 0
    fi

    echo "[$(date '+%F %T')] Downloading SRA..."

    prefetch \
        --output-directory "$RAW_DIR" \
        --max-size 50G \
        "$acc" \
        2>> "${LOG_DIR}/00_download_${acc}.log"

    echo "[$(date '+%F %T')] Converting to FASTQ..."

    fasterq-dump \
        --outdir         "$RAW_DIR" \
        --temp           "$RAW_DIR" \
        --threads        "$THREADS" \
        --split-files \
        --skip-technical \
        "${RAW_DIR}/${acc}/${acc}.sra" \
        2>> "${LOG_DIR}/00_download_${acc}.log"

    show_new_fastqs "$acc"

    # Verify output exists
    if files_ready "$acc"; then
        echo "[$(date '+%F %T')] Done: $acc"
    else
        echo "[ERROR] Missing FASTQ for $acc"
        echo "Log: ${LOG_DIR}/00_download_${acc}.log"
        exit 1
    fi
}

echo "[$(date '+%F %T')] Starting downloads"
echo "  RAW_DIR : $RAW_DIR"
echo "  THREADS : $THREADS"

download_accession "SRR7898026"
download_accession "SRR7897501"
download_accession "SRR1039508"

echo ""
echo "[$(date '+%F %T')] All downloads complete"

echo ""
echo "Files in raw/:"
find "$RAW_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort | \
    while read -r f; do echo "    $(du -sh "$f" | cut -f1)  $(basename "$f")"; done

echo ""
du -sh "$RAW_DIR" | awk '{print "Total size: " $1}'

echo ""
echo "Next: bash scripts/01_download_reference.sh"