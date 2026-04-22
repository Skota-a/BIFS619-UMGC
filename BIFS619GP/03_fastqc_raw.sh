#!/usr/bin/env bash
# =============================================================================
# 03_fastqc_raw.sh
#
# Runs FastQC on raw FASTQ files before trimming. This gives the starting
# quality profile of the data so you can compare before vs after cleanup.
#
# Main things worth checking:
#   - Per-base quality: Illumina reads often drop near the 3' end.
#   - Adapter content: indicates leftover library adapters.
#   - GC content: human RNA-seq is usually centered near expected range.
#   - Duplication: can be high in RNA-seq because abundant genes dominate.
#
# Later, compare these reports with trimmed-read FastQC results to show whether
# filtering actually improved the data.
#
# Outputs in results/03_fastqc_raw/:
#   One .html report and one extracted _fastqc/ folder per FASTQ file.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

RAW_DIR="${PROJECT_DIR}/raw"
OUT_DIR="${PROJECT_DIR}/results/03_fastqc_raw"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Collect compressed FASTQ files in raw/.
# If your files are uncompressed, change *.fastq.gz to *.fastq
mapfile -t RAW_FILES < <(find "$RAW_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort)

if [[ ${#RAW_FILES[@]} -eq 0 ]]; then
    echo "[ERROR] No .fastq.gz files found in $RAW_DIR"
    echo "        Run 00_download_sra.sh first."
    exit 1
fi

echo "[$(date '+%F %T')] Running FastQC on ${#RAW_FILES[@]} file(s)"
for f in "${RAW_FILES[@]}"; do
    echo "  -> $(basename "$f")"
done

fastqc \
    --threads  "$THREADS" \
    --outdir   "$OUT_DIR" \
    --extract \
    "${RAW_FILES[@]}" \
    2>&1 | tee "${LOG_DIR}/03_fastqc_raw.log"

echo ""
echo "[$(date '+%F %T')] Done. Reports saved to $OUT_DIR"
echo "  Open any .html file in a browser to inspect results."
echo "  MultiQC in step 07 will combine them."
