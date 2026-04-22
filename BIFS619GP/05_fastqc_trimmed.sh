#!/usr/bin/env bash
# =============================================================================
# 05_fastqc_trimmed.sh
#
# Runs FastQC again on reads after fastp trimming.
#
# This is the post-cleaning QC check and should be compared with step 03
# (raw reads). Running QC before and after trimming is standard practice and
# gives clear evidence that filtering improved read quality.
#
# Common improvements:
#   - Adapter content reduced or removed
#   - Better quality at the 3' end
#   - Read lengths become more variable after trimming
#
# Step 07 (MultiQC) will combine both raw and trimmed FastQC results into one
# summary report.
#
# Outputs in results/05_fastqc_trimmed/:
#   One .html + one _fastqc/ folder per trimmed FASTQ
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

TRIM_DIR="${PROJECT_DIR}/results/04_fastp"
OUT_DIR="${PROJECT_DIR}/results/05_fastqc_trimmed"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Collect trimmed FASTQ files only
# This avoids unrelated files from paired-end workflows
mapfile -t TRIMMED < <(find "$TRIM_DIR" -maxdepth 1 -name "*.trimmed.fastq.gz" | sort)

if [[ ${#TRIMMED[@]} -eq 0 ]]; then
    echo "[ERROR] No trimmed FASTQ files found in $TRIM_DIR"
    echo "        Run 04_fastp_trim.sh first."
    exit 1
fi

echo "[$(date '+%F %T')] Running FastQC on ${#TRIMMED[@]} trimmed file(s):"
for f in "${TRIMMED[@]}"; do
    echo "  -> $(basename "$f")"
done

fastqc \
    --threads "$THREADS" \
    --outdir  "$OUT_DIR" \
    --extract \
    "${TRIMMED[@]}" \
    2>&1 | tee "${LOG_DIR}/05_fastqc_trimmed.log"

echo ""
echo "[$(date '+%F %T')] Done. Post-trim reports in $OUT_DIR"
