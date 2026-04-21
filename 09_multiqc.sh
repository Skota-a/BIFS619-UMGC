#!/usr/bin/env bash
# =============================================================================
# 09_multiqc.sh
#
# Combines QC outputs from earlier steps into one MultiQC report.
#
# Instead of checking many separate FastQC, fastp, STAR, and featureCounts
# files, MultiQC puts the main metrics in one place for easy comparison across
# samples.
#
# Typical sections include:
#   - FastQC (raw and trimmed)
#   - fastp trimming summaries
#   - STAR alignment rates
#   - featureCounts assignment stats
#   - samtools summaries
#
# This is usually the easiest report to share in lab meetings, manuscripts,
# or thesis appendices.
#
# Output:
#   results/09_multiqc/multiqc_report.html
#   results/09_multiqc/multiqc_data/
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

RESULTS_DIR="${PROJECT_DIR}/results"
OUT_DIR="${RESULTS_DIR}/09_multiqc"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

echo "[$(date '+%F %T')] Aggregating QC reports with MultiQC"
echo "  Scanning: $RESULTS_DIR"

multiqc \
    --force \
    --outdir   "$OUT_DIR" \
    --filename "multiqc_report.html" \
    --title    "RNA-seq QC Report — GRCh38/GENCODE v46 pipeline" \
    --comment  "Samples: SRR7898026 (SE), SRR7897501 (SE), SRR8492612 (PE)" \
    --dirs \
    "${RESULTS_DIR}/03_fastqc_raw" \
    "${RESULTS_DIR}/04_fastp" \
    "${RESULTS_DIR}/05_fastqc_trimmed" \
    "${RESULTS_DIR}/06_star_align" \
    "${RESULTS_DIR}/08_featurecounts" \
    2>&1 | tee "${LOG_DIR}/09_multiqc.log"

echo ""
echo "[$(date '+%F %T')] MultiQC report written to:"
echo "  ${OUT_DIR}/multiqc_report.html"
echo ""
echo "  Open it in a browser with:"
echo "    xdg-open ${OUT_DIR}/multiqc_report.html   # Linux"
echo "    open ${OUT_DIR}/multiqc_report.html        # macOS"