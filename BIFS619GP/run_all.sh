#!/usr/bin/env bash
# =============================================================================
# run_all.sh
#
# Runs the full RNA-seq workflow from raw SRA download to final count matrix
# and QC summary plots.
#
# Each step is kept as a separate script, which makes troubleshooting easier
# and lets you rerun only the stage that changed.
#
# Usage:
#   conda activate rnaseq
#   cd /path/to/rnaseq_pipeline
#   THREADS=8 ./scripts/run_all.sh
#
# Example: rerun trimming only
#   THREADS=8 bash scripts/04_fastp_trim.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
export THREADS="${THREADS:-8}"

log_stage() {
    echo ""
    echo "=================================================================="
    echo "  STAGE: $1"
    echo "  $(date '+%F %T')"
    echo "=================================================================="
}

log_stage "00 — Download SRA data"
bash "${SCRIPT_DIR}/00_download_sra.sh"

log_stage "01 — Download GRCh38 + GENCODE GTF"
bash "${SCRIPT_DIR}/01_download_reference.sh"

log_stage "02 — Build STAR genome index"
bash "${SCRIPT_DIR}/02_build_star_index.sh"

log_stage "03 — FastQC on raw reads"
bash "${SCRIPT_DIR}/03_fastqc_raw.sh"

log_stage "04 — fastp adapter + quality trimming"
bash "${SCRIPT_DIR}/04_fastp_trim.sh"

log_stage "05 — FastQC on trimmed reads"
bash "${SCRIPT_DIR}/05_fastqc_trimmed.sh"

log_stage "06 — STAR alignment to GRCh38"
bash "${SCRIPT_DIR}/06_star_align.sh"

log_stage "07 — samtools index + flagstat"
bash "${SCRIPT_DIR}/07_samtools_index.sh"

log_stage "08 — featureCounts gene quantification"
bash "${SCRIPT_DIR}/08_featurecounts.sh"

log_stage "09 — MultiQC aggregate report"
bash "${SCRIPT_DIR}/09_multiqc.sh"

log_stage "10 — Pre/post comparison plots"
python3 "${SCRIPT_DIR}/10_compare_qc.py" --project-dir "$PROJECT_DIR"

echo ""
echo "=================================================================="
echo "  PIPELINE COMPLETE  $(date '+%F %T')"
echo "=================================================================="
echo ""
echo "  Key outputs:"
echo "    results/09_multiqc/multiqc_report.html   <- main QC report"
echo "    results/10_comparison/summary_table.md   <- pre/post table"
echo "    results/10_comparison/*.png              <- comparison plots"
echo "    results/08_featurecounts/counts_matrix.txt  <- count matrix"
echo "    results/06_star_align/*/  <- BAMs + STAR logs"
