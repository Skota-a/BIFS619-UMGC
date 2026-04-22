#!/usr/bin/env bash
# =============================================================================
# 04_fastp_trim.sh
#
# Runs fastp on all three single-end samples.
#
# Goal:
#   Remove adapter sequence, trim low-quality ends, and discard poor reads
#   before alignment. This usually improves mapping and downstream counts.
#
# Main settings:
#   Q20 cutoff              = bases below this are considered low quality
#   30% low-quality bases   = read is discarded
#   minimum length 36 bp    = very short reads are removed
#   cut_right               = trims weak quality from 3' end
#   trim_poly_g             = removes NovaSeq-style polyG tails
#
# Outputs in results/04_fastp/:
#   *.trimmed.fastq.gz   cleaned reads for STAR
#   *.fastp.html         QC report
#   *.fastp.json         summary stats
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

RAW_DIR="${PROJECT_DIR}/raw"
OUT_DIR="${PROJECT_DIR}/results/04_fastp"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Single-end trimming
# fasterq-dump often names SE files as _1.fastq, so check that first.
# ---------------------------------------------------------------------------
trim_se() {
    local sample="$1"

    local in1="${RAW_DIR}/${sample}_1.fastq"
    [[ -f "$in1" ]] || in1="${RAW_DIR}/${sample}_1.fastq.gz"
    [[ -f "$in1" ]] || in1="${RAW_DIR}/${sample}.fastq"
    [[ -f "$in1" ]] || in1="${RAW_DIR}/${sample}.fastq.gz"
    [[ -f "$in1" ]] || in1="${RAW_DIR}/${sample}.sra.fastq"
    [[ -f "$in1" ]] || in1="${RAW_DIR}/${sample}.sra.fastq.gz"
    [[ -f "$in1" ]] || {
        echo "[ERROR] Cannot find raw FASTQ for ${sample}. Checked:"
        echo "  ${RAW_DIR}/${sample}_1.fastq"
        echo "  ${RAW_DIR}/${sample}_1.fastq.gz"
        echo "  ${RAW_DIR}/${sample}.fastq"
        echo "  ${RAW_DIR}/${sample}.fastq.gz"
        echo "  ${RAW_DIR}/${sample}.sra.fastq"
        echo "  ${RAW_DIR}/${sample}.sra.fastq.gz"
        exit 1
    }

    echo "[$(date '+%F %T')] Trimming: $sample  (input: $(basename "$in1"))"

    fastp \
        --in1   "$in1" \
        --out1  "${OUT_DIR}/${sample}.trimmed.fastq.gz" \
        --json  "${OUT_DIR}/${sample}.fastp.json" \
        --html  "${OUT_DIR}/${sample}.fastp.html" \
        --report_title              "fastp — ${sample}" \
        --thread                    "$THREADS" \
        --qualified_quality_phred   20 \
        --unqualified_percent_limit 30 \
        --n_base_limit              5 \
        --length_required           36 \
        --cut_right \
        --cut_right_window_size     4 \
        --cut_right_mean_quality    20 \
        --trim_poly_g \
        --overrepresentation_analysis \
        2> "${LOG_DIR}/04_fastp_${sample}.log"

    echo "[$(date '+%F %T')] Done: $sample"
}

# ---------------------------------------------------------------------------
# Process all three samples
# ---------------------------------------------------------------------------
trim_se "SRR7898026"
trim_se "SRR7897501"
trim_se "SRR1039508"

echo ""
echo "[$(date '+%F %T')] All trimming complete."
echo "  Cleaned reads in $OUT_DIR"
echo "  Next step: bash scripts/05_fastqc_trimmed.sh"
