#!/usr/bin/env bash
# =============================================================================
# 06_star_align.sh
#
# Aligns trimmed single-end reads to GRCh38 using STAR.
#
# The script automatically uses whichever index was built in step 02:
#   - chr1-only index (low RAM test mode)
#   - full genome index
#
# If chr1 mode is used, only chromosome 1 reads can align. That is expected
# and mainly useful for checking that the workflow runs correctly.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

TRIM_DIR="${PROJECT_DIR}/results/04_fastp"
OUT_DIR="${PROJECT_DIR}/results/06_star_align"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

# Auto-detect which STAR index is available
if [[ -d "${PROJECT_DIR}/refs/star_index_chr1" ]]; then
    STAR_INDEX="${PROJECT_DIR}/refs/star_index_chr1"
    echo "[$(date '+%F %T')] Using chr1-only index: $STAR_INDEX"
    echo "  Only chr1 alignments will appear."
elif [[ -d "${PROJECT_DIR}/refs/star_index" ]]; then
    STAR_INDEX="${PROJECT_DIR}/refs/star_index"
    echo "[$(date '+%F %T')] Using full genome index: $STAR_INDEX"
else
    echo "[ERROR] No STAR index found. Run 02_build_star_index.sh first."
    exit 1
fi

mkdir -p "$OUT_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Align one single-end sample
# ---------------------------------------------------------------------------
align_se() {
    local sample="$1"
    local fq="${TRIM_DIR}/${sample}.trimmed.fastq.gz"
    local sample_out="${OUT_DIR}/${sample}"

    [[ -f "$fq" ]] || { echo "[ERROR] Missing trimmed FASTQ: $fq"; exit 1; }

    mkdir -p "$sample_out"
    echo "[$(date '+%F %T')] Aligning: $sample"

    STAR \
        --runThreadN                     "$THREADS" \
        --genomeDir                      "$STAR_INDEX" \
        --readFilesIn                    "$fq" \
        --readFilesCommand               zcat \
        --outFileNamePrefix              "${sample_out}/${sample}_" \
        --outSAMtype                     BAM SortedByCoordinate \
        --outSAMattributes               NH HI AS NM MD \
        --outSAMstrandField              intronMotif \
        --outFilterIntronMotifs          RemoveNoncanonical \
        --quantMode                      GeneCounts \
        --outFilterMultimapNmax          20 \
        --alignSJoverhangMin             8 \
        --alignSJDBoverhangMin           1 \
        --outFilterMismatchNmax          999 \
        --outFilterMismatchNoverReadLmax  0.04 \
        --alignIntronMin                 20 \
        --alignIntronMax                 1000000 \
        --alignMatesGapMax               1000000 \
        --runMode                        alignReads \
        2>&1 | tee "${LOG_DIR}/06_star_${sample}.log"

    echo "[$(date '+%F %T')] Done: $sample -> $sample_out"
}

# ---------------------------------------------------------------------------
# Run all three samples
# ---------------------------------------------------------------------------
align_se "SRR7898026"
align_se "SRR7897501"
align_se "SRR1039508"

echo ""
echo "[$(date '+%F %T')] All alignments finished."
echo ""
echo "  Alignment summaries:"
for sample in SRR7898026 SRR7897501 SRR1039508; do
    log="${OUT_DIR}/${sample}/${sample}_Log.final.out"
    if [[ -f "$log" ]]; then
        echo "  ---- $sample ----"
        grep -E "Uniquely mapped|mapped to multiple|% of reads unmapped" "$log" | \
            sed 's/^[[:space:]]*/    /'
    fi
done

echo ""
echo "  Next step: bash scripts/07_samtools_index.sh"
