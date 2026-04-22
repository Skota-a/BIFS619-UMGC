#!/usr/bin/env bash
# =============================================================================
# 07_samtools_index.sh
#
# Runs a few standard post-alignment checks on STAR BAM files.
#
# 1. samtools index
#    Creates the .bai file needed for IGV and many downstream tools.
#
# 2. samtools flagstat
#    Gives a quick alignment summary (total, mapped, duplicates, etc.).
#    Useful alongside STAR Log.final.out.
#
# 3. samtools idxstats
#    Reports read counts by chromosome. Good for spotting unusual chrM signal
#    or unexpected mapping patterns.
#
# Outputs per sample in results/06_star_align/<sample>/:
#   *.Aligned.sortedByCoord.out.bam.bai
#   *.flagstat.txt
#   *.idxstats.txt
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

ALIGN_DIR="${PROJECT_DIR}/results/06_star_align"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$LOG_DIR"

# Find sorted BAM files from STAR
mapfile -t BAMS < <(find "$ALIGN_DIR" -name "*Aligned.sortedByCoord.out.bam" | sort)

if [[ ${#BAMS[@]} -eq 0 ]]; then
    echo "[ERROR] No BAM files found under $ALIGN_DIR"
    echo "        Run 06_star_align.sh first."
    exit 1
fi

echo "[$(date '+%F %T')] Processing ${#BAMS[@]} BAM file(s)"

for bam in "${BAMS[@]}"; do
    base="${bam%.bam}"
    sample=$(basename "$(dirname "$bam")")

    echo ""
    echo "[$(date '+%F %T')] $sample"

    # Build BAM index for random access
    echo "  Indexing..."
    samtools index -@ "$THREADS" "$bam" 2>> "${LOG_DIR}/07_samtools_${sample}.log"

    # Quick mapping summary
    echo "  Running flagstat..."
    samtools flagstat \
        -@ "$THREADS" \
        "$bam" \
        > "${base}.flagstat.txt" \
        2>> "${LOG_DIR}/07_samtools_${sample}.log"

    # Counts by chromosome
    echo "  Running idxstats..."
    samtools idxstats "$bam" > "${base}.idxstats.txt" 2>> "${LOG_DIR}/07_samtools_${sample}.log"

    # Print key numbers in terminal
    echo "  --- flagstat summary for $sample ---"
    grep -E "total|mapped|properly paired|singletons" "${base}.flagstat.txt" | \
        sed 's/^/    /'

    # Check mitochondrial fraction (often used as QC signal)
    chrM_reads=$(grep -E "^chrM" "${base}.idxstats.txt" | awk '{print $3}' || echo 0)
    total_mapped=$(awk 'NR==2 {print $1}' "${base}.flagstat.txt" 2>/dev/null || echo 1)
    if [[ "$total_mapped" -gt 0 && "$chrM_reads" -gt 0 ]]; then
        chrM_pct=$(awk "BEGIN {printf \"%.1f\", ($chrM_reads/$total_mapped)*100}")
        echo "    chrM reads: ${chrM_reads} (${chrM_pct}% of mapped)"
        if (( $(echo "$chrM_pct > 20" | awk '{print ($1 > 20)}') )); then
            echo "    [WARN] chrM fraction >20% — worth noting in QC"
        fi
    fi
done

echo ""
echo "[$(date '+%F %T')] All BAMs indexed. Ready for quantification (08_featurecounts.sh)"
