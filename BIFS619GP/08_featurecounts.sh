#!/usr/bin/env bash
# =============================================================================
# 08_featurecounts.sh
#
# Uses featureCounts (Subread) to assign aligned reads to genes.
# Output is the count matrix typically used for DESeq2 / edgeR.
#
# All samples here are single-end, so paired-end options are not needed.
#
# Key settings:
#   -t exon       count reads overlapping exons
#   -g gene_id    summarize counts at gene level
#   -Q 10         ignore low-confidence alignments
#   -s 0          unstranded library
#                 Change to 1 or 2 if your prep was stranded.
#
# Output in results/08_featurecounts/:
#   counts_matrix.txt
#   counts_matrix.txt.summary
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

ALIGN_DIR="${PROJECT_DIR}/results/06_star_align"
OUT_DIR="${PROJECT_DIR}/results/08_featurecounts"
GTF="${PROJECT_DIR}/refs/gencode.v46.primary_assembly.annotation.gtf"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# featureCounts needs an uncompressed GTF
if [[ ! -f "$GTF" && -f "${GTF}.gz" ]]; then
    echo "[$(date '+%F %T')] Decompressing GTF for featureCounts..."
    pigz -dkc "${GTF}.gz" > "$GTF"
fi

[[ -f "$GTF" ]] || { echo "[ERROR] GTF not found: $GTF"; exit 1; }

# Library strandedness
# Set to 1 or 2 if your kit was stranded
STRAND=0

# Collect BAM files for all samples
BAMS=()
for sample in SRR7898026 SRR7897501 SRR1039508; do
    bam=$(find "${ALIGN_DIR}/${sample}" -name "*Aligned.sortedByCoord.out.bam" 2>/dev/null | head -1)
    [[ -n "$bam" ]] || { echo "[ERROR] BAM missing for $sample — run 06_star_align.sh first"; exit 1; }
    BAMS+=("$bam")
done

echo "[$(date '+%F %T')] Running featureCounts on ${#BAMS[@]} SE samples (strandedness: $STRAND)"

featureCounts \
    -a "$GTF" \
    -o "${OUT_DIR}/counts_matrix.txt" \
    -t exon \
    -g gene_id \
    -Q 10 \
    -s "$STRAND" \
    -T "$THREADS" \
    --verbose \
    "${BAMS[@]}" \
    2>&1 | tee "${LOG_DIR}/08_featurecounts.log"

# Replace long BAM path names with sample IDs
echo "[$(date '+%F %T')] Renaming columns to sample IDs..."
python3 - << 'PYEOF'
import os, pandas as pd

out_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "results", "08_featurecounts"
)
path = os.path.join(out_dir, "counts_matrix.txt")
if not os.path.exists(path):
    print("counts_matrix.txt not found — skipping column rename")
    exit(0)

df = pd.read_csv(path, sep="\t", comment="#", index_col=0)
count_cols = [c for c in df.columns if c.endswith(".bam")]
rename = {}
for c in count_cols:
    for acc in ["SRR7898026", "SRR7897501", "SRR1039508"]:
        if acc in c:
            rename[c] = acc
df = df[count_cols].rename(columns=rename)
df.index.name = "gene_id"
df.to_csv(path, sep="\t")
print(f"  Matrix: {df.shape[0]} genes x {df.shape[1]} samples")
print(f"  Columns: {list(df.columns)}")
print(f"  Saved: {path}")
PYEOF

echo ""
echo "[$(date '+%F %T')] featureCounts complete."
echo "  Count matrix : ${OUT_DIR}/counts_matrix.txt"
echo "  Summary      : ${OUT_DIR}/counts_matrix.txt.summary"
echo ""
echo "  Next: bash scripts/09_multiqc.sh"
echo "  Or load counts_matrix.txt into R/DESeq2 for differential expression."
