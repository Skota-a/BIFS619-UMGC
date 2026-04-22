#!/usr/bin/env bash
# =============================================================================
# quantification.sh
#
# Quantifies expression from existing STAR BAM files using featureCounts,
# then generates:
#   1. raw counts matrix
#   2. CPM matrix
#   3. TPM matrix
#   4. top 20 highly expressed genes (by mean TPM)
#
# This script is built for the current project state:
#   - SRR7898026 (single-end, completed)
#   - SRR7897501 (single-end, completed)
#
# SRR1039508 is excluded here because it produced paired-end FASTQ files and
# was not carried through the current single-end alignment workflow.
#
# Outputs:
#   results/quantification/
#     counts_matrix.txt
#     counts_matrix.txt.summary
#     counts_only.tsv
#     cpm_matrix.tsv
#     tpm_matrix.tsv
#     top20_highly_expressed_meanTPM.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

ALIGN_DIR="${PROJECT_DIR}/results/06_star_align"
REFS_DIR="${PROJECT_DIR}/refs"
OUT_DIR="${PROJECT_DIR}/results/quantification"
LOG_DIR="${PROJECT_DIR}/logs"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR" "$LOG_DIR"

# Input files
GTF_GZ="${REFS_DIR}/gencode.v46.primary_assembly.annotation.gtf.gz"
GTF="${REFS_DIR}/gencode.v46.primary_assembly.annotation.gtf"

SAMPLES=(
  "SRR7898026"
  "SRR7897501"
)

# Check dependencies
command -v featureCounts >/dev/null 2>&1 || {
    echo "[ERROR] featureCounts not found. Install Subread / featureCounts first."
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo "[ERROR] python3 not found."
    exit 1
}

# Prepare uncompressed GTF if needed
if [[ ! -f "$GTF" ]]; then
    if [[ -f "$GTF_GZ" ]]; then
        echo "[$(date '+%F %T')] Decompressing GTF..."
        gzip -dc "$GTF_GZ" > "$GTF"
    else
        echo "[ERROR] GTF not found:"
        echo "  $GTF"
        echo "  $GTF_GZ"
        exit 1
    fi
fi

# Collect BAMs
BAMS=()
for sample in "${SAMPLES[@]}"; do
    bam=$(find "${ALIGN_DIR}/${sample}" -maxdepth 1 -name "*Aligned.sortedByCoord.out.bam" | head -n 1 || true)
    if [[ -z "$bam" ]]; then
        echo "[ERROR] BAM missing for $sample"
        echo "        Expected under: ${ALIGN_DIR}/${sample}"
        exit 1
    fi
    BAMS+=("$bam")
done

echo "[$(date '+%F %T')] Running featureCounts on ${#BAMS[@]} sample(s)"
printf '  %s\n' "${BAMS[@]}"

# Run featureCounts
featureCounts \
    -a "$GTF" \
    -o "${OUT_DIR}/counts_matrix.txt" \
    -t exon \
    -g gene_id \
    -Q 10 \
    -s 0 \
    -T "$THREADS" \
    "${BAMS[@]}" \
    2>&1 | tee "${LOG_DIR}/quantification.log"

# Post-process: rename columns, generate counts/CPM/TPM, top 20 genes
export OUT_DIR
export SAMPLES_CSV="$(IFS=,; echo "${SAMPLES[*]}")"

python3 << 'PYEOF'
import os
import sys
import pandas as pd
import numpy as np

out_dir = os.environ["OUT_DIR"]
samples = os.environ["SAMPLES_CSV"].split(",")

counts_path = os.path.join(out_dir, "counts_matrix.txt")
if not os.path.exists(counts_path):
    print("[ERROR] counts_matrix.txt not found", file=sys.stderr)
    sys.exit(1)

# Read featureCounts output
df = pd.read_csv(counts_path, sep="\t", comment="#")

required_cols = ["Geneid", "Length"]
for col in required_cols:
    if col not in df.columns:
        print(f"[ERROR] Missing expected column: {col}", file=sys.stderr)
        sys.exit(1)

# Identify BAM count columns
bam_cols = [c for c in df.columns if c.endswith(".bam")]
if len(bam_cols) == 0:
    print("[ERROR] No BAM count columns found in featureCounts output", file=sys.stderr)
    sys.exit(1)

# Rename BAM columns to clean sample IDs
rename = {}
for c in bam_cols:
    for sample in samples:
        if sample in c:
            rename[c] = sample

df = df.rename(columns=rename)

# Keep metadata + sample count columns
sample_cols = [s for s in samples if s in df.columns]
if len(sample_cols) == 0:
    print("[ERROR] No sample columns found after renaming", file=sys.stderr)
    sys.exit(1)

meta = df[["Geneid", "Length"]].copy()
counts = df[sample_cols].copy()

# Force numeric
for c in sample_cols:
    counts[c] = pd.to_numeric(counts[c], errors="coerce").fillna(0)

# Save counts-only table
counts_only = pd.concat([meta, counts], axis=1)
counts_only.to_csv(os.path.join(out_dir, "counts_only.tsv"), sep="\t", index=False)

# CPM calculation
lib_sizes = counts.sum(axis=0)
cpm = counts.divide(lib_sizes, axis=1) * 1_000_000
cpm_out = pd.concat([meta[["Geneid"]], cpm], axis=1)
cpm_out.to_csv(os.path.join(out_dir, "cpm_matrix.tsv"), sep="\t", index=False)

# TPM calculation
# RPK = counts / gene_length_kb
gene_length_kb = meta["Length"] / 1000.0
rpk = counts.divide(gene_length_kb, axis=0)
scaling = rpk.sum(axis=0) / 1_000_000
tpm = rpk.divide(scaling, axis=1)
tpm_out = pd.concat([meta[["Geneid", "Length"]], tpm], axis=1)
tpm_out.to_csv(os.path.join(out_dir, "tpm_matrix.tsv"), sep="\t", index=False)

# Top 20 genes by mean TPM across samples
tpm_with_mean = tpm.copy()
tpm_with_mean["mean_TPM"] = tpm.mean(axis=1)
top20 = pd.concat([meta, tpm_with_mean], axis=1).sort_values("mean_TPM", ascending=False).head(20)
top20.to_csv(os.path.join(out_dir, "top20_highly_expressed_meanTPM.tsv"), sep="\t", index=False)

print(f"Saved counts_only.tsv with shape: {counts_only.shape}")
print(f"Saved cpm_matrix.tsv with shape: {cpm_out.shape}")
print(f"Saved tpm_matrix.tsv with shape: {tpm_out.shape}")
print("Top 20 highly expressed genes written to top20_highly_expressed_meanTPM.tsv")
PYEOF

echo ""
echo "[$(date '+%F %T')] Quantification complete."
echo "  Raw counts : ${OUT_DIR}/counts_matrix.txt"
echo "  Counts only: ${OUT_DIR}/counts_only.tsv"
echo "  CPM matrix : ${OUT_DIR}/cpm_matrix.tsv"
echo "  TPM matrix : ${OUT_DIR}/tpm_matrix.tsv"
echo "  Top genes  : ${OUT_DIR}/top20_highly_expressed_meanTPM.tsv"
