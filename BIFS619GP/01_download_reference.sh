#!/usr/bin/env bash
# =============================================================================
# 01_download_reference.sh
#
# Downloads human reference genome (GRCh38) and matching gene annotation (GTF)
# from GENCODE. These are required for alignment and gene counting downstream.
#
# Genome = sequence
# GTF = gene/exon coordinates
#
# Keep versions consistent across the project — mixing releases will break results.
#
# GENCODE is used here because it's standard in most RNA-seq workflows
# (ENCODE, GTEx, many publications). Ensembl is equivalent if you prefer it.
#
# Primary assembly includes standard chromosomes + basic scaffolds.
# It excludes alternate haplotypes to avoid mapping complications.
#
# Outputs:
#   GRCh38.primary_assembly.genome.fa.gz
#   gencode.v46.primary_assembly.annotation.gtf.gz
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

REFS_DIR="${PROJECT_DIR}/refs"
LOG_DIR="${PROJECT_DIR}/logs"

mkdir -p "$REFS_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Reference version (change here if needed)
# Don’t mix GENCODE and Ensembl versions in the same project.
# ---------------------------------------------------------------------------
GENCODE_RELEASE="46"
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

GENOME_FILE="GRCh38.primary_assembly.genome.fa.gz"
GTF_FILE="gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"

GENOME_URL="${GENCODE_BASE}/${GENOME_FILE}"
GTF_URL="${GENCODE_BASE}/${GTF_FILE}"


download_if_missing() {
    local url="$1"
    local dest="${REFS_DIR}/$(basename "$url")"

    # Skip if already downloaded
    if [[ -f "$dest" ]]; then
        echo "[$(date '+%F %T')] Exists: $(basename "$dest")"
    else
        echo "[$(date '+%F %T')] Downloading: $(basename "$dest")"

        wget \
            --continue \
            --show-progress \
            --output-document "$dest" \
            "$url" \
            2>> "${LOG_DIR}/01_download_reference.log"

        echo "[$(date '+%F %T')] Saved"
    fi
}

echo "[$(date '+%F %T')] Downloading GRCh38 + GENCODE v${GENCODE_RELEASE}"

download_if_missing "$GENOME_URL"
download_if_missing "$GTF_URL"

# Basic check to avoid silent failed downloads
for f in "${REFS_DIR}/${GENOME_FILE}" "${REFS_DIR}/${GTF_FILE}"; do
    size=$(stat -c%s "$f" 2>/dev/null || echo 0)

    if [[ "$size" -lt 1000000 ]]; then
        echo "[ERROR] File too small: $f"
        exit 1
    fi

    echo "[$(date '+%F %T')] OK: $(du -sh "$f" | cut -f1)"
done

echo ""
echo "Reference ready in $REFS_DIR"
echo "Next: 02_build_star_index.sh"
