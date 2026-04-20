#!/usr/bin/env bash
# =============================================================================
# 02_build_star_index.sh
#
# Builds a STAR genome index for alignment.
#
# Default mode uses chromosome 1 only so the pipeline can run on laptops or
# small servers (~8 GB RAM). Full GRCh38 indexing needs ~30–40 GB RAM and will
# fail on most local machines. Changed to low RAM version as UMGC virtual OS 
#couldn't handle full genome version and not enough space to download.
#
# This chr1-only setup is just for testing the pipeline end-to-end. Chr1 is
# large enough (~8% of genome) to validate that everything works correctly.
#
# When running on HPC, set:
#   USE_CHR1_ONLY=false
# to build the full genome index.
#
# Output:
#   refs/star_index_chr1/   (default, low RAM)
#   refs/star_index/        (full genome)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

REFS_DIR="${PROJECT_DIR}/refs"
LOG_DIR="${PROJECT_DIR}/logs"

THREADS="${THREADS:-8}"
READ_LENGTH="${READ_LENGTH:-100}"
OVERHANG=$((READ_LENGTH - 1))

# Use chr1 only unless explicitly disabled
USE_CHR1_ONLY="${USE_CHR1_ONLY:-true}"

mkdir -p "$LOG_DIR"

# ---------------------------------------------------------------------------
# Input reference files (from step 01)
# ---------------------------------------------------------------------------
FULL_FA="${REFS_DIR}/GRCh38.primary_assembly.genome.fa.gz"
FULL_GTF="${REFS_DIR}/gencode.v46.primary_assembly.annotation.gtf.gz"

# Basic check so we fail early if step 01 wasn't run
[[ -f "$FULL_FA"  ]] || { echo "[ERROR] Missing FASTA: run 01_download_reference.sh"; exit 1; }
[[ -f "$FULL_GTF" ]] || { echo "[ERROR] Missing GTF: run 01_download_reference.sh"; exit 1; }

# ---------------------------------------------------------------------------
# chr1-only mode (default, low RAM)
# ---------------------------------------------------------------------------
if [[ "$USE_CHR1_ONLY" == "true" ]]; then

    CHR1_FA="${REFS_DIR}/chr1.fa"
    CHR1_GTF="${REFS_DIR}/chr1.gtf"
    STAR_INDEX="${REFS_DIR}/star_index_chr1"

    echo "[$(date '+%F %T')] Mode: chr1-only (low RAM test)"

    # Pull chr1 sequence from full genome if not already done
    if [[ ! -f "$CHR1_FA" ]]; then
        echo "[$(date '+%F %T')] Extracting chr1 FASTA..."
        pigz -dkc "$FULL_FA" | \
            awk '/^>chr1$/{p=1} /^>chr[^ ]/{if(!/^>chr1$/)p=0} p' \
            > "$CHR1_FA"
    fi

    # Pull chr1 annotations from GTF
    if [[ ! -f "$CHR1_GTF" ]]; then
        echo "[$(date '+%F %T')] Extracting chr1 GTF..."
        pigz -dkc "$FULL_GTF" | \
            awk '$1 == "chr1" || /^#/' \
            > "$CHR1_GTF"
    fi

    mkdir -p "$STAR_INDEX"

    echo "[$(date '+%F %T')] Building STAR index (chr1)"

    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$CHR1_FA" \
        --sjdbGTFfile "$CHR1_GTF" \
        --sjdbOverhang "$OVERHANG" \
        --runThreadN "$THREADS" \
        --limitGenomeGenerateRAM 7000000000 \
        --genomeSAindexNbases 11 \
        2>&1 | tee "${LOG_DIR}/02_star_index_chr1.log"

# ---------------------------------------------------------------------------
# full genome mode (HPC / high RAM)
# ---------------------------------------------------------------------------
else

    GENOME_RAM_GB="${GENOME_RAM_GB:-44}"
    STAR_INDEX="${REFS_DIR}/star_index"

    GENOME_FA_UNZIPPED="${REFS_DIR}/GRCh38.primary_assembly.genome.fa"
    GTF_UNZIPPED="${REFS_DIR}/gencode.v46.primary_assembly.annotation.gtf"

    # Decompress only if needed (STAR works faster on uncompressed files)
    if [[ ! -f "$GENOME_FA_UNZIPPED" ]]; then
        echo "[$(date '+%F %T')] Decompressing genome..."
        pigz -dkc "$FULL_FA" > "$GENOME_FA_UNZIPPED"
    fi

    if [[ ! -f "$GTF_UNZIPPED" ]]; then
        echo "[$(date '+%F %T')] Decompressing GTF..."
        pigz -dkc "$FULL_GTF" > "$GTF_UNZIPPED"
    fi

    mkdir -p "$STAR_INDEX"

    echo "[$(date '+%F %T')] Building full STAR index"
    echo "  RAM limit: ${GENOME_RAM_GB} GB"

    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME_FA_UNZIPPED" \
        --sjdbGTFfile "$GTF_UNZIPPED" \
        --sjdbOverhang "$OVERHANG" \
        --runThreadN "$THREADS" \
        --limitGenomeGenerateRAM $((GENOME_RAM_GB * 1024 * 1024 * 1024)) \
        --genomeSAindexNbases 14 \
        2>&1 | tee "${LOG_DIR}/02_star_index.log"

fi

echo ""
echo "[$(date '+%F %T')] STAR index complete: $STAR_INDEX"
echo "Index size: $(du -sh "$STAR_INDEX" | cut -f1)"
echo "Next: 03_fastqc_raw.sh → 04_fastp_trim.sh → 06_star_align.sh"