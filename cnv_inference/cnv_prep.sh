#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Generate BAFExtract SNP files from RNA-seq BAM files for CaSpER
#
# Requirements:
#   - samtools
#   - BAFExtract
#
# Download the hg38 reference files from:
# https://github.com/akdess/BAFExtract
#
# Required files:
#   baf_reference/
#     ├── hg38.list
#     └── hg38/
#         ├── 1.bin
#         ├── ...
#         └── chr_ids.txt
#
# Input:
#   bam/*.markdup.sorted.bam
#
# Output:
#   baf_pileups/<sample>/
#   baf_output/<sample>.snp
#
# BAFExtract parameters
#   50  = minimum mapping quality (MAPQ)
#   20  = minimum coverage
#   4   = minimum alternative allele count
#   0.1 = minimum minor allele frequency
###############################################################################

BAFEXTRACT="/Users/annaschwager/Documents/installed_tools/BAFExtract/bin/BAFExtract"


PROJECT_DIR="/Users/annaschwager/Documents/projects/MCL/revision_analysis/rnaseq/input"

BAM_DIR="${PROJECT_DIR}/bam"
REF_DIR="${PROJECT_DIR}/baf_reference/hg38"
REF_LIST="${PROJECT_DIR}/baf_reference/hg38.list"
PILEUP_DIR="${PROJECT_DIR}/baf_pileups"
OUTPUT_DIR="${PROJECT_DIR}/baf_output"

mkdir -p "$PILEUP_DIR"
mkdir -p "$OUTPUT_DIR"

for bam in "$BAM_DIR"/*.markdup.sorted.bam
do
    sample=$(basename "$bam" .markdup.sorted.bam)

    echo "========================================="
    echo "Processing $sample"
    echo "========================================="

    # Skip completed samples
    if [[ -s "$OUTPUT_DIR/${sample}.snp" ]]; then
        echo "Output already exists. Skipping."
        continue
    fi

    mkdir -p "$PILEUP_DIR/$sample"

    echo "Generating pileup..."

    samtools view "$bam" | \
    "$BAFEXTRACT" \
        -generate_compressed_pileup_per_SAM \
        stdin \
        "$REF_LIST" \
        "$PILEUP_DIR/$sample" \
        50 \
        0

    echo "Generating SNP file..."

    "$BAFEXTRACT" \
        -get_SNVs_per_pileup \
        "$REF_LIST" \
        "$PILEUP_DIR/$sample" \
        "$REF_DIR" \
        20 \
        4 \
        0.1 \
        "$OUTPUT_DIR/${sample}.snp"

    echo "Finished $sample"
done

echo
echo "All samples processed."