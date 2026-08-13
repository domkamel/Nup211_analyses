#!/bin/bash
#
# STAR alignment, 48 h Nup211 depletion dataset
# Corresponds to GEO series GSE216106 (Kamel et al. 2024)
#
# Used STAR version=STAR_2.5.4b
#
# Reconstructed from contemporaneous notes (RNA-seq_workflow.txt).
# Reference: ASM294v2, 75 bp paired-end reads (sjdbOverhang = 75 - 1).
#
# Alignment only. Read counting was done in R with summarizeOverlaps
# (GenomicAlignments) against features retrieved from Ensembl Fungi via
# makeTxDbFromBiomart. See the accompanying R script.
#
# Requires: STAR

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REF_DIR=~/pombe-downloads
GENOME_DIR=~/genome
FASTQ_DIR=~/fastq
ALIGN_DIR=~/alignment
THREADS=8

GENOME_FA="${REF_DIR}/GCF_000002945.1_ASM294v2_genomic.fna"
GTF="${REF_DIR}/GCF_000002945.1_ASM294v2_genomic.gtf"

N_SAMPLES=6

mkdir -p "${ALIGN_DIR}"

# ---------------------------------------------------------------------------
# Step 1: Build the STAR genome index (run once)
# ---------------------------------------------------------------------------
#
#   --sjdbOverhang 74         read length (75) minus 1
#   --genomeSAindexNbases 11  scaled down for the small S. pombe genome

if [ ! -d "${GENOME_DIR}" ]; then

    mkdir -p "${GENOME_DIR}"

    STAR --runMode genomeGenerate \
         --genomeDir "${GENOME_DIR}" \
         --genomeFastaFiles "${GENOME_FA}" \
         --sjdbGTFfile "${GTF}" \
         --sjdbOverhang 74 \
         --genomeSAindexNbases 11 \
         --runThreadN "${THREADS}"

fi

# ---------------------------------------------------------------------------
# Step 2: Align each sample
# ---------------------------------------------------------------------------
#
# Each sample was sequenced across two lanes (L001, L002). STAR takes a
# comma-separated list of R1 files, then a space, then the matching list of
# R2 files, so lanes are merged during alignment.
#
# Reads are uncompressed .fastq; add --readFilesCommand zcat for .fastq.gz.

for sample_num in $(seq 1 "${N_SAMPLES}"); do

    sample="${sample_num}_S${sample_num}"

    echo "Aligning ${sample} ..."

    STAR --genomeDir "${GENOME_DIR}" \
         --readFilesIn "${FASTQ_DIR}/${sample}_L001_R1_001.fastq,${FASTQ_DIR}/${sample}_L002_R1_001.fastq" \
                       "${FASTQ_DIR}/${sample}_L001_R2_001.fastq,${FASTQ_DIR}/${sample}_L002_R2_001.fastq" \
         --runThreadN "${THREADS}" \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM 30000000000 \
         --outFileNamePrefix "${ALIGN_DIR}/${sample}"

    echo "Finished ${sample}"

done

# ---------------------------------------------------------------------------
# Step 3: Counting and differential expression
# ---------------------------------------------------------------------------
#
# Sorted BAM files were read into R with Rsamtools::BamFileList and counted
# with GenomicAlignments::summarizeOverlaps, then analyzed with DESeq2.
