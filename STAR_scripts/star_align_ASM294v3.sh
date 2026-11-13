#!/bin/bash
#
# STAR alignment + gene-level counts, S. pombe ASM294v3 (current pipeline)
#
# Environment:  conda activate rnaseq   (see envs/rnaseq.yml)
#   STAR 2.7.11b, samtools 1.24, htseq 2.1.2
#
# Reference: NCBI RefSeq ASM294v3 / GCF_000002945.2 -- the same assembly used for
# the salmon index, so STAR and salmon results are directly comparable.
#
# Outputs per sample in ${OUT_DIR}/Sample_<n>/:
#   Aligned.sortedByCoord.out.bam(.bai)   sorted, indexed alignments
#   ReadsPerGene.out.tab                  STAR gene counts (col 2 unstranded,
#                                         col 3 forward, col 4 reverse)
#   htseq_counts.tsv                      htseq-count gene counts
#   Log.final.out                         mapping summary
#
# Written for macOS /bin/bash 3.2 (no bash 4 features).

set -euo pipefail

# ---------------------------------------------------------------------------
# Settings -- edit these for each dataset
# ---------------------------------------------------------------------------

REQUIRED_STAR="2.7.11b"

REF_DIR=~/pombe-downloads
GENOME_FA="${REF_DIR}/GCF_000002945.2_ASM294v3_genomic.fna"
GTF="${REF_DIR}/GCF_000002945.2_ASM294v3_genomic.gtf"
INDEX_DIR=~/reference/star_index_ASM294v3_STAR${REQUIRED_STAR}

FASTQ_DIR=~/Documents/Zhong-HZ-17633_2024_12_12
OUT_DIR=~/star_output
THREADS=4
BAM_SORT_RAM=8000000000          # 8 GB for BAM sorting; lower on a 16 GB Mac if needed

SAMPLES="1 2 3 10 11 12"
LANES="L004 L005"
# FASTQ names are built as ${FASTQ_DIR}/<sample>_S<sample>_<lane>_R1_001.fastq.gz

# Library strandedness for htseq-count (-s): yes | no | reverse
# Salmon inferred ISR for the 16 h libraries = reverse-stranded (TruSeq-type).
# Use the matching ReadsPerGene column: reverse -> col 4, yes -> col 3, no -> col 2.
STRANDEDNESS=reverse

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

FOUND_STAR=$(STAR --version 2>/dev/null || true)
if [ "${FOUND_STAR}" != "${REQUIRED_STAR}" ]; then
    echo "ERROR: needs STAR ${REQUIRED_STAR}, found '${FOUND_STAR:-none}'. Run: conda activate rnaseq" >&2
    exit 1
fi
for tool in samtools htseq-count; do
    command -v "${tool}" >/dev/null || { echo "ERROR: ${tool} not found. Run: conda activate rnaseq" >&2; exit 1; }
done
for f in "${GENOME_FA}" "${GTF}"; do
    [ -s "${f}" ] || { echo "ERROR: missing reference file ${f}" >&2; exit 1; }
done

fastq() {   # fastq <sample> <lane> <R1|R2>
    echo "${FASTQ_DIR}/${1}_S${1}_${2}_${3}_001.fastq.gz"
}

first_sample=${SAMPLES%% *}
first_lane=${LANES%% *}
first_fq=$(fastq "${first_sample}" "${first_lane}" R1)
[ -s "${first_fq}" ] || { echo "ERROR: FASTQ not found: ${first_fq}  (check FASTQ_DIR)" >&2; exit 1; }

# Read length from the first read; sjdbOverhang = read length - 1
READ_LENGTH=$( (gunzip -c "${first_fq}" || true) | awk 'NR==2 {print length($0); exit}')
SJDB_OVERHANG=$((READ_LENGTH - 1))
echo "Read length ${READ_LENGTH} bp -> sjdbOverhang ${SJDB_OVERHANG}"

# ---------------------------------------------------------------------------
# Step 1: Genome index (built once; reused if present)
# ---------------------------------------------------------------------------
#
#   --genomeSAindexNbases 10   min(14, log2(12.6 Mb)/2 - 1) for the small S. pombe genome
#
# The index is tied to the STAR version and to sjdbOverhang. If either changes,
# delete INDEX_DIR and rerun.

if [ ! -s "${INDEX_DIR}/SA" ]; then
    echo "Building STAR index in ${INDEX_DIR} ..."
    mkdir -p "${INDEX_DIR}"
    STAR --runMode genomeGenerate \
         --genomeDir "${INDEX_DIR}" \
         --genomeFastaFiles "${GENOME_FA}" \
         --sjdbGTFfile "${GTF}" \
         --sjdbOverhang "${SJDB_OVERHANG}" \
         --genomeSAindexNbases 10 \
         --runThreadN "${THREADS}" \
         --outFileNamePrefix "${INDEX_DIR}/"
else
    built_overhang=$(awk '$1=="sjdbOverhang" {print $2}' "${INDEX_DIR}/genomeParameters.txt" 2>/dev/null || true)
    if [ -n "${built_overhang}" ] && [ "${built_overhang}" != "${SJDB_OVERHANG}" ]; then
        echo "WARNING: index built with sjdbOverhang ${built_overhang}, reads need ${SJDB_OVERHANG}." >&2
        echo "         Alignment will still run; rebuild the index for best splice-junction detection." >&2
    fi
fi

# ---------------------------------------------------------------------------
# Step 2: Align, index, count
# ---------------------------------------------------------------------------
#
# Lanes are merged during alignment: STAR takes comma-separated R1 files,
# then a space, then the matching comma-separated R2 files.
# --readFilesCommand "gunzip -c": macOS zcat does not read .gz files.

mkdir -p "${OUT_DIR}"

for sample in ${SAMPLES}; do

    r1=""; r2=""
    for lane in ${LANES}; do
        r1="${r1:+${r1},}$(fastq "${sample}" "${lane}" R1)"
        r2="${r2:+${r2},}$(fastq "${sample}" "${lane}" R2)"
    done

    out="${OUT_DIR}/Sample_${sample}"
    mkdir -p "${out}"
    echo "Aligning sample ${sample} ..."

    STAR --genomeDir "${INDEX_DIR}" \
         --readFilesIn "${r1}" "${r2}" \
         --readFilesCommand gunzip -c \
         --runThreadN "${THREADS}" \
         --outSAMtype BAM SortedByCoordinate \
         --limitBAMsortRAM "${BAM_SORT_RAM}" \
         --quantMode GeneCounts \
         --outFileNamePrefix "${out}/"

    samtools index "${out}/Aligned.sortedByCoord.out.bam"

    echo "Counting sample ${sample} with htseq-count (-s ${STRANDEDNESS}) ..."
    htseq-count \
        -f bam \
        -r pos \
        -s "${STRANDEDNESS}" \
        -t exon \
        -i gene_id \
        --additional-attr gene \
        "${out}/Aligned.sortedByCoord.out.bam" \
        "${GTF}" \
        > "${out}/htseq_counts.tsv"

    echo "Sample ${sample} complete."
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

echo ""
echo "Uniquely mapped reads:"
for d in "${OUT_DIR}"/Sample_*; do
    pct=$(awk -F'|' '/Uniquely mapped reads %/ {gsub(/[ \t]/, "", $2); print $2}' "${d}/Log.final.out")
    printf "  %-12s %s\n" "$(basename "${d}")" "${pct}"
done
