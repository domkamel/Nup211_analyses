#!/bin/bash
#
# Salmon quantification of S. pombe RNA-seq
# Batch:  Zhong-HZ-17633_2024_12_12
# Run on: 21 December 2024, Salmon v1.10.3
#
# Requires the index built by build_salmon_index.sh.
#
# Reads from lanes L004 and L005 are passed together to a single call per sample. 
# Lanes are merged during quantification rather than by concatenating FASTQ files beforehand.
#
#   -l A                 infer library type (inferred: ISR, stranded PE)
#   --validateMappings   selective alignment
#
# Used: salmon 1.10.3
#   conda env create -f envs/salmon1.10.yml     (first time only)
#   conda activate salmon1.10

set -euo pipefail

# ---------------------------------------------------------------------------
# Version check: the index and all 16 h outputs were made with salmon 1.10.3.
# An index built by one salmon version cannot be read by another.
# ---------------------------------------------------------------------------

REQUIRED_SALMON="1.10.3"
FOUND_SALMON=$(salmon --version 2>/dev/null | awk '{print $2}' || true)
if [ "${FOUND_SALMON}" != "${REQUIRED_SALMON}" ]; then
    echo "ERROR: this script needs salmon ${REQUIRED_SALMON}, found '${FOUND_SALMON:-none}'." >&2
    echo "       Run: conda activate salmon1.10   (or: conda env create -f envs/salmon1.10.yml)" >&2
    exit 1
fi

INDEX_DIR=~/reference/salmon_index
FASTQ_DIR=~/Documents/Zhong-HZ-17633_2024_12_12
OUTPUT_DIR=~/salmon_output
THREADS=4

mkdir -p "${OUTPUT_DIR}"

#for sample in {1..18}; do   --> This is for using all samples, but we only have 1,2,3,10,11,12 in this dataset.
for sample in 1 2 3 10 11 12; do

    echo "Quantifying sample ${sample} ..."

    salmon quant \
        -i "${INDEX_DIR}" \
        -l A \
        -1 "${FASTQ_DIR}/${sample}_S${sample}_L004_R1_001.fastq.gz" \
           "${FASTQ_DIR}/${sample}_S${sample}_L005_R1_001.fastq.gz" \
        -2 "${FASTQ_DIR}/${sample}_S${sample}_L004_R2_001.fastq.gz" \
           "${FASTQ_DIR}/${sample}_S${sample}_L005_R2_001.fastq.gz" \
        --validateMappings \
        -p "${THREADS}" \
        -o "${OUTPUT_DIR}/Sample_${sample}"

    echo "Sample ${sample} complete. Results in ${OUTPUT_DIR}/Sample_${sample}"

done

# ---------------------------------------------------------------------------
# Mapping rates
# ---------------------------------------------------------------------------

echo ""
echo "Mapping rates:"
for d in "${OUTPUT_DIR}"/Sample_*; do
    rate=$(grep percent_mapped "${d}/aux_info/meta_info.json" | tr -d ' ,"' | cut -d: -f2)
    printf "  %-12s %s%%\n" "$(basename "${d}")" "${rate}"
done
