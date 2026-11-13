#!/bin/bash
#
# Build the decoy-aware Salmon index for S. pombe (ASM294v3)
#
# Only need to run once. The resulting index at ~/reference/salmon_index was used for the
# Zhong-HZ-17633_2024_12_12 quantification (21 December 2024, Salmon v1.10.3).
#
# Reference files, NCBI RefSeq assembly ASM294v3 / GCF_000002945.2:
#   https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/Schizosaccharomyces_pombe/latest_assembly_versions/GCF_000002945.2_ASM294v3/
#     GCF_000002945.2_ASM294v3_rna.fna.gz       transcripts file (NM_/NR_ accession numbers)
#     GCF_000002945.2_ASM294v3_genomic.fna.gz   genome file, used to generate decoys file
#
# Verified against aux_info/meta_info.json of the original run:
#   12390 valid targets, 4 decoy targets, k-mer size 31
#   index_seq_hash e5ddfa1dea4cffb5fba35639b8f02b67d98baf123e090ed7f563d061346fd29c
#
# Used salmon version 1.10.3
#   conda env create -f envs/salmon1.10.yml     (first time only)
#   conda activate salmon1.10

set -euo pipefail # exit on error (-e), unset var (-u), and failing pipeline (-o pipefail)

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

REF_DIR=~/pombe-downloads
INDEX_DIR=~/reference/salmon_index
THREADS=4

GENOME_FASTA="${REF_DIR}/GCF_000002945.2_ASM294v3_genomic.fna.gz"
TRANSCRIPTS_FASTA="${REF_DIR}/GCF_000002945.2_ASM294v3_rna.fna.gz"
DECOYS="${REF_DIR}/decoys.txt"
COMBINED_FASTA="${REF_DIR}/trimmed_combined.fna"

# -----------------------------------------------------------------------------------
# Combined reference FASTA
# -----------------------------------------------------------------------------------
#
# Salmon requires transcript sequences FIRST and decoy (genome) sequences LAST. 
# decoys.txt names the genome records so salmon can exclude them from quantification.

gunzip -c "${GENOME_FASTA}" \
    | grep ">" \
    | cut -d " " -f 1 \
    | sed 's/>//g' \
    > "${DECOYS}"

gunzip -c "${TRANSCRIPTS_FASTA}" "${GENOME_FASTA}" > "${COMBINED_FASTA}"

echo "Combined FASTA records: $(grep -c '>' "${COMBINED_FASTA}")   expected 12394"
echo "Decoy sequences listed: $(wc -l < "${DECOYS}")   expected 4"
echo "Last records in FASTA (should be the 4 genome sequences):"
grep ">" "${COMBINED_FASTA}" | tail -4

# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------

salmon index \
    -t "${COMBINED_FASTA}" \
    -d "${DECOYS}" \
    -i "${INDEX_DIR}" \
    -k 31 \
    --threads "${THREADS}"

echo ""
echo "Index written to ${INDEX_DIR}"
echo "After the first quant run, confirm index_seq_hash in"
echo "  <output>/aux_info/meta_info.json"
echo "matches e5ddfa1dea4cffb5fba35639b8f02b67d98baf123e090ed7f563d061346fd29c"
