#!/usr/bin/env bash
set -euo pipefail

DEST_DIR="test_data"
TMP_CLONE="$(mktemp -d)"

echo ">>> Cloning nf-core/test-datasets (rnaseq branch, shallow)..."
git clone --depth 1 --single-branch --branch rnaseq \
  https://github.com/nf-core/test-datasets.git "$TMP_CLONE"

echo ">>> Copying reference files..."
mkdir -p "${DEST_DIR}/reference"
cp "${TMP_CLONE}/reference/genome.fa" "${DEST_DIR}/reference/genome.fa"
cp "${TMP_CLONE}/reference/genes.gtf" "${DEST_DIR}/reference/genes.gtf"

echo ">>> Copying paired-end read files (GSE110004 subset)..."
mkdir -p "${DEST_DIR}/reads"
READS_SRC="${TMP_CLONE}/testdata/GSE110004"

cp "${READS_SRC}/SRR6357070_1.fastq.gz" "${DEST_DIR}/reads/ctrl_rep1_R1.fastq.gz"
cp "${READS_SRC}/SRR6357070_2.fastq.gz" "${DEST_DIR}/reads/ctrl_rep1_R2.fastq.gz"
cp "${READS_SRC}/SRR6357072_1.fastq.gz" "${DEST_DIR}/reads/ctrl_rep2_R1.fastq.gz"
cp "${READS_SRC}/SRR6357072_2.fastq.gz" "${DEST_DIR}/reads/ctrl_rep2_R2.fastq.gz"
cp "${READS_SRC}/SRR6357076_1.fastq.gz" "${DEST_DIR}/reads/stress_rep1_R1.fastq.gz"
cp "${READS_SRC}/SRR6357076_2.fastq.gz" "${DEST_DIR}/reads/stress_rep1_R2.fastq.gz"
cp "${READS_SRC}/SRR6357071_1.fastq.gz" "${DEST_DIR}/reads/stress_rep2_R1.fastq.gz"
cp "${READS_SRC}/SRR6357071_2.fastq.gz" "${DEST_DIR}/reads/stress_rep2_R2.fastq.gz"

echo ">>> Cleaning up..."
rm -rf "$TMP_CLONE"

echo ">>> Verifying expected files exist..."
for f in ctrl_rep1_R1 ctrl_rep1_R2 ctrl_rep2_R1 ctrl_rep2_R2 \
         stress_rep1_R1 stress_rep1_R2 stress_rep2_R1 stress_rep2_R2; do
  if [[ ! -s "${DEST_DIR}/reads/${f}.fastq.gz" ]]; then
    echo "ERROR: ${DEST_DIR}/reads/${f}.fastq.gz is missing or empty!" >&2
    exit 1
  fi
done
for f in genome.fa genes.gtf; do
  if [[ ! -s "${DEST_DIR}/reference/${f}" ]]; then
    echo "ERROR: ${DEST_DIR}/reference/${f} is missing or empty!" >&2
    exit 1
  fi
done

echo ">>> All expected files present. Test data ready in ${DEST_DIR}/"
