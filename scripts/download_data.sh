#!/usr/bin/env bash
# Download GEO GSE179430 processed ChEC-seq coverage and the sacCer3 (S288C R64) genome.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p GSE179430/raw genome
if [ ! -f GSE179430/GSE179430_RAW.tar ]; then
  curl -L -o GSE179430/GSE179430_RAW.tar "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179430/suppl/GSE179430_RAW.tar"
fi
echo "917b9f0d952802b5497d748fb5ccaa5a  GSE179430/GSE179430_RAW.tar" | md5sum -c -
tar -xf GSE179430/GSE179430_RAW.tar -C GSE179430/raw
echo "$(ls GSE179430/raw | wc -l) sample files extracted (expect 348)"
if [ ! -f genome/sacCer3.fa ]; then
  curl -L "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz" | gunzip -c > genome/sacCer3.fa
fi
command -v samtools >/dev/null && samtools faidx genome/sacCer3.fa || echo "samtools not found; genome/sacCer3.fa.fai from the repo is valid for this file"
