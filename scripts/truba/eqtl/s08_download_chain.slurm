#!/bin/bash
#SBATCH -J eqtl_dl_chain
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:15:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# UCSC hg38 -> hg19 chain (for lifting Zeng 2026 / SingleBrain hg38 positions to GRCh37)
set -uo pipefail
mkdir -p /arf/scratch/avural/kif13a_eqtl/raw/ucsc && cd /arf/scratch/avural/kif13a_eqtl/raw/ucsc
curl -sS -L --retry 5 -o hg38ToHg19.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
ls -la; md5sum hg38ToHg19.over.chain.gz
