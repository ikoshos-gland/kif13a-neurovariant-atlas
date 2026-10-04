#!/bin/bash
#SBATCH -J eqtl_sbrain
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 7
#SBATCH --time=03:55:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# SingleBrain (Zenodo 14908182) full associations, streamed and filtered to target genes.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl/scripts
for ct in "$@"; do python -u s06_singlebrain_stream.py $ct & done
wait
# PD GWAS (Nalls 2019, GCST009325 harmonised) for the SingleBrain TOX3 positive control
G=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised
cd /arf/scratch/avural/kif13a_eqtl/raw/gwas
[ -s GCST009325.h.tsv.gz ] || curl -sS -L --retry 5 -o GCST009325.h.tsv.gz.part $G/GCST009325.h.tsv.gz && mv GCST009325.h.tsv.gz.part GCST009325.h.tsv.gz
curl -sS -L -o GCST009325.h.tsv.gz-meta.yaml $G/GCST009325.h.tsv.gz-meta.yaml
curl -sS -L -o GCST009325.md5sum.txt $G/md5sum.txt
md5sum GCST009325.h.tsv.gz
echo done
