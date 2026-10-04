#!/bin/bash
#SBATCH -J eqtl_dl_bell
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=01:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Bellenguez 2022 AD GWAS (GCST90027158, GRCh38): the AD GWAS cited by Zeng 2026, used for the
# Zeng oligodendrocyte positive control (CR1 / ADAMTS4 / APH1B x AD).
set -uo pipefail
B=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158
cd /arf/scratch/avural/kif13a_eqtl/raw/gwas
for f in GCST90027158_buildGRCh38.tsv.gz GCST90027158_buildGRCh38.tsv.gz-meta.yaml README.txt md5sum.txt; do
  o=$f; [ "$f" = README.txt ] && o=GCST90027158.README.txt; [ "$f" = md5sum.txt ] && o=GCST90027158.md5sum.txt
  [ -s $o ] || curl -sS -L --retry 5 -o $o.part $B/$f && mv $o.part $o
done
md5sum GCST90027158_buildGRCh38.tsv.gz; cat GCST90027158.md5sum.txt
