#!/bin/bash
#SBATCH -J eqtl_extract2
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# GWAS windows (incl. Bellenguez AD for the Zeng positive control) + Zeng 2026 parquet extraction.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s10_extract.py gwas &
python -u s12_zeng_extract.py &
wait
echo done
