#!/bin/bash
#SBATCH -J eqtl_big40
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 3
#SBATCH --time=02:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
mkdir -p /arf/scratch/avural/kif13a_eqtl/work
cd /arf/scratch/avural/kif13a_eqtl/scripts
for p in 1769 1925 1475; do python -u s05_big40_stream.py $p & done
wait
echo done
