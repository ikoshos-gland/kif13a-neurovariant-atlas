#!/bin/bash
#SBATCH -J eqtl_extract
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=01:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Window extraction from the raw downloads (no internet needed; light, fits debug).
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s10_extract.py "$@"
