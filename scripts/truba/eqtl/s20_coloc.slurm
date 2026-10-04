#!/bin/bash
#SBATCH -J eqtl_coloc
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=02:00:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# usage: sbatch s20_coloc.slurm build posctrl main lookup cond export   (then s30 R, s31 compare)
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s20_coloc.py "$@"
