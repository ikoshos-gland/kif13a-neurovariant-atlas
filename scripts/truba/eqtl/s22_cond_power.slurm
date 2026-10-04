#!/bin/bash
#SBATCH -J eqtl_cond2
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Second round: scale-consistent conditional coloc + LD diagnostics, sdY sensitivity, power.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s20_coloc.py sdy cond &
python -u s20_coloc.py power > ../logs/power2_${SLURM_JOB_ID}.log 2>&1 &
wait
echo done
