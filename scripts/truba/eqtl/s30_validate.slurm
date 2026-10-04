#!/bin/bash
#SBATCH -J eqtl_rval
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=00:45:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# R coloc/susieR reference run, then the Python-vs-R comparison.
set -uo pipefail
source /etc/profile.d/modules.sh
module load apps/R/4.5.2-gcc-14.1.0
cd /arf/scratch/avural/kif13a_eqtl/scripts
Rscript s30_validate.R
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
EX_N=$(grep -m1 "^example N:" ../logs/eqtl_rval_${SLURM_JOB_ID}.out | awk '{print $3}')
EX_SDY=$(grep -m1 "^example N:" ../logs/eqtl_rval_${SLURM_JOB_ID}.out | awk '{print $6}')
echo "EX_N=$EX_N EX_SDY=$EX_SDY"
EX_N=$EX_N EX_SDY=$EX_SDY python -u s31_compare.py
