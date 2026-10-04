#!/bin/bash
#SBATCH -J eqtl_pipe
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# 1000G window -> build test, positive controls, main coloc, lookups, conditional, R export;
# then (in the background of the same job) the power simulation, then the R reference run.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s10_extract.py kg
python -u s20_coloc.py build posctrl
python -u s20_coloc.py power > ../logs/power_${SLURM_JOB_ID}.log 2>&1 &
python -u s20_coloc.py main lookup cond export
source /etc/profile.d/modules.sh
module load apps/R/4.5.2-gcc-14.1.0
Rscript s30_validate.R
EX_N=$(grep -m1 "^example N:" ../logs/eqtl_pipe_${SLURM_JOB_ID}.out | awk '{print $3}')
EX_SDY=$(grep -m1 "^example N:" ../logs/eqtl_pipe_${SLURM_JOB_ID}.out | awk '{print $6}')
echo "EX_N=$EX_N EX_SDY=$EX_SDY"
EX_N=$EX_N EX_SDY=$EX_SDY python -u s31_compare.py
wait
echo "pipeline done"
