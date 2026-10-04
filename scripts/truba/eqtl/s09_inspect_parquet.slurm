#!/bin/bash
#SBATCH -J eqtl_pqinspect
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --time=00:15:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
python -u /arf/scratch/avural/kif13a_eqtl/scripts/s09_inspect_parquet.py /arf/scratch/avural/kif13a_eqtl/raw/synapse/meta-Oligo-class.parquet
