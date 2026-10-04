#!/bin/bash
#SBATCH --job-name=k13_an
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_insan_yas/logs/an_%x_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_insan_yas/logs/an_%x_%j.err
# Donor-level statistics on pseudobulk outputs. Variable: KEYS (colon separated dataset keys), SCRIPT (default analyze.py)
set -euo pipefail
BASE=/arf/scratch/avural/kif13a_insan_yas
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd "$BASE/scripts"
echo "Job $SLURM_JOB_ID on $(hostname), start $(date), KEYS=${KEYS:-}"
python -u "${SCRIPT:-analyze.py}" --base "$BASE" ${KEYS:+--keys "$KEYS"}
echo "End $(date)"
