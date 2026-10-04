#!/bin/bash
#SBATCH --job-name=k13_pb
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_insan_yas/logs/pb_%x_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_insan_yas/logs/pb_%x_%j.err
# Pseudobulk one dataset. Variables: KEY (file stem in raw/), GROUP (obs columns), WORKERS, EXTRA (e.g. --inspect-only)
set -euo pipefail
BASE=/arf/scratch/avural/kif13a_insan_yas
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd "$BASE"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "Job $SLURM_JOB_ID on $(hostname), start $(date), KEY=$KEY GROUP=$GROUP"
python -u scripts/pseudobulk.py --h5ad "raw/$KEY.h5ad" --out "out/$KEY" --group "$GROUP" --workers "${WORKERS:-6}" ${EXTRA:-}
echo "End $(date)"
