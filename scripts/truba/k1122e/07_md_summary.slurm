#!/bin/bash
#SBATCH --job-name=k13_summary_v3
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_v3_summary_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_v3_summary_%j.err

set -euo pipefail
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_md_summary_XXXXXX)
tar -xzf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK"
export PYTHONPATH="$LOCAL_WORK/site:$BASE/mcp/atlas_server/src"
export ATLAS_PROJECT_ROOT="$BASE"
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
cd "$BASE"
"$PY" scripts/k1122e_md_v3/analyze.py --run-root "$BASE/data/interim/k1122e_md_production_v3" --output-root "$BASE/data/interim/k1122e_md_results_v3"
# Preserve the finished scientific artifacts off temporary scratch as one file.
RELEASE=/arf/home/avural/kif13a_k1122e_20261004_release
mkdir -p "$RELEASE"
test ! -e "$RELEASE/md_v3_${SLURM_JOB_ID}.tar.gz"
tar -czf "$RELEASE/md_v3_${SLURM_JOB_ID}.tar.gz" scripts/k1122e_md_v3 data/interim/k1122e_md_production_v3 data/interim/k1122e_md_results_v3 data/interim/k1122e_md_opc_pilot_v3 data/interim/k1122e_rna_region data/interim/k1122e_energy data/raw/PROVENANCE.jsonl logs/md_v3*.out logs/md_v3*.err logs/opc_pilot*.out logs/opc_pilot*.err
sha256sum "$RELEASE/md_v3_${SLURM_JOB_ID}.tar.gz" > "$RELEASE/md_v3_${SLURM_JOB_ID}.tar.gz.sha256"
echo "Replica comparison and archive complete: $RELEASE/md_v3_${SLURM_JOB_ID}.tar.gz"
