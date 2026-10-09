#!/bin/bash
#SBATCH --job-name=k13_md_smoke
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_smoke_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_smoke_%j.err

set -euo pipefail
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_md_smoke_XXXXXX)
tar -xzf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK"
export PYTHONPATH="$LOCAL_WORK/site:$BASE/mcp/atlas_server/src"
export ATLAS_PROJECT_ROOT="$BASE"
export OPENMM_CPU_THREADS=8
export OMP_NUM_THREADS=8
cd "$BASE"
OUT="$BASE/data/interim/k1122e_md_smoke_v3"
for INDEX in 0 3; do
    "$PY" scripts/k1122e_md/build.py --task-index "$INDEX" --output-root "$OUT" --input-root "$BASE" --platform CPU
done
for TASK in primary_WT_r01 primary_K1122E_r01; do
    "$PY" scripts/k1122e_md/run.py --task-dir "$OUT/$TASK" --platform CPU --cpu-threads 8 --smoke
    "$PY" scripts/k1122e_md/analyze.py --task-dir "$OUT/$TASK" --allow-smoke
done
echo "CPU engine smoke passed; these short traces are excluded from scientific inference."
