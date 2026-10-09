#!/bin/bash
#SBATCH --job-name=k13_md_v3
#SBATCH --partition=akya-cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=1-12:00:00
#SBATCH --array=0-29%4
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_v3_%A_%a.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/md_v3_%A_%a.err

set -euo pipefail
source /etc/profile.d/modules.sh
module load lib/cuda/12.4
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_md_array_XXXXXX)
tar -xzf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK"
tar -xzf "$BASE/runtime/md_cuda12_plugins.tar.gz" -C "$LOCAL_WORK"
cp "$LOCAL_WORK/gpu/OpenMM.libs/lib/plugins/"*.so "$LOCAL_WORK/site/OpenMM.libs/lib/plugins/"
export PYTHONPATH="$LOCAL_WORK/site:$BASE/mcp/atlas_server/src"
export ATLAS_PROJECT_ROOT="$BASE"
export OPENMM_CPU_THREADS=10
export OMP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=1
cd "$BASE"
CODE="$BASE/scripts/k1122e_md_v3"
OUT="$BASE/data/interim/k1122e_md_production_v3"
# Submit with afterok on the complete pilot array, never a smoke job.
for PILOT in primary_WT_r01 primary_K1122E_r01 padded_WT_r01; do
    test -f "$BASE/data/interim/k1122e_md_opc_pilot_v3/$PILOT/pilot_validation.json"
done
TASK=$($PY -c 'import sys; sys.path.insert(0,"scripts/k1122e_md_v3"); from common import task_at,DEFAULT_CONFIG; import os; print(task_at(DEFAULT_CONFIG,int(os.environ["SLURM_ARRAY_TASK_ID"]))[1]["task_id"])')
if [ ! -f "$OUT/$TASK/build.json" ]; then
    "$PY" "$CODE/build.py" --task-index "$SLURM_ARRAY_TASK_ID" --output-root "$OUT" --input-root "$BASE" --platform CUDA
fi
if [ -f "$OUT/$TASK/run.json" ]; then
    "$PY" "$CODE/run.py" --task-dir "$OUT/$TASK" --platform CUDA --resume
else
    "$PY" "$CODE/run.py" --task-dir "$OUT/$TASK" --platform CUDA
fi
"$PY" "$CODE/analyze.py" --task-dir "$OUT/$TASK"
echo "Scientific MD task completed: $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID $TASK"
