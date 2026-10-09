#!/bin/bash
#SBATCH --job-name=k13_opc_pilot
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --array=0,3,15%3
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/opc_pilot_%A_%a.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/opc_pilot_%A_%a.err

set -euo pipefail
source /etc/profile.d/modules.sh
module load lib/cuda/12.4
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_cube_pilot_XXXXXX)
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
OUT="$BASE/data/interim/k1122e_md_opc_pilot_v3"
if [ "$SLURM_ARRAY_TASK_ID" -eq 0 ]; then
    "$PY" "$CODE/validate_caps.py" --input-root "$BASE" --platform CUDA --output "$OUT/cap_validation.json"
fi
TASK=$($PY -c 'import sys; sys.path.insert(0,"scripts/k1122e_md_v3"); from common import task_at,DEFAULT_CONFIG; import os; print(task_at(DEFAULT_CONFIG,int(os.environ["SLURM_ARRAY_TASK_ID"]))[1]["task_id"])')
"$PY" "$CODE/build.py" --task-index "$SLURM_ARRAY_TASK_ID" --output-root "$OUT" --input-root "$BASE" --platform CUDA
# Full 1.2 ns equilibration followed by 0.1 ns production, excluded from inference.
"$PY" "$CODE/run.py" --task-dir "$OUT/$TASK" --platform CUDA --duration-ns 0.1
"$PY" "$CODE/verify_pilot.py" --task-dir "$OUT/$TASK"
echo "ff19SB/OPC pilot passed: $TASK; excluded from scientific inference."
