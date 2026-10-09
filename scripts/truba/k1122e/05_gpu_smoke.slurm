#!/bin/bash
#SBATCH --job-name=k13_gpu_smoke
#SBATCH --partition=akya-cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --gres=gpu:1
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/gpu_smoke_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/gpu_smoke_%j.err

set -euo pipefail
source /etc/profile.d/modules.sh
module load lib/cuda/12.4
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_gpu_smoke_XXXXXX)
tar -xzf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK"
tar -xzf "$BASE/runtime/md_cuda12_plugins.tar.gz" -C "$LOCAL_WORK"
cp "$LOCAL_WORK/gpu/OpenMM.libs/lib/plugins/"*.so "$LOCAL_WORK/site/OpenMM.libs/lib/plugins/"
export PYTHONPATH="$LOCAL_WORK/site:$BASE/mcp/atlas_server/src"
export ATLAS_PROJECT_ROOT="$BASE"
export OPENMM_CPU_THREADS=10
export OMP_NUM_THREADS=10
cd "$BASE"
nvidia-smi
"$PY" -m openmm.testInstallation
OUT="$BASE/data/interim/k1122e_md_gpu_smoke"
for INDEX in 0 3; do
    "$PY" scripts/k1122e_md/build.py --task-index "$INDEX" --output-root "$OUT" --input-root "$BASE" --platform CUDA
done
for TASK in primary_WT_r01 primary_K1122E_r01; do
    "$PY" scripts/k1122e_md/run.py --task-dir "$OUT/$TASK" --platform CUDA --smoke
    "$PY" scripts/k1122e_md/analyze.py --task-dir "$OUT/$TASK" --allow-smoke
done
echo "GPU engine smoke passed; short traces are excluded from scientific inference."
