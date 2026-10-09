#!/bin/bash
#SBATCH --job-name=k13_preflight
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/preflight_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/preflight_%j.err

set -euo pipefail
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_preflight_XXXXXX)
export TMPDIR="$LOCAL_WORK"
export PIP_CACHE_DIR="$LOCAL_WORK/pip-cache"
export OPENMM_CPU_THREADS=8
tar -xzf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK"
# Add GPU plugins in a separate node-local tree, preserving the CPU bundle.
"$PY" -m pip install --no-deps --target "$LOCAL_WORK/gpu" 'openmm-cuda-12==8.3.1'
tar -czf "$BASE/runtime/md_cuda12_plugins.tar.gz" -C "$LOCAL_WORK" gpu
sha256sum "$BASE/runtime/md_cuda12_plugins.tar.gz" > "$BASE/runtime/md_cuda12_plugins.tar.gz.sha256"
export PYTHONPATH="$LOCAL_WORK/site:$BASE/mcp/atlas_server/src"
export ATLAS_PROJECT_ROOT="$BASE"
cd "$BASE"
"$PY" scripts/k1122e_md/build.py --task-index 0 --output-root "$BASE/data/interim/k1122e_md_preflight" --input-root "$BASE" --platform CPU
"$PY" scripts/k1122e_md/build.py --task-index 3 --output-root "$BASE/data/interim/k1122e_md_preflight" --input-root "$BASE" --platform CPU

# Prevent CPython 3.11 site packages leaking into system CPython 3.12 containers.
unset PYTHONPATH
apptainer exec --cleanenv --no-home /arf/sw/containers/truba-ai/pytorch-25.10-py3.sif python -c 'import sys,torch,importlib.util; print(sys.version,torch.__version__); print({m:bool(importlib.util.find_spec(m)) for m in ["omegaconf","Bio","pandas","httpx","yaml"]})'
apptainer exec --cleanenv --no-home /arf/sw/containers/truba-ai/tensorflow-25.02-tf2-py3.sif python -c 'import sys,tensorflow as tf,importlib.util; print(sys.version,tf.__version__); print({m:bool(importlib.util.find_spec(m)) for m in ["keras","tf_keras","h5py","httpx","yaml"]})'
curl --fail --location --retry 3 --output "$BASE/data/raw/alphafold/AF-Q9H1H9-F1-model_v6.pdb" https://alphafold.ebi.ac.uk/files/AF-Q9H1H9-F1-model_v6.pdb
sha256sum "$BASE/data/raw/alphafold/AF-Q9H1H9-F1-model_v6.pdb" > "$BASE/data/raw/alphafold/AF-Q9H1H9-F1-model_v6.pdb.sha256"
echo "Preflight completed, job $SLURM_JOB_ID"
