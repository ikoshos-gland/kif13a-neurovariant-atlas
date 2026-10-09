#!/bin/bash
#SBATCH --job-name=k13_env
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/env_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/env_%j.err

# Build package trees exclusively in compute-node /tmp. Preserve one archive.
set -euo pipefail
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
PY=/arf/home/avural/miniconda3/envs/allenbrain/bin/python
LOCAL_WORK=$(mktemp -d /tmp/k13_env_XXXXXX)
mkdir -p "$LOCAL_WORK/site" "$BASE/runtime" "$BASE/tools"
export TMPDIR="$LOCAL_WORK"
export PIP_CACHE_DIR="$LOCAL_WORK/pip-cache"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
"$PY" -m pip install --target "$LOCAL_WORK/site" 'openmm==8.3.1' 'mdtraj==1.11.0' 'numpy<2.3' 'scipy<1.17' biopython PyYAML httpx omegaconf
export PYTHONPATH="$LOCAL_WORK/site"
"$PY" -m pip install --no-deps --target "$LOCAL_WORK/site" https://github.com/openmm/pdbfixer/archive/refs/tags/v1.12.tar.gz
"$PY" -c 'import openmm,pdbfixer,mdtraj,sys; print(sys.version); print("openmm",openmm.__version__,"mdtraj",mdtraj.__version__); print("platforms",[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])'
"$PY" -m pip freeze --path "$LOCAL_WORK/site" > "$BASE/runtime/packages.txt"
tar -czf "$BASE/runtime/md_py311_site.tar.gz" -C "$LOCAL_WORK" site
sha256sum "$BASE/runtime/md_py311_site.tar.gz" > "$BASE/runtime/md_py311_site.tar.gz.sha256"

# Pinned official ThermoMPNN source and pretrained weights, no environment edits.
curl --fail --location --retry 3 --output "$BASE/tools/ThermoMPNN_370f76e.tar.gz" https://codeload.github.com/Kuhlman-Lab/ThermoMPNN/tar.gz/370f76ec62bd929f7425e311d8df04a0d094990f
sha256sum "$BASE/tools/ThermoMPNN_370f76e.tar.gz" > "$BASE/tools/ThermoMPNN_370f76e.tar.gz.sha256"

# Inspect framework compatibility inside existing system containers on this node.
apptainer exec --no-home /arf/sw/containers/truba-ai/pytorch-25.10-py3.sif python -c 'import sys,torch,importlib.util; print(sys.version,torch.__version__); print({m:bool(importlib.util.find_spec(m)) for m in ["pytorch_lightning","lightning","omegaconf","Bio","pandas","sklearn"]})'
apptainer exec --no-home /arf/sw/containers/truba-ai/tensorflow-25.02-tf2-py3.sif python -c 'import sys,tensorflow as tf,importlib.util; print(sys.version,tf.__version__); print({m:bool(importlib.util.find_spec(m)) for m in ["keras","tf_keras","h5py"]})'
echo "Environment preparation completed, job $SLURM_JOB_ID"
