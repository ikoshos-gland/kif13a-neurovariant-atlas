#!/bin/bash
#SBATCH -J eqtl_syn_dl
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=03:50:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Zeng ... Roussos 2026 Synapse downloads + pyarrow (pip --target, project-local) for parquet.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
W=/arf/scratch/avural/kif13a_eqtl
cd $W/scripts
if [ ! -d $W/pylib/pyarrow ]; then
  pip install --no-cache-dir --target $W/pylib pyarrow 2>&1 | tail -2
fi
PYTHONPATH=$W/pylib python -c "import pyarrow; print('pyarrow', pyarrow.__version__)"
python -u s03_synapse_download.py small 3
echo "small exit=$?"
python -u s03_synapse_download.py class 3
echo "class exit=$?"
ls -la $W/raw/synapse
