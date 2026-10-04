#!/bin/bash
#SBATCH -J eqtl_syn_list
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Lists the Synapse tree of Zeng ... Roussos 2026. The token is read inside python only.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl/scripts
python -u s00_synapse_list.py
echo "exit=$?"
