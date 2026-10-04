#!/bin/bash
#SBATCH -J eqtl_syn_retry
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=03:55:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# One retry after the user replaced the token: scope check (names/scopes/state only), then the
# small result tables, then the 8 class-level parquet files. Stops at the first 401/403.
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl/scripts
S=/arf/scratch/avural/kif13a_eqtl/raw/synapse
[ -f $S/AUTH_FAIL.txt ] && mv $S/AUTH_FAIL.txt $S/AUTH_FAIL_first_token.txt
[ -f $S/token_scopes.json ] && mv $S/token_scopes.json $S/token_scopes_first_token.json
python -u s04c_synapse_token_scopes.py | grep -v userName
python -u s03_synapse_download.py small 3
rc=$?
echo "small exit=$rc"
if [ $rc -ne 0 ]; then echo "stopping: small-table download failed"; exit $rc; fi
python -u s03_synapse_download.py class 4
echo "class exit=$?"
ls -la $S
