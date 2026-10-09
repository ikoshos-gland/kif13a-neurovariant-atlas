#!/bin/bash
#SBATCH --job-name=k13_inference
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=/arf/scratch/avural/kif13a_k1122e_20261004/logs/inference_%A_%a.out
#SBATCH --error=/arf/scratch/avural/kif13a_k1122e_20261004/logs/inference_%A_%a.err

set -euo pipefail
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
LOCAL_WORK=$(mktemp -d /tmp/k13_infer_XXXXXX)
export TMPDIR="$LOCAL_WORK"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS="$SLURM_CPUS_PER_TASK"
cd "$BASE"
BUNDLE="$BASE/data/interim/k1122e_rna_region/KIF13A_splice_bundle_20261004.json"
TORCH=/arf/sw/containers/truba-ai/pytorch-25.10-py3.sif
TF=/arf/sw/containers/truba-ai/tensorflow-25.02-tf2-py3.sif
if [ "$SLURM_ARRAY_TASK_ID" = 0 ]; then
    for DISTANCE in 50 500; do
        apptainer exec --cleanenv --no-home --env "PYTHONPATH=$BASE/mcp/atlas_server/src,ATLAS_PROJECT_ROOT=$BASE,TF_CPP_MIN_LOG_LEVEL=3,CUDA_VISIBLE_DEVICES=-1" "$TF" python scripts/k1122e_rna_region/infer_splicing.py --model spliceai --bundle "$BUNDLE" --model-dir "$BASE/data/raw/k1122e_models/spliceai" --distance "$DISTANCE" --date 20261004 --output-base "$BASE/data/interim" --device cpu --threads 8
    done
elif [ "$SLURM_ARRAY_TASK_ID" = 1 ]; then
    for DISTANCE in 50 500; do
        apptainer exec --cleanenv --no-home --env "PYTHONPATH=$BASE/mcp/atlas_server/src,ATLAS_PROJECT_ROOT=$BASE,OMP_NUM_THREADS=8" "$TORCH" python scripts/k1122e_rna_region/infer_splicing.py --model pangolin --bundle "$BUNDLE" --model-dir "$BASE/data/raw/k1122e_models/pangolin" --pangolin-model-code "$BASE/data/raw/k1122e_models/pangolin/model.py" --distance "$DISTANCE" --date 20261004 --output-base "$BASE/data/interim" --device cpu --threads 8
    done
else
    tar -xzf "$BASE/tools/ThermoMPNN_370f76e.tar.gz" -C "$LOCAL_WORK"
    # Missing inference dependencies go only to node-local /tmp, no /arf installs.
    apptainer exec --cleanenv --no-home --bind "$LOCAL_WORK:/run/k13" "$TORCH" python -m pip install --no-deps --target /run/k13/deps 'biopython==1.85' 'omegaconf==2.3.0' 'antlr4-python3-runtime==4.9.3' 'pandas==2.2.3'
    apptainer exec --cleanenv --no-home --bind "$LOCAL_WORK:/run/k13" --env "PYTHONPATH=/run/k13/deps:$BASE/mcp/atlas_server/src,ATLAS_PROJECT_ROOT=$BASE,OMP_NUM_THREADS=8,OPENBLAS_NUM_THREADS=8" "$TORCH" python scripts/k1122e_energy/thermompnn.py --source /run/k13/ThermoMPNN-370f76ec62bd929f7425e311d8df04a0d094990f --pdb "$BASE/data/raw/alphafold/AF-Q9H1H9-F1-model_v6.pdb" --sequence "$BASE/data/raw/uniprot/Q9H1H9_20260915.json" --output "$BASE/data/interim/k1122e_energy"
fi
echo "Inference completed: $SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID"
