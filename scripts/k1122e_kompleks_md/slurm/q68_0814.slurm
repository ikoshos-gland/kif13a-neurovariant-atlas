#!/bin/bash
# 2026-10-14: Rab10 Q68 ve gama-fosfat suyu, 12 kompleks kosu (10-100 ns), sonra 6'ya 6 ozet. Yalniz okur: runs/, runs_solvent/.
#SBATCH --account=avural
#SBATCH -p debug
#SBATCH -J k13_q68
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
set -uo pipefail
ROOT=/arf/scratch/avural/so3lr/complex_md; OUT=$ROOT/results_q68
SIF=/arf/home/avural/containers/mdana-openmm831-pyscf.sif
export OMP_NUM_THREADS=4
H=/tmp/${USER}_home_${SLURM_JOB_ID}; mkdir -p $H $OUT; trap 'rm -rf $H' EXIT
PY() { apptainer exec -B /arf --home $H $SIF python "$@"; }
cd $ROOT/scripts
N=0
for RUN in $ROOT/runs/*/; do
  ID=$(basename $RUN)
  ( PY q68_analiz.py run --run-dir $RUN --system-dir $ROOT/systems/$ID --out $OUT/$ID.json > $OUT/$ID.log 2>&1; echo "$ID rc=$? $(grep '^{' $OUT/$ID.log | tail -1)" ) &
  N=$((N+1)); (( N % 4 == 0 )) && wait
done
wait
PY q68_analiz.py ozet --in-dir $OUT --out $OUT/ozet
echo BITTI $(date -Is)
