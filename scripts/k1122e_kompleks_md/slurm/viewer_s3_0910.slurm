#!/bin/bash
# 2026-10-09: video icin s3 uclusunun kareleri (WT, K1122R, K1122E paired; ayni baslangic yapisi).
# Yalniz okur; cikti complex_md/viewer_0910 altina yazilir, uretim klasorlerine dokunulmaz.
#SBATCH --account=avural
#SBATCH -p debug
#SBATCH -J k13_viewer_s3
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
set -uo pipefail
SIF=/arf/home/avural/containers/mdana-openmm831-pyscf.sif
H=/tmp/${USER}_home_${SLURM_JOB_ID}; mkdir -p $H
trap 'rm -rf $H' EXIT
export OMP_NUM_THREADS=8
ROOT=/arf/scratch/avural/so3lr/complex_md; OUT=$ROOT/viewer_0910
mkdir -p $OUT
cd $ROOT/scripts
apptainer exec -B /arf --home $H $SIF python export_viewer2.py \
  --run WT_wt_s3=$ROOT/runs/WT_wt_s3 \
  --run K1122R_paired_wt_s3=$ROOT/runs/K1122R_paired_wt_s3 \
  --run K1122E_paired_wt_s3=$ROOT/runs/K1122E_paired_wt_s3 \
  --stride3d 20 --stride-chart 10 --out $OUT/viewer_s3_frames.json
echo "rc=$? $(date -Is)"
ls -la $OUT
