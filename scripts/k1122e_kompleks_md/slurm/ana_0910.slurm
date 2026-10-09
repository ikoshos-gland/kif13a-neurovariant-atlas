#!/bin/bash
# 2026-10-09: son analiz.
# A) kompleks MD: tam analizi olmayan bitmis kosular (s2'lerle ayni bayraklar: --gbsa --solvent) -> results_0910
# B) monomer v3: biten gorevlerin ek olcumleri (dssp_monomer.py), ayri klasorde; uretim klasorlerine yazilmaz.
#SBATCH --account=avural
#SBATCH -p debug
#SBATCH -J k13_ana_0910
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
set -uo pipefail
SIF=/arf/home/avural/containers/mdana-openmm831-pyscf.sif
H=/tmp/${USER}_home_${SLURM_JOB_ID}; mkdir -p $H
trap 'rm -rf $H' EXIT
PY() { apptainer exec -B /arf --home $H $SIF python "$@"; }

# ---------- B) monomer v3 (hafif, arka planda sirayla) ----------
BASE=/arf/scratch/avural/kif13a_k1122e_20261004
MOUT=$BASE/data/interim/k1122e_md_full_v3_0910
EK=$BASE/scripts/onizleme_1810
mkdir -p $MOUT
(
  for T in $BASE/data/interim/k1122e_md_production_v3/*/; do
    ID=$(basename $T)
    grep -q '"completed": true' $T/run.json || continue
    mkdir -p $MOUT/$ID
    for f in build.json run.json protein.pdb trajectory.dcd; do ln -sf $T/$f $MOUT/$ID/$f; done
    cp $T/analysis.json $T/alpha_residue_occupancies.tsv $T/contact_occupancies.tsv $T/metrics_timeseries.tsv $MOUT/$ID/
    OMP_NUM_THREADS=2 PY $EK/dssp_monomer.py $MOUT/$ID > $MOUT/$ID.ek.log 2>&1
    echo "MON $ID rc=$?"
  done
) &

# ---------- A) kompleks MD ----------
ROOT=/arf/scratch/avural/so3lr/complex_md; RES=$ROOT/results_0910
export OMP_NUM_THREADS=7 OPENMM_CPU_THREADS=7
mkdir -p $RES
cd $ROOT/scripts
PIDS=()
for RUN in $ROOT/runs/*/; do
  ID=$(basename $RUN)
  [ -f $ROOT/results_0810/$ID/analysis.json ] && continue
  grep -q '"completed": true' $RUN/run.json || { echo "CX $ID bitmemis, atlandi"; continue; }
  mkdir -p $RES/$ID
  ( PY analyze_complex.py run --run-dir $RUN --system-dir $ROOT/systems/$ID --out-dir $RES/$ID --gbsa --solvent > $RES/$ID.log 2>&1
    echo "CX $ID tam rc=$?" ) &
  PIDS+=($!)
  if (( ${#PIDS[@]} == 4 )); then wait "${PIDS[@]}"; PIDS=(); fi
done
wait
echo BITTI $(date -Is)
