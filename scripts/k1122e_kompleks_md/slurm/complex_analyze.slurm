#!/bin/bash
# Endpoints of every finished complex run + replica aggregation (CPU only, mdana image). WRITTEN 2026-10-05, NOT SUBMITTED.
#   sbatch /arf/scratch/avural/so3lr/complex_md/slurm/complex_analyze.slurm            # runs/ -> results/
#   sbatch --export=ALL,RUNS=runs_pilot,RES=results_pilot,ALLOW_SMOKE=1,DISCARD_NS=0 .../complex_analyze.slurm   # pilot
# Read-only on runs/ and runs_solvent/: per-run outputs go to $RES/<id>/, the aggregate to $RES/aggregate/.
#SBATCH --account=avural
#SBATCH -p debug
#SBATCH -J k13_cx_ana
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/so3lr/logs/%x_%j.err

set -uo pipefail
ROOT=/arf/scratch/avural/so3lr/complex_md
RUNS=$ROOT/${RUNS:-runs}
RES=$ROOT/${RES:-results}
SIF=/arf/home/avural/containers/mdana-openmm831-pyscf.sif
SMOKE=""; [[ -n "${ALLOW_SMOKE:-}" ]] && SMOKE="--allow-smoke"
DISCARD=""; [[ -n "${DISCARD_NS:-}" ]] && DISCARD="--discard-ns $DISCARD_NS"   # pilots: DISCARD_NS=0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENMM_CPU_THREADS=$SLURM_CPUS_PER_TASK
HOMEDIR=/tmp/${USER}_home_${SLURM_JOB_ID}
mkdir -p "$HOMEDIR" "$RES"
trap 'rm -rf "$HOMEDIR"' EXIT
PY() { apptainer exec -B /arf --home "$HOMEDIR" "$SIF" python "$@"; }
cd "$ROOT/scripts" || exit 1
PY analyze_complex.py selftest --out-dir "$RES/selftest" --system-dir "$ROOT/systems/WT_wt_s5" > "$RES/selftest.log" 2>&1 \
  || { echo "analysis selftest FAILED"; tail -40 "$RES/selftest.log"; exit 1; }
FAIL=0
for RUN in "$RUNS"/*/; do
  ID=$(basename "$RUN")
  grep -q '"completed": true' "$RUN/run.json" 2>/dev/null || { echo "skip unfinished $ID"; continue; }
  mkdir -p "$RES/runs/$ID/analysis"
  cp "$RUN/run.json" "$RUN/trajectory_atoms.json" "$RES/runs/$ID/"
  PY analyze_complex.py run --run-dir "$RUN" --system-dir "$ROOT/systems/$ID" --out-dir "$RES/runs/$ID/analysis" \
     --gbsa --solvent $SMOKE $DISCARD > "$RES/runs/$ID/analysis.log" 2>&1 || { echo "FAILED $ID"; FAIL=1; continue; }
  echo "analysed $ID"
done
PY analyze_complex.py aggregate --runs-root "$RES/runs" --out-dir "$RES/aggregate" $SMOKE
echo "done fail=$FAIL $(date -Is)"
exit $FAIL
