#!/bin/bash
# Rab10(GTP, Mg):KIF13A complex MD, one run per GPU (OpenMM 8.3.1 CUDA, mdana image), via bin/multigpu_run.sh.
# WRITTEN 2026-10-05, NOT SUBMITTED. Task files come from scripts/make_complex_tasks.py (complex_md/tasks/).
#
#   sbatch /arf/scratch/avural/so3lr/complex_md/slurm/complex_md_multigpu.slurm TASKFILE
#   sbatch -p barbun-cuda -c 40 /arf/scratch/avural/so3lr/complex_md/slurm/complex_md_multigpu.slurm TASKFILE   # 2x P100
#
# Resubmittable: a finished run has DONE (skipped); an unfinished one resumes from its last 1 ns checkpoint
# (run_complex.py --auto). The runs pause cleanly (exit 3, no DONE) about MARGIN_S seconds before the job end
# taken from squeue; after a hard kill the DCD/state tail beyond the checkpoint is cut back on resume.
# Size (estimate, not measured on TRUBA): 69-84k particles, about 140-175 ns/day per V100 scaled from the
# measured 207 ns/day at 58k and 127 ns/day at 95k -> about 14-18 h per 100 ns run; two runs per 2-GPU job.
#SBATCH --account=avural
#SBATCH -p akya-cuda
#SBATCH -J k13_cx_md
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=20
#SBATCH --mem=48G
#SBATCH --time=3-00:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/so3lr/logs/%x_%j.err

set -uo pipefail
TASKS=${1:?usage: complex_md_multigpu.slurm TASKFILE}
ROOT=/arf/scratch/avural/so3lr
SIF=/arf/home/avural/containers/mdana-openmm831-pyscf.sif
test -f "$SIF" || { echo "missing image $SIF"; exit 1; }
test -f "$TASKS" || { echo "missing task file $TASKS"; exit 1; }
[[ -f "$(dirname "$TASKS")/STOP" ]] && { echo "STOP file present, nothing started"; exit 0; }

END=$(squeue -h -j "$SLURM_JOB_ID" -o %e 2>/dev/null)
END_EPOCH=$(date -d "$END" +%s 2>/dev/null || echo 0)
if (( END_EPOCH > 0 )); then
  export COMPLEX_DEADLINE_EPOCH=$(( END_EPOCH - ${MARGIN_S:-1200} ))
fi
HOMEDIR=/tmp/${USER}_home_${SLURM_JOB_ID}
mkdir -p "$HOMEDIR"
trap 'rm -rf "$HOMEDIR"' EXIT
echo "host=$(hostname) job=$SLURM_JOB_ID gpus=${CUDA_VISIBLE_DEVICES:-} end=$END deadline_epoch=${COMPLEX_DEADLINE_EPOCH:-none} $(date -Is)"
nvidia-smi --query-gpu=index,name,driver_version,memory.total --format=csv,noheader
GPULOG=${TASKS%.txt}_gpu_${SLURM_JOB_ID}.csv
nvidia-smi --query-gpu=timestamp,index,utilization.gpu,memory.used --format=csv -l 300 > "$GPULOG" &
SMI=$!
bash "$ROOT/bin/multigpu_run.sh" "$TASKS" 1
RC=$?
kill "$SMI" 2>/dev/null
echo "=== status $(date -Is)"
grep -v '^#' "$TASKS" | while read -r DIR _; do
  [[ -z "$DIR" ]] && continue
  st=running; [[ -f "$DIR/DONE" ]] && st=DONE
  last=$(tail -n 1 "$DIR/run_complex.log" 2>/dev/null | cut -c1-200)
  echo "$(basename "$DIR") $st exit=$(cat "$DIR/EXIT" 2>/dev/null) :: $last"
done
echo "multigpu rc=$RC (non-zero also when runs paused before the deadline; resubmit the same task file)"
exit $RC
