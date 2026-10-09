#!/bin/bash
# Build the AmberTools 26 topology image for the Rab10:KIF13A complex arm on a debug node with internet
# (same pattern as slurm/build_image.slurm, which built mdana-openmm831-pyscf.sif). WRITTEN 2026-10-05, NOT SUBMITTED.
#   sbatch /arf/scratch/avural/so3lr/complex_md/slurm/build_ambertools.slurm
# Needs complex_md/container/{ambertools.def, ambertools-lock-linux-64.txt} on scratch.
# Only needed to rebuild topologies on TRUBA; the 12 systems were built locally with the same lock.
#SBATCH --account=avural
#SBATCH -p debug
#SBATCH -J build_ambertools
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o /arf/scratch/avural/so3lr/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/so3lr/logs/%x_%j.err

set -uo pipefail
ROOT=/arf/scratch/avural/so3lr/complex_md
NAME=ambertools
SIF=/arf/home/avural/containers/ambertools26-openmm831.sif
export APPTAINER_TMPDIR=/tmp/${USER}_aptmp_${SLURM_JOB_ID}
export APPTAINER_CACHEDIR=/tmp/${USER}_apcache_${SLURM_JOB_ID}
mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR" "$(dirname "$SIF")"
trap 'rm -rf "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"' EXIT

echo "host=$(hostname) job=$SLURM_JOB_ID date=$(date -Is)"
curl -sI https://conda.anaconda.org | head -1 || { echo "FATAL: no internet"; exit 1; }
cd "$ROOT/container" || exit 1
test -f "$NAME.def" && test -f ambertools-lock-linux-64.txt || { echo "missing $NAME.def or lock file"; exit 1; }
apptainer build --force --fakeroot "$SIF.partial" "$NAME.def"
RC=$?
echo "BUILD_RC=$RC"
if [ $RC -eq 0 ]; then
  mv -f "$SIF.partial" "$SIF"
  ls -la "$SIF"
  sha256sum "$SIF" > "$SIF.sha256"
  apptainer exec "$SIF" sh -c 'cat /opt/*-freeze.txt' > "$ROOT/container/${NAME}_freeze_${SLURM_JOB_ID}.txt"
  apptainer test "$SIF"
  RC=$?
else
  rm -f "$SIF.partial"
fi
echo "DONE rc=$RC $(date -Is)"
exit $RC
