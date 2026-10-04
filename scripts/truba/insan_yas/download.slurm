#!/bin/bash
#SBATCH --job-name=k13_dl
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_insan_yas/logs/dl_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_insan_yas/logs/dl_%j.err
# Downloads CELLxGENE Discover h5ad files (public, CC BY 4.0) listed in datasets.tsv.
# Internet exists only on debug nodes, so downloads run here, never on the login node.
set -uo pipefail
BASE=/arf/scratch/avural/kif13a_insan_yas
mkdir -p "$BASE/raw"
cd "$BASE/raw"
echo "Job $SLURM_JOB_ID on $(hostname), start $(date)"
KEYS="${KEYS:-all}"; KEYS="${KEYS//:/,}"  # sbatch --export splits on commas, so pass keys separated by colons
tail -n +2 "$BASE/scripts/datasets.tsv" | while IFS=$'\t' read -r key did dvid title size; do
  if [ "$KEYS" != "all" ] && ! echo ",$KEYS," | grep -q ",$key,"; then continue; fi
  f="$key.h5ad"
  url="https://datasets.cellxgene.cziscience.com/$dvid.h5ad"
  echo "== $key  $url  (expected ~${size} GB)"
  for try in 1 2 3; do
    curl -sS -L -C - --retry 5 --retry-delay 10 -o "$f" "$url" && break
    echo "retry $try for $key"
  done
  remote=$(curl -sS -I -L "$url" | grep -i '^content-length' | tail -1 | tr -dc '0-9')
  local_sz=$(stat -c %s "$f")
  echo "$key remote=$remote local=$local_sz $( [ "$remote" = "$local_sz" ] && echo OK || echo MISMATCH )"
  echo -e "$key\t$dvid\t$url\t$local_sz\t$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$BASE/raw/DOWNLOADS.tsv"
done
df -h "$BASE" | tail -1
du -sh "$BASE/raw"
echo "End $(date)"
