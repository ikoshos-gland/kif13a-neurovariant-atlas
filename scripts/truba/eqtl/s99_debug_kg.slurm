#!/bin/bash
#SBATCH -J eqtl_dbg
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=00:10:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd /arf/scratch/avural/kif13a_eqtl
python - <<'PY'
import numpy as np, gzip
z=np.load("work/kg_eur_chr6.npz")
print(z["rsid"][:5], z["rsid"].dtype, z["pos"][:5], z["pos"].min(), z["pos"].max(), z["af"][:5])
print(z["G"].shape, z["G"][:3,:10])
with gzip.open("raw/kg/chr6.vcf.gz","rt") as fh:
    n=0
    for line in fh:
        if line.startswith("#"): continue
        print(line[:120]); n+=1
        if n>2: break
PY
