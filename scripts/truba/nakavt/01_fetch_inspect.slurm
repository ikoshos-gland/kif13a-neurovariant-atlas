#!/bin/bash
#SBATCH --job-name=k13_nakavt_fetch
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_nakavt/logs/fetch_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_nakavt/logs/fetch_%j.err

# C ekseni (nakavt), adim 1: indirme (debug = internet) + duz metne cevirme.
set -euo pipefail
BASE=/arf/scratch/avural/kif13a_nakavt
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd "$BASE"
echo "Job $SLURM_JOB_ID on $(hostname), start $(date)"

# Excel ve PDF okumak icin iki kucuk saf-Python paketi, yalniz bu klasore (conda agacina degil)
python - <<'EOF' || pip install --quiet --no-deps --target "$BASE/pylib" openpyxl et-xmlfile pypdf
import sys; sys.path.insert(0, "/arf/scratch/avural/kif13a_nakavt/pylib")
import openpyxl, pypdf; print("openpyxl", openpyxl.__version__, "pypdf", pypdf.__version__)
EOF

python -u scripts/fetch_sources.py files
python -u scripts/fetch_sources.py ensembl ENST00000259711 ENST00000498729 ENST00000524189
python -u scripts/fetch_sources.py gnomad KIF13A:ENST00000259711 KIF1A:ENST00000498729 KIF13B:ENST00000524189
python -u scripts/inspect_sources.py
echo "End $(date)"
