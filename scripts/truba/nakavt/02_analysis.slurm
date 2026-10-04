#!/bin/bash
#SBATCH --job-name=k13_nakavt_calc
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/avural/kif13a_nakavt/logs/calc_%j.out
#SBATCH --error=/arf/scratch/avural/kif13a_nakavt/logs/calc_%j.err

# C ekseni (nakavt), adim 2: sinonim dahil gnomAD cekimi (internet) + butun hesaplar.
# Girdi: adim 1'in raw/ ve interim/sheets/ ciktilari, interim/project/ altina kopyalanmis proje tablolari.
set -euo pipefail
BASE=/arf/scratch/avural/kif13a_nakavt
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
cd "$BASE"
echo "Job $SLURM_JOB_ID on $(hostname), start $(date)"
python -u scripts/fetch_sources.py gnomad_all KIF13A:ENST00000259711 KIF1A:ENST00000498729 KIF13B:ENST00000524189
python -u scripts/nakavt_analysis.py
echo "End $(date)"
