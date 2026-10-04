#!/bin/bash
#SBATCH -J eqtl_dl_r2
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=12G
#SBATCH --time=03:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Round 2 downloads (debug = internet):
#  * 1000 Genomes chr6 VCF: resume the interrupted transfer (curl -C -), then verify gzip integrity
#  * Zeng 2026 EN / IN class parquet: resume with HTTP Range (syn_common.download), md5-checked
#  * SingleBrain narrow extraction by HTTP Range bisection (s06_singlebrain_stream.py)
set -uo pipefail
source /arf/home/avural/miniconda3/etc/profile.d/conda.sh
conda activate allenbrain
W=/arf/scratch/avural/kif13a_eqtl
cd $W/scripts
(
  cd $W/raw/kg
  K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
  for i in 1 2 3 4 5 6; do
    curl -sS -L -C - --connect-timeout 60 -o chr6.vcf.gz.part $K && break
    echo "kg curl attempt $i failed, resuming"; sleep 20
  done
  s=$(stat -c %s chr6.vcf.gz.part)
  echo "kg size $s (expected 959144250)"
  if [ "$s" = "959144250" ]; then mv chr6.vcf.gz.part chr6.vcf.gz; fi
) > $W/logs/kg_resume_${SLURM_JOB_ID}.log 2>&1 &
python -u s03_synapse_download.py meta-EN-class.parquet,meta-IN-class.parquet 2 > $W/logs/zeng_resume_${SLURM_JOB_ID}.log 2>&1 &
for ct in OD OPC Ast Ext IN MG End; do
  python -u s06_singlebrain_stream.py $ct > $W/logs/sbrain_${ct}_${SLURM_JOB_ID}.log 2>&1 &
  sleep 2
done
wait
echo "== all done"; tail -3 $W/logs/kg_resume_${SLURM_JOB_ID}.log; tail -3 $W/logs/zeng_resume_${SLURM_JOB_ID}.log
for ct in OD OPC Ast Ext IN MG End; do tail -2 $W/logs/sbrain_${ct}_${SLURM_JOB_ID}.log; done
