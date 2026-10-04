#!/bin/bash
#SBATCH -J eqtl_dl_open
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=03:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Open-access downloads (debug partition is the only one with internet).
# Bryois 2022 Zenodo (chr6 + chr1 for the positive control), GWAS Catalog sumstats,
# BIG40 33k sumstats for three IDPs, 1000 Genomes phase 3 chr6 (GRCh37) for LD.
set -uo pipefail
W=/arf/scratch/avural/kif13a_eqtl
R=$W/raw
mkdir -p $R/bryois $R/gwas $R/big40 $R/kg
cd $R
UA="kif13a-atlas/1.0 (academic; TRUBA)"
dl() {  # url dest
  local u="$1" d="$2"
  if [ -s "$d" ]; then echo "have $d"; return 0; fi
  curl -sS -L --retry 5 --retry-delay 10 --connect-timeout 60 -A "$UA" -o "$d.part" "$u" && mv "$d.part" "$d"
  echo "$(date +%T) $d $(stat -c %s "$d" 2>/dev/null) bytes  <- $u"
}

echo "== Zenodo record (Bryois 2022)"
dl "https://zenodo.org/api/records/7276971" bryois/zenodo_7276971.json
dl "https://zenodo.org/api/records/5543734" bryois/zenodo_5543734.json
ZB="https://zenodo.org/records/7276971/files"
for ct in Astrocytes Endothelial.cells Excitatory.neurons Inhibitory.neurons Microglia Oligodendrocytes "OPCs...COPs" Pericytes pb; do
  dl "$ZB/${ct}.6.gz?download=1" "bryois/${ct}.6.gz"
done
for ct in Oligodendrocytes "OPCs...COPs"; do
  dl "$ZB/${ct}.1.gz?download=1" "bryois/${ct}.1.gz"
done
dl "$ZB/snp_pos.txt.gz?download=1" bryois/snp_pos.txt.gz

echo "== GWAS Catalog"
G=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
dl $G/GCST90860001-GCST90861000/GCST90860937/GCST90860937.tsv.gz gwas/GCST90860937.tsv.gz
dl $G/GCST90860001-GCST90861000/GCST90860937/GCST90860937.tsv.gz-meta.yaml gwas/GCST90860937.tsv.gz-meta.yaml
dl $G/GCST90860001-GCST90861000/GCST90860937/md5sum.txt gwas/GCST90860937.md5sum.txt
dl $G/GCST90095001-GCST90096000/GCST90095129/GCST90095129_buildGRCh37.tsv.gz gwas/GCST90095129_buildGRCh37.tsv.gz
dl $G/GCST90095001-GCST90096000/GCST90095129/GCST90095129_buildGRCh37.tsv.gz-meta.yaml gwas/GCST90095129_buildGRCh37.tsv.gz-meta.yaml
dl $G/GCST90095001-GCST90096000/GCST90095129/md5sum.txt gwas/GCST90095129.md5sum.txt
dl $G/GCST90012001-GCST90013000/GCST90012877/GCST90012877_buildGRCh37.tsv.gz gwas/GCST90012877_buildGRCh37.tsv.gz
dl $G/GCST90012001-GCST90013000/GCST90012877/GCST90012877_buildGRCh37.tsv.gz-meta.yaml gwas/GCST90012877_buildGRCh37.tsv.gz-meta.yaml
dl $G/GCST90012001-GCST90013000/GCST90012877/md5sum.txt gwas/GCST90012877.md5sum.txt
dl $G/GCST90012001-GCST90013000/GCST90012877/README.txt gwas/GCST90012877.README.txt
(cd gwas && for f in *.md5sum.txt; do cat "$f"; echo; done > all_md5.txt; md5sum GCST90860937.tsv.gz GCST90095129_buildGRCh37.tsv.gz GCST90012877_buildGRCh37.tsv.gz > computed_md5.txt)

echo "== BIG40 33k (GRCh37): L3 / ICVF / FA left anterior corona radiata"
for p in 1769 1925 1475; do
  dl "https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k/${p}.txt.gz" big40/${p}.txt.gz
done

echo "== 1000 Genomes phase 3 chr6 (GRCh37)"
K=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
dl $K/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz kg/chr6.vcf.gz
dl $K/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi kg/chr6.vcf.gz.tbi
dl $K/integrated_call_samples_v3.20130502.ALL.panel kg/panel.txt

echo "== done $(date)"
ls -la $R/*
