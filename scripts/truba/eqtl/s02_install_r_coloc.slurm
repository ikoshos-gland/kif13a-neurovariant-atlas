#!/bin/bash
#SBATCH -J eqtl_r_coloc
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=01:30:00
#SBATCH -o /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.out
#SBATCH -e /arf/scratch/avural/kif13a_eqtl/logs/%x_%j.err
# Installs the reference R implementation of coloc (Wallace lab, CRAN) into a project-local
# library on scratch. Used ONLY to validate the Python coloc.abf / SuSiE re-implementation on
# identical inputs. Needs internet -> debug partition.
set -uo pipefail
source /etc/profile.d/modules.sh
module load apps/R/4.5.2-gcc-14.1.0
export R_LIBS_USER=/arf/scratch/avural/kif13a_eqtl/Rlib
mkdir -p $R_LIBS_USER
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
for (p in c("coloc")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib)
library(coloc); cat("coloc version", as.character(packageVersion("coloc")), "\n")
cat("susieR version", as.character(packageVersion("susieR")), "\n")
'
echo "exit=$?"
du -sh $R_LIBS_USER; find $R_LIBS_USER | wc -l
