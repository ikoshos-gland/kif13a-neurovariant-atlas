"""Compare the Python coloc.abf with CRAN coloc on identical inputs (written by s30_validate.R)."""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coloc_py as cp  # noqa: E402

RES = "/arf/scratch/avural/kif13a_eqtl/results"
RV = f"{RES}/r_validation"
rows = []

# 1) coloc_test_data D1 x D2 (both quant; sdY given in the R object, printed in the R log)
d1 = pd.read_csv(f"{RV}/example_D1.tsv", sep="\t")
d2 = pd.read_csv(f"{RV}/example_D2.tsv", sep="\t")
R = pd.read_csv(f"{RV}/example_D1D2_R_summary.tsv", sep="\t").iloc[0]
sdY = float(os.environ.get("EX_SDY", "1"))
l1 = cp.labf_beta(d1.beta, d1.varbeta, sdY=sdY)
l2 = cp.labf_beta(d2.beta, d2.varbeta, sdY=sdY)
r = cp.coloc_abf(l1, l2)
for h in range(5):
    rows.append({"case": "coloc_test_data D1xD2 (beta)", "quantity": f"PP.H{h}",
                 "R_coloc": float(R[f"PP.H{h}.abf"]), "python": float(r[f"PP.H{h}"])})
Rp = pd.read_csv(f"{RV}/example_D1D2_pvalue_R_summary.tsv", sep="\t").iloc[0]
N = float(os.environ.get("EX_N", "500"))
from scipy import stats
p1 = 2 * stats.norm.sf(np.abs(d1.beta / np.sqrt(d1.varbeta)))
p2 = 2 * stats.norm.sf(np.abs(d2.beta / np.sqrt(d2.varbeta)))
rp = cp.coloc_abf(cp.labf_p(p1, d1.MAF, N, sdY=sdY), cp.labf_p(p2, d2.MAF, N, sdY=sdY))
for h in range(5):
    rows.append({"case": "coloc_test_data D1xD2 (p-values+MAF+N)", "quantity": f"PP.H{h}",
                 "R_coloc": float(Rp[f"PP.H{h}.abf"]), "python": float(rp[f"PP.H{h}"])})

# 2) real KIF13A inputs
rr = pd.read_csv(f"{RV}/real_coloc_abf_R.tsv", sep="\t")
for _, x in rr.iterrows():
    d = pd.read_csv(f"{RES}/r_input/{x['input']}.tsv", sep="\t")
    if d.beta.isna().all():
        lg = cp.labf_p(d.p, d.maf_ref, np.median(d.n), sdY=1)
    else:
        lg = cp.labf_beta(d.beta, d.se ** 2, sdY=1)
    le = cp.labf_beta(d.beta_e, d.se_e ** 2, sdY=1)
    r = cp.coloc_abf(lg, le)
    for h in (3, 4):
        rows.append({"case": x["input"], "quantity": f"PP.H{h}", "R_coloc": float(x[f"PP.H{h}.abf"]),
                     "python": float(r[f"PP.H{h}"])})
out = pd.DataFrame(rows)
out["abs_diff"] = (out.R_coloc - out.python).abs()
out.to_csv(f"{RES}/validation_python_vs_R.tsv", sep="\t", index=False)
print(out.to_string())
print("max abs diff:", out.abs_diff.max())
