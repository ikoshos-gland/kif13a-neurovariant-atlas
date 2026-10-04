"""Leak-free consensus check: mean of within-comparison score ranks, no fitted weights.

Uses the scores already saved by evaluate_calibration_sets.py. Label-trained models (REVEL, MetaRNN)
are excluded from the consensus and reported separately for contrast.

    mcp\\atlas_server\\.venv\\Scripts\\python.exe scripts\\consensus_rank_mean.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from atlas_mcp import core  # noqa: E402
from evaluate_calibration_sets import MODELS, auroc, grouped_bootstrap_ci, sens_at_spec  # noqa: E402

CONSENSUS = ["am_pathogenicity", "esm1b", "primateai", "cadd_phred"]


def main() -> int:
    path = sorted((core.INTERIM / "clinvar_curation").glob("KIF1A_calibration_scores_*.tsv"))[-1]
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))
    for r in rows:
        for k in MODELS:
            r[k] = float(r[k]) if r[k] not in ("", None) else None

    pos = [r for r in rows if r["set"] == "missense_positive_core" and r["domain"] in ("motor", "neck_linker")]
    neg = [r for r in rows if r["set"] == "missense_negative_motor_population"]
    pool = pos + neg
    usable = [k for k in CONSENSUS if all(r[k] is not None for r in pool)]
    print(f"consensus over {usable}  n+={len(pos)} n-={len(neg)}")

    # rank within the pooled comparison so each model contributes on the same scale
    ranks = {}
    for k in usable:
        sign = MODELS[k][1]
        values = np.array([sign * r[k] for r in pool], float)
        order = np.argsort(values, kind="mergesort")
        rk = np.empty(values.size, float)
        rk[order] = np.arange(values.size)
        ranks[k] = rk / (values.size - 1)
    mean_rank = np.mean([ranks[k] for k in usable], axis=0)
    for r, m in zip(pool, mean_rank):
        r["consensus"] = float(m)

    p = [(r["cv_group"], r["consensus"]) for r in pos]
    n = [(r["cv_group"], r["consensus"]) for r in neg]
    a = auroc(np.array([s for _, s in p]), np.array([s for _, s in n]))
    ci = grouped_bootstrap_ci(p, n)
    sens, thr = sens_at_spec(np.array([s for _, s in p]), np.array([s for _, s in n]))
    print(f"  rank-mean consensus  AUROC {a:.3f}  95% CI {ci[0]:.3f}-{ci[1]:.3f}  sens@90%spec {sens:.2f}")
    best = max(usable, key=lambda k: auroc(np.array([MODELS[k][1] * r[k] for r in pos]),
                                           np.array([MODELS[k][1] * r[k] for r in neg])))
    print(f"  best single non-label-trained model: {MODELS[best][0]} "
          f"{auroc(np.array([MODELS[best][1] * r[best] for r in pos]), np.array([MODELS[best][1] * r[best] for r in neg])):.3f}")
    worst = sorted(pos, key=lambda r: r["consensus"])[:6]
    print("  hardest positives (lowest consensus):", [(r["protein_change"], round(r["consensus"], 2)) for r in worst])
    top_neg = sorted(neg, key=lambda r: -r["consensus"])[:6]
    print("  hardest negatives (highest consensus):", [(r["protein_change"], round(r["consensus"], 2), r["gnomad_af_max"][:7]) for r in top_neg])
    core.save_tsv("clinvar_curation", f"KIF1A_consensus_rank_mean_{core.today()}.tsv",
                  [{k: r.get(k) for k in ("set", "protein_change", "cv_group", "label", "consensus", *usable)} for r in pool],
                  base=core.INTERIM, source="derived: rank-mean consensus of non-label-trained models",
                  urls=[], tool="scripts/consensus_rank_mean.py", params={"models": usable})
    return 0


if __name__ == "__main__":
    sys.exit(main())
