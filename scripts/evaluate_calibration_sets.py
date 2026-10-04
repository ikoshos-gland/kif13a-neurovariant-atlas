"""Measure variant-effect models on the curated KIF1A calibration sets.

Primary comparison: ClinVar core positives inside the motor domain vs gnomAD motor-domain population
controls (AC >= 3). Secondary: the same models on the ClinVar benign set (motor-free, so it also measures
how much the domain confound inflates the numbers) and on the held-out PP3-dependent positives.

Confidence intervals come from a residue-grouped bootstrap: residues, not variants, are resampled, so
several substitutions of the same residue never land in different folds of the same draw.

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\evaluate_calibration_sets.py
"""

from __future__ import annotations

import asyncio
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, variants  # noqa: E402

core.load_dotenv()
CUR = core.INTERIM / "clinvar_curation"
VEP_REGION = "https://rest.ensembl.org/vep/human/region"
RNG = np.random.default_rng(20260915)
N_BOOT = 2000

# model key in the VEP-derived row -> (label, +1 if higher means more damaging else -1)
MODELS = {
    "am_pathogenicity": ("AlphaMissense", 1), "cadd_phred": ("CADD v1.7", 1), "revel": ("REVEL", 1),
    "esm1b": ("ESM1b (dbNSFP)", -1), "primateai": ("PrimateAI", 1), "metarnn": ("MetaRNN", 1),
}
LABEL_TRAINED = {"revel", "metarnn"}  # trained on clinical labels: report but keep out of the consensus


def auroc(pos: np.ndarray, neg: np.ndarray) -> float | None:
    if pos.size == 0 or neg.size == 0:
        return None
    order = np.argsort(np.concatenate([pos, neg]), kind="mergesort")
    ranks = np.empty(order.size, dtype=float)
    ranks[order] = np.arange(1, order.size + 1)
    # average ranks within ties
    values = np.concatenate([pos, neg])[order]
    i = 0
    while i < values.size:
        j = i
        while j + 1 < values.size and values[j + 1] == values[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return float((ranks[:pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size))


def grouped_bootstrap_ci(pos_rows: list[tuple[str, float]], neg_rows: list[tuple[str, float]]) -> tuple[float, float] | None:
    """pos/neg rows are (group_id, score). Resample groups with replacement."""
    if not pos_rows or not neg_rows:
        return None
    pos_groups: dict[str, list[float]] = defaultdict(list)
    neg_groups: dict[str, list[float]] = defaultdict(list)
    for g, s in pos_rows:
        pos_groups[g].append(s)
    for g, s in neg_rows:
        neg_groups[g].append(s)
    pk, nk = list(pos_groups), list(neg_groups)
    draws = []
    for _ in range(N_BOOT):
        p = np.concatenate([pos_groups[pk[i]] for i in RNG.integers(0, len(pk), len(pk))])
        n = np.concatenate([neg_groups[nk[i]] for i in RNG.integers(0, len(nk), len(nk))])
        a = auroc(np.asarray(p, float), np.asarray(n, float))
        if a is not None:
            draws.append(a)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def sens_at_spec(pos: np.ndarray, neg: np.ndarray, spec: float = 0.90) -> tuple[float, float] | None:
    """Sensitivity at a fixed specificity, plus the score threshold that achieves it."""
    if pos.size == 0 or neg.size == 0:
        return None
    threshold = float(np.quantile(neg, spec))
    return float((pos >= threshold).mean()), threshold


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    def rank(x):
        order = np.argsort(x, kind="mergesort")
        r = np.empty(x.size, float)
        r[order] = np.arange(x.size)
        return r
    ra, rb = rank(a), rank(b)
    return float(np.corrcoef(ra, rb)[0, 1])


async def annotate(keys: list[str]) -> dict[str, dict]:
    lines = sorted((f"{k.split('-')[0]} {k.split('-')[1]} . {k.split('-')[2]} {k.split('-')[3]}" for k in keys),
                   key=variants.vcf_sort_key)
    out: dict[str, dict] = {}
    raw_all = []
    for batch in core.chunks(lines, 200):
        raw = await core.post_json(VEP_REGION, {"variants": [b + " . . ." for b in batch]}, params=variants.VEP_PARAMS,
                                   headers={"Content-Type": "application/json"}, timeout=600)
        raw_all += raw
        for v in raw:
            chrom, pos, _dot, ref, alt = v["input"].split()[:5]
            out[f"{chrom}-{pos}-{ref}-{alt}"] = variants._flatten_vep(v)
        print(f"  VEP {len(out)}/{len(lines)}", flush=True)
    core.save_json("vep", f"KIF1A_calibration_scores_vep_{core.today()}.raw.json", raw_all, compress=True,
                   source="Ensembl VEP REST", urls=[VEP_REGION], tool="scripts/evaluate_calibration_sets.py",
                   params=variants.VEP_PARAMS)
    return out


async def main() -> int:
    sets = list(csv.DictReader(open(sorted(CUR.glob("KIF1A_calibration_sets_*.tsv"))[-1], encoding="utf-8"), delimiter="\t"))
    missense = [r for r in sets if r["set"] != "splice_positive" and r["set"] != "reference_lof" and r["variant_key"]]
    scores = await annotate(sorted({r["variant_key"] for r in missense}))

    rows = []
    for r in missense:
        v = scores.get(r["variant_key"], {})
        row = {**{k: r[k] for k in ("set", "label", "source", "id", "protein_change", "protein_pos", "cv_group",
                                    "domain", "residue_identical_in_kif13a", "functional_literature", "gnomad_af_max")}}
        for key in MODELS:
            raw = v.get(key)
            try:
                row[key] = float(raw)
            except (TypeError, ValueError):
                row[key] = None
        row["vep_consequence"] = v.get("most_severe_consequence")
        rows.append(row)
    prov = dict(source="derived: Ensembl VEP scores on curated KIF1A calibration sets", urls=[VEP_REGION],
                tool="scripts/evaluate_calibration_sets.py")
    print(core.save_tsv("clinvar_curation", f"KIF1A_calibration_scores_{core.today()}.tsv", rows, base=core.INTERIM, **prov))

    def pick(set_name: str, motor_only: bool = False) -> list[dict]:
        out = [r for r in rows if r["set"] == set_name]
        return [r for r in out if r["domain"] in ("motor", "neck_linker")] if motor_only else out

    comparisons = {
        "PRIMARY motor: core positives vs gnomAD AC>=3": (pick("missense_positive_core", True), pick("missense_negative_motor_population")),
        "motor: core positives vs gnomAD AC>=2 (incl. AC2)": (pick("missense_positive_core", True),
                                                              pick("missense_negative_motor_population") + pick("missense_negative_motor_population_ac2")),
        "motor: functional-literature positives vs gnomAD AC>=3": ([r for r in pick("missense_positive_core", True) if r["functional_literature"]],
                                                                   pick("missense_negative_motor_population")),
        "motor: held-out PP3-dependent positives vs gnomAD AC>=3": (pick("heldout_pp3_dependent_positive", True),
                                                                    pick("missense_negative_motor_population")),
        "WHOLE PROTEIN (confounded): core positives vs ClinVar B/LB": (pick("missense_positive_core"),
                                                                       pick("missense_negative_clinvar_outside_motor")),
    }

    report: dict[str, dict] = {}
    for name, (pos_rows, neg_rows) in comparisons.items():
        print(f"\n{name}  n+={len(pos_rows)} n-={len(neg_rows)}")
        report[name] = {"n_pos": len(pos_rows), "n_neg": len(neg_rows), "models": {}}
        for key, (label, sign) in MODELS.items():
            p = [(r["cv_group"], sign * r[key]) for r in pos_rows if r[key] is not None]
            n = [(r["cv_group"], sign * r[key]) for r in neg_rows if r[key] is not None]
            a = auroc(np.array([s for _, s in p]), np.array([s for _, s in n]))
            if a is None:
                continue
            ci = grouped_bootstrap_ci(p, n)
            sens = sens_at_spec(np.array([s for _, s in p]), np.array([s for _, s in n]))
            flag = " (label-trained)" if key in LABEL_TRAINED else ""
            print(f"  {label:<16}{flag:<16} AUROC {a:.3f}  95% CI {ci[0]:.3f}-{ci[1]:.3f}"
                  f"  sens@90%spec {sens[0]:.2f} (esik {sign * sens[1]:.3g})  n+={len(p)} n-={len(n)}")
            report[name]["models"][label] = {"auroc": round(a, 3), "ci95": [round(ci[0], 3), round(ci[1], 3)],
                                             "sens_at_90_spec": round(sens[0], 3), "threshold": round(sign * sens[1], 4),
                                             "n_pos": len(p), "n_neg": len(n), "label_trained": key in LABEL_TRAINED}
        # baselines that use no variant-level model at all
        for base_name, getter in (("paralog residue identical in KIF13A", lambda r: 1.0 if r["residue_identical_in_kif13a"] == "True" else 0.0),
                                  ("is motor domain", lambda r: 1.0 if r["domain"] in ("motor", "neck_linker") else 0.0)):
            a = auroc(np.array([getter(r) for r in pos_rows]), np.array([getter(r) for r in neg_rows]))
            if a is not None:
                print(f"  baseline {base_name:<38} AUROC {a:.3f}")
                report[name]["models"][f"baseline: {base_name}"] = {"auroc": round(a, 3)}

    # model agreement on the primary comparison
    prim = pick("missense_positive_core", True) + pick("missense_negative_motor_population")
    print("\nmodel agreement (Spearman) on the primary comparison")
    keys = [k for k in MODELS if sum(1 for r in prim if r[k] is not None) > 100]
    print("       " + "".join(f"{MODELS[k][0][:9]:>11}" for k in keys))
    for k1 in keys:
        line = f"  {MODELS[k1][0][:9]:<9}"
        for k2 in keys:
            both = [(r[k1] * MODELS[k1][1], r[k2] * MODELS[k2][1]) for r in prim if r[k1] is not None and r[k2] is not None]
            line += f"{spearman(np.array([a for a, _ in both]), np.array([b for _, b in both])):>11.2f}"
        print(line)
    print("\nscore coverage on the primary comparison:",
          {MODELS[k][0]: f"{sum(1 for r in prim if r[k] is not None)}/{len(prim)}" for k in MODELS})
    core.save_json("clinvar_curation", f"KIF1A_calibration_metrics_{core.today()}.json", report, base=core.INTERIM, **prov)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
