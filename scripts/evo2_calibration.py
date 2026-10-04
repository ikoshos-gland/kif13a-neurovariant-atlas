"""Put Evo 2 40B through the same KIF1A exam as the other models.

Scores the primary comparison (ClinVar core positives in the motor domain vs gnomAD motor controls with
AC >= 3) in the gene-sense orientation, which is the reverse complement of the plus strand for KIF1A.
Score = summed change in per-base log-probability from the variant base onward (the hosted NIM shifts
upstream bases slightly, so those are excluded; see docs/MCP_KURULUM.md).

Cost: 2 NIM calls per variant. Partial results are written after every variant, so a rerun resumes.

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\evo2_calibration.py
"""

from __future__ import annotations

import asyncio
import csv
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from atlas_mcp import core, genome_models as g2  # noqa: E402
from evaluate_calibration_sets import auroc, grouped_bootstrap_ci, sens_at_spec  # noqa: E402

core.load_dotenv()
W, HALF = 8192, 4096
COMP = str.maketrans("ACGTN", "TGCAN")
REGION = "https://rest.ensembl.org/sequence/region/human"
OUT = core.INTERIM / "evo2" / f"KIF1A_calibration_evo2_w{W}.tsv"
SETS = {"missense_positive_core": 1, "missense_negative_motor_population": 0}


_region: tuple[str, int, str] | None = None  # (chrom, start, sequence): one fetch covers every window


async def load_region(chrom: str, lo: int, hi: int) -> None:
    """Fetch one span covering all variant windows, so each variant costs no Ensembl call."""
    global _region
    start, end = lo - HALF, hi + HALF
    url = f"{REGION}/{chrom}:{start}..{end}:1"
    seq = (await core.request("GET", url, headers={"Accept": "text/plain"}, timeout=300)).text.strip().upper()
    _region = (chrom, start, seq)
    print(f"region {chrom}:{start}-{end} fetched, {len(seq)} bp", flush=True)


async def score_variant(chrom: str, pos: int, ref: str, alt: str) -> dict:
    start = pos - HALF + 1
    if _region and _region[0] == chrom and _region[1] <= start and start + W - 1 <= _region[1] + len(_region[2]) - 1:
        offset = start - _region[1]
        plus = _region[2][offset:offset + W]
    else:
        url = f"{REGION}/{chrom}:{start}..{start + W - 1}:1"
        plus = (await core.request("GET", url, headers={"Accept": "text/plain"})).text.strip().upper()
    if plus[HALF - 1] != ref:
        raise ValueError(f"reference base at {chrom}:{pos} is {plus[HALF - 1]}, not {ref}")
    # gene-sense orientation (KIF1A is on the minus strand)
    seq = plus[::-1].translate(COMP)
    idx = W - HALF
    a = alt.translate(COMP)
    alt_seq = seq[:idx] + a + seq[idx + 1:]
    lp_ref = g2.token_logprobs(seq, await g2.evo2_logits(seq))
    lp_alt = g2.token_logprobs(alt_seq, await g2.evo2_logits(alt_seq))
    d = lp_alt - lp_ref
    return {"evo2_delta": round(float(d[idx:].sum()), 4), "evo2_site_term": round(float(d[idx]), 4),
            "evo2_downstream": round(float(d[idx + 1:].sum()), 4),
            "evo2_upstream_noise": round(float(d[:idx].sum()), 4)}


async def main() -> int:
    sets_file = sorted((core.INTERIM / "clinvar_curation").glob("KIF1A_calibration_sets_*.tsv"))[-1]
    rows = [r for r in csv.DictReader(open(sets_file, encoding="utf-8"), delimiter="\t")
            if r["set"] in SETS and (r["set"] != "missense_positive_core" or r["domain"] in ("motor", "neck_linker"))]
    done: dict[str, dict] = {}
    if OUT.exists():
        done = {r["variant_key"]: r for r in csv.DictReader(open(OUT, encoding="utf-8"), delimiter="\t")}
        print(f"resuming: {len(done)} variants already scored")

    results = list(done.values())
    todo = [r for r in rows if r["variant_key"] not in done]
    positions = [int(r["variant_key"].split("-")[1]) for r in todo]
    if positions:
        await load_region(todo[0]["variant_key"].split("-")[0], min(positions), max(positions))
    t0, n_start = time.monotonic(), len(results)
    limit = asyncio.Semaphore(4)  # a few NIM calls in flight; the hosted endpoint handles this fine

    async def one(r: dict) -> dict | None:
        chrom, pos, ref, alt = r["variant_key"].split("-")
        async with limit:
            try:
                s = await score_variant(chrom, int(pos), ref, alt)
            except Exception as exc:
                print(f"  FAIL {r['protein_change']}: {type(exc).__name__}: {exc}"[:200], flush=True)
                return None
        return {"variant_key": r["variant_key"], "set": r["set"], "label": SETS[r["set"]],
                "protein_change": r["protein_change"], "protein_pos": r["protein_pos"], "cv_group": r["cv_group"],
                "am_pathogenicity": r["am_pathogenicity"], **s}

    for chunk in core.chunks(todo, 8):
        results += [x for x in await asyncio.gather(*(one(r) for r in chunk)) if x]
        OUT.parent.mkdir(parents=True, exist_ok=True)
        with open(OUT, "w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
        rate = (time.monotonic() - t0) / max(1, len(results) - n_start)
        print(f"  {len(results)}/{len(rows)} scored, {rate:.1f}s per variant, "
              f"{(len(rows) - len(results)) * rate / 60:.0f} min left", flush=True)

    core.save_tsv("evo2", OUT.name, results, base=core.INTERIM, source="NVIDIA NIM arc/evo2-40b forward (gene-sense)",
                  urls=[g2.NIM_FORWARD, REGION], tool="scripts/evo2_calibration.py", params={"window": W})

    pos_rows = [r for r in results if int(r["label"]) == 1]
    neg_rows = [r for r in results if int(r["label"]) == 0]
    print(f"\nEvo 2 40B, motor domain: n+={len(pos_rows)} n-={len(neg_rows)}")
    for name, key, sign in (("Evo 2 delta (site+downstream)", "evo2_delta", -1),
                            ("Evo 2 site term only", "evo2_site_term", -1),
                            ("AlphaMissense (same variants)", "am_pathogenicity", 1)):
        p = [(r["cv_group"], sign * float(r[key])) for r in pos_rows if r[key] not in ("", None)]
        n = [(r["cv_group"], sign * float(r[key])) for r in neg_rows if r[key] not in ("", None)]
        a = auroc(np.array([s for _, s in p]), np.array([s for _, s in n]))
        ci = grouped_bootstrap_ci(p, n)
        sens, thr = sens_at_spec(np.array([s for _, s in p]), np.array([s for _, s in n]))
        print(f"  {name:<32} AUROC {a:.3f}  95% CI {ci[0]:.3f}-{ci[1]:.3f}  sens@90%spec {sens:.2f} (esik {sign * thr:.3g})")

    def rank(x):
        order = np.argsort(x, kind="mergesort")
        out = np.empty(x.size)
        out[order] = np.arange(x.size)
        return out
    both = [(-float(r["evo2_delta"]), float(r["am_pathogenicity"])) for r in results if r["am_pathogenicity"] not in ("", None)]
    e, am = np.array([x for x, _ in both]), np.array([y for _, y in both])
    print(f"  Evo 2 ve AlphaMissense uyumu (Spearman): {np.corrcoef(rank(e), rank(am))[0, 1]:.2f}")
    worst = sorted(pos_rows, key=lambda r: float(r["evo2_delta"]))[:5]
    best_neg = sorted(neg_rows, key=lambda r: float(r["evo2_delta"]))[:5]
    print("  Evo 2'ye gore en bozucu pozitifler:", [(r["protein_change"], r["evo2_delta"]) for r in worst])
    print("  Evo 2'nin en bozucu bulduğu negatifler:", [(r["protein_change"], r["evo2_delta"]) for r in best_neg])
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
