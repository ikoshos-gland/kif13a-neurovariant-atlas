"""Rank the KIF13A missense variants seen in gnomAD, using the thresholds calibrated on KIF1A.

Evidence per variant:
  - AlphaMissense, ESM1b, CADD, PrimateAI from one Ensembl VEP run (same source as the calibration)
  - the gene-calibrated AlphaMissense cutoff 0.968 (90 percent specificity against gnomAD motor controls
    in KIF1A; the tool's own generic cutoff is 0.564)
  - paralog evidence: the equivalent KIF1A residue, whether it is identical, whether KIF1A carries a
    curated pathogenic missense there, and whether gnomAD tolerates missense there in KIF1A
  - gnomAD frequency (a dominant de novo disorder allele should be absent or a singleton)

Tiers (transparent rules, no fitted weights):
  A  AM >= cutoff, identical paralog residue, KIF1A pathogenic variant at that residue
  B  AM >= cutoff, identical paralog residue in the motor domain, no KIF1A disease variant there
  C  AM >= cutoff, elsewhere (tail, non-identical residue or unmapped)
  D  AM <  cutoff but KIF1A pathogenic at an identical residue (model and paralog disagree)
  -  everything else

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\build_kif13a_candidates.py
"""

from __future__ import annotations

import asyncio
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, variants  # noqa: E402

core.load_dotenv()
VEP_REGION = "https://rest.ensembl.org/vep/human/region"
AM_CUTOFF = 0.968          # calibrated on KIF1A, 90 percent specificity
AM_GENERIC = 0.564         # AlphaMissense default likely-pathogenic cutoff, kept for comparison
AF_DOMINANT_MAX = 1e-4     # above this, a fully penetrant dominant allele is implausible

# KIF13A domains, MANE numbering = UniProt Q9H1H9 (verified 2026-09-15)
DOMAINS = [(1, 4, "N_term"), (5, 352, "motor"), (353, 380, "neck_linker"), (381, 468, "CC1_region"),
           (469, 519, "FHA"), (520, 1805, "tail")]


def domain(pos: int) -> str:
    return next((name for a, b, name in DOMAINS if a <= pos <= b), "tail")


def latest(pattern: str, base: Path = core.RAW) -> Path:
    return sorted(base.glob(pattern))[-1]


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


async def vep_scores(keys: list[str]) -> dict[str, dict]:
    cache = core.INTERIM / "vep" / f"KIF13A_missense_scores_{core.today()}.tsv"
    if cache.exists():
        print(f"using cached {core.rel(cache)}")
        return {r["variant_key"]: r for r in table(cache)}
    lines = sorted((f"{k.split('-')[0]} {k.split('-')[1]} . {k.split('-')[2]} {k.split('-')[3]} . . ." for k in keys),
                   key=variants.vcf_sort_key)
    rows: dict[str, dict] = {}
    for i, batch in enumerate(core.chunks(lines, 200)):
        raw = await core.post_json(VEP_REGION, {"variants": batch}, params=variants.VEP_PARAMS,
                                   headers={"Content-Type": "application/json"}, timeout=900)
        for v in raw:
            chrom, pos, _dot, ref, alt = v["input"].split()[:5]
            f = variants._flatten_vep(v)
            rows[f"{chrom}-{pos}-{ref}-{alt}"] = {
                "variant_key": f"{chrom}-{pos}-{ref}-{alt}", "hgvsp": f.get("hgvsp"), "protein_pos": f.get("protein_pos"),
                "aa_ref": f.get("aa_ref"), "aa_alt": f.get("aa_alt"), "am_pathogenicity": f.get("am_pathogenicity"),
                "am_class": f.get("am_class"), "cadd_phred": f.get("cadd_phred"), "esm1b": f.get("esm1b"),
                "primateai": f.get("primateai"), "revel": f.get("revel"), "consequence": f.get("most_severe_consequence"),
            }
        print(f"  VEP batch {i + 1}: {len(rows)}/{len(lines)}", flush=True)
    core.save_tsv("vep", cache.name, list(rows.values()), base=core.INTERIM, source="Ensembl VEP REST",
                  urls=[VEP_REGION], tool="scripts/build_kif13a_candidates.py", params=variants.VEP_PARAMS)
    return rows


async def main() -> int:
    gn = [r for r in table(latest("gnomad/KIF13A_ENST00000259711_gnomad_r4_all_*.tsv")) if r["consequence"] == "missense_variant"]
    print(f"KIF13A gnomAD missense variants: {len(gn)}")
    scores = await vep_scores([r["variant_id"] for r in gn])

    # paralog side: invert the KIF1A -> KIF13A map, then mark KIF1A evidence per residue
    inv: dict[int, dict] = {}
    for r in table(core.INTERIM / "paralog_maps" / "KIF1A_to_KIF13A_pairwise_MANE.tsv"):
        inv[int(r["KIF13A_pos"])] = {"kif1a_pos": int(r["KIF1A_pos"]), "kif1a_aa": r["KIF1A_aa"], "kif13a_aa": r["KIF13A_aa"]}
    sets = table(latest("clinvar_curation/KIF1A_calibration_sets_*.tsv", core.INTERIM))
    kif1a_path: dict[int, list[str]] = defaultdict(list)
    kif1a_tolerated: dict[int, list[str]] = defaultdict(list)
    for r in sets:
        if not r["protein_pos"]:
            continue
        pos = int(r["protein_pos"])
        if r["set"] == "missense_positive_core":
            kif1a_path[pos].append(r["protein_change"] + (" [fonksiyonel]" if r["functional_literature"] else ""))
        elif r["set"].startswith("missense_negative_motor_population"):
            kif1a_tolerated[pos].append(r["protein_change"])

    rows = []
    for g in gn:
        s = scores.get(g["variant_id"], {})
        pos = int(s.get("protein_pos") or g["protein_pos"] or 0) or None
        if not pos:
            continue
        am = float(s["am_pathogenicity"]) if s.get("am_pathogenicity") not in (None, "") else None
        af = float(g["joint_af"]) if g["joint_af"] else None
        p = inv.get(pos, {})
        identical = bool(p) and p.get("kif1a_aa") == p.get("kif13a_aa")
        path_here = kif1a_path.get(p.get("kif1a_pos", -1), [])
        tol_here = kif1a_tolerated.get(p.get("kif1a_pos", -1), [])
        dom = domain(pos)
        over = am is not None and am >= AM_CUTOFF
        if over and identical and path_here:
            tier = "A"
        elif over and identical and dom in ("motor", "neck_linker"):
            tier = "B"
        elif over:
            tier = "C"
        elif identical and path_here:
            tier = "D"
        else:
            tier = "-"
        rows.append({
            "tier": tier, "variant_id": g["variant_id"], "hgvsc": g["hgvsc"], "hgvsp": g["hgvsp"] or s.get("hgvsp"),
            "protein_change": f"{s.get('aa_ref') or g['aa_ref']}{pos}{s.get('aa_alt') or g['aa_alt']}",
            "protein_pos": pos, "domain": dom,
            "am_pathogenicity": am, "am_class": s.get("am_class"), "am_over_calibrated_cutoff": over,
            "am_over_generic_cutoff": am is not None and am >= AM_GENERIC,
            "esm1b": s.get("esm1b"), "cadd_phred": s.get("cadd_phred"), "primateai": s.get("primateai"),
            "gnomad_joint_ac": g["joint_ac"], "gnomad_joint_af": af, "gnomad_hom": g["joint_hom"],
            "gnomad_filters": g["joint_filters"],
            # AC0 means every carrier genotype failed quality filters: the variant is not confidently observed
            "confidently_observed": bool(g["joint_ac"]) and int(g["joint_ac"]) > 0 and "AC0" not in (g["joint_filters"] or ""),
            "af_too_common_for_dominant": af is not None and af >= AF_DOMINANT_MAX,
            "kif1a_pos": p.get("kif1a_pos"), "kif1a_aa": p.get("kif1a_aa"), "residue_identical": identical,
            "kif1a_pathogenic_at_residue": ";".join(path_here), "kif1a_tolerated_at_residue": ";".join(tol_here),
        })

    rows.sort(key=lambda r: ("ABCD-".index(r["tier"]), -(r["am_pathogenicity"] or 0)))
    out = core.save_tsv("", f"KIF13A_candidates_{core.today()}.tsv", rows, base=ROOT / "data" / "processed",
                        source="derived: gnomAD v4.1 + Ensembl VEP + KIF1A calibration sets + paralog map",
                        urls=[VEP_REGION], tool="scripts/build_kif13a_candidates.py",
                        params={"am_cutoff": AM_CUTOFF, "af_max": AF_DOMINANT_MAX})
    print(out)

    conf = [r for r in rows if r["confidently_observed"]]
    print(f"\nmissense rows: {len(rows)} | confidently observed (AC >= 1, no AC0 filter): {len(conf)}"
          f" | AC0-filtered (no reliable carrier): {len(rows) - len(conf)}")
    print("tiers (all):", Counter(r["tier"] for r in rows).most_common())
    print("tiers (confidently observed only):", Counter(r["tier"] for r in conf).most_common())
    print("AlphaMissense over calibrated cutoff 0.968:", sum(1 for r in rows if r["am_over_calibrated_cutoff"]),
          "| over generic cutoff 0.564:", sum(1 for r in rows if r["am_over_generic_cutoff"]))
    print("by domain (tier A+B+C, confident):", Counter(r["domain"] for r in conf if r["tier"] in "ABC").most_common())
    print("\nTIER A, confidently observed (model + paralog agree)")
    for r in [r for r in conf if r["tier"] == "A"]:
        print(f"  {r['protein_change']:<9} {r['hgvsc']:<14} AM {r['am_pathogenicity']:.3f}  AC={r['gnomad_joint_ac']:<5}"
              f" AF={r['gnomad_joint_af']:.2g}  KIF1A {r['kif1a_aa']}{r['kif1a_pos']}: {r['kif1a_pathogenic_at_residue']}")
    print("\nTIER A, AC0-filtered (not confidently seen in gnomAD; treat as theoretical)")
    for r in [r for r in rows if r["tier"] == "A" and not r["confidently_observed"]]:
        print(f"  {r['protein_change']:<9} {r['hgvsc']:<14} AM {r['am_pathogenicity']:.3f}  AC={r['gnomad_joint_ac']:<5}"
              f" AF={r['gnomad_joint_af']:.2g}  KIF1A {r['kif1a_aa']}{r['kif1a_pos']}: {r['kif1a_pathogenic_at_residue']}")
    print("\nTIER B, confidently observed, first 15 by AlphaMissense")
    for r in [r for r in conf if r["tier"] == "B"][:15]:
        print(f"  {r['protein_change']:<9} {r['hgvsc']:<14} AM {r['am_pathogenicity']:.3f}  AC={r['gnomad_joint_ac']:<5}"
              f" KIF1A {r['kif1a_aa']}{r['kif1a_pos']}" + (f"  (KIF1A'da tolere: {r['kif1a_tolerated_at_residue']})" if r["kif1a_tolerated_at_residue"] else ""))
    print("\nTIER D, confidently observed (paralog says pathogenic residue, AlphaMissense below cutoff)")
    for r in [r for r in conf if r["tier"] == "D"][:10]:
        print(f"  {r['protein_change']:<9} AM {(r['am_pathogenicity'] or 0):.3f}  KIF1A {r['kif1a_aa']}{r['kif1a_pos']}: {r['kif1a_pathogenic_at_residue']}")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
