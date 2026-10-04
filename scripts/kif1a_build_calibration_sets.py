"""Build the KIF1A calibration sets used to evaluate and weight variant-effect models before transfer to KIF13A.

Inputs (run the two curation scripts first):
  data/interim/clinvar_curation/KIF1A_clinvar_curated_*.tsv        (kif1a_clinvar_curation.py)
  data/interim/clinvar_curation/KIF1A_gnomad_motor_controls_*.tsv  (kif1a_curation_extras.py)
  data/interim/clinvar_curation/KIF1A_literature_mentions_*.tsv    (kif1a_curation_extras.py)
  data/interim/paralog_maps/KIF1A_to_KIF13A_pairwise_MANE.tsv

Output: data/interim/clinvar_curation/KIF1A_calibration_sets_<date>.tsv, one row per variant per set.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

CUR = core.INTERIM / "clinvar_curation"
MOTOR_END = 365
POPULATION_MIN_AC = 3  # AC >= 3 main negative set; AC 2 kept as a sensitivity set

SET_NOTES = {
    "missense_positive_core": "ClinVar P/LP missense, >=1 star, not PP3-dependent, gnomAD AF < 1e-5 (plus literature rescue)",
    "missense_negative_motor_population": f"gnomAD v4.1 motor-domain missense with AC >= {POPULATION_MIN_AC}, not in ClinVar P/LP/conflict; "
                                          "circular for AlphaMissense (population-frequency training)",
    "missense_negative_motor_population_ac2": "same, AC == 2 (sensitivity only)",
    "missense_negative_clinvar_outside_motor": "ClinVar B/LB missense, >=1 star; none lie in the motor domain",
    "splice_positive": "ClinVar P/LP canonical splice donor/acceptor (Evo 2 / SpliceAI calibration)",
    "heldout_pp3_dependent_positive": "single-submitter LP citing in-silico evidence: evaluate separately to expose circularity",
    "sensitivity_conflict_P_vs_VUS": "conflicting P/LP vs VUS: sensitivity analysis only",
    "reference_lof": "P/LP nonsense/frameshift: HSN2C recessive / haploinsufficiency reference, not missense labels",
}


def latest(pattern: str) -> list[dict]:
    path = sorted(CUR.glob(pattern))[-1]
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def main() -> int:
    curated = latest("KIF1A_clinvar_curated_*.tsv")
    controls = latest("KIF1A_gnomad_motor_controls_*.tsv")
    literature = {r["variant"]: r["papers"] for r in latest("KIF1A_literature_mentions_*.tsv")}
    pmap = {int(r["KIF1A_pos"]): r for r in csv.DictReader(
        open(core.INTERIM / "paralog_maps" / "KIF1A_to_KIF13A_pairwise_MANE.tsv", encoding="utf-8"), delimiter="\t")}

    def base(pos, ref, alt):
        m = pmap.get(int(pos)) if pos else None
        key = f"{ref}{pos}{alt}" if pos else ""
        return {
            "protein_change": key, "protein_pos": pos, "cv_group": f"res{pos}" if pos else "",
            "kif13a_pos": m["KIF13A_pos"] if m else "", "kif13a_ref": m["KIF13A_aa"] if m else "",
            "residue_identical_in_kif13a": bool(m) and m["KIF13A_aa"] == ref,
            "functional_literature": literature.get(key, ""),
        }

    rows = []
    for r in curated:
        tier, key = r["tier"], f"{r['aa_ref']}{r['protein_pos']}{r['aa_alt']}"
        rescued = tier == "exclude" and r["reason"] == "no assertion criteria" and r["group"] == "missense" and key in literature
        set_name = {
            "positive_gold": "missense_positive_core", "positive_silver": "missense_positive_core",
            "negative_other": "missense_negative_clinvar_outside_motor", "negative_frequency": "missense_negative_clinvar_outside_motor",
            "positive_splice": "splice_positive", "positive_circular": "heldout_pp3_dependent_positive",
            "conflict_P_vs_VUS": "sensitivity_conflict_P_vs_VUS", "lof_reference": "reference_lof",
        }.get(tier, "missense_positive_core" if rescued else None)
        if not set_name:
            continue
        rows.append({
            "set": set_name, "label": 1 if set_name in ("missense_positive_core", "splice_positive", "heldout_pp3_dependent_positive") else
            0 if set_name.startswith("missense_negative") else "",
            "source": "ClinVar", "id": r["accession"], "hgvsc": r["hgvsc"], "variant_key": r["variant_key"], "domain": r["domain"],
            **base(r["protein_pos"], r["aa_ref"], r["aa_alt"]),
            "classification": r["classification"], "stars": r["stars"], "tier": "literature_rescue" if rescued else tier,
            "gnomad_af_max": r["gnomad_af_max"], "am_pathogenicity": r["am_pathogenicity"], "cadd_phred": r["cadd_phred"],
            "spliceai_ds_max": r["spliceai_ds_max"],
        })
    for c in controls:
        ac = int(c["joint_ac"] or 0)
        set_name = ("missense_negative_motor_population" if ac >= POPULATION_MIN_AC else
                    "missense_negative_motor_population_ac2" if ac == 2 else None)
        if not set_name:
            continue
        rows.append({
            "set": set_name, "label": 0, "source": "gnomAD v4.1", "id": c["variant_id"], "hgvsc": c["hgvsc"],
            "variant_key": c["variant_id"], "domain": "motor" if int(c["protein_pos"]) <= 354 else "neck_linker",
            **base(c["protein_pos"], c["aa_ref"], c["aa_alt"]),
            "classification": f"population AC={ac} hom={c['joint_hom']}", "stars": "", "tier": "",
            "gnomad_af_max": c["joint_af"], "am_pathogenicity": c["am"], "cadd_phred": c["cadd_phred"], "spliceai_ds_max": "",
        })

    out = core.save_tsv("clinvar_curation", f"KIF1A_calibration_sets_{core.today()}.tsv", rows, base=core.INTERIM,
                        source="derived: curated ClinVar + gnomAD motor controls + literature + paralog map",
                        urls=[], tool="scripts/kif1a_build_calibration_sets.py", params={"population_min_ac": POPULATION_MIN_AC})
    print(out)
    for s, n in Counter(r["set"] for r in rows).most_common():
        sub = [r for r in rows if r["set"] == s]
        motor = sum(1 for r in sub if r["domain"] in ("motor", "neck_linker"))
        ident = sum(1 for r in sub if r["residue_identical_in_kif13a"])
        func = sum(1 for r in sub if r["functional_literature"])
        print(f"  {n:>4}  {s:<42} motor={motor:<4} identical_residue_in_KIF13A={ident:<4} functional_lit={func:<3} {SET_NOTES[s][:70]}")
    core_pos = [r for r in rows if r["set"] == "missense_positive_core"]
    print("  core positives: unique residues", len({r["protein_pos"] for r in core_pos}),
          "| residues identical in KIF13A", len({r["protein_pos"] for r in core_pos if r["residue_identical_in_kif13a"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
