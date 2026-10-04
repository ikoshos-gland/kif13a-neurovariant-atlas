"""Follow-up to kif1a_clinvar_curation.py: the ClinVar benign set has no motor-domain variants, so
domain location alone separates the classes. This script looks for within-motor negatives and for an
independent (literature, functional) positive subset.

  1. inspect how the ClinVar missense negatives were classified (frequency vs in-silico wording, AF)
  2. gnomAD v4.1 motor-domain missense seen in >= 2 / 3 / 5 alleles as population tolerance controls,
     and AUROC of AlphaMissense / CADD / REVEL against them (within the motor domain only)
  3. variants named in open-access functional KIF1A papers (Europe PMC full text), validated against
     the MANE protein sequence, intersected with the ClinVar positives
  4. MaveDB search for KIF1A deep mutational scans

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\kif1a_curation_extras.py
"""

from __future__ import annotations

import asyncio
import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from atlas_mcp import core  # noqa: E402
from kif1a_clinvar_curation import auroc  # noqa: E402

core.load_dotenv()
MOTOR_END = 365  # motor domain + neck linker, MANE numbering (identical to canonical up to 393)
AA3 = {v: k for k, v in core.AA3.items()}

# Verified in Europe PMC on 2026-09-15 (open access full text)
PAPERS = {
    "PMC8054982": "Boyle 2021 HGG Adv, PMID 33880452 (severity + motility classes)",
    "PMC9371658": "Anazawa 2022 PNAS, PMID 35917346 (dominant-negative KAND)",
    "PMC6744892": "Chiba 2019 PNAS, PMID 31455732 (SPG variants hyperactivate KIF1A)",
    "PMC13254124": "Shatarupa 2026 Nat Commun, PMID 41912506 (R350 salt bridge)",
}


def latest(pattern: str, base: Path = core.RAW) -> Path:
    return sorted(base.glob(pattern))[-1]


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


async def main() -> int:
    curated = table(latest("clinvar_curation/KIF1A_clinvar_curated_*.tsv", core.INTERIM))
    seq = "".join((core.RAW / "sequences" / "KIF1A_ENST00000498729.fasta").read_text().splitlines()[1:])

    # ---------------------------------------------------------------- 1. ClinVar negatives
    negs = [r for r in curated if r["tier"].startswith("negative") and r["group"] == "missense"]
    def af_bin(a):
        if a in ("", None):
            return "absent"
        a = float(a)
        return "<1e-5" if a < 1e-5 else "1e-5..1e-4" if a < 1e-4 else "1e-4..1e-3" if a < 1e-3 else ">=1e-3"
    print("1. ClinVar missense negatives:", len(negs))
    print("   AF bins:", Counter(af_bin(r["gnomad_af_max"]) for r in negs))
    print("   wording (frequency, in_silico):", Counter((r["text_frequency"], r["text_in_silico"]) for r in negs))
    print("   classification:", Counter(r["classification"] for r in negs), "| stars:", Counter(r["stars"] for r in negs))

    # ---------------------------------------------------------------- 2. gnomAD motor controls
    gn = table(latest("gnomad/KIF1A_ENST00000498729_gnomad_r4_coding_*.tsv"))
    am = {row["protein_variant"]: float(row["am_pathogenicity"])
          for row in csv.DictReader(open(latest("alphamissense/AF-Q12756-F1-aa-substitutions_*.csv"), encoding="utf-8"))}
    clinvar_labelled = {(r["protein_pos"], r["aa_alt"]) for r in curated if r["class"] in ("P/LP", "conflict")}
    pos = [r for r in curated if r["tier"] in ("positive_gold", "positive_silver") and r["protein_pos"] and int(r["protein_pos"]) <= MOTOR_END]

    def am_score(p, ref, alt):
        return am.get(f"{ref}{p}{alt}")

    controls = []
    for r in gn:
        if r["consequence"] != "missense_variant" or not r["protein_pos"] or int(r["protein_pos"]) > MOTOR_END:
            continue
        if (r["protein_pos"], r["aa_alt"]) in clinvar_labelled:
            continue
        if seq[int(r["protein_pos"]) - 1] != r["aa_ref"]:
            raise RuntimeError(f"MANE sequence mismatch at {r['hgvsp']}")
        controls.append({**r, "am": am_score(int(r["protein_pos"]), r["aa_ref"], r["aa_alt"])})
    print(f"\n2. gnomAD motor-domain (1-{MOTOR_END}) missense not labelled in ClinVar: {len(controls)}")
    print("   by allele count:", Counter("AC1" if int(r["joint_ac"] or 0) == 1 else "AC2" if int(r["joint_ac"] or 0) == 2
                                        else "AC3-4" if int(r["joint_ac"] or 0) < 5 else "AC>=5" for r in controls))
    p_am = [am_score(int(r["protein_pos"]), r["aa_ref"], r["aa_alt"]) for r in pos]
    p_am = [x for x in p_am if x is not None]
    p_cadd = [float(r["cadd_phred"]) for r in pos if r["cadd_phred"]]
    p_revel = [float(r["revel"]) for r in pos if r["revel"]]
    for min_ac in (1, 2, 3, 5):
        c = [r for r in controls if int(r["joint_ac"] or 0) >= min_ac]
        print(f"   within motor: positives n={len(pos)} vs gnomAD AC>={min_ac} n={len(c)}:"
              f"  AlphaMissense {auroc(p_am, [r['am'] for r in c if r['am'] is not None])}"
              f"  CADD {auroc(p_cadd, [float(r['cadd_phred']) for r in c if r['cadd_phred']])}"
              f"  REVEL {auroc(p_revel, [float(r['revel']) for r in c if r['revel']])}")
    strong = sorted((r for r in controls if int(r["joint_ac"] or 0) >= 2), key=lambda r: -(r["am"] or 0))[:8]
    print("   highest-AlphaMissense tolerated controls (AC>=2):",
          [(r["hgvsp"], r["joint_ac"], r["am"]) for r in strong])
    core.save_tsv("clinvar_curation", f"KIF1A_gnomad_motor_controls_{core.today()}.tsv", controls, base=core.INTERIM,
                  source="derived: gnomAD v4.1 MANE transcript variants + AlphaFold DB AlphaMissense",
                  urls=[], tool="scripts/kif1a_curation_extras.py", params={"motor_end": MOTOR_END})

    # ---------------------------------------------------------------- 3. literature full text
    print("\n3. functional papers (Europe PMC full text)")
    mention_rx = re.compile(r"\b(?:p\.)?(?:([ACDEFGHIKLMNPQRSTVWY])(\d{1,4})([ACDEFGHIKLMNPQRSTVWY])|"
                            r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val)(\d{1,4})"
                            r"(Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val))\b")
    pos_keys = {f"{r['aa_ref']}{r['protein_pos']}{r['aa_alt']}": r for r in curated if r["class"] == "P/LP" and r["protein_pos"]}
    literature: dict[str, set[str]] = {}
    for pmcid, label in PAPERS.items():
        url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
        try:
            xml = (await core.request("GET", url, headers={"Accept": "application/xml"}, timeout=120)).text
        except Exception as exc:
            print(f"   {label}: full text not available ({str(exc)[:60]})")
            continue
        text = re.sub(r"<[^>]+>", " ", xml)
        found = set()
        for m in mention_rx.finditer(text):
            if m.group(1):
                ref, p, alt = m.group(1), int(m.group(2)), m.group(3)
            else:
                ref, p, alt = core.AA3[m.group(4)], int(m.group(5)), core.AA3[m.group(6)]
            if 1 <= p <= len(seq) and seq[p - 1] == ref and ref != alt:
                found.add(f"{ref}{p}{alt}")
        for v in found:
            literature.setdefault(v, set()).add(label.split(",")[0])
        in_clinvar = sorted(v for v in found if v in pos_keys)
        print(f"   {label}: {len(found)} KIF1A-consistent variant names, {len(in_clinvar)} are ClinVar P/LP")
    lit_pos = [r for r in curated if r["tier"].startswith("positive") and f"{r['aa_ref']}{r['protein_pos']}{r['aa_alt']}" in literature]
    print("   ClinVar positives also characterised in these papers:", len(lit_pos),
          "| by tier:", Counter(r["tier"] for r in lit_pos))
    not_in_clinvar = sorted((v for v in literature if v not in pos_keys and int(re.sub(r"\D", "", v)) <= MOTOR_END),
                            key=lambda v: int(re.sub(r"\D", "", v)))
    print("   motor variants named in papers but not ClinVar P/LP:", not_in_clinvar[:30])
    core.save_tsv("clinvar_curation", f"KIF1A_literature_mentions_{core.today()}.tsv",
                  [{"variant": v, "papers": sorted(s), "clinvar_tier": (next((r["tier"] for r in curated
                    if f"{r['aa_ref']}{r['protein_pos']}{r['aa_alt']}" == v), ""))} for v, s in sorted(literature.items())],
                  base=core.INTERIM, source="derived: Europe PMC full text regex, validated against MANE sequence",
                  urls=list(PAPERS), tool="scripts/kif1a_curation_extras.py")

    # ---------------------------------------------------------------- 4. MaveDB
    try:
        hits = await core.post_json("https://api.mavedb.org/api/v1/score-sets/search", {"text": "KIF1A"},
                                    headers={"Content-Type": "application/json"}, timeout=120)
        items = hits if isinstance(hits, list) else hits.get("scoreSets", hits.get("items", []))
        print(f"\n4. MaveDB score sets matching KIF1A: {len(items)}",
              [(i.get("urn"), i.get("title", "")[:60]) for i in items[:5]])
    except Exception as exc:
        print("\n4. MaveDB search failed:", str(exc)[:120])
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
