"""Curate KIF1A ClinVar records into calibration sets for the KIF13A atlas.

For every candidate record (all P/LP; missense, splice and in-frame B/LB; conflicting missense/splice):
  1. fetch the VCV XML (submitters, method, origin, de novo, ACMG criteria in comments, first ClinVar date)
  2. annotate SNVs with Ensembl VEP on MANE (protein position, AlphaMissense, CADD, REVEL, SpliceAI, gnomAD AF)
  3. assign a tier and an exclusion reason (circularity, population frequency, review status, domain)
  4. report label counts and how much AUROC is driven by domain location instead of the variant itself

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\kif1a_clinvar_curation.py
"""

from __future__ import annotations

import asyncio
import csv
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, variants  # noqa: E402

core.load_dotenv()
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
VEP_REGION = "https://rest.ensembl.org/vep/human/region"
OUT = "clinvar_curation"

P_LP = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
B_LB = {"benign", "likely benign", "benign/likely benign"}

# KIF1A domains in MANE numbering (UniProt Q12756 canonical + isoform offsets from config/targets.yaml)
DOMAINS = [
    (1, 354, "motor"), (355, 365, "neck_linker"), (366, 383, "neck_coil_CC1"), (384, 524, "CC_region"),
    (525, 581, "FHA"), (582, 1675, "stalk"), (1676, 1774, "PH"), (1775, 1791, "C_terminal"),
]

TEXT_FLAGS = {
    "in_silico": r"\bPP3\b|\bBP4\b|in[- ]silico|computational|REVEL|CADD|SIFT|PolyPhen|AlphaMissense|algorithm|predict(ed|ion) (tool|to be)",
    "functional": r"\bPS3\b|\bBS3\b|functional (stud|assay|evidence|analys)|motility|in vitro|experimental stud|single[- ]molecule|axonal transport",
    "frequency": r"\bBA1\b|\bBS1\b|\bBS2\b|\bPM2\b|gnomAD|population database|allele frequency|ExAC|1000 Genomes",
    "de_novo": r"\bPS2\b|\bPM6\b|de novo",
    "same_residue": r"\bPS1\b|\bPM5\b|same amino acid change|different (amino acid|missense) change",
    "segregation": r"\bPP1\b|\bBS4\b|segregat",
}


def cls(c: str) -> str:
    c = (c or "").lower()
    if c in P_LP:
        return "P/LP"
    if c in B_LB:
        return "B/LB"
    if c == "uncertain significance":
        return "VUS"
    return "conflict" if "conflicting" in c else "other"


def csq_group(mc: str) -> str:
    if "missense variant" in mc:
        return "missense"
    if "splice donor" in mc or "splice acceptor" in mc:
        return "splice_canonical"
    if "nonsense" in mc or "frameshift" in mc:
        return "lof"
    if "inframe" in mc:
        return "inframe"
    if "synonymous" in mc:
        return "synonymous"
    return "intron" if "intron" in mc else "other"


def domain(pos: int | None) -> str | None:
    if not pos:
        return None
    return next((name for a, b, name in DOMAINS if a <= pos <= b), None)


def auroc(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """Mann-Whitney AUROC with ties counted as 0.5; None if either class is empty."""
    if not pos_scores or not neg_scores:
        return None
    wins = sum((p > n) + 0.5 * (p == n) for p in pos_scores for n in neg_scores)
    return round(wins / (len(pos_scores) * len(neg_scores)), 3)


# --------------------------------------------------------------------------- ClinVar XML


async def fetch_vcv(ids: list[str]) -> list[ET.Element]:
    key = core.secret("NCBI_API_KEY")
    archives: list[ET.Element] = []
    for i, batch in enumerate(core.chunks(ids, 50)):
        cache = core.RAW / "clinvar" / "vcv_xml" / f"KIF1A_vcv_{core.today()}_{i:03d}.xml"
        if cache.exists():
            content = cache.read_bytes()
        else:
            params = {"db": "clinvar", "id": ",".join(batch), "rettype": "vcv", "is_variationid": "", **({"api_key": key} if key else {})}
            content = (await core.request("GET", EUTILS, params=params, headers={"Accept": "application/xml"}, timeout=300)).content
            core.save_bytes("clinvar/vcv_xml", cache.name, content, source="NCBI E-utilities ClinVar VCV XML",
                            urls=[EUTILS], tool="scripts/kif1a_clinvar_curation.py", params={"n_ids": len(batch)})
        archives += ET.fromstring(content).findall("VariationArchive")
    return archives


def parse_vcv(va: ET.Element) -> dict:
    germ = va.find(".//Classifications/GermlineClassification")
    methods, origins, submitters, comments, pmids = set(), set(), [], [], set()
    for ca in va.findall(".//ClinicalAssertionList/ClinicalAssertion"):
        acc = ca.find("ClinVarAccession")
        submitters.append(acc.get("SubmitterName") if acc is not None else None)
        methods.update(m.text for m in ca.findall(".//ObservedIn/Method/MethodType") if m.text)
        origins.update(o.text for o in ca.findall(".//ObservedIn/Sample/Origin") if o.text)
        comments += [c.text for c in ca.findall("./Classification/Comment") if c.text]
        pmids.update(x.text for x in ca.findall(".//Citation/ID[@Source='PubMed']") if x.text)
    text = " ".join(comments)
    flags = {name: bool(re.search(rx, text, re.I)) for name, rx in TEXT_FLAGS.items()}
    flags["de_novo"] = flags["de_novo"] or any("de novo" in o for o in origins)
    return {
        "variation_id": va.get("VariationID"),
        "date_created": germ.get("DateCreated") if germ is not None else None,
        "n_submitters": int(germ.get("NumberOfSubmitters", 0)) if germ is not None else 0,
        "explanation": (germ.findtext("Explanation") or "") if germ is not None else "",
        "submitters": sorted({s for s in submitters if s}),
        "methods": sorted(methods), "origins": sorted(origins), "n_pmids": len(pmids),
        "literature_only": bool(methods) and methods <= {"literature only"},
        **{f"text_{k}": v for k, v in flags.items()},
    }


# --------------------------------------------------------------------------- tiers


def assign_tier(r: dict) -> tuple[str, str]:
    c, g, stars = r["class"], r["group"], int(r["stars"] or 0)
    af = r["gnomad_af_max"]
    if r["variant_type"] != "single nucleotide variant" and g in ("missense", "splice_canonical"):
        return "exclude", "not an SNV"
    if c == "P/LP":
        if g == "lof":
            return "lof_reference", "LoF: recessive HSN2C / haploinsufficiency question, not a missense calibration label"
        if g == "inframe":
            return "inframe_reference", "in-frame indel: not scorable by AlphaMissense/ESM"
        if stars == 0:
            return "exclude", "no assertion criteria"
        if g == "splice_canonical":
            return "positive_splice", "canonical splice P/LP: Evo 2 / SpliceAI positive"
        if g != "missense":
            return "exclude", f"consequence {g}"
        if r["protein_pos"] is None:
            return "exclude", "no MANE protein position"
        if af is not None and af >= 1e-5:
            return "positive_review", f"gnomAD AF {af:.2g} too common for a dominant de novo disorder (recessive allele or misclassified)"
        independent = r["text_de_novo"] or r["text_functional"] or r["n_submitters"] >= 2 or r["classification"].lower() != "likely pathogenic"
        if r["classification"].lower() == "likely pathogenic" and r["n_submitters"] < 2 and r["text_in_silico"] and not (r["text_de_novo"] or r["text_functional"]):
            return "positive_circular", "single-submitter LP whose comment cites in-silico evidence (PP3 may be decisive)"
        return ("positive_gold" if (stars >= 2 or r["text_de_novo"] or r["text_functional"]) else "positive_silver",
                "independent evidence (>=2 stars, de novo or functional)" if independent else "P/LP")
    if c == "B/LB":
        if g != "missense":
            return "negative_noncoding_reference", f"{g} B/LB"
        if stars == 0:
            return "exclude", "no assertion criteria"
        if (af is not None and af >= 1e-4) or (r["text_frequency"] and not r["text_in_silico"]):
            return "negative_frequency", "benign by population frequency: circular for AlphaMissense (trained on population frequency)"
        if r["text_in_silico"] and not r["text_frequency"]:
            return "negative_circular", "benign largely by BP4 / in-silico evidence"
        return "negative_other", "benign by other evidence (e.g. BS2, segregation, low-frequency but observed)"
    if c == "conflict" and g in ("missense", "splice_canonical"):
        e = r["explanation"].lower()  # e.g. "Pathogenic (2); Uncertain significance (1)"
        has_p, has_b = "pathogenic" in e, "benign" in e
        if has_p and not has_b:
            return "conflict_P_vs_VUS", "sensitivity analysis only"
        if has_b and not has_p:
            return "conflict_B_vs_VUS", "sensitivity analysis only"
        return "exclude", "conflict spans P and B"
    return "exclude", f"{c} / {g}"


# --------------------------------------------------------------------------- main


async def main() -> int:
    tsv = sorted((core.RAW / "clinvar").glob("KIF1A_clinvar_*.tsv"))[-1]
    rows = list(csv.DictReader(open(tsv, encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["class"], r["group"] = cls(r["classification"]), csq_group(r["molecular_consequence"] or "")
    cand = [r for r in rows if r["class"] == "P/LP"
            or (r["class"] == "B/LB" and r["group"] in ("missense", "splice_canonical", "inframe"))
            or (r["class"] == "conflict" and r["group"] in ("missense", "splice_canonical"))]
    print(f"source {core.rel(tsv)}: {len(rows)} records, {len(cand)} candidates", flush=True)

    parsed = {p["variation_id"]: p for p in map(parse_vcv, await fetch_vcv([r["uid"] for r in cand]))}
    print(f"VCV XML parsed: {len(parsed)}", flush=True)

    vep_in, key_of = [], {}
    for r in cand:
        if r["spdi"] and r["ref"] and r["alt"] and len(r["ref"]) == 1 and len(r["alt"]) == 1:
            _acc, p0, ref, alt = r["spdi"].split(":")
            line = f"2 {int(p0) + 1} . {ref} {alt} . . ."
            vep_in.append(line)
            key_of[r["uid"]] = f"2-{int(p0) + 1}-{ref}-{alt}"
    vep_in.sort(key=variants.vcf_sort_key)  # VEP rejects unsorted VCF-style input
    vep_rows: dict[str, dict] = {}
    raw_all = []
    for batch in core.chunks(vep_in, 200):
        raw = await core.post_json(VEP_REGION, {"variants": batch}, params=variants.VEP_PARAMS,
                                   headers={"Content-Type": "application/json"}, timeout=600)
        raw_all += raw
        for v in raw:
            f = variants._flatten_vep(v)
            chrom, pos, _dot, ref, alt = v["input"].split()[:5]
            vep_rows[f"{chrom}-{pos}-{ref}-{alt}"] = f
        print(f"VEP annotated {len(vep_rows)}/{len(vep_in)}", flush=True)
    core.save_json("vep", f"KIF1A_clinvar_candidates_vep_{core.today()}.raw.json", raw_all, compress=True,
                   source="Ensembl VEP REST", urls=[VEP_REGION], tool="scripts/kif1a_clinvar_curation.py", params=variants.VEP_PARAMS)

    curated = []
    for r in cand:
        p = parsed.get(r["uid"], {})
        v = vep_rows.get(key_of.get(r["uid"], ""), {})
        afs = [x for x in (v.get("gnomad_exomes_af"), v.get("gnomad_genomes_af")) if isinstance(x, (int, float))]
        pos = v.get("protein_pos") if v.get("mane_select") else None
        row = {
            "variation_id": r["uid"], "accession": r["accession"], "hgvsc": r["hgvsc"], "hgvsp_vep": v.get("hgvsp"),
            "protein_pos": pos, "aa_ref": v.get("aa_ref"), "aa_alt": v.get("aa_alt"), "domain": domain(pos),
            "classification": r["classification"], "class": r["class"], "group": r["group"], "stars": r["stars"],
            "variant_type": r["variant_type"], "variant_key": key_of.get(r["uid"]),
            "gnomad_af_max": max(afs) if afs else None,
            "am_pathogenicity": v.get("am_pathogenicity"), "cadd_phred": v.get("cadd_phred"), "revel": v.get("revel"),
            "spliceai_ds_max": v.get("spliceai_ds_max"), "last_evaluated": r["last_evaluated"],
            **{k: p.get(k) for k in ("date_created", "n_submitters", "explanation", "submitters", "methods", "origins",
                                     "n_pmids", "literature_only", *[f"text_{t}" for t in TEXT_FLAGS])},
        }
        row["tier"], row["reason"] = assign_tier(row)
        curated.append(row)

    prov = dict(source="derived: ClinVar VCV XML + Ensembl VEP", urls=[EUTILS, VEP_REGION], tool="scripts/kif1a_clinvar_curation.py")
    print(core.save_tsv(OUT, f"KIF1A_clinvar_curated_{core.today()}.tsv", curated, base=core.INTERIM, **prov))

    # ------------------------------------------------------------------ report
    print("\nTIERS")
    for t, n in Counter(r["tier"] for r in curated).most_common():
        print(f"  {n:>4}  {t}")
    print("\nexclusion / review reasons")
    for t, n in Counter((r["tier"], r["reason"]) for r in curated if r["tier"] in ("exclude", "positive_review", "positive_circular", "negative_circular")).most_common(12):
        print(f"  {n:>4}  {t[0]}: {t[1]}")

    pos_core = [r for r in curated if r["tier"] in ("positive_gold", "positive_silver")]
    neg_all = [r for r in curated if r["tier"] in ("negative_frequency", "negative_other")]
    uniq = {(r["protein_pos"], r["aa_alt"]) for r in pos_core}
    print(f"\nmissense positives (gold+silver): {len(pos_core)} records, {len(uniq)} unique protein changes, "
          f"{len({r['protein_pos'] for r in pos_core})} positions")
    print("  by domain:", Counter(r["domain"] for r in pos_core).most_common())
    print("  recurrent positions:", Counter(r["protein_pos"] for r in pos_core).most_common(8))
    print(f"missense negatives: {len(neg_all)}  by domain:", Counter(r["domain"] for r in neg_all).most_common())
    print("  negative tiers in motor domain:", Counter(r["tier"] for r in neg_all if r["domain"] == "motor"))

    def scores(rs, key):
        return [float(r[key]) for r in rs if r[key] not in (None, "")]

    print("\nAUROC (positives gold+silver vs negatives)")
    for label, negs in (("all negatives", neg_all), ("negatives excl. frequency-based", [r for r in neg_all if r["tier"] == "negative_other"])):
        in_motor = [1.0 if r["domain"] in ("motor", "neck_linker") else 0.0 for r in pos_core]
        out_motor = [1.0 if r["domain"] in ("motor", "neck_linker") else 0.0 for r in negs]
        print(f"  [{label}] n+={len(pos_core)} n-={len(negs)}  AlphaMissense {auroc(scores(pos_core,'am_pathogenicity'), scores(negs,'am_pathogenicity'))}"
              f"  CADD {auroc(scores(pos_core,'cadd_phred'), scores(negs,'cadd_phred'))}  REVEL {auroc(scores(pos_core,'revel'), scores(negs,'revel'))}"
              f"  'is motor domain' baseline {auroc(in_motor, out_motor)}")
    mp = [r for r in pos_core if r["domain"] == "motor"]
    mn = [r for r in neg_all if r["domain"] == "motor"]
    print(f"  [motor domain only] n+={len(mp)} n-={len(mn)}  AlphaMissense {auroc(scores(mp,'am_pathogenicity'), scores(mn,'am_pathogenicity'))}"
          f"  CADD {auroc(scores(mp,'cadd_phred'), scores(mn,'cadd_phred'))}  REVEL {auroc(scores(mp,'revel'), scores(mn,'revel'))}")

    print("\nevidence flags among missense P/LP candidates")
    mis_p = [r for r in curated if r["class"] == "P/LP" and r["group"] == "missense"]
    for k in ("text_in_silico", "text_functional", "text_de_novo", "text_frequency", "text_same_residue", "literature_only"):
        print(f"  {k:<20} {sum(1 for r in mis_p if r[k])}/{len(mis_p)}")
    print("  top submitters:", Counter(s for r in mis_p for s in (r["submitters"] or [])).most_common(6))
    print("  first ClinVar year:", sorted(Counter((r["date_created"] or "")[:4] for r in mis_p).items()))
    summary = {"tiers": Counter(r["tier"] for r in curated), "n_candidates": len(cand)}
    core.save_json(OUT, f"KIF1A_clinvar_curation_summary_{core.today()}.json", summary, base=core.INTERIM, **prov)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
