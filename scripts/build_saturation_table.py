"""Saturation lookup for KIF13A missense: all 1805 x 19 = 34295 possible amino acid substitutions.

Why this exists. `scripts/build_kif13a_candidates.py` ranks the 2671 missense variants that gnomAD has
actually seen. That is a population catalogue, not a disease catalogue: 99 of the 111 curated KIF1A
pathogenic missense variants (89 percent) are absent from gnomAD, and so is the only KIF13A patient
variant with a phenotype (p.Lys1122Glu, VCV002570672). A variant that causes a dominant paediatric
disorder is *expected* to be missing from a population database, so a gnomAD-anchored list can never
contain the candidate we most want. This table therefore enumerates the whole substitution space and
annotates every cell, with "is it in gnomAD" demoted from a filter to just another column.

Per substitution:
  - AlphaMissense score and class (per-protein table for Q9H1H9, canonical == MANE, no offset)
    plus the gene-calibrated cutoff 0.968 from docs/kalibrasyon_olcumleri.md (generic cutoff 0.564)
  - whether one nucleotide change of the MANE CDS (ENST00000259711) can produce it, the codon change
    and the implied c. notation
  - gnomAD v4.1 status: observed with a confident carrier, AC0-filtered, or never seen
  - ClinVar record if any (germline classification, review stars, traits)
  - paralog evidence: equivalent KIF1A residue, identity, curated KIF1A pathogenic / tolerated variants
  - structure (motor domain and neck linker only, positions 1-380): class, nucleotide and tubulin
    distances, relative SASA, pLDDT. Blank elsewhere. The tail was never assessed structurally.
  - the tier rule of build_kif13a_candidates.py, reused (constants and domain map are imported from it,
    and the run verifies that this table reproduces every tier in the existing candidate file)

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\build_saturation_table.py
"""

from __future__ import annotations

import asyncio
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from atlas_mcp import core  # noqa: E402
import build_kif13a_candidates as base  # noqa: E402  (domain map, cutoffs, tier rule live there)

PROTEIN_LEN = 1805
TRANSCRIPT = "ENST00000259711"
ENSEMBL_CDS = f"https://rest.ensembl.org/sequence/id/{TRANSCRIPT}"
AA_ALL = "ACDEFGHIKLMNPQRSTVWY"

CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def read_fasta(path: Path) -> str:
    return "".join(l.strip() for l in path.read_text(encoding="utf-8").splitlines() if not l.startswith(">"))


def tier_of(am: float | None, identical: bool, path_here: list[str], dom: str) -> str:
    """The tier rule of scripts/build_kif13a_candidates.py, transcribed unchanged.

    A  AM >= cutoff, identical paralog residue, KIF1A pathogenic variant at that residue
    B  AM >= cutoff, identical paralog residue in the motor domain / neck linker
    C  AM >= cutoff, elsewhere
    D  AM <  cutoff but KIF1A pathogenic at an identical residue
    -  everything else
    The run below proves this reproduces every tier of data/processed/KIF13A_candidates_*.tsv.
    """
    over = am is not None and am >= base.AM_CUTOFF
    if over and identical and path_here:
        return "A"
    if over and identical and dom in ("motor", "neck_linker"):
        return "B"
    if over:
        return "C"
    if identical and path_here:
        return "D"
    return "-"


async def mane_cds(protein: str) -> str:
    """CDS of the MANE transcript, fetched once and cached under data/raw/sequences."""
    cache = core.RAW / "sequences" / f"KIF13A_{TRANSCRIPT}_cds.fasta"
    if cache.exists():
        print(f"using cached {core.rel(cache)}")
        cds = read_fasta(cache)
    else:
        text = (await core.request("GET", ENSEMBL_CDS, params={"type": "cds"},
                                   headers={"Accept": "text/x-fasta"})).text
        core.save_bytes("sequences", cache.name, text.encode("utf-8"), source="Ensembl REST",
                        urls=[ENSEMBL_CDS], tool="scripts/build_saturation_table.py",
                        params={"transcript": TRANSCRIPT, "type": "cds"})
        print(f"fetched {core.rel(cache)}")
        cds = "".join(l.strip() for l in text.splitlines() if not l.startswith(">"))
    cds = cds.upper()
    if len(cds) != 3 * (PROTEIN_LEN + 1):
        raise SystemExit(f"CDS length {len(cds)} != {3 * (PROTEIN_LEN + 1)} (1805 aa + stop)")
    translated = "".join(CODONS[cds[i:i + 3]] for i in range(0, len(cds) - 3, 3))
    if translated != protein:
        first = next(i for i, (a, b) in enumerate(zip(translated, protein)) if a != b)
        raise SystemExit(f"CDS does not translate to the MANE protein, first mismatch at residue {first + 1}")
    if CODONS[cds[-3:]] != "*":
        raise SystemExit("CDS does not end in a stop codon")
    print(f"CDS ok: {len(cds)} nt, translates to the 1805 aa MANE protein, ends in {cds[-3:]}")
    return cds


def snv_paths(codon: str) -> dict[str, list[tuple[int, str, str]]]:
    """alt aa -> [(codon offset 0-2, reference base, alternate base), ...] reachable by one SNV."""
    out: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for i in range(3):
        for b in "ACGT":
            if b == codon[i]:
                continue
            aa = CODONS[codon[:i] + b + codon[i + 1:]]
            if aa not in ("*", CODONS[codon]):
                out[aa].append((i, codon[i], b))
    return out


INK, MUTED, PALE, MID, ACCENT = "#0b0b0b", "#52514e", "#d8d7d0", "#7d7c75", "#d1495b"


def funnel_figure(stages: list[tuple[str, str, int, int]], kif1a: tuple[int, int]) -> str:
    """One panel: how the substitution space narrows, and how much of each stage gnomAD can even see.

    stages: (label, sublabel, n, n confidently observed in gnomAD), widest first.
    kif1a: (curated KIF1A pathogenic missense, how many of them gnomAD has seen) as the reference bar.
    """
    import math

    W, H, X0, BAR = 1160, 600, 300, 560
    top = math.log10(stages[0][2])
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
         'font-family="Arial, Helvetica, sans-serif">', f'<rect width="{W}" height="{H}" fill="#ffffff"/>',
         '<defs><pattern id="none" width="6" height="6" patternTransform="rotate(45)" '
         f'patternUnits="userSpaceOnUse"><rect width="6" height="6" fill="#ffffff"/>'
         f'<line x1="0" y1="0" x2="0" y2="6" stroke="{PALE}" stroke-width="3"/></pattern></defs>',
         f'<text x="40" y="34" font-size="15" fill="{INK}">KIF13A missense: butun degisim uzayi ve gnomAD kor noktasi</text>',
         f'<text x="40" y="56" font-size="11" fill="{MUTED}">Dolu kisim gnomAD v4.1&#8217;de guvenilir tasiyicisi olan degisimler; tarali kisim gnomAD&#8217;in hic gormedigi ya da AC0 ile eledigi degisimler.</text>',
         f'<text x="40" y="72" font-size="11" fill="{MUTED}">Genislik log10 olcekli. Eski aday listesi yalnizca dolu dilimlerden secilebiliyordu.</text>']
    y = 110
    for label, sub, n, seen in stages:
        w = BAR * math.log10(n) / top
        wseen = w * seen / n
        o.append(f'<text x="{X0 - 14}" y="{y + 15}" font-size="12" fill="{INK}" text-anchor="end">{label}</text>')
        o.append(f'<text x="{X0 - 14}" y="{y + 30}" font-size="10" fill="{MUTED}" text-anchor="end">{sub}</text>')
        o.append(f'<rect x="{X0}" y="{y}" width="{w:.1f}" height="26" fill="url(#none)" stroke="{MID}" stroke-width="0.8"/>')
        o.append(f'<rect x="{X0}" y="{y}" width="{max(wseen, 1.2):.1f}" height="26" fill="{MID}"/>')
        o.append(f'<text x="{X0 + w + 12}" y="{y + 18}" font-size="12" fill="{INK}">{n}</text>')
        o.append(f'<text x="{X0 + w + 12}" y="{y + 32}" font-size="10" fill="{ACCENT}">{n - seen} tanesi gnomAD&#8217;de guvenilir degil'
                 f' (%{100 * (n - seen) / n:.0f})</text>')
        y += 66
    y += 18
    o.append(f'<line x1="40" y1="{y - 22}" x2="{W - 40}" y2="{y - 22}" stroke="{PALE}" stroke-width="1"/>')
    n1, seen1 = kif1a
    w = BAR * math.log10(n1) / top
    o.append(f'<text x="{X0 - 14}" y="{y + 15}" font-size="12" fill="{ACCENT}" text-anchor="end">Referans: KIF1A patojenik missense</text>')
    o.append(f'<text x="{X0 - 14}" y="{y + 30}" font-size="10" fill="{MUTED}" text-anchor="end">kurate cekirdek set, klinik olarak kanitlanmis</text>')
    o.append(f'<rect x="{X0}" y="{y}" width="{w:.1f}" height="26" fill="url(#none)" stroke="{ACCENT}" stroke-width="0.8"/>')
    o.append(f'<rect x="{X0}" y="{y}" width="{w * seen1 / n1:.1f}" height="26" fill="{ACCENT}"/>')
    o.append(f'<text x="{X0 + w + 12}" y="{y + 18}" font-size="12" fill="{INK}">{n1}</text>')
    o.append(f'<text x="{X0 + w + 12}" y="{y + 32}" font-size="10" fill="{ACCENT}">{n1 - seen1} tanesi gnomAD&#8217;de yok'
             f' (%{100 * (n1 - seen1) / n1:.0f})</text>')
    o.append(f'<text x="40" y="{H - 44}" font-size="11" fill="{INK}">Okunusu: bir hastalik varyantinin gnomAD&#8217;de bulunmamasi beklenen seydir. KIF1A&#8217;da kurate patojenik missense&#8217;lerin %89&#8217;u gnomAD&#8217;de yok.</text>')
    o.append(f'<text x="40" y="{H - 26}" font-size="11" fill="{MUTED}">Tek fenotipli KIF13A hasta varyanti p.Lys1122Glu (VCV002570672, AlphaMissense 0,8044) de gnomAD&#8217;de yok, bu yuzden eski aday listesine hic giremedi.</text>')
    o.append("</svg>")
    return "\n".join(o)


async def main() -> int:
    protein = read_fasta(core.RAW / "sequences" / f"KIF13A_{TRANSCRIPT}.fasta").upper()
    if len(protein) != PROTEIN_LEN:
        raise SystemExit(f"protein length {len(protein)} != {PROTEIN_LEN}")
    cds = await mane_cds(protein)

    # ---------------------------------------------------------------- AlphaMissense (Q9H1H9 canonical == MANE)
    am_path = sorted((core.RAW / "alphamissense").glob("AF-Q9H1H9-F1-aa-substitutions_*.csv"))[-1]
    am: dict[str, dict] = {}
    for r in csv.DictReader(open(am_path, encoding="utf-8")):
        am[r["protein_variant"]] = r
    print(f"AlphaMissense rows: {len(am)} from {core.rel(am_path)}")
    bad = [v for v in am if protein[int(v[1:-1]) - 1] != v[0]]
    if bad:
        raise SystemExit(f"AlphaMissense reference residue disagrees with the MANE protein, e.g. {bad[:5]}")

    # ---------------------------------------------------------------- gnomAD v4.1, MANE transcript query
    gn: dict[str, dict] = {}
    for g in table(base.latest("gnomad/KIF13A_ENST00000259711_gnomad_r4_all_*.tsv")):
        if g["consequence"] != "missense_variant" or not g["protein_pos"]:
            continue
        key = f"{g['aa_ref']}{g['protein_pos']}{g['aa_alt']}"
        ac = int(g["joint_ac"] or 0)
        confident = ac > 0 and "AC0" not in (g["joint_filters"] or "")
        prev = gn.get(key)
        # one protein change can arise from several codon-level variants: keep the best-supported one
        if prev is None or (confident, ac) > (prev["_confident"], int(prev["joint_ac"] or 0)):
            gn[key] = {**g, "_confident": confident}
    print(f"gnomAD missense protein changes: {len(gn)}")

    # ---------------------------------------------------------------- ClinVar (germline_classification)
    cv: dict[str, list[dict]] = defaultdict(list)
    for c in table(base.latest("clinvar/KIF13A_clinvar_*.tsv")):
        if "missense" not in (c["molecular_consequence"] or "") or not c["protein_pos"]:
            continue
        cv[f"{c['aa_ref']}{c['protein_pos']}{c['aa_alt']}"].append(c)
    print(f"ClinVar missense protein changes: {len(cv)}")

    # ---------------------------------------------------------------- paralog transfer from KIF1A
    inv: dict[int, dict] = {}
    for r in table(core.INTERIM / "paralog_maps" / "KIF1A_to_KIF13A_pairwise_MANE.tsv"):
        inv[int(r["KIF13A_pos"])] = {"kif1a_pos": int(r["KIF1A_pos"]), "kif1a_aa": r["KIF1A_aa"],
                                     "kif13a_aa": r["KIF13A_aa"]}
    kif1a_path: dict[int, list[str]] = defaultdict(list)
    kif1a_tol: dict[int, list[str]] = defaultdict(list)
    for r in table(base.latest("clinvar_curation/KIF1A_calibration_sets_*.tsv", core.INTERIM)):
        if not r["protein_pos"]:
            continue
        pos = int(r["protein_pos"])
        if r["set"] == "missense_positive_core":
            kif1a_path[pos].append(r["protein_change"] + (" [fonksiyonel]" if r["functional_literature"] else ""))
        elif r["set"].startswith("missense_negative_motor_population"):
            kif1a_tol[pos].append(r["protein_change"])
    print(f"paralog map covers {len(inv)} KIF13A residues; KIF1A pathogenic residues: {len(kif1a_path)}")

    # ---------------------------------------------------------------- structure axis (positions 1-380 only)
    st: dict[int, dict] = {}
    for r in table(base.latest("structure/KIF13A_motor_residue_structure_*.tsv", core.INTERIM)):
        st[int(r["protein_pos"])] = r
    print(f"structure rows: {len(st)} (positions {min(st)}-{max(st)}); everything above is left blank")

    # ---------------------------------------------------------------- the 34295 rows
    rows: list[dict] = []
    for pos in range(1, PROTEIN_LEN + 1):
        ref = protein[pos - 1]
        codon = cds[3 * (pos - 1):3 * pos]
        paths = snv_paths(codon)
        dom = base.domain(pos)
        p = inv.get(pos, {})
        identical = bool(p) and p.get("kif1a_aa") == p.get("kif13a_aa")
        path_here = kif1a_path.get(p.get("kif1a_pos", -1), [])
        tol_here = kif1a_tol.get(p.get("kif1a_pos", -1), [])
        # the sharpest form of paralog transfer: KIF1A is pathogenic with the *same* replacement residue
        path_alts = {m.group(1) for c in path_here for m in [re.match(r"^[A-Z]\d+([A-Z])", c)] if m}
        s = st.get(pos)
        for alt in AA_ALL:
            if alt == ref:
                continue
            change = f"{ref}{pos}{alt}"
            a = am.get(change, {})
            score = float(a["am_pathogenicity"]) if a.get("am_pathogenicity") else None
            hits = paths.get(alt, [])
            g = gn.get(change)
            records = cv.get(change, [])
            rows.append({
                "protein_change": change, "protein_pos": pos, "aa_ref": ref, "aa_alt": alt, "domain": dom,
                "tier": tier_of(score, identical, path_here, dom),
                "am_pathogenicity": score, "am_class": a.get("am_class"),
                "am_over_calibrated_cutoff": score is not None and score >= base.AM_CUTOFF,
                "am_over_generic_cutoff": score is not None and score >= base.AM_GENERIC,
                "snv_accessible": bool(hits), "n_snv_paths": len(hits), "codon_ref": codon,
                "codon_change": ";".join(f"{codon}>{codon[:i] + b + codon[i + 1:]}" for i, _r, b in hits),
                "hgvsc_snv": ";".join(f"c.{3 * (pos - 1) + i + 1}{r0}>{b}" for i, r0, b in hits),
                "gnomad_status": ("observed_confident" if g and g["_confident"]
                                  else "ac0_filtered" if g else "not_observed"),
                "gnomad_variant_id": (g or {}).get("variant_id"), "gnomad_hgvsc": (g or {}).get("hgvsc"),
                "gnomad_joint_ac": (g or {}).get("joint_ac"),
                "gnomad_joint_af": float(g["joint_af"]) if g and g["joint_af"] else None,
                "gnomad_hom": (g or {}).get("joint_hom"), "gnomad_filters": (g or {}).get("joint_filters"),
                "confidently_observed": bool(g and g["_confident"]),
                "af_too_common_for_dominant": bool(g and g["joint_af"]
                                                   and float(g["joint_af"]) >= base.AF_DOMINANT_MAX),
                "clinvar_accession": ";".join(c["accession"] for c in records),
                "clinvar_classification": ";".join(c["classification"] for c in records),
                "clinvar_review_stars": max((int(c["stars"] or 0) for c in records), default=None),
                "clinvar_traits": ";".join(dict.fromkeys(t for c in records for t in (c["traits"] or "").split(";") if t)),
                "kif1a_pos": p.get("kif1a_pos"), "kif1a_aa": p.get("kif1a_aa"), "residue_identical": identical,
                "kif1a_pathogenic_at_residue": ";".join(path_here),
                "kif1a_same_substitution_pathogenic": identical and alt in path_alts,
                "kif1a_tolerated_at_residue": ";".join(tol_here),
                "structure_assessed": s is not None,
                "structural_class": (s or {}).get("structural_class"),
                "dist_nucleotide_A": (s or {}).get("dist_nucleotide_A"),
                "dist_tubulin_A": (s or {}).get("dist_tubulin_A"),
                "nearest_tubulin": (s or {}).get("nearest_tubulin"),
                "rel_sasa": (s or {}).get("rel_sasa"), "plddt": (s or {}).get("plddt"),
                # same definition as scripts/structure_axis.py: confident model and a surviving superposition pair
                "structure_evidence_ok": bool(s and s["plddt"] and float(s["plddt"]) >= 70
                                              and s["superpose_pair"] == "True"),
            })

    expected = PROTEIN_LEN * 19
    if len(rows) != expected:
        raise SystemExit(f"{len(rows)} rows, expected {expected}")
    missing_am = [r["protein_change"] for r in rows if r["am_pathogenicity"] is None]
    print(f"\nrows: {len(rows)} | without an AlphaMissense score: {len(missing_am)}")

    out = core.save_tsv("", f"KIF13A_saturation_lookup_{core.today()}.tsv", rows,
                        base=ROOT / "data" / "processed",
                        source="derived: AlphaMissense Q9H1H9 + Ensembl MANE CDS + gnomAD v4.1 + ClinVar"
                               " + KIF1A paralog map and calibration sets + AlphaFold/8UTN structure axis",
                        urls=[ENSEMBL_CDS], tool="scripts/build_saturation_table.py",
                        params={"transcript": TRANSCRIPT, "protein_len": PROTEIN_LEN,
                                "am_cutoff": base.AM_CUTOFF, "am_generic": base.AM_GENERIC,
                                "af_max": base.AF_DOMINANT_MAX})
    print(out)

    # ---------------------------------------------------------------- sanity checks
    by_change = {r["protein_change"]: r for r in rows}
    checks: list[tuple[str, bool, str]] = []

    k = by_change["K1122E"]
    checks.append(("K1122E present with AlphaMissense 0.8044",
                   k["am_pathogenicity"] == 0.8044, f"AM={k['am_pathogenicity']}"))
    checks.append(("K1122E has no paralog mapping", k["kif1a_pos"] is None,
                   f"kif1a_pos={k['kif1a_pos']}"))
    checks.append(("K1122E absent from gnomAD", k["gnomad_status"] == "not_observed", k["gnomad_status"]))
    checks.append(("K1122E carries ClinVar VCV002570672",
                   "VCV002570672" in k["clinvar_accession"], k["clinvar_accession"]))
    checks.append(("K1122E reachable as c.3364A>G", k["hgvsc_snv"] == "c.3364A>G", k["hgvsc_snv"]))
    checks.append(("K1122E not structurally assessed (tail)", not k["structure_assessed"], ""))

    old = table(base.latest("KIF13A_candidates_*.tsv", ROOT / "data" / "processed"))
    mism = [(o["protein_change"], o["tier"], by_change[o["protein_change"]]["tier"]) for o in old
            if o["protein_change"] in by_change and by_change[o["protein_change"]]["tier"] != o["tier"]]
    absent = [o["protein_change"] for o in old if o["protein_change"] not in by_change]
    old_a = [o["protein_change"] for o in old if o["tier"] == "A"]
    a_kept = sum(1 for c in old_a if by_change.get(c, {}).get("tier") == "A")
    checks.append((f"all {len(old)} existing candidate rows reproduce their tier",
                   not mism and not absent, f"{len(mism)} mismatches, {len(absent)} not found"))
    checks.append((f"all {len(old_a)} existing tier A keep tier A", a_kept == len(old_a),
                   f"{a_kept}/{len(old_a)}"))
    obs = {r["protein_change"] for r in rows if r["gnomad_status"] != "not_observed"}
    checks.append(("every gnomAD missense protein change is in the table",
                   len(obs) == len(gn), f"{len(obs)} vs {len(gn)}"))

    print("\nsanity checks")
    for name, ok, detail in checks:
        print(f"  [{'ok ' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if mism[:5]:
        print("  first tier mismatches:", mism[:5])

    # ---------------------------------------------------------------- headline numbers
    over = [r for r in rows if r["am_over_calibrated_cutoff"]]
    over_snv = [r for r in over if r["snv_accessible"]]
    print(f"\nAlphaMissense >= {base.AM_CUTOFF} (gene-calibrated): {len(over)} of {len(rows)}"
          f" | of those reachable by one nucleotide change: {len(over_snv)}")
    print(f"AlphaMissense >= {base.AM_GENERIC} (generic): {sum(1 for r in rows if r['am_over_generic_cutoff'])}")
    print(f"reachable by one nucleotide change overall: {sum(1 for r in rows if r['snv_accessible'])}")
    print("tiers (all):", Counter(r["tier"] for r in rows).most_common())
    print("tiers (SNV-reachable only):", Counter(r["tier"] for r in rows if r["snv_accessible"]).most_common())
    print("gnomAD status:", Counter(r["gnomad_status"] for r in rows).most_common())
    for t in "ABD":
        sub = [r for r in rows if r["tier"] == t and r["snv_accessible"]]
        seen = [r for r in sub if r["gnomad_status"] == "observed_confident"]
        print(f"  tier {t}: {len(sub)} reachable, {len(seen)} confidently in gnomAD,"
              f" {len(sub) - len(seen)} never confidently seen")

    exact = [r for r in rows if r["kif1a_same_substitution_pathogenic"] and r["snv_accessible"]]
    exact_a = [r for r in exact if r["tier"] == "A"]
    exact_seen = [r for r in exact if r["gnomad_status"] == "observed_confident"]
    print(f"\nreachable substitutions whose exact twin is pathogenic in KIF1A: {len(exact)}"
          f" (of which tier A: {len(exact_a)})"
          f" | confidently in gnomAD: {len(exact_seen)} | never confidently seen: {len(exact) - len(exact_seen)}")
    print("  never confidently seen, by AlphaMissense:")
    seen_names = {r["protein_change"] for r in exact_seen}
    for r in sorted((r for r in exact if r["protein_change"] not in seen_names),
                    key=lambda r: -(r["am_pathogenicity"] or 0)):
        print(f"    {r['protein_change']:<8} AM {r['am_pathogenicity']:.4f}  {r['gnomad_status']:<18}"
              f" {r['structural_class'] or '-':<42} pLDDT {r['plddt'] or '-':<6}"
              f" KIF1A {r['kif1a_pathogenic_at_residue']}  {r['hgvsc_snv']}")

    fig = ROOT / "figures" / "KIF13A_missense_saturation_blindspot.svg"
    fig.write_text(funnel_figure(
        [("Butun olasi degisimler", "1805 aa x 19", len(rows), len({r["protein_change"] for r in rows if r["gnomad_status"] == "observed_confident"})),
         ("Tek nukleotidle erisilebilir", "MANE CDS ENST00000259711", len(snv_reach := [r for r in rows if r["snv_accessible"]]),
          sum(1 for r in snv_reach if r["gnomad_status"] == "observed_confident")),
         (f"AlphaMissense >= {base.AM_CUTOFF}", "gen-kalibre esik, erisilebilir olanlar", len(over_snv),
          sum(1 for r in over_snv if r["gnomad_status"] == "observed_confident")),
         ("Katman A", "model + ozdes paralog kalinti + KIF1A patojenik", len([r for r in rows if r["tier"] == "A" and r["snv_accessible"]]),
          sum(1 for r in rows if r["tier"] == "A" and r["snv_accessible"] and r["gnomad_status"] == "observed_confident")),
         ("Katman A + birebir paralog esi", "ayni kalinti, ayni yerine gecen aminoasit", len(exact_a),
          sum(1 for r in exact_a if r["gnomad_status"] == "observed_confident"))],
        (111, 12)), encoding="utf-8")
    print(f"\n{core.rel(fig)}")

    summary = {
        "rows": len(rows), "protein_len": PROTEIN_LEN, "transcript": TRANSCRIPT,
        "kif1a_exact_substitution_match_snv": len(exact),
        "kif1a_exact_substitution_match_snv_tier_a": len(exact_a),
        "kif1a_exact_substitution_match_snv_observed": len(exact_seen),
        "am_cutoff_calibrated": base.AM_CUTOFF, "am_cutoff_generic": base.AM_GENERIC,
        "am_over_calibrated": len(over), "am_over_calibrated_snv_reachable": len(over_snv),
        "am_over_generic": sum(1 for r in rows if r["am_over_generic_cutoff"]),
        "snv_reachable": sum(1 for r in rows if r["snv_accessible"]),
        "gnomad_status": dict(Counter(r["gnomad_status"] for r in rows)),
        "tiers_all": dict(Counter(r["tier"] for r in rows)),
        "tiers_snv_reachable": dict(Counter(r["tier"] for r in rows if r["snv_accessible"])),
        "checks": [{"check": n, "passed": ok, "detail": d} for n, ok, d in checks],
    }
    core.save_json("", f"KIF13A_saturation_summary_{core.today()}.json", summary,
                   base=ROOT / "data" / "processed", source="derived",
                   tool="scripts/build_saturation_table.py")
    return 0 if all(ok for _n, ok, _d in checks) else 1


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
