"""Regional missense constraint: is the KIF13A motor domain constrained once measured locally?

The project's strongest counter-argument is gene-level: KIF13A missense Z = 2.79, observed/expected
missense = 0.861. But that is an average over 1805 amino acids of which only 352 are the motor domain.
A tolerant 1285-residue tail can swamp a constrained motor domain in the average. This script measures
missense constraint *inside* the motor domain.

Method (Samocha 2014 framework, PMID 25086666, calibrated on the gene's own synonymous variants):

  1. MANE exon structure and the plus-strand genomic sequence from Ensembl REST, turned into a CDS
     where every coding base knows its genomic neighbours (so the trinucleotide context at an
     exon-exon junction is the genomic one, not a spliced-together fake).
  2. Every possible single-nucleotide change in the CDS is enumerated and classified as synonymous,
     missense or nonsense from the standard genetic code.
  3. Each possible change gets a trinucleotide-context mutation rate from the Samocha 2014 table.
  4. The rate is turned into an expected *probability of being seen in gnomAD* by a calibration fitted
     on SYNONYMOUS sites only, pooled across the whole gene panel: possible synonymous changes are
     binned by mutation rate and the observed fraction per bin is the expected probability. This is
     what makes sequencing depth, sample size and demographic history cancel out, and it handles the
     saturation of high-rate CpG contexts in a 730k-exome callset, which a linear rate model does not.
  5. Observed = gnomAD v4.1 exome SNVs on the MANE transcript, PASS, AF <= 0.1 %.
  6. observed/expected with a 90 % Poisson confidence interval, per domain and in sliding windows.

Sanity checks the script runs and prints (see docs/bolgesel_kisitlilik.md):
  - KIF1A motor domain must come out strongly constrained (positive control)
  - synonymous observed/expected must be near 1.0 in every window (calibration check)
  - the gene-level KIF13A number must reproduce gnomAD's published 1901/2209 = 0.861
  - the motor-domain observed missense count is reported so a wide interval cannot be sold as a point

Everything downloaded is cached under data/raw with a provenance sidecar; nothing is re-downloaded.

    python scripts\\regional_constraint.py
    python scripts\\regional_constraint.py --refresh --no-figure

NOTE on the interpreter: the project venv (mcp\\atlas_server\\.venv) is currently refused by Windows
Application Control ("Uygulama Denetimi ilkesi bu dosyayi engelledi") because uv's downloaded CPython
is blocked, so this script was run with C:\\Users\\pc\\anaconda3\\python.exe. It only needs
httpx, pyyaml and numpy, so either interpreter works once the venv is usable again.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL = "https://rest.ensembl.org"
# Samocha et al. 2014, Nat Genet 46:944-950, PMID 25086666, DOI 10.1038/ng.3050.
# The 192-row trinucleotide table (64 contexts x 3 alternates) redistributed by denovonear (MIT).
RATES_URL = "https://raw.githubusercontent.com/jeremymcrae/denovonear/master/src/denovonear/data/rates.txt"

# KIF13A domains in MANE numbering, the project's authoritative table (config/targets.yaml, HANDOFF.md)
KIF13A_DOMAINS = [(5, 352, "motor"), (353, 380, "boyun baglayici"), (381, 468, "CC1"),
                  (469, 519, "FHA"), (520, 1805, "kuyruk")]

# The comparison frame. KIF13A first, then the paralogs the brief names, then the wider kinesin panel
# so that "motor domains are constrained in general" can be separated from "this motor domain".
# AUTOSOMES ONLY. KIF4A sits on chrX and a first run with it in the panel poisoned the calibration:
# hemizygous males give the X a different allele number, its synonymous observed/expected came out at
# 0.56, and that dragged every other gene's expectation down by several per cent.
PANEL = ["KIF13A", "KIF1A", "KIF13B", "KIF1B", "KIF1C",
         "KIF5A", "KIF5C", "KIF3A", "KIF21A", "KIF16B", "KIF14"]
AUTOSOMES = {str(i) for i in range(1, 23)}

MAX_AF = 0.001          # gnomAD's constraint pipeline counts variants below 0.1 %
N_RATE_BINS = 24        # calibration bins over the mutation rate
WINDOW, STEP = 100, 10  # sliding window in codons
FINE_WINDOW, FINE_STEP = 50, 5

CODON_SNV = re.compile(r"^c\.(\d+)([ACGT])>([ACGT])$")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

GENETIC_CODE = {}
for _i, _b1 in enumerate("TCAG"):
    for _b2 in "TCAG":
        for _b3 in "TCAG":
            GENETIC_CODE[_b1 + _b2 + _b3] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[
                len(GENETIC_CODE)]


# --------------------------------------------------------------------------- small helpers


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def cached(subdir: str, pattern: str, base: Path = core.RAW) -> Path | None:
    folder = base / subdir
    if not folder.exists():
        return None
    hits = sorted(p for p in folder.glob(pattern) if not p.name.endswith(".prov.json"))
    return hits[-1] if hits else None


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


async def gnomad(query: str) -> dict:
    data = await core.post_json(GNOMAD_API, {"query": query}, timeout=600)
    if data.get("errors"):
        raise RuntimeError(f"gnomAD GraphQL error: {str(data['errors'])[:400]}")
    return data["data"]


# --------------------------------------------------------------------------- 1. mutation rates


async def mutation_rates(refresh: bool, symmetrise: bool = True) -> dict[tuple[str, str], float]:
    """Samocha 2014 trinucleotide mutation rate table: (context, alt base) -> per-generation rate.

    Germline mutation is strand symmetric, so a context and its reverse complement should carry the
    same rate; the published table's pairs differ by up to ~13 % because each was estimated
    separately from finite intergenic 1000 Genomes data. Averaging the two removes that estimation
    noise and, more importantly, removes the arbitrary choice of which strand to read the gene on
    (KIF13A, KIF1A and KIF13B are all minus strand, KIF1B and KIF1C plus strand).
    """
    hit = None if refresh else cached("mutation_rates", "samocha2014_trinucleotide_rates_*.txt")
    if hit:
        text = hit.read_text(encoding="utf-8")
        print(f"  using cached {core.rel(hit)}")
    else:
        raw = await core.get_bytes(RATES_URL)
        rec = core.save_bytes("mutation_rates", f"samocha2014_trinucleotide_rates_{core.today()}.txt", raw,
                              source=("Samocha et al. 2014 Nat Genet 46:944-950 (PMID 25086666, "
                                      "doi 10.1038/ng.3050) trinucleotide mutation rate table, "
                                      "redistributed in denovonear (MIT)"),
                              urls=[RATES_URL], tool="scripts/regional_constraint.py", params={})
        print(f"  downloaded {rec['path']}")
        text = raw.decode("utf-8")

    rates: dict[tuple[str, str], float] = {}
    for row in csv.DictReader(text.splitlines(), delimiter=" " if " " in text.splitlines()[0] else "\t"):
        frm, to = row["from"].strip().upper(), row["to"].strip().upper()
        if len(frm) != 3 or len(to) != 3 or frm[0] != to[0] or frm[2] != to[2] or frm[1] == to[1]:
            raise RuntimeError(f"unexpected rate row {row}")
        rates[(frm, to[1])] = float(row["mu_snp"])
    if len(rates) != 192:
        raise RuntimeError(f"expected 192 rate rows, got {len(rates)}")
    devs = [abs(rates[(c, a)] - rates[(revcomp(c), a.translate(COMP))]) / rates[(c, a)]
            for c, a in rates]
    print(f"  192 context rates; strand-symmetry deviation median {np.median(devs):.3f}, "
          f"max {max(devs):.3f}")
    if symmetrise:
        rates = {(c, a): 0.5 * (rates[(c, a)] + rates[(revcomp(c), a.translate(COMP))])
                 for c, a in rates}
        print("  rates symmetrised over reverse-complement pairs")
    return rates


# --------------------------------------------------------------------------- 2. sequence and structure


async def transcript_structure(enst: str, refresh: bool) -> dict:
    hit = None if refresh else cached("ensembl", f"{enst}_structure_*.json")
    if hit:
        raw = json.loads(hit.read_text(encoding="utf-8"))
    else:
        url = f"{ENSEMBL}/lookup/id/{enst}?expand=1;utr=1"
        raw = await core.get_json(url, headers={"Content-Type": "application/json"})
        core.save_json("ensembl", f"{enst}_structure_{core.today()}.json", raw,
                       source="Ensembl REST lookup (expand, utr)", urls=[url],
                       tool="scripts/regional_constraint.py", params={"transcript": enst})
    tr = raw["Translation"]
    cds_lo, cds_hi = min(tr["start"], tr["end"]), max(tr["start"], tr["end"])
    strand = raw["strand"]
    exons = sorted(raw["Exon"], key=lambda e: min(e["start"], e["end"]), reverse=strand < 0)
    coding = []
    for i, ex in enumerate(exons, start=1):
        lo, hi = min(ex["start"], ex["end"]), max(ex["start"], ex["end"])
        c_lo, c_hi = max(lo, cds_lo), min(hi, cds_hi)
        if c_hi >= c_lo:
            coding.append({"exon_number": i, "genomic_start": c_lo, "genomic_end": c_hi,
                           "coding_length": c_hi - c_lo + 1})
    return {"transcript": f"{raw['id']}.{raw['version']}", "chrom": raw["seq_region_name"],
            "strand": strand, "protein_id": tr["id"], "protein_length": tr["length"],
            "gene_start": min(raw["start"], raw["end"]), "gene_end": max(raw["start"], raw["end"]),
            "n_exons": len(exons), "coding_exons": coding,
            "cds_length": sum(e["coding_length"] for e in coding)}


async def genomic_sequence(chrom: str, start: int, stop: int, tag: str, refresh: bool) -> str:
    """Plus-strand genomic sequence for the gene span, padded, so every CDS base has real neighbours."""
    hit = None if refresh else cached("ensembl", f"{tag}_genomic_*.txt")
    if hit:
        return hit.read_text(encoding="utf-8").strip().upper()
    url = f"{ENSEMBL}/sequence/region/human/{chrom}:{start}..{stop}:1"
    raw = await core.get_bytes(url, headers={"Content-Type": "text/plain", "Accept": "text/plain"})
    core.save_bytes("ensembl", f"{tag}_genomic_{core.today()}.txt", raw,
                    source="Ensembl REST sequence/region, GRCh38 plus strand", urls=[url],
                    tool="scripts/regional_constraint.py",
                    params={"region": f"{chrom}:{start}-{stop}", "strand": 1})
    return raw.decode("ascii").strip().upper()


def build_cds(struct: dict, plus_seq: str, seq_start: int) -> tuple[str, list[int], list[str]]:
    """CDS in transcript orientation, plus each base's genomic position and trinucleotide context.

    The context is taken from the genome, so a codon split across an exon-exon junction gets its real
    intronic neighbour rather than the next exon's first base.
    """
    def plus(pos: int) -> str:
        return plus_seq[pos - seq_start]

    minus = struct["strand"] < 0
    bases, gpos = [], []
    for ex in struct["coding_exons"]:
        rng = range(ex["genomic_end"], ex["genomic_start"] - 1, -1) if minus else \
            range(ex["genomic_start"], ex["genomic_end"] + 1)
        for p in rng:
            gpos.append(p)
            bases.append(plus(p).translate(COMP) if minus else plus(p))
    contexts = []
    for p in gpos:
        tri = plus(p - 1) + plus(p) + plus(p + 1)
        contexts.append(revcomp(tri) if minus else tri)
    return "".join(bases), gpos, contexts


# --------------------------------------------------------------------------- 3. gnomAD


async def gnomad_variants(symbol: str, enst: str, refresh: bool) -> list[dict]:
    """Every coding SNV on the MANE transcript, with exome and joint allele counts and filters."""
    hit = None if refresh else cached("gnomad", f"{symbol}_{enst}_gnomad_r4_codingsnv_*.tsv")
    if hit:
        print(f"    using cached {core.rel(hit)}")
        return table(hit)
    query = """
    { transcript(transcript_id: "%s", reference_genome: GRCh38) {
        transcript_id transcript_version
        variants(dataset: gnomad_r4) {
          variant_id pos ref alt hgvsc hgvsp consequence flags
          exome { ac an filters } genome { ac an filters } joint { ac an filters }
        } } }""" % enst
    tx = (await gnomad(query))["transcript"]
    rows = []
    for v in tx["variants"]:
        m = CODON_SNV.match((v.get("hgvsc") or "").split(":")[-1].strip())
        if not m:
            continue  # intronic, UTR, indel, or a multi-base HGVS we cannot place on a codon
        ex, ge, jo = v.get("exome") or {}, v.get("genome") or {}, v.get("joint") or {}
        rows.append({
            "variant_id": v["variant_id"], "hgvsc": v["hgvsc"], "hgvsp": v.get("hgvsp"),
            "consequence": v.get("consequence"), "flags": ";".join(v.get("flags") or []),
            "cds_pos": int(m.group(1)), "cds_ref": m.group(2), "cds_alt": m.group(3),
            "exome_ac": ex.get("ac") or 0, "exome_an": ex.get("an") or 0,
            "exome_filters": ";".join(ex.get("filters") or []),
            "genome_ac": ge.get("ac") or 0, "genome_an": ge.get("an") or 0,
            "genome_filters": ";".join(ge.get("filters") or []),
            "joint_ac": jo.get("ac") or 0, "joint_an": jo.get("an") or 0,
            "joint_filters": ";".join(jo.get("filters") or []),
        })
    rec = core.save_tsv("gnomad", f"{symbol}_{enst}_gnomad_r4_codingsnv_{core.today()}.tsv", rows,
                        source="gnomAD v4.1 GraphQL transcript query, coding SNVs only",
                        urls=[GNOMAD_API], tool="scripts/regional_constraint.py",
                        params={"gene": symbol, "transcript": enst, "dataset": "gnomad_r4"})
    print(f"    downloaded {rec['path']} ({len(rows)} coding SNVs)")
    return table(cached("gnomad", f"{symbol}_{enst}_gnomad_r4_codingsnv_*.tsv"))


async def gnomad_coverage(symbol: str, struct: dict, refresh: bool) -> dict[int, float]:
    """Per-base gnomAD v4.1 exome coverage over the coding exons (the region query is not binned)."""
    hit = None if refresh else cached("gnomad", f"{symbol}_exome_coverage_*.tsv")
    if hit:
        print(f"    using cached {core.rel(hit)}")
        return {int(r["pos"]): float(r["over_20"] or 0) for r in table(hit)}
    rows = []
    for ex in struct["coding_exons"]:
        q = """{ region(chrom: "%s", start: %d, stop: %d, reference_genome: GRCh38) {
                 coverage(dataset: gnomad_r4) { exome { pos mean median over_20 } } } }""" % (
            struct["chrom"], ex["genomic_start"], ex["genomic_end"])
        cov = ((await gnomad(q))["region"]["coverage"] or {}).get("exome") or []
        for c in cov:
            if ex["genomic_start"] <= c["pos"] <= ex["genomic_end"]:
                rows.append({"pos": c["pos"], "mean": c["mean"], "median": c["median"],
                             "over_20": c["over_20"], "exon_number": ex["exon_number"]})
    rec = core.save_tsv("gnomad", f"{symbol}_exome_coverage_{core.today()}.tsv", rows,
                        source="gnomAD v4.1 GraphQL region coverage (exome), per coding exon",
                        urls=[GNOMAD_API], tool="scripts/regional_constraint.py",
                        params={"gene": symbol, "dataset": "gnomad_r4"})
    print(f"    downloaded {rec['path']} ({len(rows)} bases)")
    return {r["pos"]: r["over_20"] for r in rows}


CONSTRAINT_Q = """{ gene(gene_symbol: "%s", reference_genome: %s) {
  gene_id symbol mane_select_transcript { ensembl_id ensembl_version }
  gnomad_constraint { obs_mis exp_mis oe_mis oe_mis_lower oe_mis_upper obs_syn exp_syn oe_syn
                      obs_lof exp_lof oe_lof oe_lof_upper mis_z syn_z lof_z pli flags } } }"""


async def published_constraint(symbols: list[str], refresh: bool) -> dict:
    """gnomAD's own gene-level constraint. GRCh38 returns v4.1, GRCh37 returns v2.1.1 (project trap 3)."""
    hit = None if refresh else cached("gnomad", "regional_constraint_reference_*.json")
    if hit:
        out = json.loads(hit.read_text(encoding="utf-8"))
        if all(s in out.get("GRCh38", {}) for s in symbols):
            print(f"  using cached {core.rel(hit)}")
            return out
    out = {"GRCh38": {}, "GRCh37": {}}
    for build in ("GRCh38", "GRCh37"):
        for s in symbols:
            try:
                g = (await gnomad(CONSTRAINT_Q % (s, build)))["gene"]
            except Exception as exc:  # a gene may be absent from one build's canonical set
                out[build][s] = {"error": str(exc)[:120]}
                continue
            c = dict(g.get("gnomad_constraint") or {})
            c["mane_select_transcript"] = (g.get("mane_select_transcript") or {}).get("ensembl_id")
            c["mane_select_version"] = (g.get("mane_select_transcript") or {}).get("ensembl_version")
            out[build][s] = c
    rec = core.save_json("gnomad", f"regional_constraint_reference_{core.today()}.json", out,
                         source=("gnomAD GraphQL gene constraint; GRCh38 serves v4.1 and GRCh37 "
                                 "serves v2.1.1 constraint from the same endpoint"),
                         urls=[GNOMAD_API], tool="scripts/regional_constraint.py",
                         params={"genes": symbols})
    print(f"  downloaded {rec['path']}")
    return out


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TITLE_P = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
TITLE_TX = re.compile(r"(N[MP]_\d+\.\d+)")


async def clinvar_pathogenic_missense(gene: str, refresh: bool) -> list[dict]:
    """ClinVar germline P/LP missense records, with the protein position parsed from the title.

    Uses the Properties syntax, not `pathogenic[clinsig]`: Entrez silently rewrites the latter to
    an all-fields search and returns an inflated count (project trap 14).
    """
    hit = None if refresh else cached("clinvar", f"{gene}_clinvar_plp_missense_*.tsv")
    if hit:
        return table(hit)
    key = core.secret("NCBI_API_KEY")
    common = {"db": "clinvar", "retmode": "json", **({"api_key": key} if key else {})}
    term = (f'{gene}[gene] AND "missense variant"[molecular consequence] AND '
            f'("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])')
    search = await core.get_json(f"{EUTILS}/esearch.fcgi",
                                 params={**common, "term": term, "retmax": "1000"})
    ids = search["esearchresult"]["idlist"]
    rows = []
    for block in core.chunks(ids, 100):
        summ = await core.get_json(f"{EUTILS}/esummary.fcgi", params={**common, "id": ",".join(block)})
        res = summ["result"]
        for uid in res["uids"]:
            r = res[uid]
            title = r.get("title") or ""
            germ = r.get("germline_classification") or {}
            m, t = TITLE_P.search(title), TITLE_TX.search(title)
            rows.append({"gene": gene, "accession": r.get("accession"), "title": title,
                         "transcript_in_title": t.group(1) if t else "",
                         "protein_pos": int(m.group(2)) if m else "",
                         "classification": germ.get("description") or "",
                         "review_status": germ.get("review_status") or ""})
    core.save_tsv("clinvar", f"{gene}_clinvar_plp_missense_{core.today()}.tsv", rows,
                  source="NCBI E-utilities clinvar esearch+esummary, germline P/LP missense",
                  urls=[f"{EUTILS}/esearch.fcgi"], tool="scripts/regional_constraint.py",
                  params={"term": term})
    return rows


async def protein_domains(symbol: str, ensp: str, refresh: bool) -> list[dict]:
    """Pfam domains in the translation's own numbering, so no UniProt isoform offset can creep in."""
    hit = None if refresh else cached("ensembl", f"{symbol}_{ensp}_pfam_*.json")
    if hit:
        return json.loads(hit.read_text(encoding="utf-8"))
    url = f"{ENSEMBL}/overlap/translation/{ensp}?type=Pfam"
    raw = await core.get_json(url, headers={"Content-Type": "application/json"})
    core.save_json("ensembl", f"{symbol}_{ensp}_pfam_{core.today()}.json", raw,
                   source="Ensembl REST overlap/translation (Pfam), translation numbering",
                   urls=[url], tool="scripts/regional_constraint.py",
                   params={"gene": symbol, "translation": ensp})
    return raw


# --------------------------------------------------------------------------- 4. enumeration


def enumerate_snvs(gene: str, cds: str, gpos: list[int], contexts: list[str],
                   rates: dict, protein_length: int) -> list[dict]:
    """Every possible single-nucleotide change in the CDS, classified and rate-tagged."""
    rows = []
    n = len(cds)
    for i in range(n):                      # 0-based CDS index
        codon_no = i // 3 + 1
        if codon_no > protein_length:       # the terminal stop codon: not a missense/synonymous site
            continue
        frame = i % 3
        codon = cds[i - frame: i - frame + 3]
        if len(codon) != 3 or codon not in GENETIC_CODE:
            continue
        aa_ref = GENETIC_CODE[codon]
        ctx = contexts[i]
        if "N" in ctx or cds[i] != ctx[1]:
            raise RuntimeError(f"{gene}: context {ctx} does not match CDS base {cds[i]} at {i + 1}")
        for alt in "ACGT":
            if alt == cds[i]:
                continue
            new = codon[:frame] + alt + codon[frame + 1:]
            aa_alt = GENETIC_CODE[new]
            if codon_no == 1:
                klass = "start_lost" if aa_alt != aa_ref else "synonymous"
            elif aa_alt == aa_ref:
                klass = "synonymous"
            elif aa_alt == "*":
                klass = "nonsense"
            else:
                klass = "missense"
            rows.append({"gene": gene, "cds_pos": i + 1, "genomic_pos": gpos[i], "codon": codon_no,
                         "frame": frame, "ref": cds[i], "alt": alt, "context": ctx,
                         "aa_ref": aa_ref, "aa_alt": aa_alt, "klass": klass,
                         "mu": rates[(ctx, alt)],
                         "cpg": ctx[1:] == "CG" or ctx[:2] == "CG"})
    return rows


def attach_observed(gene: str, possible: list[dict], variants: list[dict], cds: str,
                    max_af: float) -> tuple[list[dict], dict]:
    """Mark which possible changes gnomAD has actually seen. Returns the rows and a QC dict."""
    index = {(int(r["cds_pos"]), r["alt"]): r for r in possible}
    for r in possible:
        r["observed"] = 0
        r["exome_ac"] = 0
        r["exome_af"] = 0.0
    qc = Counter()
    for v in variants:
        pos, ref, alt = int(v["cds_pos"]), v["cds_ref"], v["cds_alt"]
        if not (1 <= pos <= len(cds)):
            qc["outside_cds"] += 1
            continue
        if cds[pos - 1] != ref:
            qc["ref_mismatch"] += 1
            continue
        qc["ref_match"] += 1
        ac, an = int(v["exome_ac"] or 0), int(v["exome_an"] or 0)
        af = ac / an if an else 0.0
        filters = v["exome_filters"]
        if ac < 1:
            qc["exome_ac0"] += 1
            continue
        if filters:
            qc["exome_filtered"] += 1
            continue
        if af > max_af:
            qc["above_af_cutoff"] += 1
            continue
        row = index.get((pos, alt))
        if row is None:
            qc["not_in_enumeration"] += 1   # stop codon positions, which we do not enumerate
            continue
        row["observed"] = 1
        row["exome_ac"] = ac
        row["exome_af"] = af
        qc["counted"] += 1
    return possible, dict(qc)


# --------------------------------------------------------------------------- 5. calibration


def fit_calibration(all_rows: list[dict], n_bins: int) -> dict:
    """P(seen in gnomAD | mutation rate), fitted on synonymous sites pooled over the gene panel.

    Binning by rate rather than assuming expected = c * mu is what lets the model survive the
    saturation of high-rate CpG transitions in a 730k-exome callset.
    """
    syn = [r for r in all_rows if r["klass"] == "synonymous"]
    mus = np.array([r["mu"] for r in syn])
    obs = np.array([r["observed"] for r in syn], dtype=float)
    edges = np.unique(np.quantile(mus, np.linspace(0, 1, n_bins + 1)))
    edges[0], edges[-1] = -np.inf, np.inf
    idx = np.clip(np.searchsorted(edges, mus, side="right") - 1, 0, len(edges) - 2)
    bins = []
    for b in range(len(edges) - 1):
        sel = idx == b
        if not sel.any():
            continue
        bins.append({"bin": b, "mu_lo": float(mus[sel].min()), "mu_hi": float(mus[sel].max()),
                     "n_possible_syn": int(sel.sum()), "n_observed_syn": int(obs[sel].sum()),
                     "p_observed": float(obs[sel].mean())})
    lookup = {b["bin"]: b["p_observed"] for b in bins}
    return {"edges": edges.tolist(), "bins": bins, "lookup": lookup,
            "n_syn_possible": len(syn), "n_syn_observed": int(obs.sum())}


def apply_calibration(rows: list[dict], cal: dict) -> None:
    edges = np.array(cal["edges"])
    mus = np.array([r["mu"] for r in rows])
    idx = np.clip(np.searchsorted(edges, mus, side="right") - 1, 0, len(edges) - 2)
    for r, b in zip(rows, idx):
        r["p_exp"] = cal["lookup"].get(int(b), 0.0)


# --------------------------------------------------------------------------- 6. statistics


_OE_GRID = np.arange(0.0, 3.0005, 0.001)


def poisson_oe_ci(obs: int, exp: float, ci: float = 0.90) -> tuple[float, float]:
    """gnomAD's LOEUF interval, applied to missense: the Poisson likelihood over a grid of o/e."""
    if exp <= 0:
        return (float("nan"), float("nan"))
    lam = exp * _OE_GRID
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = -lam + obs * np.log(np.where(lam > 0, lam, 1e-300))
    w = np.exp(ll - ll.max())
    w /= w.sum()
    c = np.cumsum(w)
    lo = float(_OE_GRID[np.searchsorted(c, (1 - ci) / 2)])
    hi = float(_OE_GRID[min(np.searchsorted(c, 1 - (1 - ci) / 2), len(_OE_GRID) - 1)])
    return lo, hi


def rate_ratio_test(obs_a: int, exp_a: float, obs_b: int, exp_b: float) -> tuple[float, float]:
    """Is region A depleted relative to region B? Conditional binomial test on the counts.

    Given obs_a + obs_b variants, the share falling in A is Binomial with p = exp_a/(exp_a+exp_b)
    under the null that both regions carry the same observed/expected. This is scale free, so it
    survives any error in the absolute calibration as long as the error hits both regions equally.
    """
    n = obs_a + obs_b
    if n == 0 or exp_a <= 0 or exp_b <= 0:
        return (float("nan"), float("nan"))
    p0 = exp_a / (exp_a + exp_b)
    ratio = ((obs_a / exp_a) / (obs_b / exp_b)) if obs_b else float("inf")
    mean, sd = n * p0, math.sqrt(n * p0 * (1 - p0))
    if sd == 0:
        return (ratio, float("nan"))
    z = (abs(obs_a - mean) - 0.5) / sd
    return (ratio, max(0.0, min(1.0, math.erfc(max(z, 0.0) / math.sqrt(2)))))


def poisson_p_depletion(obs: int, exp: float) -> float:
    """One-sided P(X <= obs) under Poisson(exp): how surprising the depletion is."""
    if exp <= 0:
        return float("nan")
    if exp > 400:  # the exact sum underflows; the normal approximation is exact enough this far out
        z = (obs + 0.5 - exp) / math.sqrt(exp)
        return 0.5 * math.erfc(-z / math.sqrt(2))
    total, term = 0.0, math.exp(-exp)
    for k in range(0, obs + 1):
        if k:
            term *= exp / k
        total += term
    return min(1.0, total)


def region_stats(rows: list[dict], lo: int, hi: int, label: str, gene: str) -> dict:
    sel = [r for r in rows if lo <= r["codon"] <= hi]
    out = {"gene": gene, "region": label, "aa_start": lo, "aa_end": hi, "n_aa": hi - lo + 1}
    for klass, tag in (("missense", "mis"), ("synonymous", "syn"), ("nonsense", "non")):
        k = [r for r in sel if r["klass"] == klass]
        obs = sum(r["observed"] for r in k)
        exp = sum(r["p_exp"] for r in k)
        oe = obs / exp if exp > 0 else float("nan")
        lo_ci, hi_ci = poisson_oe_ci(obs, exp)
        out[f"possible_{tag}"] = len(k)
        out[f"obs_{tag}"] = obs
        out[f"exp_{tag}"] = round(exp, 2)
        out[f"oe_{tag}"] = round(oe, 4) if exp > 0 else ""
        out[f"oe_{tag}_lower"] = round(lo_ci, 4)
        out[f"oe_{tag}_upper"] = round(hi_ci, 4)
    # missense expectation re-scaled by this region's own synonymous o/e: immune to any local
    # calibration bias (coverage, CpG methylation) that hits synonymous and missense alike
    syn_oe = out["oe_syn"] if out["oe_syn"] != "" else 0
    if syn_oe:
        adj_exp = out["exp_mis"] * float(syn_oe)
        out["oe_mis_syn_normalised"] = round(out["obs_mis"] / adj_exp, 4) if adj_exp > 0 else ""
        lo_a, hi_a = poisson_oe_ci(out["obs_mis"], adj_exp)
        out["oe_mis_syn_normalised_lower"] = round(lo_a, 4)
        out["oe_mis_syn_normalised_upper"] = round(hi_a, 4)
    else:
        out["oe_mis_syn_normalised"] = ""
        out["oe_mis_syn_normalised_lower"] = ""
        out["oe_mis_syn_normalised_upper"] = ""
    out["p_missense_depletion"] = f"{poisson_p_depletion(out['obs_mis'], out['exp_mis']):.3g}"
    return out


def sliding(rows: list[dict], protein_length: int, gene: str, window: int, step: int) -> list[dict]:
    out = []
    for start in range(1, protein_length - window + 2, step):
        out.append(region_stats(rows, start, start + window - 1,
                                f"w{start}-{start + window - 1}", gene))
    return out


def grade(oe_upper: float) -> str:
    """The verdict scale the brief asks for, keyed on the upper bound so a wide interval cannot pass."""
    if oe_upper >= 0.95:
        return "hayir, tolere ediyor"
    if oe_upper >= 0.80:
        return "zayif kisitli"
    if oe_upper >= 0.60:
        return "orta kisitli"
    return "guclu kisitli"


# --------------------------------------------------------------------------- 7. figure

FIG_W, FIG_H = 1180, 1120
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
C_MIS, C_SYN, C_BAND, C_DOM = "#c2362f", "#3f7fb8", "#f0c9c5", "#cfe3f7"


def build_figure(profiles: dict, struct_rows: list[dict], motor_rows: list[dict], gene_level: dict,
                 path: Path) -> None:
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{FIG_H}" '
         f'viewBox="0 0 {FIG_W} {FIG_H}" font-family="Arial, Helvetica, sans-serif">',
         f'<rect width="{FIG_W}" height="{FIG_H}" fill="#ffffff"/>',
         f'<text x="40" y="28" font-size="15" fill="{INK}">KIF13A bolgesel missense kisitliligi: '
         f'gen duzeyindeki 0,86 ortalamasi motor domaini gizliyor mu</text>',
         f'<text x="40" y="46" font-size="11" fill="{MUTED}">gnomAD v4.1 ekzom, MANE transkripti; '
         f'beklenti Samocha 2014 trinukleotid hizlarindan, panelin sinonim varyantlariyla kalibre '
         f'edildi. Bant = %90 Poisson araligi.</text>']

    L, R = 70, 250
    plotw = FIG_W - L - R
    panels = [("KIF13A", 80, 200), ("KIF1A", 320, 200)]

    for gene, top, h in panels:
        prof = profiles[gene]
        length = max(p["aa_end"] for p in prof)
        ymin, ymax = 0.0, 1.6

        def gx(aa):
            return L + plotw * (aa - 1) / (length - 1)

        def gy(v):
            return top + h - h * (min(max(v, ymin), ymax) - ymin) / (ymax - ymin)

        o.append(f'<text x="{L}" y="{top - 12}" font-size="12" fill="{INK}">'
                 f'{"A" if gene == "KIF13A" else "B"}. {gene} '
                 f'{"(hedef)" if gene == "KIF13A" else "(pozitif kontrol)"} kayan pencere, '
                 f'{WINDOW} kodon, {STEP} kodon adim</text>')
        o.append(f'<rect x="{L}" y="{top}" width="{plotw}" height="{h}" fill="#fcfcfa" stroke="{GRID}"/>')
        # domain shading, KIF13A only (the project has an authoritative domain table for it)
        if gene == "KIF13A":
            for lo, hi, name in KIF13A_DOMAINS[:4]:
                o.append(f'<rect x="{gx(lo):.1f}" y="{top}" width="{max(1.0, gx(hi) - gx(lo)):.1f}" '
                         f'height="{h}" fill="{C_DOM}" opacity="0.45"/>')
            o.append(f'<text x="{gx(178):.1f}" y="{top + 12}" font-size="10" fill="#2f6ea8" '
                     f'text-anchor="middle">motor 5-352</text>')
        else:
            o.append(f'<rect x="{gx(5):.1f}" y="{top}" width="{gx(354) - gx(5):.1f}" height="{h}" '
                     f'fill="{C_DOM}" opacity="0.45"/>')
            o.append(f'<text x="{gx(180):.1f}" y="{top + 12}" font-size="10" fill="#2f6ea8" '
                     f'text-anchor="middle">motor 5-354</text>')
        for t in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
            o.append(f'<line x1="{L}" y1="{gy(t):.1f}" x2="{L + plotw}" y2="{gy(t):.1f}" '
                     f'stroke="{GRID}" stroke-dasharray="{"0" if t == 1.0 else "2,3"}"/>')
            o.append(f'<text x="{L - 8}" y="{gy(t) + 3.5:.1f}" font-size="10" fill="{MUTED}" '
                     f'text-anchor="end">{t:.2f}</text>')
        # 90 % band for missense
        band = ([f'{gx(p["aa_start"] + WINDOW / 2):.1f},{gy(p["oe_mis_upper"]):.1f}' for p in prof]
                + [f'{gx(p["aa_start"] + WINDOW / 2):.1f},{gy(p["oe_mis_lower"]):.1f}'
                   for p in reversed(prof)])
        o.append(f'<polygon points="{" ".join(band)}" fill="{C_BAND}" opacity="0.75"/>')
        for key, col, wdt in (("oe_syn", C_SYN, 1.3), ("oe_mis", C_MIS, 2.0)):
            pts = " ".join(f'{gx(p["aa_start"] + WINDOW / 2):.1f},{gy(float(p[key])):.1f}'
                           for p in prof if p[key] != "")
            o.append(f'<polyline points="{pts}" fill="none" stroke="{col}" stroke-width="{wdt}"/>')
        gl = gene_level[gene]
        o.append(f'<line x1="{L}" y1="{gy(gl["oe_mis_published"]):.1f}" x2="{L + plotw}" '
                 f'y2="{gy(gl["oe_mis_published"]):.1f}" stroke="#8a8a84" stroke-width="1" '
                 f'stroke-dasharray="6,4"/>')
        o.append(f'<text x="{L + plotw - 4}" y="{gy(gl["oe_mis_published"]) - 5:.1f}" font-size="10" '
                 f'fill="{MUTED}" text-anchor="end">gen duzeyi gnomAD o/e = '
                 f'{gl["oe_mis_published"]:.3f}</text>')
        for tick in (1, length // 2, length):
            o.append(f'<text x="{gx(tick):.1f}" y="{top + h + 15}" font-size="10" fill="{MUTED}" '
                     f'text-anchor="middle">{tick}</text>')
        o.append(f'<text x="{L + plotw / 2:.1f}" y="{top + h + 30}" font-size="10" fill="{MUTED}" '
                 f'text-anchor="middle">aminoasit konumu ({length} aa)</text>')

    # legend
    lx, ly = L + plotw + 24, 96
    for label, col in (("missense gozlenen/beklenen", C_MIS), ("sinonim (kalibrasyon kontrolu)", C_SYN),
                       ("%90 Poisson araligi", C_BAND)):
        o.append(f'<rect x="{lx}" y="{ly - 8}" width="16" height="4" fill="{col}"/>')
        o.append(f'<text x="{lx + 22}" y="{ly}" font-size="10" fill="{MUTED}">{label}</text>')
        ly += 16

    # panel C: motor-domain missense o/e across the kinesin panel
    top, h = 600, 190
    o.append(f'<text x="{L}" y="{top - 12}" font-size="12" fill="{INK}">C. Motor domain missense '
             f'gozlenen/beklenen, kinezin paneli (nokta = tahmin, cizgi = %90 aralik). '
             f'"Motor domainler zaten kisitlidir" itirazinin cerceve olcumu.</text>')
    vals = sorted(motor_rows, key=lambda r: r["oe_mis"])
    xmin, xmax = 0.0, 1.35

    def cx(v):
        return L + plotw * (min(max(v, xmin), xmax) - xmin) / (xmax - xmin)

    for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2):
        o.append(f'<line x1="{cx(t):.1f}" y1="{top}" x2="{cx(t):.1f}" y2="{top + h}" '
                 f'stroke="{GRID}" stroke-dasharray="{"0" if t == 1.0 else "2,3"}"/>')
        o.append(f'<text x="{cx(t):.1f}" y="{top + h + 14}" font-size="10" fill="{MUTED}" '
                 f'text-anchor="middle">{t:.1f}</text>')
    stepy = h / (len(vals) + 1)
    for i, r in enumerate(vals):
        yy = top + stepy * (i + 1)
        is_t = r["gene"] in ("KIF13A", "KIF1A")
        col = C_MIS if r["gene"] == "KIF13A" else ("#2a78d6" if r["gene"] == "KIF1A" else MUTED)
        o.append(f'<line x1="{cx(r["oe_mis_lower"]):.1f}" y1="{yy:.1f}" '
                 f'x2="{cx(r["oe_mis_upper"]):.1f}" y2="{yy:.1f}" stroke="{col}" stroke-width="1.4"/>')
        o.append(f'<circle cx="{cx(r["oe_mis"]):.1f}" cy="{yy:.1f}" r="{4 if is_t else 3}" fill="{col}"/>')
        o.append(f'<text x="{L - 8}" y="{yy + 3.5:.1f}" font-size="10" fill="{INK if is_t else MUTED}" '
                 f'text-anchor="end">{r["gene"]}</text>')
        o.append(f'<text x="{cx(r["oe_mis_upper"]) + 8:.1f}" y="{yy + 3.5:.1f}" font-size="9" '
                 f'fill="{MUTED}">{r["oe_mis"]:.2f} ({r["oe_mis_lower"]:.2f}-{r["oe_mis_upper"]:.2f}), '
                 f'{r["obs_mis"]} gozlem, ClinVar P/LP '
                 f'{r.get("clinvar_plp_missense_in_motor", "?")}</text>')
    o.append(f'<text x="{L + plotw / 2:.1f}" y="{top + h + 30}" font-size="10" fill="{MUTED}" '
             f'text-anchor="middle">motor domain missense gozlenen/beklenen '
             f'(kucuk = kisitlilik daha guclu); "ClinVar P/LP" = o motor domainde siniflandirilmis '
             f'patojenik/olasi patojenik missense sayisi</text>')

    # panel D: the same structural classes in KIF13A and in KIF1A. This is the answer to
    # "motor domains are compact and buried, so of course they look constrained".
    top, h = 860, 200
    o.append(f'<text x="{L}" y="{top - 26}" font-size="12" fill="{INK}">D. Motor domainin ICINDE '
             f'yapisal sinifa gore kisitlilik (1-380). Ayni siniflar superpozisyon eslesmesiyle '
             f'KIF1A\'ya tasindi.</text>')
    o.append(f'<text x="{L}" y="{top - 11}" font-size="10" fill="{MUTED}">Gradyan ayni yonde, ama '
             f'KIF13A her sinifta 2-4 kat daha sig. KIF13A\'nin motor yuzeyi hic kisitli degil; '
             f'KIF1A\'ninki kisitli.</text>')
    order = ["yuzey", "gomulu cekirdek", "mikrotubul arayuzu", "ATP cebi"]
    byk = {(r["gene"], r["structural_class"]): r for r in struct_rows}
    xmin, xmax = 0.0, 1.25
    LD = 190                      # panel D needs room for the class names on the left
    cw = FIG_W - LD - 260

    def dx(v):
        return LD + cw * (min(max(v, xmin), xmax) - xmin) / (xmax - xmin)

    for t in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25):
        o.append(f'<line x1="{dx(t):.1f}" y1="{top}" x2="{dx(t):.1f}" y2="{top + h}" stroke="{GRID}" '
                 f'stroke-dasharray="{"0" if t == 1.0 else "2,3"}"/>')
        o.append(f'<text x="{dx(t):.1f}" y="{top + h + 14}" font-size="10" fill="{MUTED}" '
                 f'text-anchor="middle">{t:.2f}</text>')
    rowh = h / (len(order) + 0.5)
    for i, cname in enumerate(order):
        yy = top + rowh * (i + 0.6)
        r13, r1a = byk.get(("KIF13A", cname)), byk.get(("KIF1A", cname))
        o.append(f'<text x="{LD - 10}" y="{yy + 4:.1f}" font-size="10" fill="{INK}" '
                 f'text-anchor="end">{cname}</text>')
        if r13:
            o.append(f'<text x="{LD - 10}" y="{yy + 18:.1f}" font-size="9" fill="{MUTED}" '
                     f'text-anchor="end">{r13["n_residues"]} kalinti, katman A '
                     f'{r13["tier_a_candidates"]}</text>')
        for j, (rr, col, lbl) in enumerate(((r13, C_MIS, "KIF13A"), (r1a, "#2a78d6", "KIF1A"))):
            if not rr or rr["oe_mis"] == "":
                continue
            yb = yy + j * 17
            o.append(f'<line x1="{dx(rr["oe_mis_lower"]):.1f}" y1="{yb:.1f}" '
                     f'x2="{dx(rr["oe_mis_upper"]):.1f}" y2="{yb:.1f}" stroke="{col}" '
                     f'stroke-width="1.4"/>')
            o.append(f'<circle cx="{dx(rr["oe_mis"]):.1f}" cy="{yb:.1f}" r="4" fill="{col}"/>')
            o.append(f'<text x="{dx(rr["oe_mis_upper"]) + 8:.1f}" y="{yb + 3.5:.1f}" font-size="9" '
                     f'fill="{col}">{lbl} {rr["oe_mis"]:.3f} '
                     f'({rr["oe_mis_lower"]:.2f}-{rr["oe_mis_upper"]:.2f}), {rr["obs_mis"]} gozlem</text>')
    o.append(f'<text x="{LD + cw / 2:.1f}" y="{top + h + 30}" font-size="10" fill="{MUTED}" '
             f'text-anchor="middle">missense gozlenen/beklenen; yapisal siniflar '
             f'docs/yapi_ekseni.md kalinti tablosundan</text>')

    o.append(f'<text x="40" y="{FIG_H - 12}" font-size="10" fill="{MUTED}">Kaynak: gnomAD v4.1 '
             f'(GraphQL transkript sorgusu, GRCh38), Ensembl REST MANE ekzon yapisi ve GRCh38 dizisi, '
             f'Samocha 2014 trinukleotid mutasyon hizi tablosu (PMID 25086666). '
             f'Olcum tarihi {core.today()}. scripts/regional_constraint.py</text>')
    o.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(o), encoding="utf-8")
    print(f"  wrote {core.rel(path)} ({path.stat().st_size} bytes)")


# --------------------------------------------------------------------------- main


async def resolve_panel(symbols: list[str], refresh: bool) -> dict[str, str]:
    """MANE transcript per gene: config/targets.yaml where we verified it, gnomAD otherwise."""
    out = {}
    known = core.load_targets()["genes"]
    for s in symbols:
        if s in known:
            out[s] = core.unversioned(known[s]["mane_transcript"])
    missing = [s for s in symbols if s not in out]
    if missing:
        hit = None if refresh else cached("gnomad", "panel_mane_transcripts_*.json")
        extra = json.loads(hit.read_text(encoding="utf-8")) if hit else {}
        for s in missing:
            if s not in extra:
                g = (await gnomad(CONSTRAINT_Q % (s, "GRCh38")))["gene"]
                extra[s] = (g.get("mane_select_transcript") or {}).get("ensembl_id")
        core.save_json("gnomad", f"panel_mane_transcripts_{core.today()}.json", extra,
                       source="gnomAD GraphQL gene query, mane_select_transcript", urls=[GNOMAD_API],
                       tool="scripts/regional_constraint.py", params={"genes": missing})
        out.update({s: extra[s] for s in missing if extra.get(s)})
    return out


async def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=",".join(PANEL))
    ap.add_argument("--coverage-genes", default="KIF13A,KIF1A")
    ap.add_argument("--max-af", type=float, default=MAX_AF)
    ap.add_argument("--bins", type=int, default=N_RATE_BINS)
    ap.add_argument("--raw-rates", action="store_true",
                    help="do not symmetrise the mutation rate table over reverse-complement pairs")
    ap.add_argument("--tag", default="", help="suffix for the output files of a sensitivity run")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    core.load_dotenv()
    stamp = core.today() + (f"_{args.tag}" if args.tag else "")
    panel = [g.strip() for g in args.panel.split(",") if g.strip()]
    if "KIF13A" not in panel or "KIF1A" not in panel:
        raise SystemExit("the panel must contain KIF13A (target) and KIF1A (positive control)")
    cov_genes = {g.strip() for g in args.coverage_genes.split(",") if g.strip()}
    outdir = "regional_constraint"

    print("== mutation rate table")
    rates = await mutation_rates(args.refresh, symmetrise=not args.raw_rates)

    print("\n== gnomAD published constraint (GRCh38 = v4.1, GRCh37 = v2.1.1)")
    published = await published_constraint(panel, args.refresh)

    print("\n== resolving MANE transcripts")
    tx = await resolve_panel(panel, args.refresh)
    print("  " + ", ".join(f"{g}:{t}" for g, t in tx.items()))

    genes: dict[str, dict] = {}
    all_rows: list[dict] = []
    for gene in panel:
        enst = tx.get(gene)
        if not enst:
            print(f"\n{gene}: no MANE transcript, skipped")
            continue
        print(f"\n== {gene} ({enst})")
        struct = await transcript_structure(enst, args.refresh)
        if struct["chrom"] not in AUTOSOMES:
            raise RuntimeError(f"{gene} is on chr{struct['chrom']}: sex chromosomes have a different "
                               f"allele number and break the shared synonymous calibration")
        known = core.load_targets()["genes"].get(gene)
        if known and struct["protein_length"] != known["protein_length"]:
            raise RuntimeError(f"{gene}: Ensembl {struct['protein_length']} aa vs targets.yaml "
                               f"{known['protein_length']} aa")
        pad = 5
        seq = await genomic_sequence(struct["chrom"], struct["gene_start"] - pad,
                                     struct["gene_end"] + pad, f"{gene}_{enst}", args.refresh)
        cds, gpos, contexts = build_cds(struct, seq, struct["gene_start"] - pad)
        protein = "".join(GENETIC_CODE[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))
        if not (cds.startswith("ATG") and protein.endswith("*")
                and len(protein) == struct["protein_length"] + 1
                and "*" not in protein[:-1]):
            raise RuntimeError(f"{gene}: CDS does not translate cleanly "
                               f"({len(cds)} nt -> {len(protein)} aa, start {cds[:3]})")
        print(f"  CDS {len(cds)} nt -> {struct['protein_length']} aa, clean ORF, "
              f"{len(struct['coding_exons'])} coding exons, strand {struct['strand']}")

        variants = await gnomad_variants(gene, enst, args.refresh)
        possible = enumerate_snvs(gene, cds, gpos, contexts, rates, struct["protein_length"])
        possible, qc = attach_observed(gene, possible, variants, cds, args.max_af)
        print(f"  {len(possible)} possible SNVs, {qc.get('counted', 0)} observed and counted; QC {qc}")
        if qc.get("ref_mismatch"):
            raise RuntimeError(f"{gene}: {qc['ref_mismatch']} gnomAD HGVS reference bases disagree "
                               f"with the Ensembl CDS")

        cov = await gnomad_coverage(gene, struct, args.refresh) if gene in cov_genes else {}
        if cov:
            vals = np.array([cov.get(p, 0.0) for p in gpos], dtype=float)
            print(f"  exome coverage over 20x: median {np.median(vals):.3f}, "
                  f"{(vals < 0.8).sum()} of {len(vals)} CDS bases below 0.80")
            for r in possible:
                r["cov_over_20"] = cov.get(r["genomic_pos"], "")

        genes[gene] = {"struct": struct, "rows": possible, "qc": qc, "cds": cds,
                       "coverage": cov, "protein_length": struct["protein_length"]}
        all_rows.extend(possible)

    print(f"\n== calibration on {sum(1 for r in all_rows if r['klass'] == 'synonymous')} possible "
          f"synonymous changes pooled over {len(genes)} genes")
    cal = fit_calibration(all_rows, args.bins)
    apply_calibration(all_rows, cal)
    print(f"  {cal['n_syn_observed']} / {cal['n_syn_possible']} synonymous seen "
          f"({cal['n_syn_observed'] / cal['n_syn_possible']:.3f} overall)")
    for b in cal["bins"]:
        print(f"    mu {b['mu_lo']:.3e}-{b['mu_hi']:.3e}  n={b['n_possible_syn']:6d}  "
              f"P(gorulme)={b['p_observed']:.3f}")
    core.save_tsv(outdir, f"calibration_bins_{stamp}.tsv", cal["bins"], base=core.INTERIM,
                  source="derived: proportion of possible synonymous SNVs observed in gnomAD v4.1, "
                         "binned by Samocha 2014 trinucleotide mutation rate",
                  urls=[GNOMAD_API, RATES_URL], tool="scripts/regional_constraint.py",
                  params={"genes": list(genes), "n_bins": args.bins, "max_af": args.max_af})

    # ---------------- gene level, the reproduction check
    gene_level = {}
    for gene, g in genes.items():
        whole = region_stats(g["rows"], 1, g["protein_length"], "gen butunu", gene)
        pub = published["GRCh38"].get(gene, {})
        gene_level[gene] = {
            **whole,
            "obs_mis_published": pub.get("obs_mis"), "exp_mis_published": pub.get("exp_mis"),
            "oe_mis_published": pub.get("oe_mis"), "obs_syn_published": pub.get("obs_syn"),
            "exp_syn_published": pub.get("exp_syn"), "oe_syn_published": pub.get("oe_syn"),
            "mis_z_published": pub.get("mis_z"),
            "obs_mis_delta_pct": (round(100 * (whole["obs_mis"] - pub["obs_mis"]) / pub["obs_mis"], 1)
                                  if pub.get("obs_mis") else ""),
            "exp_mis_delta_pct": (round(100 * (whole["exp_mis"] - pub["exp_mis"]) / pub["exp_mis"], 1)
                                  if pub.get("exp_mis") else ""),
            "oe_mis_delta_pct": (round(100 * (whole["oe_mis"] - pub["oe_mis"]) / pub["oe_mis"], 1)
                                 if pub.get("oe_mis") else ""),
            "published_version": "gnomAD v4.1 (GraphQL, GRCh38)",
            "v2_oe_mis": (published["GRCh37"].get(gene) or {}).get("oe_mis"),
            "v2_mis_z": (published["GRCh37"].get(gene) or {}).get("mis_z"),
            "v2_transcript": (published["GRCh37"].get(gene) or {}).get("mane_select_transcript"),
        }
    print("\n== gene level: ours vs gnomAD v4.1 published")
    for gene, r in gene_level.items():
        print(f"  {gene:7s} obs_mis {r['obs_mis']:5d} (gnomAD {r['obs_mis_published']})  "
              f"exp_mis {r['exp_mis']:8.1f} (gnomAD {r['exp_mis_published']:.1f})  "
              f"o/e {r['oe_mis']:.3f} (gnomAD {r['oe_mis_published']:.3f})  "
              f"syn o/e {r['oe_syn']:.3f} (gnomAD {r['oe_syn_published']:.3f})")
    core.save_tsv(outdir, f"gene_level_reproduction_{stamp}.tsv", list(gene_level.values()),
                  base=core.INTERIM, source="derived: this pipeline vs gnomAD v4.1 published constraint",
                  urls=[GNOMAD_API], tool="scripts/regional_constraint.py", params={"max_af": args.max_af})

    # Our calibration forces the PANEL's synonymous observed/expected to 1.0, while gnomAD's
    # genome-wide model leaves these genes at a mean synonymous o/e of about 0.95. Our expectations
    # are therefore a few per cent lower than gnomAD's and our o/e a few per cent higher, uniformly.
    # This scalar puts every regional number back on gnomAD's published scale so the motor-domain
    # value can be read next to the familiar gene-level 0.861 without an apples-to-oranges jump.
    scale = {}
    for gene, r in gene_level.items():
        scale[gene] = (r["oe_mis_published"] / float(r["oe_mis"])
                       if r.get("oe_mis_published") and r["oe_mis"] else 1.0)
    print("  gnomAD olcegine cevirme katsayilari: "
          + ", ".join(f"{g} x{s:.3f}" for g, s in scale.items()))

    def anchor(rows: list[dict]) -> list[dict]:
        for r in rows:
            s = scale.get(r["gene"], 1.0)
            for key in ("oe_mis", "oe_mis_lower", "oe_mis_upper"):
                v = r.get(key)
                r[key + "_gnomad_scaled"] = round(float(v) * s, 4) if v not in ("", None) else ""
            r["gnomad_scale_factor"] = round(s, 4)
        return rows

    # ---------------- the 52 tier A candidates, so the verdict can be applied to them directly
    cand_path = ROOT / "data" / "processed" / "KIF13A_candidates_20260915.tsv"
    tierA = [r for r in table(cand_path) if r.get("tier") == "A"] if cand_path.exists() else []
    tierA_pos = sorted(int(r["protein_pos"]) for r in tierA)
    sfile_early = cached("structure", "KIF13A_motor_residue_structure_*.tsv", base=core.INTERIM)
    cls_early = ({int(r["protein_pos"]): r["structural_class"] for r in table(sfile_early)}
                 if sfile_early else {})
    tierA_struct = [cls_early.get(p, "(bilinmiyor)") for p in tierA_pos]
    print(f"\n== katman A adaylari: {len(tierA)} varyant, konum araligi "
          f"{tierA_pos[0] if tierA_pos else '-'}-{tierA_pos[-1] if tierA_pos else '-'}, "
          f"yapisal sinif dagilimi {dict(Counter(tierA_struct).most_common())}")

    # ---------------- per domain
    per_domain: dict[str, list[dict]] = {}
    kif13a_domains = []
    for lo, hi, name in KIF13A_DOMAINS:
        kif13a_domains.append(region_stats(genes["KIF13A"]["rows"], lo, hi, name, "KIF13A"))
    # every KIF13A domain tested against the tail, the within-gene contrast that needs no calibration
    tail_ref = kif13a_domains[-1]
    for r in kif13a_domains:
        if r["region"] == "kuyruk":
            r["ratio_vs_tail"], r["p_vs_tail"] = 1.0, ""
            continue
        ratio, p = rate_ratio_test(r["obs_mis"], r["exp_mis"], tail_ref["obs_mis"], tail_ref["exp_mis"])
        r["ratio_vs_tail"], r["p_vs_tail"] = round(ratio, 4), f"{p:.3g}"
    anchor(kif13a_domains)
    per_domain["KIF13A"] = kif13a_domains
    print("\n== KIF13A per domain")
    for r in kif13a_domains:
        print(f"  {r['region']:16s} {r['aa_start']:5d}-{r['aa_end']:<5d} "
              f"missense {r['obs_mis']:5d}/{r['exp_mis']:8.1f} o/e {r['oe_mis']:.3f} "
              f"[{r['oe_mis_lower']:.3f}-{r['oe_mis_upper']:.3f}]  "
              f"gnomAD olceginde {r['oe_mis_gnomad_scaled']:.3f}  "
              f"sinonim o/e {r['oe_syn']:.3f}  kuyruga karsi p={r['p_vs_tail']}  "
              f"-> {grade(r['oe_mis_upper'])}")

    # ---------------- inside the motor domain: sub-regions and structural classes
    sub_rows = []
    for lo, hi, label in [(5, 180, "motor N yarisi 5-180"), (181, 352, "motor C yarisi 181-352"),
                          (270, 360, "MT arayuzu / anahtar bolgesi 270-360"),
                          (1, 80, "ilk 80 kodon (5' CpG adasi bolgesi)"),
                          # the region MTR and CCR independently point at as KIF13A's only depleted
                          # stretch; measured here on gnomAD v4.1 so it can be read on the same scale
                          (353, 519, "boyun+CC1+FHA 353-519 (MTR/CCR'nin isaret ettigi)"),
                          (381, 519, "CC1+FHA 381-519")]:
        r = region_stats(genes["KIF13A"]["rows"], lo, hi, label, "KIF13A")
        ratio, p = rate_ratio_test(r["obs_mis"], r["exp_mis"], tail_ref["obs_mis"], tail_ref["exp_mis"])
        r["ratio_vs_tail"], r["p_vs_tail"] = round(ratio, 4), f"{p:.3g}"
        sub_rows.append(r)
    anchor(sub_rows)
    print("\n== KIF13A motor domain, sub-regions")
    for r in sub_rows:
        print(f"  {r['region']:38s} missense {r['obs_mis']:4d}/{r['exp_mis']:7.1f} "
              f"o/e {r['oe_mis']:.3f} [{r['oe_mis_lower']:.3f}-{r['oe_mis_upper']:.3f}]  "
              f"sinonim o/e {r['oe_syn']:.3f}  sinonim-normalize {r['oe_mis_syn_normalised']}  "
              f"kuyruga karsi p={r['p_vs_tail']}")

    # constraint per structural class, joined to the structure axis residue table (residues 1-380)
    # The same classes are carried over to KIF1A through the superposition's residue pairing, so the
    # generic objection ("buried residues are constrained in any globular domain") can be measured
    # rather than argued: if KIF13A's ATP pocket is constrained only as much as burial explains, the
    # KIF1A column will show the same ratio between classes, only deeper.
    struct_rows = []
    sfile = cached("structure", "KIF13A_motor_residue_structure_*.tsv", base=core.INTERIM)
    if sfile:
        sres = table(sfile)
        maps = {"KIF13A": {int(r["protein_pos"]): r["structural_class"] for r in sres},
                "KIF1A": {int(r["kif1a_pos"]): r["structural_class"] for r in sres
                          if (r.get("kif1a_pos") or "").strip().isdigit()}}
        print(f"\n== 1-380, constraint by structural class ({core.rel(sfile)}); "
              f"KIF1A siniflari ayni dosyanin superpozisyon eslesmesinden tasindi "
              f"({len(maps['KIF1A'])} kalinti)")
        for gene in ("KIF13A", "KIF1A"):
            cls = maps[gene]
            by_class = defaultdict(list)
            for r in genes[gene]["rows"]:
                c = cls.get(r["codon"])
                if c:
                    by_class[c].append(r)
            for c, rws in sorted(by_class.items(), key=lambda kv: -len(kv[1])):
                obs = sum(x["observed"] for x in rws if x["klass"] == "missense")
                exp = sum(x["p_exp"] for x in rws if x["klass"] == "missense")
                sobs = sum(x["observed"] for x in rws if x["klass"] == "synonymous")
                sexp = sum(x["p_exp"] for x in rws if x["klass"] == "synonymous")
                lo_ci, hi_ci = poisson_oe_ci(obs, exp)
                syn_oe = sobs / sexp if sexp else 0
                row = {"gene": gene, "structural_class": c,
                       "n_residues": sum(1 for v in cls.values() if v == c),
                       "obs_mis": obs, "exp_mis": round(exp, 2),
                       "oe_mis": round(obs / exp, 4) if exp else "",
                       "oe_mis_lower": round(lo_ci, 4), "oe_mis_upper": round(hi_ci, 4),
                       "obs_syn": sobs, "exp_syn": round(sexp, 2),
                       "oe_syn": round(syn_oe, 4) if sexp else "",
                       "oe_mis_syn_normalised": (round(obs / (exp * syn_oe), 4)
                                                 if exp and syn_oe else ""),
                       "p_missense_depletion": f"{poisson_p_depletion(obs, exp):.3g}",
                       "tier_a_candidates": (sum(1 for x in tierA_struct if x == c)
                                             if gene == "KIF13A" else "")}
                anchor([row])
                struct_rows.append(row)
                print(f"  {gene:7s} {c:42s} {row['n_residues']:4d} kalinti  "
                      f"missense {obs:4d}/{exp:7.1f} o/e {row['oe_mis']:.3f} "
                      f"[{lo_ci:.3f}-{hi_ci:.3f}]  sinonim {sobs}/{sexp:.1f}={row['oe_syn']}  "
                      f"sin-norm {row['oe_mis_syn_normalised']}  katman A {row['tier_a_candidates']}")
    core.save_tsv(outdir, f"KIF13A_motor_subregions_{stamp}.tsv", sub_rows + struct_rows,
                  base=core.INTERIM,
                  source="derived: regional missense constraint inside the KIF13A motor domain, "
                         "by sub-region and by structure-axis structural class",
                  urls=[GNOMAD_API], tool="scripts/regional_constraint.py", params={})

    # ---------------- the kinesin frame: motor domain of every panel gene
    motor_rows = []
    motor_bounds = {}
    for gene, g in genes.items():
        if gene == "KIF13A":
            lo, hi = 5, 352
            src = "proje domain tablosu (config/targets.yaml, HANDOFF.md)"
        elif gene == "KIF1A":
            lo, hi = 5, 354
            src = "proje domain tablosu"
        else:
            pfam = await protein_domains(gene, g["struct"]["protein_id"], args.refresh)
            kin = [d for d in pfam if "kinesin" in (d.get("description") or d.get("id") or "").lower()
                   or d.get("id") == "PF00225"]
            if not kin:
                continue
            k = max(kin, key=lambda d: d["end"] - d["start"])
            lo, hi = int(k["start"]), int(k["end"])
            src = f"Ensembl Pfam {k.get('id')} ceviri numaralandirmasinda"
        motor_bounds[gene] = (lo, hi, src)
        r = region_stats(g["rows"], lo, hi, "motor", gene)
        r["boundary_source"] = src
        r["protein_length"] = g["protein_length"]
        tail = region_stats(g["rows"], hi + 1, g["protein_length"], "motor sonrasi", gene)
        r["tail_oe_mis"] = tail["oe_mis"]
        r["tail_obs_mis"] = tail["obs_mis"]
        r["tail_exp_mis"] = tail["exp_mis"]
        r["tail_oe_syn"] = tail["oe_syn"]
        r["motor_minus_tail"] = round(float(r["oe_mis"]) - float(tail["oe_mis"]), 4)
        ratio, p = rate_ratio_test(r["obs_mis"], r["exp_mis"], tail["obs_mis"], tail["exp_mis"])
        r["ratio_motor_over_tail"], r["p_motor_vs_tail"] = round(ratio, 4), f"{p:.3g}"
        r["grade"] = grade(r["oe_mis_upper"])
        motor_rows.append(r)
    anchor(motor_rows)

    # Does the constraint measure track clinical reality? Count the ClinVar pathogenic missense
    # records that actually sit in each motor domain. This is what separates "motor domains are
    # constrained in general" from "the constrained ones are the disease genes".
    print("\n== ClinVar P/LP missense inside each motor domain")
    for r in motor_rows:
        recs = await clinvar_pathogenic_missense(r["gene"], args.refresh)
        placed = [x for x in recs if str(x.get("protein_pos") or "").isdigit()]
        inmotor = [x for x in placed if r["aa_start"] <= int(x["protein_pos"]) <= r["aa_end"]]
        mane_rs = (core.load_targets()["genes"].get(r["gene"]) or {}).get("refseq_transcript")
        r["clinvar_plp_missense_total"] = len(recs)
        r["clinvar_plp_missense_placed"] = len(placed)
        r["clinvar_plp_missense_in_motor"] = len(inmotor)
        r["clinvar_title_transcripts"] = ";".join(
            f"{t}:{n}" for t, n in Counter(x["transcript_in_title"] for x in placed).most_common(3))
        r["clinvar_title_matches_mane"] = bool(mane_rs and any(
            x["transcript_in_title"] == mane_rs for x in placed))
        print(f"  {r['gene']:7s} {len(recs):4d} P/LP missense, {len(inmotor):4d} motor domainde "
              f"({r['clinvar_title_transcripts']})")

    print("\n== kinesin frame: motor domain vs everything after it")
    for r in sorted(motor_rows, key=lambda x: x["oe_mis"]):
        print(f"  {r['gene']:7s} motor {r['aa_start']:4d}-{r['aa_end']:<4d} "
              f"o/e {r['oe_mis']:.3f} [{r['oe_mis_lower']:.3f}-{r['oe_mis_upper']:.3f}] "
              f"(gnomAD olceginde {r['oe_mis_gnomad_scaled']:.3f})  n={r['obs_mis']:4d}  "
              f"kuyruk o/e {r['tail_oe_mis']:.3f}  motor/kuyruk {r['ratio_motor_over_tail']:.3f} "
              f"p={r['p_motor_vs_tail']}  -> {r['grade']}")
    core.save_tsv(outdir, f"motor_domain_panel_{stamp}.tsv", motor_rows, base=core.INTERIM,
                  source="derived: regional missense constraint, kinesin motor domains",
                  urls=[GNOMAD_API, ENSEMBL, RATES_URL], tool="scripts/regional_constraint.py",
                  params={"panel": list(genes), "max_af": args.max_af})

    # per-domain rows for every panel gene (motor / after motor), plus KIF13A's fine domains
    domain_rows = list(kif13a_domains)
    for gene, (lo, hi, src) in motor_bounds.items():
        if gene == "KIF13A":
            continue
        domain_rows.append({**region_stats(genes[gene]["rows"], lo, hi, "motor", gene),
                            "boundary_source": src})
        domain_rows.append({**region_stats(genes[gene]["rows"], hi + 1,
                                           genes[gene]["protein_length"], "motor sonrasi", gene),
                            "boundary_source": src})
    anchor(domain_rows)
    core.save_tsv(outdir, f"per_domain_{stamp}.tsv", domain_rows, base=core.INTERIM,
                  source="derived: regional missense constraint per domain",
                  urls=[GNOMAD_API, ENSEMBL, RATES_URL], tool="scripts/regional_constraint.py",
                  params={"kif13a_domains": KIF13A_DOMAINS, "max_af": args.max_af})

    # ---------------- sliding windows
    profiles, fine = {}, {}
    for gene, g in genes.items():
        profiles[gene] = anchor(sliding(g["rows"], g["protein_length"], gene, WINDOW, STEP))
        fine[gene] = anchor(sliding(g["rows"], g["protein_length"], gene, FINE_WINDOW, FINE_STEP))
    core.save_tsv(outdir, f"sliding_window_{WINDOW}aa_{stamp}.tsv",
                  [r for gene in profiles for r in profiles[gene]], base=core.INTERIM,
                  source="derived: sliding-window regional missense constraint",
                  urls=[GNOMAD_API], tool="scripts/regional_constraint.py",
                  params={"window": WINDOW, "step": STEP})
    core.save_tsv(outdir, f"sliding_window_{FINE_WINDOW}aa_{stamp}.tsv",
                  [r for gene in fine for r in fine[gene]], base=core.INTERIM,
                  source="derived: sliding-window regional missense constraint, fine",
                  urls=[GNOMAD_API], tool="scripts/regional_constraint.py",
                  params={"window": FINE_WINDOW, "step": FINE_STEP})

    # ---------------- sanity checks
    k13, k1a = gene_level["KIF13A"], gene_level["KIF1A"]
    k13_motor = next(r for r in motor_rows if r["gene"] == "KIF13A")
    k1a_motor = next(r for r in motor_rows if r["gene"] == "KIF1A")
    syn_windows = [float(p["oe_syn"]) for p in profiles["KIF13A"] if p["oe_syn"] != ""]
    syn_windows_1a = [float(p["oe_syn"]) for p in profiles["KIF1A"] if p["oe_syn"] != ""]

    def syn_dispersion(gene: str) -> dict:
        """Is the spread of synonymous o/e across NON-overlapping windows just Poisson noise?

        A 100-codon window expects only about 45 synonymous variants, so a raw band like
        0.80-1.25 is the wrong test: Poisson noise alone puts roughly a third of windows outside
        it. The right question is whether the dispersion exceeds Poisson, which is what a
        chi-square goodness of fit on disjoint windows answers.
        """
        g = genes[gene]
        chi, df, worst = 0.0, 0, (None, 1.0)
        for start in range(1, g["protein_length"] - WINDOW + 2, WINDOW):
            r = region_stats(g["rows"], start, start + WINDOW - 1, "x", gene)
            if r["exp_syn"] <= 5:
                continue
            chi += (r["obs_syn"] - r["exp_syn"]) ** 2 / r["exp_syn"]
            df += 1
            if abs(float(r["oe_syn"]) - 1) > abs(worst[1] - 1):
                worst = (f"{start}-{start + WINDOW - 1}", float(r["oe_syn"]))
        return {"chi2": round(chi, 1), "df": df, "chi2_over_df": round(chi / df, 3) if df else None,
                "worst_window": worst[0], "worst_oe_syn": round(worst[1], 3)}

    disp13, disp1a = syn_dispersion("KIF13A"), syn_dispersion("KIF1A")
    checks = [
        {"check": "KIF1A motor domaini pozitif kontrol olarak kisitli cikiyor mu",
         "result": f"o/e {k1a_motor['oe_mis']:.3f} "
                   f"[{k1a_motor['oe_mis_lower']:.3f}-{k1a_motor['oe_mis_upper']:.3f}], "
                   f"{k1a_motor['obs_mis']} gozlem -> {k1a_motor['grade']}",
         "pass": bool(k1a_motor["oe_mis_upper"] < 0.60)},
        {"check": "KIF13A gen duzeyi gnomAD v4.1 yayimlanan o/e missense degerini yeniden uretiyor mu",
         "result": f"bizim {k13['oe_mis']:.3f} vs gnomAD {k13['oe_mis_published']:.3f} "
                   f"({k13['oe_mis_delta_pct']:+.1f} %); gozlenen {k13['obs_mis']} vs "
                   f"{k13['obs_mis_published']} ({k13['obs_mis_delta_pct']:+.1f} %); "
                   f"beklenen {k13['exp_mis']:.1f} vs {k13['exp_mis_published']:.1f} "
                   f"({k13['exp_mis_delta_pct']:+.1f} %)",
         "pass": bool(abs(float(k13["oe_mis_delta_pct"])) < 10)},
        {"check": "sinonim o/e her pencerede 1,0 civarinda mi (KIF13A)",
         "result": f"medyan {np.median(syn_windows):.3f}, aralik "
                   f"{min(syn_windows):.3f}-{max(syn_windows):.3f}; ortusmeyen pencerelerde "
                   f"ki-kare/sd = {disp13['chi2_over_df']} ({disp13['df']} pencere), en sapan "
                   f"pencere {disp13['worst_window']} o/e {disp13['worst_oe_syn']}",
         "pass": bool(0.9 <= np.median(syn_windows) <= 1.1 and (disp13["chi2_over_df"] or 9) <= 1.6)},
        {"check": "sinonim o/e her pencerede 1,0 civarinda mi (KIF1A)",
         "result": f"medyan {np.median(syn_windows_1a):.3f}, aralik "
                   f"{min(syn_windows_1a):.3f}-{max(syn_windows_1a):.3f}; ki-kare/sd = "
                   f"{disp1a['chi2_over_df']} ({disp1a['df']} pencere), en sapan pencere "
                   f"{disp1a['worst_window']} o/e {disp1a['worst_oe_syn']}",
         "pass": bool(0.9 <= np.median(syn_windows_1a) <= 1.2
                      and (disp1a["chi2_over_df"] or 9) <= 1.6)},
        {"check": "KIF13A motor domaininde kac gozlenen missense var (aralik genisligi)",
         "result": f"{k13_motor['obs_mis']} gozlenen, {k13_motor['exp_mis']:.1f} beklenen; "
                   f"%90 aralik genisligi {k13_motor['oe_mis_upper'] - k13_motor['oe_mis_lower']:.3f}",
         "pass": bool(k13_motor["obs_mis"] >= 100)},
        {"check": "KIF13A motor domaininde sinonim o/e (yerel kalibrasyon kontrolu)",
         "result": f"{k13_motor['oe_syn']:.3f} "
                   f"[{k13_motor['oe_syn_lower']:.3f}-{k13_motor['oe_syn_upper']:.3f}], "
                   f"{k13_motor['obs_syn']} gozlenen / {k13_motor['exp_syn']:.1f} beklenen",
         "pass": bool(0.75 <= float(k13_motor["oe_syn"]) <= 1.3)},
        {"check": "KIF13A motor domaini kendi kuyruguna gore anlamli olcude tukenmis mi "
                  "(kalibrasyondan bagimsiz ic karsilastirma)",
         "result": f"motor/kuyruk oran {k13_motor['ratio_motor_over_tail']:.3f}, "
                   f"p={k13_motor['p_motor_vs_tail']}; ayni karsilastirma KIF1A'da "
                   f"{k1a_motor['ratio_motor_over_tail']:.3f}, p={k1a_motor['p_motor_vs_tail']}",
         "pass": bool(float(k13_motor["p_motor_vs_tail"]) < 0.05)},
        {"check": "kisitlilik olcumu klinik gerceklikle uyusuyor mu (motor domaindeki ClinVar P/LP "
                  "missense sayisi, kisitlilik sirasina karsi)",
         "result": "; ".join(f"{r['gene']} o/e {r['oe_mis']:.2f} -> "
                             f"{r['clinvar_plp_missense_in_motor']} P/LP"
                             for r in sorted(motor_rows, key=lambda x: x["oe_mis"])),
         "pass": True},
        {"check": "nonsense o/e motor domainde LoF eksenini tekrar veriyor mu (ic tutarlilik)",
         "result": f"KIF13A motor nonsense o/e {k13_motor['oe_non']:.3f} "
                   f"({k13_motor['obs_non']}/{k13_motor['exp_non']:.1f}); "
                   f"KIF1A motor {k1a_motor['oe_non']:.3f} "
                   f"({k1a_motor['obs_non']}/{k1a_motor['exp_non']:.1f})",
         "pass": True},
    ]
    print("\n== sanity checks")
    for c in checks:
        print(f"  [{'GECTI' if c['pass'] else 'KALDI'}] {c['check']}\n        {c['result']}")
    core.save_tsv(outdir, f"sanity_checks_{stamp}.tsv", checks, base=core.INTERIM,
                  source="derived", urls=[], tool="scripts/regional_constraint.py", params={})

    # ---------------- candidate-list impact
    summary = {
        "measured_utc_date": stamp,
        "method": {
            "mutation_rates": "Samocha 2014 trinucleotide table (PMID 25086666)",
            "rates_url": RATES_URL,
            "observed": f"gnomAD v4.1 exome SNVs, PASS, AC>=1, AF<={args.max_af}",
            "calibration": f"P(observed|mu) from {cal['n_syn_possible']} possible synonymous SNVs "
                           f"pooled over {len(genes)} genes, {len(cal['bins'])} rate bins",
            "interval": "90 % Poisson likelihood interval on o/e (the gnomAD LOEUF construction)",
        },
        "kif13a": {
            "gene_level": gene_level["KIF13A"],
            "per_domain": kif13a_domains,
            "motor": k13_motor,
            "verdict": grade(k13_motor["oe_mis_upper"]),
        },
        "kif1a_positive_control": {"motor": k1a_motor, "verdict": grade(k1a_motor["oe_mis_upper"]),
                                   "passed": bool(k1a_motor["oe_mis_upper"] < 0.60)},
        "kinesin_frame": sorted(
            [{"gene": r["gene"], "motor_aa": f"{r['aa_start']}-{r['aa_end']}",
              "motor_oe": r["oe_mis"], "motor_lo": r["oe_mis_lower"], "motor_hi": r["oe_mis_upper"],
              "motor_oe_gnomad_scaled": r["oe_mis_gnomad_scaled"], "obs": r["obs_mis"],
              "tail_oe": r["tail_oe_mis"], "diff": r["motor_minus_tail"],
              "ratio_motor_over_tail": r["ratio_motor_over_tail"], "p": r["p_motor_vs_tail"],
              "clinvar_plp_missense_in_motor": r["clinvar_plp_missense_in_motor"],
              "clinvar_plp_missense_total": r["clinvar_plp_missense_total"],
              "grade": r["grade"], "bounds": r["boundary_source"]}
             for r in motor_rows], key=lambda x: x["motor_oe"]),
        "sanity_checks": checks,
        "motor_subregions": sub_rows,
        "structural_classes": struct_rows,
        "tier_a": {"n": len(tierA), "position_range": [tierA_pos[0], tierA_pos[-1]] if tierA_pos else [],
                   "in_motor_5_352": sum(1 for p in tierA_pos if 5 <= p <= 352),
                   "in_270_360": sum(1 for p in tierA_pos if 270 <= p <= 360),
                   "in_motor_n_half_5_180": sum(1 for p in tierA_pos if 5 <= p <= 180),
                   "structural_classes": dict(Counter(tierA_struct).most_common())},
        "calibration_bins": cal["bins"],
    }
    rec = core.save_json(outdir, f"regional_constraint_summary_{stamp}.json", summary,
                         base=core.INTERIM, source="derived", urls=[GNOMAD_API, ENSEMBL, RATES_URL],
                         tool="scripts/regional_constraint.py", params={"panel": panel})
    print(f"\n  wrote {rec['path']}")

    # the full per-SNV table for KIF13A, so the windows can be recomputed without the network
    core.save_tsv(outdir, f"KIF13A_possible_snv_constraint_{stamp}.tsv", genes["KIF13A"]["rows"],
                  base=core.INTERIM,
                  source="derived: every possible KIF13A CDS SNV with context, rate, calibrated "
                         "expectation and whether gnomAD v4.1 has seen it",
                  urls=[GNOMAD_API, ENSEMBL, RATES_URL], tool="scripts/regional_constraint.py",
                  params={"transcript": tx["KIF13A"], "max_af": args.max_af})

    if not args.no_figure:
        print("\n== figure")
        build_figure(profiles, struct_rows, motor_rows, gene_level,
                     ROOT / "figures" / "KIF13A_bolgesel_kisitlilik.svg")

    print(f"\nVERDICT  KIF13A motor 5-352 missense o/e = {k13_motor['oe_mis']:.3f} "
          f"[{k13_motor['oe_mis_lower']:.3f}-{k13_motor['oe_mis_upper']:.3f}] "
          f"-> {grade(k13_motor['oe_mis_upper'])}   "
          f"(kuyruk 520-1805 o/e = {kif13a_domains[-1]['oe_mis']:.3f}, "
          f"gen butunu {gene_level['KIF13A']['oe_mis']:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
