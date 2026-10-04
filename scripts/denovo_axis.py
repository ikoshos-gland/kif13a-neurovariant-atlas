"""The de novo axis: is KIF13A hit by de novo variants in neurodevelopmental trios?

KIF13A has no established disease. gnomAD already says it tolerates loss of function better than
any ClinGen haploinsufficiency-3 neurodevelopmental gene (observed/expected LoF 0.372 versus
0.037-0.287), so severe fully penetrant haploinsufficiency is unlikely. What survives is mild or
incompletely penetrant dose sensitivity, and the test that settles that is de novo variant burden in
large trio cohorts: if losing one copy of KIF13A causes a developmental disorder, KIF13A must be hit
by de novo loss-of-function variants more often than chance.

This script is written so that it can return a clean zero, and a clean zero is the result.

  1. Kaplanis et al. 2020 Nature (PMID 33057194), 31,058 developmental disorder trios
     (DDD + GeneDx + Radboud). Supplementary Table 1 is the full list of 45,221 de novo variants;
     Supplementary Table 2 is the per-gene count plus the Samocha-framework denovoWEST statistics.
     Its coordinates are GRCh37; the script proves that rather than assuming it, lifts them to
     GRCh38 and re-annotates them on the MANE transcript so every protein position is MANE 1-1805.
  2. Gene4Denovo2 (PMID 41024706), the aggregated de novo database, queried live through the same
     JSON API its own web interface uses, variant level and gene level, for all five genes.
  3. denovo-db v1.6.1 (PMID 27907889), non-SSC portion. The SSC portion is deliberately NOT
     downloaded: it carries a stated use restriction.
  4. gnomAD v2.1.1 per-gene trinucleotide mutation rates (mu_lof, mu_mis, mu_syn), which are the
     Samocha framework rates, so that expected = 2 * n_trios * mu can be computed per gene and the
     observed count can be turned into a Poisson p value instead of being reported naked.
  5. SFARI Gene, live, as an independent "is this gene on anybody's list" check.

Controls run through the identical pipeline: ANKRD11 and CHD8 (established haploinsufficiency,
must be strongly enriched), KIF1A (disease gene, missense mechanism, must show de novo missense but
not de novo LoF), KIF13B (paralog, expected null).

    python scripts\\denovo_axis.py
    python scripts\\denovo_axis.py --refresh

Nothing is re-downloaded if it is already in data/raw (see data/raw/PROVENANCE.jsonl).
Outputs land in data/interim/denovo/.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import re
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

import httpx  # noqa: E402

from atlas_mcp import core  # noqa: E402

# --------------------------------------------------------------------------- constants

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
ENSEMBL = "https://rest.ensembl.org"
G4D2_API = "https://java.genemed.tech/admin-api/backend/data/hg38/{table}/geneSymbol/page"
G4D2_UI = "https://genemed.tech/gene4denovo2/#/genemed/search?reference=hg38&queryby=geneSymbol&terms={gene}"
G4D2_VARIANT_TABLE = "gene4denovo2_all"
G4D2_GENE_TABLE = "gene4denovo2_Candidate_Gene2"
DENOVODB_TSV = "https://denovo-db.gs.washington.edu/denovo-db.non-ssc-samples.variants.tsv.gz"
DENOVODB_SSC = "https://denovo-db.gs.washington.edu/denovo-db.ssc-samples.variants.tsv.gz"
GNOMAD_RATES = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/"
    "release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
)
SFARI_GENE = "https://gene.sfari.org/database/human-gene/{gene}"
KAPLANIS_PMCID = "PMC7116826"
KAPLANIS_SUPPL = f"{EPMC}/{KAPLANIS_PMCID}/supplementaryFiles"
KAPLANIS_PMID = "33057194"
KAPLANIS_TRIOS = 31058  # verified live from the abstract, not assumed

BROWSER_UA = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124 Safari/537.36"
)

GENES = ["KIF13A", "KIF1A", "KIF13B", "ANKRD11", "CHD8"]

# KIF13A domains, MANE numbering (identical table to scripts/lof_axis.py)
DOMAINS = [(5, 352, "motor"), (353, 380, "boyun"), (381, 468, "CC1"), (469, 519, "FHA"), (520, 1805, "kuyruk")]
KIF13A_LEN = 1805
MOTOR_START, MOTOR_END = 5, 352

# Kaplanis / VEP consequence terms that count as protein-truncating
LOF_TERMS = {"frameshift_variant", "stop_gained", "splice_acceptor_variant", "splice_donor_variant"}
MIS_TERMS = {"missense_variant"}
SYN_TERMS = {"synonymous_variant"}

# Gene4Denovo2 ExonicFunc / Func labels
G4D2_LOF = {"frameshift deletion", "frameshift insertion", "stopgain", "stoploss",
            "frameshift block substitution"}
G4D2_MIS = {"nonsynonymous SNV"}
G4D2_SYN = {"synonymous SNV"}

# denovo-db FunctionClass labels
DDB_LOF = {"frameshift", "stop-gained", "splice-donor", "splice-acceptor", "start-lost",
           "frameshift-near-splice", "stop-gained-near-splice"}
DDB_MIS = {"missense", "missense-near-splice"}
DDB_SYN = {"synonymous", "synonymous-near-splice"}

# The three known human KIF13A observations, to be looked for in every cohort
KNOWN = {
    "Lam2021_mosaic_frameshift": {"hgvsc": "c.2357dupA", "protein_pos": 786, "pmid": "33526817"},
    "K1122E_tail_VUS": {"hgvsc": "c.3364A>G", "protein_pos": 1122, "vcv": "VCV002570672"},
    "del_253kb": {"vcv": "VCV000442594", "note": "253 kb 6p22.3 deletion, KIF13A only"},
}

OUT = core.INTERIM / "denovo"


# --------------------------------------------------------------------------- small helpers


def http_get(url: str, *, binary: bool = False, headers: dict | None = None, timeout: float = 300.0):
    hdrs = {"User-Agent": BROWSER_UA, "Accept": "*/*"}
    hdrs.update(headers or {})
    with httpx.Client(timeout=timeout, follow_redirects=True) as c:
        r = c.get(url, headers=hdrs)
        r.raise_for_status()
        return r.content if binary else r.text


def http_post_json(url: str, body: dict, *, timeout: float = 300.0) -> dict:
    hdrs = {
        "User-Agent": BROWSER_UA,
        "Content-Type": "application/json",
        "Accept": "application/json",
        "Origin": "https://genemed.tech",
        "Referer": "https://genemed.tech/",
    }
    with httpx.Client(timeout=timeout, follow_redirects=True) as c:
        r = c.post(url, json=body, headers=hdrs)
        r.raise_for_status()
        return r.json()


def cached(subdir: str, pattern: str, base: Path = core.RAW) -> Path | None:
    folder = base / subdir
    if not folder.exists():
        return None
    hits = sorted(p for p in folder.glob(pattern) if not p.name.endswith(".prov.json"))
    return hits[-1] if hits else None


def domain_of(pos: int | None) -> str:
    if not pos:
        return ""
    for lo, hi, name in DOMAINS:
        if lo <= pos <= hi:
            return name
    return "domain disi"


def poisson_p_ge(k: int, lam: float) -> float:
    """P(X >= k) for X ~ Poisson(lam), summed over the UPPER tail in log space.

    The naive 1 - CDF form underflows to exactly 0 for the positive controls (ANKRD11 has 111
    observed against 0.31 expected), which would print a fake p = 0. Summing the upper tail in
    log space keeps those p values as real numbers down to the double-precision floor.
    """
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    logs: list[float] = []
    i = k
    logterm = -lam + i * math.log(lam) - math.lgamma(i + 1)
    logs.append(logterm)
    peak = logterm
    while True:
        i += 1
        logterm += math.log(lam) - math.log(i)
        peak = max(peak, logterm)
        logs.append(logterm)
        if i > lam and logterm < peak - 80:
            break
        if i > k + 100000:
            break
    mx = max(logs)
    return math.exp(mx) * sum(math.exp(l - mx) for l in logs)


def poisson_upper_95(k: int) -> float:
    """Upper limit of the exact 95% Poisson confidence interval for k observed events.

    Uses the chi-square relation lambda_hi = 0.5 * chi2_{0.975, 2(k+1)}. For k = 0 this is
    the familiar 3.689 for a two-sided 95% interval; the commonly quoted 'rule of three'
    (3.0) is the one-sided 95% bound. We report the one-sided bound, which is the relevant
    one for 'how much enrichment can we still not rule out'.
    """
    # one-sided upper 95%: solve P(X <= k | lambda) = 0.05
    lo, hi = 0.0, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        # P(X <= k | mid)
        term = math.exp(-mid)
        cum = term
        for i in range(1, k + 1):
            term *= mid / i
            cum += term
        if cum > 0.05:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def poisson_lower_95(k: int) -> float:
    """Lower limit of the one-sided 95% Poisson interval for k observed events (0 for k = 0)."""
    if k <= 0:
        return 0.0
    lo, hi = 0.0, float(k) * 10 + 20
    for _ in range(200):
        mid = (lo + hi) / 2
        # P(X >= k | mid) = 0.05  ->  solve
        if poisson_p_ge(k, mid) < 0.05:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def tr(x: float, digits: int = 3) -> str:
    """Turkish decimal comma."""
    if x is None:
        return ""
    if isinstance(x, str):
        return x
    s = f"{x:.{digits}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{digits}f}"
    return s.replace(".", ",")


# --------------------------------------------------------------------------- 1. mutation rates


def mutation_rates(refresh: bool) -> dict[str, dict]:
    """gnomAD v2.1.1 per-gene trinucleotide mutation rates (the Samocha framework rates)."""
    print("\n[1] gnomAD v2.1.1 per-gene mutation rates (mu_lof / mu_mis / mu_syn)")
    hit = None if refresh else cached("denovo", "gnomad_v211_lof_metrics_by_gene_*.txt.bgz")
    if hit:
        raw = hit.read_bytes()
        print(f"  using cached {core.rel(hit)}")
    else:
        raw = http_get(GNOMAD_RATES, binary=True)
        rec = core.save_bytes(
            "denovo", f"gnomad_v211_lof_metrics_by_gene_{core.today()}.txt.bgz", raw,
            source="gnomAD v2.1.1 constraint, per-gene LoF metrics (carries mu_syn/mu_mis/mu_lof)",
            urls=[GNOMAD_RATES], tool="scripts/denovo_axis.py",
        )
        print(f"  downloaded {rec['path']} ({rec['bytes']} bytes)")

    fh = gzip.open(io.BytesIO(raw), "rt", encoding="utf-8", errors="replace")
    hdr = fh.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    out: dict[str, dict] = {}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        g = p[idx["gene"]]
        if g not in GENES:
            continue
        rec = {
            "gene": g,
            "transcript": p[idx["transcript"]],
            "mu_lof": float(p[idx["mu_lof"]]),
            "mu_mis": float(p[idx["mu_mis"]]),
            "mu_syn": float(p[idx["mu_syn"]]),
            "cds_length": int(p[idx["cds_length"]]),
        }
        # keep the row whose transcript matches our MANE transcript when we have one
        prev = out.get(g)
        if prev is None:
            out[g] = rec
        else:
            try:
                mane = core.unversioned(core.target(g)["mane_transcript"])
            except Exception:
                mane = None
            if mane and rec["transcript"] == mane:
                out[g] = rec
    for g in GENES:
        r = out.get(g)
        if not r:
            print(f"  {g}: ABSENT from the constraint file")
            continue
        aa = r["cds_length"] / 3
        print(f"  {g:8s} {r['transcript']}  CDS {r['cds_length']} nt = {aa:.1f} aa  "
              f"mu_lof {r['mu_lof']:.4g}  mu_mis {r['mu_mis']:.4g}")
    return out


def sanity_check_lengths(rates: dict[str, dict]) -> list[str]:
    """CDS length is an independent check on the protein length trap (1805 vs 1770)."""
    notes = []
    for g in ("KIF13A", "KIF1A", "KIF13B"):
        r = rates.get(g)
        if not r:
            continue
        t = core.target(g)
        expected_cds = (t["protein_length"] + 1) * 3  # +1 for the stop codon
        ok = r["cds_length"] == expected_cds or r["cds_length"] == t["protein_length"] * 3
        aa = r["cds_length"] / 3
        notes.append(f"{g}: CDS {r['cds_length']} nt -> {aa:.0f} aa, config says {t['protein_length']} aa "
                     f"[{'TUTUYOR' if ok else 'TUTMUYOR'}]")
        mane = core.unversioned(t["mane_transcript"])
        notes.append(f"{g}: rate table transcript {r['transcript']} vs MANE {mane} "
                     f"[{'TUTUYOR' if r['transcript'] == mane else 'TUTMUYOR'}]")
    return notes


# --------------------------------------------------------------------------- 2. Kaplanis 2020


def kaplanis(refresh: bool) -> tuple[list[dict], dict[str, dict], dict]:
    """Supplementary Table 1 (all DNMs) and Table 2 (per-gene counts + denovoWEST)."""
    print(f"\n[2] Kaplanis et al. 2020 Nature (PMID {KAPLANIS_PMID}), {KAPLANIS_TRIOS:,} DD trios")
    hit = None if refresh else cached("denovo", "kaplanis2020_supplementary_*.zip")
    if hit:
        blob = hit.read_bytes()
        print(f"  using cached {core.rel(hit)}")
    else:
        blob = http_get(KAPLANIS_SUPPL, binary=True)
        rec = core.save_bytes(
            "denovo", f"kaplanis2020_supplementary_{core.today()}.zip", blob,
            source=f"Europe PMC supplementary files for {KAPLANIS_PMCID} (Kaplanis 2020 Nature)",
            urls=[KAPLANIS_SUPPL], tool="scripts/denovo_axis.py",
            params={"pmid": KAPLANIS_PMID, "pmcid": KAPLANIS_PMCID, "n_trios": KAPLANIS_TRIOS},
        )
        print(f"  downloaded {rec['path']} ({rec['bytes']} bytes)")

    z = zipfile.ZipFile(io.BytesIO(blob))
    t1 = next(n for n in z.namelist() if "SuppTable1" in n)
    t2 = next(n for n in z.namelist() if "SuppTable2" in n)

    # -- Table 1: every de novo variant
    lines = z.read(t1).decode("utf-8", errors="replace").splitlines()
    hdr = lines[0].split("\t")
    dnms: list[dict] = []
    studies = Counter()
    total = 0
    for line in lines[1:]:
        p = line.split("\t")
        if len(p) < len(hdr):
            continue
        d = dict(zip(hdr, p))
        total += 1
        studies[d["study"]] += 1
        if d["symbol"] in GENES:
            dnms.append(d)
    print(f"  Supplementary Table 1: {total:,} de novo variants, studies {dict(studies)}")

    # -- Table 2: per-gene counts and denovoWEST
    import openpyxl
    wb = openpyxl.load_workbook(io.BytesIO(z.read(t2)), read_only=True)
    ws = wb[wb.sheetnames[0]]
    it = ws.iter_rows(values_only=True)
    cols = list(next(it))
    per_gene: dict[str, dict] = {}
    n_genes = 0
    for row in it:
        n_genes += 1
        if row[0] in GENES:
            per_gene[row[0]] = dict(zip(cols, row))
    wb.close()
    print(f"  Supplementary Table 2: {n_genes:,} genes, sheet '{ws.title}'")

    meta = {"n_trios": KAPLANIS_TRIOS, "n_dnm_total": total, "studies": dict(studies),
            "n_genes_table2": n_genes, "table2_sheet": ws.title,
            "table1_file": t1, "table2_file": t2}
    return dnms, per_gene, meta


def kaplanis_build(dnms: list[dict], refresh: bool) -> list[dict]:
    """Prove the build of the Kaplanis coordinates, then annotate KIF13A rows on MANE.

    Kaplanis reports no assembly in the table header. We test it: VEP each KIF13A row on both
    builds and keep the one whose consequence matches the consequence the table states.
    """
    print("\n[2b] Kaplanis coordinate build test + MANE annotation")
    kif13a = [d for d in dnms if d["symbol"] == "KIF13A"]
    if not kif13a:
        print("  no KIF13A rows to test")
        return []

    hit = None if refresh else cached("denovo", "kaplanis_KIF13A_mane_annotation_*.json")
    if hit:
        payload = json.loads(hit.read_text(encoding="utf-8"))
        print(f"  using cached {core.rel(hit)}")
        return payload["rows"]

    def vep(host: str, chrom: str, pos: int, alt: str) -> list[dict]:
        url = f"{host}/vep/human/region/{chrom}:{pos}-{pos}/{alt}?content-type=application/json&mane=1&hgvs=1"
        try:
            return json.loads(http_get(url, headers={"Content-Type": "application/json"}, timeout=120))
        except Exception:
            return []

    def kif13a_tc(payload: list[dict]) -> dict | None:
        if not payload:
            return None
        for t in payload[0].get("transcript_consequences", []):
            if t.get("transcript_id") == "ENST00000259711":
                return t
        return None

    # -- build test
    votes = Counter()
    for d in kif13a:
        for host, label in ((ENSEMBL, "GRCh38"), ("https://grch37.rest.ensembl.org", "GRCh37")):
            t = kif13a_tc(vep(host, d["chrom"], int(d["pos"]), d["alt"]))
            if t and d["consequence"] in t.get("consequence_terms", []):
                votes[label] += 1
            time.sleep(0.15)
    build = votes.most_common(1)[0][0] if votes else "bilinmiyor"
    print(f"  build test: {dict(votes)} of {len(kif13a)} rows match the stated consequence -> {build}")

    # -- annotate on MANE (lift over first when the source build is GRCh37)
    rows = []
    for d in kif13a:
        chrom, pos, alt = d["chrom"], int(d["pos"]), d["alt"]
        pos38 = pos
        if build == "GRCh37":
            url = f"{ENSEMBL}/map/human/GRCh37/{chrom}:{pos}..{pos}/GRCh38?content-type=application/json"
            m = json.loads(http_get(url, headers={"Content-Type": "application/json"}, timeout=120))
            maps = m.get("mappings") or []
            if not maps:
                print(f"  {chrom}:{pos} could not be lifted over, skipping")
                continue
            pos38 = maps[0]["mapped"]["start"]
            time.sleep(0.15)
        t = kif13a_tc(vep(ENSEMBL, chrom, pos38, alt))
        time.sleep(0.15)
        prot = t.get("protein_start") if t else None
        rows.append({
            "kaynak": "Kaplanis 2020",
            "id": d["id"],
            "study": d["study"],
            "fenotip": "gelisimsel bozukluk (DD)",
            "build_kaynak": build,
            "grch37": f"{chrom}:{pos}{d['ref']}>{alt}" if build == "GRCh37" else "",
            "grch38": f"chr{chrom}:{pos38}{d['ref']}>{alt}",
            "tablo_consequence": d["consequence"],
            "vep_consequence": ",".join(t.get("consequence_terms", [])) if t else "",
            "hgvsc_mane": t.get("hgvsc") if t else "",
            "hgvsp_mane": t.get("hgvsp") if t else "",
            "transkript": "ENST00000259711.11 (MANE, NM_022113.6)",
            "protein_pos_MANE": prot,
            "domain": domain_of(prot),
            "sinif": ("LoF" if d["consequence"] in LOF_TERMS else
                      "missense" if d["consequence"] in MIS_TERMS else
                      "sinonim" if d["consequence"] in SYN_TERMS else d["consequence"]),
        })
        print(f"    {d['id']:16s} {d['consequence']:20s} -> {rows[-1]['hgvsp_mane']} "
              f"pos {prot} ({rows[-1]['domain']})")

    core.save_json(
        "denovo", f"kaplanis_KIF13A_mane_annotation_{core.today()}.json",
        {"build": build, "votes": dict(votes), "rows": rows}, base=core.RAW,
        source="Ensembl REST assembly map + VEP, KIF13A de novo rows of Kaplanis 2020 on MANE",
        urls=[f"{ENSEMBL}/vep/human/region", f"{ENSEMBL}/map/human/GRCh37"],
        tool="scripts/denovo_axis.py",
    )
    return rows


# --------------------------------------------------------------------------- 3. Gene4Denovo2


def gene4denovo2(refresh: bool) -> tuple[dict[str, list[dict]], dict[str, list[dict]]]:
    """Live query of the Gene4Denovo2 JSON API, variant level and gene level, for all five genes."""
    print("\n[3] Gene4Denovo2 (PMID 41024706), live API")
    variants: dict[str, list[dict]] = {}
    genelevel: dict[str, list[dict]] = {}

    hit = None if refresh else cached("denovo", "gene4denovo2_pull_*.json")
    if hit:
        payload = json.loads(hit.read_text(encoding="utf-8"))
        print(f"  using cached {core.rel(hit)}")
        return payload["variants"], payload["genelevel"]

    for gene in GENES:
        for table, mode, sink in ((G4D2_VARIANT_TABLE, "variant", variants),
                                  (G4D2_GENE_TABLE, "gene", genelevel)):
            body = {"reference": "hg38", "queryby": "geneSymbol", "terms": gene,
                    "orderBy": "", "order": "", "page": 1, "pageNo": 1, "currentPage": 1,
                    "pageSize": 2000, "filterParam": "[]", "mode": mode}
            data = http_post_json(G4D2_API.format(table=table), body)["data"]
            recs = data.get("records") or []
            sink[gene] = recs
            print(f"  {gene:8s} {mode:8s} total {data.get('total'):>5}  fetched {len(recs):>5}")
            time.sleep(0.3)

    core.save_json(
        "denovo", f"gene4denovo2_pull_{core.today()}.json",
        {"variants": variants, "genelevel": genelevel}, base=core.RAW,
        source="Gene4Denovo2 backend API (the endpoint its own web UI calls), hg38, geneSymbol query",
        urls=[G4D2_API.format(table=G4D2_VARIANT_TABLE), G4D2_API.format(table=G4D2_GENE_TABLE)],
        tool="scripts/denovo_axis.py",
        params={"genes": GENES, "reference": "hg38", "tables": [G4D2_VARIANT_TABLE, G4D2_GENE_TABLE]},
    )
    return variants, genelevel


G4D2_AA_RE = re.compile(r"p\.([A-Z])(\d+)([A-Za-z*]+)")
# denovo-db writes protein changes as p.(Q237H) / p.(N1092=), with the parentheses
PROT_RE = re.compile(r"p\.\(?([A-Z])(\d+)")


def g4d2_protein_pos(rec: dict) -> tuple[int | None, str, str]:
    """Protein position on the MANE transcript.

    Gene4Denovo2 ships two annotation columns. AAChange_refGene uses whichever RefSeq transcripts
    the RefGene build carries (for KIF13A that is NM_001105567 / NM_001105568, a DIFFERENT isoform
    whose numbering is 13 residues lower). AAChange_mrefGene uses the MANE RefSeq NM_022113.
    We read mrefGene and record which transcript the number came from.
    """
    for field, label in (("AAChange_mrefGene", "mref"), ("AAChange_refGene", "ref")):
        val = rec.get(field) or ""
        if not val or val == ".":
            continue
        for chunk in val.split(","):
            parts = chunk.split(":")
            if len(parts) < 4:
                continue
            tx = parts[1]
            m = G4D2_AA_RE.search(chunk)
            if m:
                return int(m.group(2)), tx, chunk
    return None, "", ""


def g4d2_class(rec: dict) -> str:
    func = (rec.get("Func_refGene") or "").strip()
    ex = (rec.get("ExonicFunc_mrefGene") or rec.get("ExonicFunc_refGene") or "").strip()
    if func in ("splicing", "exonic;splicing"):
        return "LoF"
    if ex in G4D2_LOF:
        return "LoF"
    if ex in G4D2_MIS:
        return "missense"
    if ex in G4D2_SYN:
        return "sinonim"
    if func in ("exonic",):
        return ex or "exonic-diger"
    return f"kodlamayan ({func})" if func else "bilinmiyor"


# --------------------------------------------------------------------------- 4. denovo-db


def denovodb(refresh: bool) -> tuple[list[dict], dict]:
    print("\n[4] denovo-db v1.6.1 (PMID 27907889), non-SSC portion")
    hit = None if refresh else cached("denovo", "denovodb_nonssc_v161_*.tsv.gz")
    if hit:
        blob = hit.read_bytes()
        print(f"  using cached {core.rel(hit)}")
    else:
        blob = http_get(DENOVODB_TSV, binary=True)
        rec = core.save_bytes(
            "denovo", f"denovodb_nonssc_v161_{core.today()}.tsv.gz", blob,
            source="denovo-db v1.6.1, non-SSC samples, variant TSV",
            urls=[DENOVODB_TSV], tool="scripts/denovo_axis.py",
            params={"note": "SSC portion deliberately not downloaded: stated use restriction"},
        )
        print(f"  downloaded {rec['path']} ({rec['bytes']} bytes)")

    fh = gzip.open(io.BytesIO(blob), "rt", encoding="utf-8", errors="replace")
    version = fh.readline().strip().lstrip("#")
    hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
    rows: list[dict] = []
    total = 0
    phenos = Counter()
    studies = Counter()
    samples: dict[str, set] = defaultdict(set)
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr):
            continue
        total += 1
        d = dict(zip(hdr, p))
        studies[d["StudyName"]] += 1
        phenos[d["PrimaryPhenotype"]] += 1
        samples[d["StudyName"]].add(d["SampleID"])
        if d["Gene"] in GENES:
            rows.append(d)
    meta = {"version": version, "n_variants_total": total,
            "n_studies": len(studies), "n_samples_seen": sum(len(v) for v in samples.values()),
            "top_phenotypes": dict(phenos.most_common(12))}
    print(f"  {version}: {total:,} variants, {len(studies)} studies, "
          f"{meta['n_samples_seen']:,} distinct sample IDs")
    print(f"  five-gene rows: {len(rows)}")
    return rows, meta


# --------------------------------------------------------------------------- 5. SFARI


def sfari(refresh: bool) -> dict[str, str]:
    print("\n[5] SFARI Gene, live")
    hit = None if refresh else cached("denovo", "sfari_gene_status_*.json")
    if hit:
        out = json.loads(hit.read_text(encoding="utf-8"))
        print(f"  using cached {core.rel(hit)}")
        return out
    out = {}
    release = ""
    for gene in GENES:
        try:
            html = http_get(SFARI_GENE.format(gene=gene), timeout=120)
        except Exception as exc:
            out[gene] = f"HATA: {type(exc).__name__}"
            continue
        txt = re.sub(r"<[^>]+>", " ", re.sub(r"<script.*?</script>|<style.*?</style>", "", html,
                                             flags=re.S | re.I))
        txt = re.sub(r"\s+", " ", txt)
        m = re.search(r"Latest Release:\s*(\S+\s*\S*)", txt)
        if m:
            release = m.group(1).strip()
        if f"{gene} is not in the SFARI:Gene database" in txt:
            out[gene] = "listede degil"
        else:
            score = re.search(r"Gene Score\s*:?\s*(\d+S?)\s*([A-Za-z ,]{0,40})", txt)
            out[gene] = (f"listede (score {score.group(1)}, {score.group(2).strip()})"
                         if score else "listede (score okunamadi)")
        print(f"  {gene:8s} {out[gene]}")
        time.sleep(0.3)
    out["_release"] = release
    core.save_json("denovo", f"sfari_gene_status_{core.today()}.json", out, base=core.RAW,
                   source="SFARI Gene human gene pages", tool="scripts/denovo_axis.py",
                   urls=[SFARI_GENE.format(gene=g) for g in GENES])
    return out


# --------------------------------------------------------------------------- 6. analysis


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="re-download everything")
    args = ap.parse_args()

    core.load_dotenv()
    OUT.mkdir(parents=True, exist_ok=True)
    stamp = core.today()

    rates = mutation_rates(args.refresh)
    length_notes = sanity_check_lengths(rates)
    print("\n  sanity checks:")
    for n in length_notes:
        print("   ", n)

    kap_dnms, kap_gene, kap_meta = kaplanis(args.refresh)
    kap_kif13a_rows = kaplanis_build(kap_dnms, args.refresh)
    g4d2_var, g4d2_gene = gene4denovo2(args.refresh)
    ddb_rows, ddb_meta = denovodb(args.refresh)
    sfari_status = sfari(args.refresh)

    # ---------------------------------------------------------------- 6.1 Kaplanis per-gene counts
    print("\n[6] Kaplanis per-gene de novo counts and Poisson test")
    kap_counts: dict[str, dict] = {}
    for g in GENES:
        row = kap_gene.get(g, {})
        def n(key):
            v = row.get(key)
            try:
                return int(v)
            except (TypeError, ValueError):
                return 0
        lof = n("frameshift_variant") + n("stop_gained") + n("splice_acceptor_variant") + n("splice_donor_variant")
        mis = n("missense_variant")
        syn = n("synonymous_variant")
        r = rates.get(g, {})
        exp_lof = 2 * KAPLANIS_TRIOS * r.get("mu_lof", 0.0)
        exp_mis = 2 * KAPLANIS_TRIOS * r.get("mu_mis", 0.0)
        exp_syn = 2 * KAPLANIS_TRIOS * r.get("mu_syn", 0.0)
        kap_counts[g] = {
            "gen": g,
            "obs_LoF": lof, "exp_LoF": exp_lof,
            "p_LoF": poisson_p_ge(lof, exp_lof),
            "oran_LoF": (lof / exp_lof) if exp_lof else None,
            "obs_missense": mis, "exp_missense": exp_mis,
            "p_missense": poisson_p_ge(mis, exp_mis),
            "oran_missense": (mis / exp_mis) if exp_mis else None,
            "obs_sinonim": syn, "exp_sinonim": exp_syn,
            "p_sinonim": poisson_p_ge(syn, exp_syn),
            "frameshift": n("frameshift_variant"), "stop_gained": n("stop_gained"),
            "splice_acceptor": n("splice_acceptor_variant"), "splice_donor": n("splice_donor_variant"),
            "inframe_del": n("inframe_deletion"), "inframe_ins": n("inframe_insertion"),
            "denovoWEST_p": row.get("denovoWEST_p_full"),
            "PTV_enrichment_p": row.get("PTV_enrichment_p"),
            "s_het": row.get("s_het"),
            "kaplanis_significant": row.get("significant"),
            "ust_sinir_LoF_95": poisson_upper_95(lof),
            "ust_sinir_zenginlesme_LoF": (poisson_upper_95(lof) / exp_lof) if exp_lof else None,
        }
        c = kap_counts[g]
        print(f"  {g:8s} LoF {lof:3d} / bekl {exp_lof:6.3f} (p={c['p_LoF']:.3g})   "
              f"missense {mis:3d} / bekl {exp_mis:6.3f} (p={c['p_missense']:.3g})   "
              f"sinonim {syn:2d} / bekl {exp_syn:5.3f}   denovoWEST p={row.get('denovoWEST_p_full')}")

    core.save_tsv("denovo", f"kaplanis_per_gene_poisson_{stamp}.tsv",
                  [kap_counts[g] for g in GENES], base=core.INTERIM,
                  source=f"Kaplanis 2020 Supplementary Table 2 counts + gnomAD v2.1.1 mu rates, "
                         f"expected = 2 * {KAPLANIS_TRIOS} * mu",
                  urls=[KAPLANIS_SUPPL, GNOMAD_RATES], tool="scripts/denovo_axis.py",
                  params={"n_trios": KAPLANIS_TRIOS, "pmid": KAPLANIS_PMID})

    # ---------------------------------------------------------------- 6.2 Kaplanis raw rows
    kap_rows_out = []
    for d in kap_dnms:
        kap_rows_out.append({
            "gen": d["symbol"], "id": d["id"], "study": d["study"],
            "build": "GRCh37", "chrom": d["chrom"], "pos": d["pos"],
            "ref": d["ref"], "alt": d["alt"], "consequence": d["consequence"],
            "altprop_child": d["altprop_child"],
        })
    core.save_tsv("denovo", f"kaplanis_five_gene_dnms_{stamp}.tsv", kap_rows_out, base=core.INTERIM,
                  source="Kaplanis 2020 Supplementary Table 1, rows for the five project genes",
                  urls=[KAPLANIS_SUPPL], tool="scripts/denovo_axis.py")
    if kap_kif13a_rows:
        core.save_tsv("denovo", f"kaplanis_KIF13A_mane_rows_{stamp}.tsv", kap_kif13a_rows,
                      base=core.INTERIM,
                      source="Kaplanis 2020 KIF13A de novo rows lifted to GRCh38 and annotated on MANE",
                      urls=[KAPLANIS_SUPPL, f"{ENSEMBL}/vep/human/region"], tool="scripts/denovo_axis.py")

    # ---------------------------------------------------------------- 6.3 Gene4Denovo2 tables
    print("\n[7] Gene4Denovo2 per-gene breakdown")
    g4d2_rows_out = []
    g4d2_summary: dict[str, Counter] = {}
    for g in GENES:
        cnt: Counter = Counter()
        for rec in g4d2_var.get(g, []):
            cls = g4d2_class(rec)
            pheno = rec.get("Phenotype_Short") or ""
            cnt[(pheno, cls)] += 1
            pos, tx, chunk = g4d2_protein_pos(rec)
            g4d2_rows_out.append({
                "gen": g,
                "kohort_calismasi": rec.get("Study"),
                "PMID": rec.get("PMID"),
                "fenotip": rec.get("Phenotype"),
                "fenotip_kisa": pheno,
                "platform": rec.get("Platform"),
                "sample_id": rec.get("Sample_ID"),
                "grch38": f"{rec.get('Chr')}:{rec.get('Start')}{rec.get('Ref')}>{rec.get('Alt')}",
                "Func_refGene": rec.get("Func_refGene"),
                "ExonicFunc_mrefGene": rec.get("ExonicFunc_mrefGene"),
                "sinif": cls,
                "AAChange_mrefGene": rec.get("AAChange_mrefGene"),
                "protein_pos": pos,
                "protein_pos_transkripti": tx,
                "domain": domain_of(pos) if g == "KIF13A" else "",
                "ACMG": rec.get("ACMG"),
            })
        g4d2_summary[g] = cnt
        coding = sum(v for (ph, cl), v in cnt.items() if cl in ("LoF", "missense", "sinonim"))
        lof_case = sum(v for (ph, cl), v in cnt.items() if cl == "LoF" and ph != "Control")
        lof_ctrl = sum(v for (ph, cl), v in cnt.items() if cl == "LoF" and ph == "Control")
        mis_case = sum(v for (ph, cl), v in cnt.items() if cl == "missense" and ph != "Control")
        print(f"  {g:8s} toplam {len(g4d2_var.get(g, [])):4d} DNM, kodlayan {coding:3d}  "
              f"LoF(vaka) {lof_case}  LoF(kontrol) {lof_ctrl}  missense(vaka) {mis_case}")

    core.save_tsv("denovo", f"gene4denovo2_five_gene_dnms_{stamp}.tsv", g4d2_rows_out,
                  base=core.INTERIM,
                  source="Gene4Denovo2 variant-level rows for the five project genes (hg38)",
                  urls=[G4D2_API.format(table=G4D2_VARIANT_TABLE)], tool="scripts/denovo_axis.py")

    g4d2_gene_out = []
    for g in GENES:
        for rec in g4d2_gene.get(g, []):
            rec = dict(rec)
            rec["_gen"] = g
            g4d2_gene_out.append(rec)
    core.save_tsv("denovo", f"gene4denovo2_gene_level_{stamp}.tsv", g4d2_gene_out, base=core.INTERIM,
                  source="Gene4Denovo2 gene-level (Candidate_Gene2) summary rows",
                  urls=[G4D2_API.format(table=G4D2_GENE_TABLE)], tool="scripts/denovo_axis.py")

    # ---------------------------------------------------------------- 6.4 denovo-db
    print("\n[8] denovo-db per-gene breakdown")
    print("  NOT: denovo-db her varyanti HER transkript icin bir satir olarak yazar, yani ham satir")
    print("  sayisi varyant sayisi degildir. Asagisi (orneklem, konum, alel) uzerinden tekillestirilmistir.")
    ddb_out = []
    ddb_summary: dict[str, Counter] = {}
    ddb_probands: dict[str, int] = {}
    for d in ddb_rows:
        try:
            ddb_probands[d["StudyName"]] = max(ddb_probands.get(d["StudyName"], 0), int(d["NumProbands"]))
        except (TypeError, ValueError):
            pass
    for g in GENES:
        cnt: Counter = Counter()
        seen: set = set()
        mane_rs = core.unversioned(core.target(g)["refseq_transcript"]) if g in ("KIF13A", "KIF1A", "KIF13B") else None
        # first pass: for each unique (sample, position, allele) pick the MANE-transcript row if present
        by_var: dict[tuple, list[dict]] = defaultdict(list)
        for d in ddb_rows:
            if d["Gene"] != g:
                continue
            by_var[(d["SampleID"], d["Chr"], d["Position"], d["Variant"])].append(d)
        for key, group in by_var.items():
            pick = None
            if mane_rs:
                pick = next((x for x in group if core.unversioned(x["Transcript"] or "") == mane_rs), None)
            d = pick or group[0]
            fc = d["FunctionClass"]
            cls = ("LoF" if fc in DDB_LOF else "missense" if fc in DDB_MIS
                   else "sinonim" if fc in DDB_SYN else f"diger ({fc})")
            pheno = d["PrimaryPhenotype"]
            cnt[(pheno, cls)] += 1
            seen.add(key)
            m = PROT_RE.search(d.get("ProteinVariant") or "")
            pos = int(m.group(2)) if m else None
            on_mane = bool(mane_rs and core.unversioned(d["Transcript"] or "") == mane_rs)
            ddb_out.append({
                "gen": g, "sample_id": d["SampleID"], "study": d["StudyName"],
                "PMID": d["PubmedID"], "fenotip": pheno,
                "sekans_tipi": d["SequenceType"], "num_probands": d["NumProbands"],
                "build": "GRCh37 (denovo-db v1.6.1)",
                "pozisyon": f"chr{d['Chr']}:{d['Position']}{d['Variant']}",
                "transkript": d["Transcript"],
                "transkript_MANE_mi": "evet" if on_mane else "hayir",
                "FunctionClass": fc, "sinif": cls,
                "cDnaVariant": d["cDnaVariant"], "ProteinVariant": d["ProteinVariant"],
                "protein_pos": pos if on_mane else None,
                "protein_pos_kaynak_transkriptte": pos,
                "n_transkript_satiri": len(group),
            })
        ddb_summary[g] = cnt
        coding = sum(v for (ph, cl), v in cnt.items() if cl in ("LoF", "missense", "sinonim"))
        lof_case = sum(v for (ph, cl), v in cnt.items()
                       if cl == "LoF" and ph not in ("control", "controlNeuro"))
        raw = sum(1 for d in ddb_rows if d["Gene"] == g)
        print(f"  {g:8s} ham satir {raw:4d} -> tekil varyant {sum(cnt.values()):3d}, "
              f"kodlayan {coding:3d}  LoF(vaka) {lof_case}")

    core.save_tsv("denovo", f"denovodb_five_gene_dnms_{stamp}.tsv", ddb_out, base=core.INTERIM,
                  source="denovo-db v1.6.1 non-SSC rows for the five project genes",
                  urls=[DENOVODB_TSV], tool="scripts/denovo_axis.py")

    # ---------------------------------------------------------------- 6.5 KIF13A union + motor
    print("\n[9] KIF13A kodlayan de novo varyantlarin kaynaklar arasi birlesimi (GRCh38, MANE)")
    print("  Kaynaklar ortusuyor (ayni birey birden fazla yayinda olabilir), o yuzden birlesim")
    print("  (GRCh38 konum, ref, alt) uzerinden tekillestirilir ve hangi kaynak/orneklem destekliyorsa yazilir.")

    union: dict[tuple, dict] = {}

    def add(key, *, kaynak, study, pmid, pheno, sample, hgvsp, pos, cls, note=""):
        e = union.setdefault(key, {"grch38": key[0], "ref": key[1], "alt": key[2],
                                   "hgvsp_MANE": hgvsp, "protein_pos_MANE": pos,
                                   "domain": domain_of(pos), "sinif": cls,
                                   "kaynaklar": set(), "calismalar": set(), "pmidler": set(),
                                   "fenotipler": set(), "ornekler": set(), "notlar": set()})
        e["kaynaklar"].add(kaynak)
        if study:
            e["calismalar"].add(str(study))
        if pmid:
            e["pmidler"].add(str(pmid))
        if pheno:
            e["fenotipler"].add(str(pheno))
        if sample:
            e["ornekler"].add(str(sample))
        if note:
            e["notlar"].add(note)
        if pos and not e["protein_pos_MANE"]:
            e.update({"protein_pos_MANE": pos, "domain": domain_of(pos), "hgvsp_MANE": hgvsp})

    KIF13A_POS_RE = re.compile(r"chr?6:(\d+)")

    for r in kap_kif13a_rows:
        if r["sinif"] not in ("missense", "LoF", "sinonim"):
            continue
        m = KIF13A_POS_RE.search(r["grch38"])
        if not m:
            continue
        ref, alt = r["grch38"].split(":")[1][len(m.group(1)):].split(">")
        add((int(m.group(1)), ref, alt), kaynak="Kaplanis 2020", study=r["study"],
            pmid=KAPLANIS_PMID, pheno="gelisimsel bozukluk (DD)", sample=r["id"],
            hgvsp=r["hgvsp_mane"], pos=r["protein_pos_MANE"], cls=r["sinif"])

    for row in g4d2_rows_out:
        if row["gen"] != "KIF13A" or row["sinif"] not in ("missense", "LoF", "sinonim"):
            continue
        m = KIF13A_POS_RE.search(row["grch38"] or "")
        if not m:
            continue
        tail = (row["grch38"].split(":")[1])[len(m.group(1)):]
        ref, _, alt = tail.partition(">")
        add((int(m.group(1)), ref, alt), kaynak="Gene4Denovo2", study=row["kohort_calismasi"],
            pmid=row["PMID"], pheno=row["fenotip_kisa"], sample=row["sample_id"],
            hgvsp=row["AAChange_mrefGene"], pos=row["protein_pos"], cls=row["sinif"],
            note=f"ACMG={row.get('ACMG')}" if row.get("ACMG") else "")

    # denovo-db is GRCh37: lift its KIF13A coding rows to GRCh38 so they join the same key space
    ddb_kif13a = [r for r in ddb_out if r["gen"] == "KIF13A" and r["sinif"] in ("missense", "LoF", "sinonim")]
    lift_cache_path = cached("denovo", "denovodb_KIF13A_liftover_*.json") if not args.refresh else None
    if lift_cache_path:
        lift = json.loads(lift_cache_path.read_text(encoding="utf-8"))
        print(f"  using cached {core.rel(lift_cache_path)}")
    else:
        lift = {}
        for r in ddb_kif13a:
            m = re.search(r"chr6:(\d+)", r["pozisyon"])
            if not m:
                continue
            p37 = m.group(1)
            if p37 in lift:
                continue
            url = f"{ENSEMBL}/map/human/GRCh37/6:{p37}..{p37}/GRCh38?content-type=application/json"
            try:
                mm = json.loads(http_get(url, headers={"Content-Type": "application/json"}, timeout=120))
                maps = mm.get("mappings") or []
                lift[p37] = maps[0]["mapped"]["start"] if maps else None
            except Exception:
                lift[p37] = None
            time.sleep(0.15)
        core.save_json("denovo", f"denovodb_KIF13A_liftover_{stamp}.json", lift, base=core.RAW,
                       source="Ensembl REST GRCh37 -> GRCh38 map for denovo-db KIF13A coding rows",
                       urls=[f"{ENSEMBL}/map/human/GRCh37"], tool="scripts/denovo_axis.py")

    for r in ddb_kif13a:
        m = re.search(r"chr6:(\d+)(.*)$", r["pozisyon"])
        if not m:
            continue
        p38 = lift.get(m.group(1))
        if not p38:
            continue
        allele = m.group(2)
        ref, _, alt = allele.partition(">")
        add((int(p38), ref, alt), kaynak="denovo-db", study=r["study"], pmid=r["PMID"],
            pheno=r["fenotip"], sample=r["sample_id"], hgvsp=r["ProteinVariant"],
            pos=r["protein_pos"], cls=r["sinif"])

    union_rows = []
    for key in sorted(union, key=lambda k: union[k]["protein_pos_MANE"] or 0):
        e = union[key]
        union_rows.append({
            "grch38": f"chr6:{e['grch38']}{e['ref']}>{e['alt']}",
            "hgvsp_MANE": e["hgvsp_MANE"],
            "protein_pos_MANE": e["protein_pos_MANE"],
            "domain": e["domain"],
            "sinif": e["sinif"],
            "motor_domainde": "evet" if (e["protein_pos_MANE"] and
                                         MOTOR_START <= e["protein_pos_MANE"] <= MOTOR_END) else "hayir",
            "kaynaklar": ";".join(sorted(e["kaynaklar"])),
            "calismalar": ";".join(sorted(e["calismalar"])),
            "pmidler": ";".join(sorted(e["pmidler"])),
            "fenotipler": ";".join(sorted(e["fenotipler"])),
            "ornekler": ";".join(sorted(e["ornekler"])),
            "n_ornek": len(e["ornekler"]),
            "notlar": ";".join(sorted(e["notlar"])),
        })
    core.save_tsv("denovo", f"KIF13A_denovo_union_{stamp}.tsv", union_rows, base=core.INTERIM,
                  source="Deduplicated union of KIF13A coding de novo variants across Kaplanis 2020, "
                         "Gene4Denovo2 and denovo-db, all on GRCh38 / MANE ENST00000259711.11",
                  urls=[KAPLANIS_SUPPL, G4D2_API.format(table=G4D2_VARIANT_TABLE), DENOVODB_TSV],
                  tool="scripts/denovo_axis.py")

    for r in union_rows:
        print(f"  {r['grch38']:24s} {str(r['protein_pos_MANE']):>5} {r['domain']:8s} {r['sinif']:9s} "
              f"{'MOTOR' if r['motor_domainde'] == 'evet' else '     '} "
              f"n_ornek={r['n_ornek']} [{r['fenotipler']}] {r['kaynaklar']}")

    coding_union = [r for r in union_rows if r["sinif"] in ("missense", "LoF")]
    motor_union = [r for r in coding_union if r["motor_domainde"] == "evet"]
    all_kif13a_coding = [(r["kaynaklar"], r["hgvsp_MANE"], r["protein_pos_MANE"], r["domain"], r["sinif"])
                         for r in coding_union]
    motor_hits = [(r["kaynaklar"], r["hgvsp_MANE"], r["protein_pos_MANE"], r["sinif"]) for r in motor_union]
    motor_frac = (MOTOR_END - MOTOR_START + 1) / KIF13A_LEN
    p_motor = 1 - (1 - motor_frac) ** len(coding_union) if coding_union else None
    print(f"\n  motor alan (5-352) = {MOTOR_END - MOTOR_START + 1} kalinti = KIF13A'nin "
          f"%{motor_frac*100:.1f}'i (gorev metnindeki 352/1805 = %19,5 yaklasimi)")
    print(f"  tekillestirilmis kodlayan (missense+LoF) de novo: {len(coding_union)}, "
          f"motorda: {len(motor_union)}")
    if p_motor is not None:
        print(f"  uzunluk oranina gore en az bir motor isabeti gorme olasiligi: {p_motor:.3f} "
              f"(yani {len(motor_union)} isabet surpriz degil)")

    # ---------------------------------------------------------------- 6.6 known observations
    print("\n[10] Bilinen uc KIF13A gozlemi bu kohortlarda var mi")
    known_report = []
    for label, spec in KNOWN.items():
        want_pos = spec.get("protein_pos")
        found = []
        for src, hgvsp, pos, dom, cls in all_kif13a_coding:
            if want_pos and pos == want_pos:
                found.append(f"{src}: {hgvsp}")
        # also scan every KIF13A row (including non-coding) for the exact HGVSc
        hgvsc = spec.get("hgvsc")
        if hgvsc:
            needle = hgvsc.replace(">", "").replace("c.", "")
            for row in g4d2_rows_out:
                if row["gen"] != "KIF13A":
                    continue
                aa = (row.get("AAChange_mrefGene") or "")
                if hgvsc in aa or (needle and needle in aa.replace(">", "")):
                    found.append(f"Gene4Denovo2: {aa}")
        known_report.append({"gozlem": label, **{k: v for k, v in spec.items()},
                             "kohortlarda_bulundu": "; ".join(found) if found else "HAYIR"})
        print(f"  {label:32s} -> {'; '.join(found) if found else 'bulunamadi'}")

    core.save_tsv("denovo", f"known_observations_crosscheck_{stamp}.tsv", known_report,
                  base=core.INTERIM,
                  source="Cross-check of the three known human KIF13A observations against the "
                         "assembled de novo cohorts",
                  urls=[G4D2_API.format(table=G4D2_VARIANT_TABLE), KAPLANIS_SUPPL, DENOVODB_TSV],
                  tool="scripts/denovo_axis.py")

    # ---------------------------------------------------------------- 6.7 five-gene table
    five = []
    for g in GENES:
        c = kap_counts[g]
        g4 = g4d2_summary[g]
        dd = ddb_summary[g]
        five.append({
            "gen": g,
            "rol": {"KIF13A": "hedef", "KIF1A": "paralog, missense mekanizmasi",
                    "KIF13B": "paralog, null beklenir", "ANKRD11": "pozitif kontrol (HI=3)",
                    "CHD8": "pozitif kontrol (HI=3)"}[g],
            "kaplanis_obs_LoF": c["obs_LoF"], "kaplanis_exp_LoF": round(c["exp_LoF"], 3),
            "kaplanis_p_LoF": c["p_LoF"],
            "kaplanis_obs_mis": c["obs_missense"], "kaplanis_exp_mis": round(c["exp_missense"], 3),
            "kaplanis_p_mis": c["p_missense"],
            "kaplanis_obs_syn": c["obs_sinonim"], "kaplanis_exp_syn": round(c["exp_sinonim"], 3),
            "denovoWEST_p": c["denovoWEST_p"], "kaplanis_significant": c["kaplanis_significant"],
            "g4d2_LoF_vaka": sum(v for (ph, cl), v in g4.items() if cl == "LoF" and ph != "Control"),
            "g4d2_LoF_kontrol": sum(v for (ph, cl), v in g4.items() if cl == "LoF" and ph == "Control"),
            "g4d2_missense_vaka": sum(v for (ph, cl), v in g4.items() if cl == "missense" and ph != "Control"),
            "g4d2_toplam_DNM": len(g4d2_var.get(g, [])),
            "denovodb_LoF_vaka": sum(v for (ph, cl), v in dd.items()
                                     if cl == "LoF" and ph not in ("control", "controlNeuro")),
            "denovodb_missense_vaka": sum(v for (ph, cl), v in dd.items()
                                          if cl == "missense" and ph not in ("control", "controlNeuro")),
            "sfari": sfari_status.get(g, ""),
            "mu_lof": rates.get(g, {}).get("mu_lof"),
            "mu_mis": rates.get(g, {}).get("mu_mis"),
        })
    core.save_tsv("denovo", f"five_gene_comparison_{stamp}.tsv", five, base=core.INTERIM,
                  source="Five-gene de novo comparison across Kaplanis 2020, Gene4Denovo2, denovo-db, SFARI",
                  urls=[KAPLANIS_SUPPL, G4D2_API.format(table=G4D2_VARIANT_TABLE), DENOVODB_TSV],
                  tool="scripts/denovo_axis.py")

    print("\n[11] Bes gen karsilastirmasi (Kaplanis, 31.058 DD triosu)")
    print(f"  {'gen':9s} {'obsLoF':>7s} {'beklLoF':>8s} {'p':>10s} {'obsMis':>7s} {'beklMis':>8s} {'p':>10s} {'denovoWEST':>12s}")
    for r in five:
        print(f"  {r['gen']:9s} {r['kaplanis_obs_LoF']:7d} {r['kaplanis_exp_LoF']:8.3f} "
              f"{r['kaplanis_p_LoF']:10.3g} {r['kaplanis_obs_mis']:7d} {r['kaplanis_exp_mis']:8.3f} "
              f"{r['kaplanis_p_mis']:10.3g} {str(r['denovoWEST_p'])[:12]:>12s}")

    # ------------------------------------------------- 6.7 aggregate, calibrated on synonymous DNMs
    print("\n[11b] Gene4Denovo2 toplamı: etkin trio sayisi sinonim de novo'lardan kalibre edilir")
    print("  Gerekce: sinonim de novo varyantlar notrdur, yani sayilari yalniz maruziyeti")
    print("  (kac trio dizilendi) olcer. Bes genin sinonim toplamini ayni bes genin sinonim")
    print("  mutasyon hizi toplamina bolmek, kohortun etkin buyuklugunu verir.")
    g4_lvl: dict[str, dict] = {}
    for g in GENES:
        for rec in g4d2_gene.get(g, []):
            if rec.get("Phenotype_Short") == "ALL Disorders":
                g4_lvl[g] = rec
    syn_tot = sum(int(g4_lvl[g].get("Syn") or 0) for g in GENES if g in g4_lvl)
    mu_syn_tot = sum(rates[g]["mu_syn"] for g in GENES if g in g4_lvl)
    n_eff = syn_tot / (2 * mu_syn_tot) if mu_syn_tot else 0.0
    # Poisson uncertainty on the synonymous total propagates straight to n_eff
    n_eff_lo = poisson_lower_95(syn_tot) / (2 * mu_syn_tot) if mu_syn_tot else 0.0
    n_eff_hi = poisson_upper_95(syn_tot) / (2 * mu_syn_tot) if mu_syn_tot else 0.0
    print(f"  bes gende 'ALL Disorders' sinonim de novo toplami : {syn_tot}")
    print(f"  ayni bes genin mu_syn toplami                     : {mu_syn_tot:.6g}")
    print(f"  etkin trio sayisi (tahmin)                        : {n_eff:,.0f} "
          f"(%95 {n_eff_lo:,.0f}-{n_eff_hi:,.0f})")
    print(f"  Gene4Denovo2'nin kendi beyani                     : 130.439 birey, 96 fenotip")

    agg = []
    for g in GENES:
        rec = g4_lvl.get(g, {})
        lof = int(rec.get("LoF") or 0)
        syn = int(rec.get("Syn") or 0)
        exp_lof = 2 * n_eff * rates[g]["mu_lof"]
        agg.append({
            "gen": g,
            "g4d2_LoF_ALL_Disorders": lof,
            "g4d2_Syn_ALL_Disorders": syn,
            "beklenen_LoF_n_eff": round(exp_lof, 3),
            "zenginlesme": round(lof / exp_lof, 1) if exp_lof else None,
            "poisson_p": poisson_p_ge(lof, exp_lof),
            "g4d2_kendi_Pval": rec.get("Pval"),
            "g4d2_kendi_Qval": rec.get("Qval"),
            "g4d2_guven": rec.get("Confidence_level"),
            "ust_sinir_zenginlesme_95": round(poisson_upper_95(lof) / exp_lof, 1) if exp_lof else None,
        })
        print(f"  {g:8s} LoF {lof:4d} / bekl {exp_lof:6.3f} = {agg[-1]['zenginlesme'] or 0:6.1f}x  "
              f"p={agg[-1]['poisson_p']:.3g}  (G4D2 kendi p={rec.get('Pval')}, {rec.get('Confidence_level')})")
    core.save_tsv("denovo", f"gene4denovo2_aggregate_poisson_{stamp}.tsv", agg, base=core.INTERIM,
                  source="Gene4Denovo2 ALL Disorders LoF counts against expectation at an effective "
                         "trio count calibrated from synonymous de novo variants of the same five genes",
                  urls=[G4D2_API.format(table=G4D2_GENE_TABLE), GNOMAD_RATES],
                  tool="scripts/denovo_axis.py",
                  params={"n_eff": n_eff, "syn_total": syn_tot, "mu_syn_total": mu_syn_tot})

    # ---------------------------------------------------------------- 6.75 power / sensitivity
    print("\n[12] Duyarlilik: kohort buyudukce 0 gozlem neyi disliyor")
    mu_lof = rates["KIF13A"]["mu_lof"]
    mu_mis = rates["KIF13A"]["mu_mis"]
    sens = []
    for n_trio in (31058, 50000, 100000, 200000, 500000):
        exp_lof = 2 * n_trio * mu_lof
        upper = poisson_upper_95(0)
        sens.append({
            "n_trio (VARSAYIM)": n_trio,
            "beklenen_LoF": round(exp_lof, 4),
            "0_gozlemin_95_ust_siniri_olay": round(upper, 3),
            "dislanan_zenginlesme_kati": round(upper / exp_lof, 2) if exp_lof else None,
            "ANKRD11_benzeri_360x_icin_beklenen": round(exp_lof * 360, 1),
        })
        print(f"  {n_trio:>7,} trio -> beklenen LoF {exp_lof:6.3f}; 0 gozlem "
              f"{upper/exp_lof:6.1f}x uzeri zenginlesmeyi disliyor")
    print("  NOT: yalniz ilk satir (31.058) olculmus bir kohorttur; digerleri varsayimsal ve")
    print("  yalniz 'bu eksen ne kadar ileri goturulebilir' sorusunu yanitlamak icindir.")
    core.save_tsv("denovo", f"power_sensitivity_{stamp}.tsv", sens, base=core.INTERIM,
                  source="Poisson power sensitivity for KIF13A de novo LoF at assumed cohort sizes; "
                         "only n=31058 is a measured cohort",
                  urls=[GNOMAD_RATES], tool="scripts/denovo_axis.py",
                  params={"mu_lof": mu_lof, "mu_mis": mu_mis})

    # ---------------------------------------------------------------- 6.8 machine-readable summary
    summary = {
        "tarih_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc)
                     .isoformat(timespec="seconds"),
        "kaplanis": kap_meta,
        "denovodb": ddb_meta,
        "sfari": sfari_status,
        "rates": rates,
        "length_sanity": length_notes,
        "kaplanis_per_gene": kap_counts,
        "five_gene": five,
        "kif13a_union": union_rows,
        "motor": {"fraction_of_cds": motor_frac, "n_coding_denovo": len(all_kif13a_coding),
                  "n_in_motor": len(motor_hits), "p_at_least_one_by_length": p_motor,
                  "hits": [{"kaynak": s, "hgvsp": h, "pos": p, "sinif": c} for s, h, p, c in motor_hits]},
        "known_observations": known_report,
        "power_sensitivity": sens,
        "denovodb_probands_per_study": ddb_probands,
    }
    core.save_json("denovo", f"denovo_axis_summary_{stamp}.json", summary, base=core.INTERIM,
                   source="Machine-readable summary of the de novo axis",
                   urls=[KAPLANIS_SUPPL, G4D2_API.format(table=G4D2_VARIANT_TABLE), DENOVODB_TSV,
                         GNOMAD_RATES], tool="scripts/denovo_axis.py")

    k = kap_counts["KIF13A"]
    print("\n" + "=" * 78)
    print("KIF13A, Kaplanis 2020, 31.058 gelisimsel bozukluk triosu:")
    print(f"  gozlenen de novo LoF      : {k['obs_LoF']}")
    print(f"  beklenen de novo LoF      : {k['exp_LoF']:.3f}")
    print(f"  Poisson p (tek yonlu)     : {k['p_LoF']:.3g}")
    print(f"  0 gozlemin %95 ust siniri : {k['ust_sinir_LoF_95']:.2f} olay "
          f"= {k['ust_sinir_zenginlesme_LoF']:.1f}x zenginlesme")
    print(f"  gozlenen/beklenen missense: {k['obs_missense']} / {k['exp_missense']:.3f} "
          f"(p={k['p_missense']:.3g})")
    print("=" * 78)
    print(f"\nCiktilar: {core.rel(OUT)}")


if __name__ == "__main__":
    main()
