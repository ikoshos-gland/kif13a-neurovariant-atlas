"""Scan human paralog pairs of the KIF1A -> KIF13A shape: one sibling clinically rich, the other clinically empty.

The atlas calibrates variant-effect models on KIF1A (135 ClinVar P/LP missense) and transfers the thresholds
plus the paralog-conservation feature to KIF13A (0 P/LP missense). That is one case study. This script asks
how many other human paralog pairs have the same shape, and measures - cheaply, without running the pipeline -
how strong the transfer signal would be in each of them.

Stages (each cached; rerun a single one with --stage):
  sources  ClinVar gene_specific_summary.txt + HGNC complete set (symbol <-> ENSG <-> Entrez)
  rich     exact ClinVar P/LP missense count per gene (E-utilities esearch), prefiltered by the summary table
  pairs    Ensembl REST homology (type=paralogues) for every rich gene, then the sibling's P/LP missense count
  metrics  per-pair feasibility: MANE lengths, global identity, transfer signal, gnomAD constraint, clustering
  rank     composite ranking plus explicit weakness flags

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\paralog_pairs_scan.py --stage all

TRAP (verified 2026-09-16): `pathogenic[clinsig]` is NOT a ClinVar Entrez field. Entrez silently rewrites it
to `pathogenic[All Fields]`, which matched 456 KIF1A missense records instead of the true 135 P/LP. The only
form that reproduces the curated count is `("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])`.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

core.load_dotenv()

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENSEMBL = "https://rest.ensembl.org"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
CLINVAR_GENE_SUMMARY = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/gene_specific_summary.txt"
HGNC_COMPLETE = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"

OUT = "paralog_pairs"                 # data/interim/paralog_pairs
CACHE = core.RAW / "paralog_scan"     # gitignored downloads + per-gene caches

PLP_PROPS = '("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])'
MISSENSE = '"missense variant"[molecular consequence]'

TITLE_RE = re.compile(r"^(NM_\d+)\.\d+\(([^)]+)\):(c\.[^ ]+)(?: \((p\.[^)]+)\))?")
HGVSP_RE = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


# --------------------------------------------------------------------------- tiny cache helpers


def cache_path(*parts: str) -> Path:
    p = CACHE.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def load_cache(name: str) -> dict[str, Any]:
    p = cache_path(name)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def dump_cache(name: str, obj: Any) -> None:
    cache_path(name).write_text(json.dumps(obj), encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def latest(pattern: str) -> Path:
    hits = sorted((core.INTERIM / OUT).glob(pattern))
    if not hits:
        raise SystemExit(f"missing input {pattern} in data/interim/{OUT}; run the earlier stage first")
    return hits[-1]


def num(x: Any) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------- stage: sources


async def stage_sources(args: argparse.Namespace) -> dict[str, Any]:
    """ClinVar per-gene summary and the HGNC symbol map. Both are single files, both get provenance."""
    out = {}
    for label, url, fname in [
        ("clinvar_gene_summary", CLINVAR_GENE_SUMMARY, "gene_specific_summary.txt"),
        ("hgnc_complete_set", HGNC_COMPLETE, "hgnc_complete_set.txt"),
    ]:
        dest = cache_path(fname)
        if dest.exists() and not args.refresh:
            out[label] = {"path": core.rel(dest), "cached": True, "bytes": dest.stat().st_size}
            continue
        content = await core.get_bytes(url, timeout=600)
        out[label] = core.save_bytes("paralog_scan", fname, content, source=label, urls=[url],
                                     tool="scripts/paralog_pairs_scan.py::sources", params={})
    # the ClinVar file states its own release date on line 1
    first = cache_path("gene_specific_summary.txt").read_text(encoding="utf-8", errors="replace").split("\n", 1)[0]
    out["clinvar_release_line"] = first.strip()
    return out


def hgnc_maps() -> tuple[dict[str, dict[str, str]], dict[str, str]]:
    """symbol -> hgnc row (protein-coding, approved) and ENSG -> symbol."""
    rows = read_tsv(cache_path("hgnc_complete_set.txt"))
    by_symbol, by_ensg = {}, {}
    for r in rows:
        if r.get("status") != "Approved" or r.get("locus_group") != "protein-coding gene":
            continue
        ensg = (r.get("ensembl_gene_id") or "").strip()
        sym = r["symbol"]
        by_symbol[sym] = {"symbol": sym, "ensembl_gene_id": ensg, "entrez_id": r.get("entrez_id", ""),
                          "gene_group": r.get("gene_group", ""), "name": r.get("name", "")}
        if ensg:
            by_ensg[ensg] = sym
    return by_symbol, by_ensg


def clinvar_gene_summary() -> dict[str, dict[str, Any]]:
    """symbol -> {'plp_alleles': int, 'total_alleles': int}. '-' means none reported."""
    text = cache_path("gene_specific_summary.txt").read_text(encoding="utf-8", errors="replace")
    lines = [l for l in text.splitlines() if l and not l.startswith("#Overview")]
    reader = csv.DictReader(lines, delimiter="\t")
    out = {}
    for r in reader:
        sym = (r.get("#Symbol") or "").strip()
        if not sym:
            continue
        plp = r.get("Alleles_reported_Pathogenic_Likely_pathogenic", "-")
        out[sym] = {"plp_alleles": int(plp) if plp not in ("-", "", None) else 0,
                    "total_alleles": int(r.get("Total_alleles") or 0)}
    return out


# --------------------------------------------------------------------------- ClinVar counting


async def clinvar_plp_missense_count(symbol: str, sem: asyncio.Semaphore) -> int:
    key = core.secret("NCBI_API_KEY")
    term = f"{symbol}[gene] AND single_gene[prop] AND {PLP_PROPS} AND {MISSENSE}"
    params = {"db": "clinvar", "retmode": "json", "retmax": 0, "term": term, **({"api_key": key} if key else {})}
    async with sem:
        data = await core.get_json(f"{EUTILS}/esearch.fcgi", params=params, timeout=120)
    res = data["esearchresult"]
    if res.get("errorlist", {}).get("phrasesnotfound"):
        # an unmatched gene symbol returns count 0 with a warning; treat it as 0 but keep it visible
        return 0
    return int(res["count"])


async def count_many(symbols: list[str], cache_name: str, args: argparse.Namespace) -> dict[str, int]:
    cached: dict[str, Any] = load_cache(cache_name)
    todo = [s for s in symbols if s not in cached or args.refresh]
    if todo:
        sem = asyncio.Semaphore(args.ncbi_concurrency)
        done = 0
        t0 = time.monotonic()
        for block in core.chunks(todo, 200):
            results = await asyncio.gather(*[clinvar_plp_missense_count(s, sem) for s in block],
                                           return_exceptions=True)
            for sym, res in zip(block, results):
                if isinstance(res, Exception):
                    print(f"    ! {sym}: {type(res).__name__} {res}", flush=True)
                    continue
                cached[sym] = res
            done += len(block)
            dump_cache(cache_name, cached)
            rate = done / max(time.monotonic() - t0, 1e-9)
            print(f"    ClinVar counts {done}/{len(todo)} ({rate:.1f}/s)", flush=True)
    return {s: int(cached[s]) for s in symbols if s in cached}


# --------------------------------------------------------------------------- stage: rich


async def stage_rich(args: argparse.Namespace) -> dict[str, Any]:
    by_symbol, _ = hgnc_maps()
    summary = clinvar_gene_summary()

    # Prefilter is lossless for the >=min_rich rule: P/LP missense alleles cannot exceed total P/LP alleles.
    universe = sorted(s for s, v in summary.items()
                      if v["plp_alleles"] >= args.min_rich and s in by_symbol and by_symbol[s]["ensembl_gene_id"])
    print(f"  prefilter: {len(universe)} approved protein-coding genes with >= {args.min_rich} P/LP alleles "
          f"in the ClinVar gene summary", flush=True)

    counts = await count_many(universe, "clinvar_plp_missense_counts.json", args)
    rows = []
    for sym in universe:
        n = counts.get(sym)
        if n is None or n < args.min_rich:
            continue
        rows.append({"symbol": sym, "ensembl_gene_id": by_symbol[sym]["ensembl_gene_id"],
                     "entrez_id": by_symbol[sym]["entrez_id"], "gene_group": by_symbol[sym]["gene_group"],
                     "clinvar_plp_missense": n, "clinvar_plp_alleles_summary": summary[sym]["plp_alleles"],
                     "clinvar_total_alleles_summary": summary[sym]["total_alleles"]})
    rows.sort(key=lambda r: -r["clinvar_plp_missense"])
    saved = core.save_tsv(OUT, f"paralog_scan_rich_genes_{core.today()}.tsv", rows, base=core.INTERIM,
                          source="NCBI E-utilities ClinVar esearch + ClinVar gene_specific_summary + HGNC",
                          urls=[f"{EUTILS}/esearch.fcgi", CLINVAR_GENE_SUMMARY, HGNC_COMPLETE],
                          tool="scripts/paralog_pairs_scan.py::rich",
                          params={"min_rich": args.min_rich, "term_template": f"<GENE>[gene] AND single_gene[prop] AND {PLP_PROPS} AND {MISSENSE}"})
    return {"n_prefiltered": len(universe), "n_rich": len(rows), "saved": saved,
            "top": [(r["symbol"], r["clinvar_plp_missense"]) for r in rows[:12]]}


# --------------------------------------------------------------------------- stage: pairs


async def ensembl_paralogs(symbol: str, timeout: float, retries: int) -> list[dict[str, Any]]:
    data = await core.get_json(
        f"{ENSEMBL}/homology/symbol/human/{symbol}",
        params={"type": "paralogues", "format": "condensed", "target_species": "homo_sapiens",
                "content-type": "application/json"}, timeout=timeout, retries=retries)
    blocks = data.get("data") or []
    return blocks[0].get("homologies", []) if blocks else []


async def fetch_homology(symbols: list[str], args: argparse.Namespace) -> dict[str, Any]:
    """Worker pool, not blocked batches: Compara paralogy for a big gene family can take minutes and one
    slow gene must not hold up the other workers. Two passes, the second with a longer timeout."""
    cached: dict[str, Any] = load_cache("ensembl_paralogs.json")
    todo = [s for s in symbols if s not in cached or args.refresh]
    if not todo:
        return cached
    for attempt, (timeout, retries) in enumerate([(args.ensembl_timeout, 2), (args.ensembl_timeout * 4, 3)], 1):
        if not todo:
            break
        queue: asyncio.Queue[str] = asyncio.Queue()
        for s in todo:
            queue.put_nowait(s)
        failed: list[str] = []
        done = [0]
        total, t0 = len(todo), time.monotonic()

        async def worker() -> None:
            while True:
                try:
                    sym = queue.get_nowait()
                except asyncio.QueueEmpty:
                    return
                try:
                    homs = await ensembl_paralogs(sym, timeout, retries)
                    cached[sym] = [{"id": h.get("id"), "type": h.get("type"),
                                    "taxonomy_level": h.get("taxonomy_level")} for h in homs]
                except Exception as exc:
                    failed.append(sym)
                    print(f"    ! homology {sym}: {type(exc).__name__}", flush=True)
                done[0] += 1
                if done[0] % 25 == 0:
                    dump_cache("ensembl_paralogs.json", cached)
                    el = time.monotonic() - t0
                    print(f"    Ensembl homology pass {attempt}: {done[0]}/{total}  {el:.0f}s elapsed, "
                          f"~{el / done[0] * (total - done[0]) / 60:.0f} min left", flush=True)

        await asyncio.gather(*[worker() for _ in range(args.ensembl_concurrency)])
        dump_cache("ensembl_paralogs.json", cached)
        print(f"    pass {attempt}: {total - len(failed)}/{total} resolved in "
              f"{time.monotonic() - t0:.0f}s, {len(failed)} failed", flush=True)
        todo = failed
    if todo:
        print(f"    ! homology never resolved for {len(todo)} genes: {', '.join(sorted(todo))}", flush=True)
        for s in todo:
            cached[s] = None
        dump_cache("ensembl_paralogs.json", cached)
    return cached


async def stage_pairs(args: argparse.Namespace) -> dict[str, Any]:
    rich = read_tsv(latest("paralog_scan_rich_genes_*.tsv"))
    _, by_ensg = hgnc_maps()
    summary = clinvar_gene_summary()

    cached = await fetch_homology([r["symbol"] for r in rich], args)
    n_failed = sum(1 for r in rich if cached.get(r["symbol"]) is None)

    # candidate pairs: sibling must be clinically empty in the summary table before we spend an esearch on it
    candidates: list[dict[str, Any]] = []
    sibling_symbols: set[str] = set()
    for r in rich:
        for h in (cached.get(r["symbol"]) or []):
            sib = by_ensg.get(h["id"] or "")
            if not sib or sib == r["symbol"]:
                continue
            s = summary.get(sib, {"plp_alleles": 0, "total_alleles": 0})
            if s["plp_alleles"] > args.sibling_prefilter:
                continue
            candidates.append({"rich_gene": r["symbol"], "rich_plp_missense": int(r["clinvar_plp_missense"]),
                               "sibling_gene": sib, "sibling_ensembl_gene_id": h["id"],
                               "paralog_type": h["type"], "duplication_level": h["taxonomy_level"],
                               "sibling_plp_alleles_summary": s["plp_alleles"],
                               "sibling_total_alleles_summary": s["total_alleles"]})
            sibling_symbols.add(sib)
    print(f"  {len(candidates)} candidate pairs, {len(sibling_symbols)} distinct siblings to count", flush=True)

    sib_counts = await count_many(sorted(sibling_symbols), "clinvar_plp_missense_counts.json", args)
    qualified = []
    for c in candidates:
        n = sib_counts.get(c["sibling_gene"])
        if n is None or n > args.max_sibling:
            continue
        qualified.append({**c, "sibling_plp_missense": n,
                          "asymmetry": c["rich_plp_missense"] - n})
    qualified.sort(key=lambda r: (-r["rich_plp_missense"], r["sibling_gene"]))

    saved_all = core.save_tsv(OUT, f"paralog_scan_pairs_all_{core.today()}.tsv", candidates, base=core.INTERIM,
                              source="Ensembl REST homology (paralogues, condensed) + HGNC + ClinVar gene summary",
                              urls=[f"{ENSEMBL}/homology/symbol/human/", HGNC_COMPLETE, CLINVAR_GENE_SUMMARY],
                              tool="scripts/paralog_pairs_scan.py::pairs",
                              params={"sibling_prefilter_plp_alleles": args.sibling_prefilter})
    saved_q = core.save_tsv(OUT, f"paralog_scan_pairs_qualified_{core.today()}.tsv", qualified, base=core.INTERIM,
                            source="Ensembl REST homology + ClinVar esearch counts for both siblings",
                            urls=[f"{ENSEMBL}/homology/symbol/human/", f"{EUTILS}/esearch.fcgi"],
                            tool="scripts/paralog_pairs_scan.py::pairs",
                            params={"min_rich": args.min_rich, "max_sibling": args.max_sibling})
    return {"n_candidate_pairs": len(candidates), "n_qualified_pairs": len(qualified),
            "n_rich_genes_without_homology": n_failed,
            "n_distinct_rich": len({q["rich_gene"] for q in qualified}),
            "n_distinct_siblings": len({q["sibling_gene"] for q in qualified}),
            "saved_all": saved_all, "saved_qualified": saved_q}


# --------------------------------------------------------------------------- gene facts (gnomAD, sequence)

GENE_FIELDS = ("symbol mane_select_transcript { ensembl_id refseq_id } "
               "gnomad_constraint { oe_lof_upper oe_lof pli mis_z oe_mis oe_mis_upper exp_mis obs_mis flags }")


async def gnomad_genes_batch(symbols: list[str], args: argparse.Namespace) -> dict[str, Any]:
    """MANE Select transcript + gnomAD v4.1 constraint. GraphQL aliases let one request carry many genes,
    which is both far faster and far politer than one request per gene."""
    cache: dict[str, Any] = load_cache("gnomad_genes.json")
    todo = [s for s in sorted(set(symbols)) if s not in cache or args.refresh]
    for i, block in enumerate(core.chunks(todo, args.gnomad_batch)):
        parts = "".join(f'g{j}: gene(gene_symbol: "{s}", reference_genome: GRCh38) {{ {GENE_FIELDS} }} '
                        for j, s in enumerate(block))
        try:
            data = await core.post_json(GNOMAD_API, {"query": "{ %s }" % parts}, timeout=300)
        except Exception as exc:
            print(f"    ! gnomAD block {i}: {type(exc).__name__}", flush=True)
            continue
        payload = data.get("data") or {}
        for j, s in enumerate(block):
            cache[s] = payload.get(f"g{j}") or {}
        dump_cache("gnomad_genes.json", cache)
        print(f"    gnomAD constraint {min((i + 1) * args.gnomad_batch, len(todo))}/{len(todo)}", flush=True)
    return cache


def cached_sequence(symbol: str, enst: str) -> str | None:
    path = cache_path("sequences", f"{symbol}_{enst}.fasta")
    if not path.exists():
        return None
    return "".join(path.read_text(encoding="utf-8").splitlines()[1:])


async def mane_proteins_batch(wanted: list[tuple[str, str]], args: argparse.Namespace) -> None:
    """POST /sequence/id takes a list of transcript ids, so 50 proteins come back per request."""
    todo = [(sym, enst) for sym, enst in sorted(set(wanted)) if cached_sequence(sym, enst) is None]
    by_enst = {enst: sym for sym, enst in todo}
    ids = sorted(by_enst)
    for i, block in enumerate(core.chunks(ids, args.sequence_batch)):
        try:
            resp = await core.request("POST", f"{ENSEMBL}/sequence/id", params={"type": "protein"},
                                      json_body={"ids": block},
                                      headers={"Content-Type": "application/json", "Accept": "application/json"},
                                      timeout=300)
            entries = resp.json()
        except Exception as exc:
            print(f"    ! sequence block {i}: {type(exc).__name__}", flush=True)
            continue
        for e in entries:
            enst = e.get("query") or e.get("id")
            sym = by_enst.get(enst)
            if not sym or not e.get("seq"):
                continue
            seq = e["seq"].replace("*", "")
            cache_path("sequences", f"{sym}_{enst}.fasta").write_text(
                f">{sym}|{enst}|MANE\n{seq}\n", encoding="utf-8")
        print(f"    MANE proteins {min((i + 1) * args.sequence_batch, len(ids))}/{len(ids)}", flush=True)


async def clinvar_plp_missense_positions(symbol: str, mane_refseq: str | None) -> dict[str, Any]:
    """Protein positions of the gene's P/LP missense records, restricted to titles on the MANE RefSeq."""
    cache = load_cache("clinvar_positions.json")
    if symbol in cache:
        return cache[symbol]
    key = core.secret("NCBI_API_KEY")
    base = {"db": "clinvar", "retmode": "json", **({"api_key": key} if key else {})}
    term = f"{symbol}[gene] AND single_gene[prop] AND {PLP_PROPS} AND {MISSENSE}"
    search = await core.get_json(f"{EUTILS}/esearch.fcgi",
                                 params={**base, "term": term, "usehistory": "y", "retmax": 0}, timeout=120)
    res = search["esearchresult"]
    total = int(res["count"])
    positions: list[int] = []
    off_transcript = 0
    unparsed = 0
    for start in range(0, total, 300):
        page = await core.get_json(f"{EUTILS}/esummary.fcgi",
                                   params={**base, "query_key": res["querykey"], "WebEnv": res["webenv"],
                                           "retstart": start, "retmax": 300}, timeout=300)
        result = page.get("result", {})
        for uid in result.get("uids", []):
            m = TITLE_RE.match(result[uid].get("title", ""))
            if not m:
                unparsed += 1
                continue
            nm, prot = m.group(1), m.group(4)
            if mane_refseq and nm != mane_refseq:
                off_transcript += 1
                continue
            pm = HGVSP_RE.search(prot or "")
            if not pm:
                unparsed += 1
                continue
            positions.append(int(pm.group(2)))
    out = {"n_records": total, "positions": positions, "n_off_mane_transcript": off_transcript,
           "n_unparsed": unparsed}
    cache = load_cache("clinvar_positions.json")
    cache[symbol] = out
    dump_cache("clinvar_positions.json", cache)
    return out


async def uniprot_domains(symbol: str) -> list[dict[str, Any]]:
    cache = load_cache("uniprot_domains.json")
    if symbol in cache:
        return cache[symbol]
    data = await core.get_json(UNIPROT, params={
        "query": f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true",
        "fields": "accession,length,ft_domain,ft_region", "format": "json", "size": 1}, timeout=180)
    feats = []
    for entry in data.get("results", [])[:1]:
        for f in entry.get("features", []):
            if f.get("type") in ("Domain", "Region"):
                loc = f.get("location", {})
                feats.append({"type": f["type"], "description": f.get("description", ""),
                              "start": (loc.get("start") or {}).get("value"),
                              "end": (loc.get("end") or {}).get("value"),
                              "accession": entry.get("primaryAccession")})
    cache = load_cache("uniprot_domains.json")
    cache[symbol] = feats
    dump_cache("uniprot_domains.json", cache)
    return feats


# --------------------------------------------------------------------------- alignment + metrics


def align_pair(src_seq: str, tgt_seq: str) -> dict[int, int]:
    """Same settings as the project's paralog_map tool: global NW-Gotoh, BLOSUM62, gap -10/-0.5, free end gaps."""
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10.0, -0.5
    aligner.end_gap_score = 0.0
    aln = aligner.align(src_seq, tgt_seq)[0]
    mapping: dict[int, int] = {}
    for (s0, s1), (t0, _t1) in zip(*aln.aligned.tolist()):
        for i in range(s1 - s0):
            mapping[s0 + i + 1] = t0 + i + 1
    return mapping


def hotspot(positions: list[int], length: int) -> dict[str, Any]:
    """How tightly do the pathogenic residues sit together? Two transparent statistics, no domain model needed."""
    uniq = sorted(set(positions))
    if not uniq:
        return {"window50_start": None, "window50_end": None, "window50_frac_of_protein": None,
                "best_decile_occupancy": None}
    half = max(1, (len(uniq) + 1) // 2)
    best = min(((uniq[i + half - 1] - uniq[i], uniq[i], uniq[i + half - 1])
                for i in range(len(uniq) - half + 1)), key=lambda t: t[0])
    width = max(1, int(round(0.10 * length)))
    occ = 0
    for i, p in enumerate(uniq):
        j = i
        while j < len(uniq) and uniq[j] <= p + width:
            j += 1
        occ = max(occ, j - i)
    return {"window50_start": best[1], "window50_end": best[2],
            "window50_frac_of_protein": round((best[0] + 1) / length, 3),
            "best_decile_occupancy": round(occ / len(uniq), 3)}


def constraint_score(c: dict[str, Any]) -> float:
    """0-1 summary of 'is this gene under selection at all', from LOEUF and missense Z."""
    loeuf, misz = num(c.get("oe_lof_upper")), num(c.get("mis_z"))
    a = 1.0 if loeuf is not None and loeuf <= 0.35 else 0.7 if loeuf is not None and loeuf <= 0.6 else \
        0.35 if loeuf is not None and loeuf <= 1.0 else 0.0
    b = 1.0 if misz is not None and misz >= 3 else 0.7 if misz is not None and misz >= 2 else \
        0.35 if misz is not None and misz >= 1 else 0.0
    return round((a + b) / 2, 3)


def pair_metrics(pair: dict[str, Any], genes: dict[str, Any], positions: dict[str, Any],
                 args: argparse.Namespace) -> dict[str, Any]:
    """Pure computation on the prefetched caches: no network, so every qualified pair can be scored."""
    rich, sib = pair["rich_gene"], pair["sibling_gene"]
    row: dict[str, Any] = {"rich_gene": rich, "sibling_gene": sib,
                           "rich_plp_missense": int(pair["rich_plp_missense"]),
                           "sibling_plp_missense": int(pair["sibling_plp_missense"]),
                           "paralog_type": pair["paralog_type"], "duplication_level": pair["duplication_level"],
                           "status": "ok", "note": ""}

    g_rich, g_sib = genes.get(rich) or {}, genes.get(sib) or {}
    for tag, g in (("rich", g_rich), ("sibling", g_sib)):
        mane = (g or {}).get("mane_select_transcript") or {}
        c = (g or {}).get("gnomad_constraint") or {}
        row[f"{tag}_mane_enst"] = mane.get("ensembl_id")
        row[f"{tag}_mane_refseq"] = mane.get("refseq_id")
        row[f"{tag}_loeuf"] = num(c.get("oe_lof_upper"))
        row[f"{tag}_pli"] = num(c.get("pli"))
        row[f"{tag}_mis_z"] = num(c.get("mis_z"))
        row[f"{tag}_oe_mis_upper"] = num(c.get("oe_mis_upper"))
    if not row["rich_mane_enst"] or not row["sibling_mane_enst"]:
        row.update(status="skipped", note="no MANE Select transcript in gnomAD for one of the two genes")
        return row

    src_seq = cached_sequence(rich, row["rich_mane_enst"])
    tgt_seq = cached_sequence(sib, row["sibling_mane_enst"])
    if not src_seq or not tgt_seq:
        row.update(status="skipped", note="Ensembl did not return the MANE protein for one of the two genes")
        return row
    row["rich_len"], row["sibling_len"] = len(src_seq), len(tgt_seq)
    if max(len(src_seq), len(tgt_seq)) > args.max_protein_len:
        row.update(status="skipped", note=f"protein longer than --max-protein-len {args.max_protein_len}; "
                                          "global alignment memory blows up")
        return row

    mapping = align_pair(src_seq, tgt_seq)
    ident_pairs = sum(1 for s, t in mapping.items() if src_seq[s - 1] == tgt_seq[t - 1])
    row["aligned_pairs"] = len(mapping)
    row["identity_over_aligned"] = round(ident_pairs / len(mapping), 3) if mapping else None
    row["identity_over_shorter_protein"] = round(ident_pairs / min(len(src_seq), len(tgt_seq)), 3)

    cv = positions.get(rich)
    if not cv:
        row.update(status="skipped", note="ClinVar positions for the rich gene were not fetched")
        return row
    pos_all = [p for p in cv["positions"] if 1 <= p <= len(src_seq)]
    row["n_plp_positions_on_mane"] = len(pos_all)
    row["n_plp_off_mane_transcript"] = cv["n_off_mane_transcript"]
    uniq = sorted(set(pos_all))
    row["n_plp_unique_residues"] = len(uniq)
    if not uniq:
        row.update(status="skipped", note="no P/LP missense record could be placed on the MANE RefSeq transcript")
        return row

    aligned_uniq = [p for p in uniq if p in mapping]
    row["gap_fraction_at_plp_residues"] = round(1 - len(aligned_uniq) / len(uniq), 3)
    ident_res = [p for p in aligned_uniq if src_seq[p - 1] == tgt_seq[mapping[p] - 1]]
    row["transfer_signal_residues"] = round(len(ident_res) / len(uniq), 3)
    ident_var = sum(1 for p in pos_all if p in mapping and src_seq[p - 1] == tgt_seq[mapping[p] - 1])
    row["transfer_signal_variants"] = round(ident_var / len(pos_all), 3)
    row["n_identical_residues"] = len(ident_res)

    hs = hotspot(uniq, len(src_seq))
    row.update(hs)
    w0, w1 = hs["window50_start"], hs["window50_end"]
    # honest background: identity of the same region, over every aligned residue in it, not just the pathogenic ones
    win = [p for p in range(w0, w1 + 1) if p in mapping]
    row["background_identity_in_window"] = round(
        sum(1 for p in win if src_seq[p - 1] == tgt_seq[mapping[p] - 1]) / len(win), 3) if win else None
    if row["background_identity_in_window"] is not None:
        row["transfer_lift"] = round(row["transfer_signal_residues"] - row["background_identity_in_window"], 3)
    else:
        row["transfer_lift"] = None
    row["sibling_constraint_score"] = constraint_score((g_sib or {}).get("gnomad_constraint") or {})
    row["rich_constraint_score"] = constraint_score((g_rich or {}).get("gnomad_constraint") or {})
    row["rich_hotspot_domains"] = ""
    return row


async def stage_metrics(args: argparse.Namespace) -> dict[str, Any]:
    pairs = read_tsv(latest("paralog_scan_pairs_qualified_*.tsv"))
    if args.top > 0:
        wanted = {tuple(p.split("/")) for p in (args.include or [])}
        forced = [p for p in pairs if (p["rich_gene"], p["sibling_gene"]) in wanted]
        # one sibling can pair with several rich genes; keep its richest partner so a single big family
        # cannot flood the shortlist
        seen: set[str] = set()
        ordered: list[dict[str, Any]] = []
        for p in sorted(pairs, key=lambda r: -int(r["rich_plp_missense"])):
            if p["sibling_gene"] in seen:
                continue
            seen.add(p["sibling_gene"])
            ordered.append(p)
        chosen = forced + [p for p in ordered if p not in forced][: args.top]
    else:
        chosen = pairs

    symbols = sorted({p["rich_gene"] for p in chosen} | {p["sibling_gene"] for p in chosen})
    print(f"  scoring {len(chosen)} pairs over {len(symbols)} distinct genes", flush=True)
    genes = await gnomad_genes_batch(symbols, args)

    wanted_seqs = []
    for s in symbols:
        enst = ((genes.get(s) or {}).get("mane_select_transcript") or {}).get("ensembl_id")
        if enst:
            wanted_seqs.append((s, enst))
    print(f"  {len(wanted_seqs)}/{len(symbols)} genes have a MANE Select transcript in gnomAD", flush=True)
    await mane_proteins_batch(wanted_seqs, args)

    rich_genes = sorted({p["rich_gene"] for p in chosen})
    positions: dict[str, Any] = {}
    for i, sym in enumerate(rich_genes, 1):
        refseq = ((genes.get(sym) or {}).get("mane_select_transcript") or {}).get("refseq_id")
        try:
            positions[sym] = await clinvar_plp_missense_positions(sym, refseq)
        except Exception as exc:
            print(f"    ! ClinVar positions {sym}: {type(exc).__name__}", flush=True)
        if i % 50 == 0:
            print(f"    ClinVar positions {i}/{len(rich_genes)}", flush=True)

    rows = []
    for p in chosen:
        try:
            rows.append(pair_metrics(p, genes, positions, args))
        except Exception as exc:  # keep the scan going, record the failure
            rows.append({"rich_gene": p["rich_gene"], "sibling_gene": p["sibling_gene"], "status": "error",
                         "note": f"{type(exc).__name__}: {exc}"})
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"  scored {len(ok)}/{len(rows)} pairs "
          f"({Counter(r['note'][:40] for r in rows if r['status'] != 'ok').most_common(4)})", flush=True)

    if args.domains and ok:
        top = sorted(ok, key=feasibility, reverse=True)[: args.domain_top]
        for i, r in enumerate(top, 1):
            try:
                doms = await uniprot_domains(r["rich_gene"])
            except Exception:
                continue
            w0, w1 = r["window50_start"], r["window50_end"]
            hits = [d for d in doms if d["start"] and d["end"] and not (d["end"] < w0 or d["start"] > w1)]
            r["rich_hotspot_domains"] = "; ".join(
                f"{d['description'] or d['type']} {d['start']}-{d['end']}" for d in hits[:4])
            if i % 20 == 0:
                print(f"    UniProt domains {i}/{len(top)}", flush=True)

    saved = core.save_tsv(OUT, f"paralog_scan_metrics_{core.today()}.tsv", rows, base=core.INTERIM,
                          source="gnomAD v4.1 GraphQL + Ensembl MANE proteins + ClinVar esummary + UniProt features",
                          urls=[GNOMAD_API, f"{ENSEMBL}/sequence/id/", f"{EUTILS}/esummary.fcgi", UNIPROT],
                          tool="scripts/paralog_pairs_scan.py::metrics",
                          params={"top": args.top, "max_protein_len": args.max_protein_len,
                                  "aligner": "global NW-Gotoh BLOSUM62 gap -10/-0.5 free end gaps"})
    return {"n_scored": sum(1 for r in rows if r["status"] == "ok"), "n_rows": len(rows), "saved": saved}


# --------------------------------------------------------------------------- stage: rank

def val(r: dict[str, Any], key: str, default: float) -> float:
    """Missing or unparseable becomes `default`. Never write `x or default`: a legitimate 0.0 is falsy and
    would be silently replaced, which is exactly how a perfect-coverage pair gets punished."""
    v = r.get(key)
    return default if v is None or v == "" else float(v)


WEAK_RULES = [
    ("low_identity", lambda r: val(r, "background_identity_in_window", 0.0) < 0.30),
    ("region_absent_in_sibling", lambda r: val(r, "gap_fraction_at_plp_residues", 1.0) > 0.25),
    ("sibling_unconstrained", lambda r: val(r, "sibling_constraint_score", 0.0) < 0.35),
    ("no_conservation_lift", lambda r: val(r, "transfer_lift", -1.0) <= 0.05),
    ("thin_labels", lambda r: val(r, "n_plp_unique_residues", 0.0) < 25),
]


def feasibility(r: dict[str, Any]) -> float:
    """Explicit, arbitrary-but-stated weights. Every component is in the table so a reader can re-rank."""
    transfer = val(r, "transfer_signal_residues", 0.0)
    lift = max(0.0, val(r, "transfer_lift", 0.0))
    cover = 1.0 - val(r, "gap_fraction_at_plp_residues", 1.0)
    sib = val(r, "sibling_constraint_score", 0.0)
    labels = min(val(r, "n_plp_unique_residues", 0.0) / 100.0, 1.0)
    return round(0.30 * transfer + 0.25 * min(lift / 0.30, 1.0) + 0.15 * cover + 0.20 * sib + 0.10 * labels, 4)


def stage_rank(args: argparse.Namespace) -> dict[str, Any]:
    rows = read_tsv(latest("paralog_scan_metrics_*.tsv"))
    out = []
    for r in rows:
        if r.get("status") != "ok":
            continue
        rec = {k: (num(v) if k not in ("rich_gene", "sibling_gene", "paralog_type", "duplication_level",
                                       "rich_hotspot_domains", "rich_mane_enst", "sibling_mane_enst",
                                       "rich_mane_refseq", "sibling_mane_refseq", "status", "note")
                   else v) for k, v in r.items()}
        rec["weak_flags"] = ";".join(name for name, test in WEAK_RULES if test(rec))
        rec["n_weak_flags"] = len(rec["weak_flags"].split(";")) if rec["weak_flags"] else 0
        rec["feasibility"] = feasibility(rec)
        out.append(rec)
    out.sort(key=lambda r: -r["feasibility"])
    for i, r in enumerate(out, 1):
        r["rank"] = i
    cols = ["rank", "feasibility", "rich_gene", "sibling_gene", "rich_plp_missense", "sibling_plp_missense",
            "rich_len", "sibling_len", "identity_over_aligned", "n_plp_unique_residues",
            "transfer_signal_residues", "background_identity_in_window", "transfer_lift",
            "gap_fraction_at_plp_residues", "window50_start", "window50_end", "window50_frac_of_protein",
            "best_decile_occupancy", "rich_loeuf", "rich_mis_z", "sibling_loeuf", "sibling_pli",
            "sibling_mis_z", "sibling_constraint_score", "weak_flags", "n_weak_flags",
            "duplication_level", "rich_hotspot_domains"]
    ordered = [{c: r.get(c) for c in cols} | {k: v for k, v in r.items() if k not in cols} for r in out]
    saved = core.save_tsv(OUT, f"paralog_scan_ranked_{core.today()}.tsv", ordered, base=core.INTERIM,
                          source="derived: paralog_scan_metrics ranked by the feasibility score in this script",
                          urls=[], tool="scripts/paralog_pairs_scan.py::rank",
                          params={"weights": "0.30 transfer + 0.25 min(lift/0.30,1) + 0.15 coverage "
                                             "+ 0.20 sibling_constraint + 0.10 min(residues/100,1)"})
    print(f"\n  {'#':>2} {'pair':<20} {'feas':>5} {'n+':>4} {'ident':>6} {'transf':>7} {'bg':>5} {'lift':>6} "
          f"{'gap':>5} {'sibLOEUF':>8}  flags")
    for r in ordered[:20]:
        print(f"  {r['rank']:>2} {r['rich_gene'] + '/' + r['sibling_gene']:<20} {r['feasibility']:>5} "
              f"{int(r['rich_plp_missense']):>4} {r['identity_over_aligned']:>6} "
              f"{r['transfer_signal_residues']:>7} {r['background_identity_in_window']:>5} "
              f"{r['transfer_lift']:>6} {r['gap_fraction_at_plp_residues']:>5} "
              f"{r['sibling_loeuf'] if r['sibling_loeuf'] is not None else float('nan'):>8.2f}  {r['weak_flags']}")
    return {"n_ranked": len(ordered), "saved": saved,
            "n_clean": sum(1 for r in ordered if r["n_weak_flags"] == 0)}


# --------------------------------------------------------------------------- main


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", default="all", choices=["all", "sources", "rich", "pairs", "metrics", "rank"])
    p.add_argument("--min-rich", type=int, default=30, help="minimum ClinVar P/LP missense in the rich gene")
    p.add_argument("--max-sibling", type=int, default=2, help="maximum ClinVar P/LP missense in the sibling")
    p.add_argument("--sibling-prefilter", type=int, default=8,
                   help="skip siblings with more than this many P/LP alleles in the ClinVar gene summary")
    p.add_argument("--top", type=int, default=0,
                   help="pairs to score in the metrics stage; 0 (default) scores every qualified pair, "
                        "which is affordable because the per-pair work is one local alignment")
    p.add_argument("--include", nargs="*", default=["KIF1A/KIF13A"],
                   help="pairs to score regardless of rank when --top > 0, as RICH/SIBLING")
    p.add_argument("--max-protein-len", type=int, default=6000)
    p.add_argument("--domains", action="store_true", default=True)
    p.add_argument("--no-domains", dest="domains", action="store_false")
    p.add_argument("--domain-top", type=int, default=500,
                   help="fetch UniProt domain names only for this many best-scoring pairs")
    p.add_argument("--gnomad-batch", type=int, default=25, help="genes per gnomAD GraphQL request")
    p.add_argument("--sequence-batch", type=int, default=50, help="transcripts per Ensembl POST /sequence/id")
    p.add_argument("--ncbi-concurrency", type=int, default=6, help="NCBI allows 10 requests/s with an API key")
    p.add_argument("--ensembl-concurrency", type=int, default=12, help="Ensembl allows 15 requests/s")
    p.add_argument("--ensembl-timeout", type=float, default=60.0,
                   help="first-pass homology timeout; slow gene families are retried once with 4x this")
    p.add_argument("--refresh", action="store_true", help="ignore caches and re-download")
    return p


async def run(args: argparse.Namespace) -> int:
    started = time.monotonic()
    print(f"paralog_pairs_scan  stage={args.stage}  min_rich={args.min_rich}  max_sibling={args.max_sibling}")
    print(f"  project root {core.ROOT}")
    print(f"  NCBI_API_KEY loaded: {bool(core.secret('NCBI_API_KEY'))} (10 req/s instead of 3)")
    print(f"  cache {core.rel(CACHE)}  output data/interim/{OUT}\n")

    order = ["sources", "rich", "pairs", "metrics", "rank"]
    stages = order if args.stage == "all" else [args.stage]
    for name in stages:
        print(f"[{name}]", flush=True)
        if name == "sources":
            res = await stage_sources(args)
        elif name == "rich":
            res = await stage_rich(args)
        elif name == "pairs":
            res = await stage_pairs(args)
        elif name == "metrics":
            res = await stage_metrics(args)
        else:
            res = stage_rank(args)
        print(json.dumps(res, indent=1, default=str), "\n", flush=True)

    print(f"done in {time.monotonic() - started:.0f}s")
    print("provenance: data/raw/PROVENANCE.jsonl plus a .prov.json next to every table written")
    return 0


def main() -> int:
    return asyncio.run(run(build_parser().parse_args()))


if __name__ == "__main__":
    sys.exit(main())
