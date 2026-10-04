"""C ekseni (nakavt): kaynak dosyalari ve gnomAD/Ensembl cekimleri.

TRUBA debug bolumunde (internet olan tek bolum) SLURM isi olarak calisir.
Her indirme icin raw/MANIFEST.jsonl satiri yazar: url, utc zaman, bayt, sha256.
Yalniz standart kutuphane kullanir (urllib), boylece conda ortamina paket kurmak gerekmez.

Kullanim:
  python fetch_sources.py files           # ek tablolar, G2P, gnomAD kisitlama, popEVE
  python fetch_sources.py gnomad GEN:ENST [GEN:ENST ...]   # gnomAD v4.1 populasyon bazli LoF
  python fetch_sources.py gnomad_all GEN:ENST [...]        # ayni, sinonim dahil butun varyantlar
  python fetch_sources.py ensembl ENST [ENST ...]          # MANE ekzon yapisi
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = Path("/arf/scratch/avural/kif13a_nakavt")
RAW = BASE / "raw"
MANIFEST = RAW / "MANIFEST.jsonl"
UA = "kif13a-neurovariant-atlas/0.1 (academic research; TRUBA)"

SPR = "https://static-content.springer.com/esm/art%3A10.1038%2F{doi}/MediaObjects/{fn}"
FILES = [
    # Genes & Health, Kim 2026 Nat Genet, PMID 41896352, PMC13083238
    ("gh_supp_tables", SPR.format(doi="s41588-026-02553-7", fn="41588_2026_2553_MOESM4_ESM.xlsx")),
    ("gh_supp_info", SPR.format(doi="s41588-026-02553-7", fn="41588_2026_2553_MOESM1_ESM.pdf")),
    # RGC-ME, Sun 2024 Nature, PMID 38768635, PMC11254753
    ("rgc_supp_info", SPR.format(doi="s41586-024-07556-0", fn="41586_2024_7556_MOESM1_ESM.docx")),
    ("rgc_st1", SPR.format(doi="s41586-024-07556-0", fn="41586_2024_7556_MOESM3_ESM.xlsx")),
    ("rgc_st2", SPR.format(doi="s41586-024-07556-0", fn="41586_2024_7556_MOESM4_ESM.xlsx")),
    ("rgc_st6", SPR.format(doi="s41586-024-07556-0", fn="41586_2024_7556_MOESM8_ESM.xlsx")),
    # Chundru 2024 Nat Genet, PMID 39313616, PMC11525179
    ("chundru_supp_tables", SPR.format(doi="s41588-024-01910-8", fn="41588_2024_1910_MOESM4_ESM.xlsx")),
    ("chundru_supp_data", SPR.format(doi="s41588-024-01910-8", fn="41588_2024_1910_MOESM5_ESM.zip")),
    ("chundru_supp_info", SPR.format(doi="s41588-024-01910-8", fn="41588_2024_1910_MOESM1_ESM.pdf")),
    # popEVE, Orenbuch 2025 Nat Genet, PMID 41286104, PMC12695638
    ("popeve_supp_tables", SPR.format(doi="s41588-025-02400-1", fn="41588_2025_2400_MOESM3_ESM.xlsx")),
    ("popeve_KIF13A_NP_071396.4", "https://data.evemodel.org/popeve/v1.1/mutations/NP_071396.4.csv"),
    ("popeve_KIF1A_NP_004312.2", "https://data.evemodel.org/popeve/v1.1/mutations/NP_004312.2.csv"),
    # DDG2P (Gene2Phenotype, aylik surum)
    ("ddg2p", "https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2026_09_28/DDG2P_2026-09-28.csv.gz"),
    # gnomAD v4.1 gen kisitlamasi
    ("gnomad_v41_constraint", "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"),
]

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL = "https://rest.ensembl.org"


def now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def http(url: str, data: bytes | None = None, headers: dict | None = None, tries: int = 6,
         timeout: int = 600) -> bytes:
    hdr = {"User-Agent": UA}
    hdr.update(headers or {})
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=hdr)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:  # noqa: PERF203
            last = e
            code = getattr(e, "code", None)
            if code in (400, 401, 403, 404):
                raise
            wait = 5 * (i + 1)
            print(f"    retry {i + 1}/{tries} after {wait}s: {e}", flush=True)
            time.sleep(wait)
    raise RuntimeError(f"failed {url}: {last}")


def record(path: Path, url: str, body: bytes, extra: dict | None = None) -> dict:
    path.write_bytes(body)
    rec = {"file": str(path.relative_to(BASE)), "url": url, "retrieved_utc": now(),
           "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()}
    rec.update(extra or {})
    with MANIFEST.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"  {rec['file']}  {rec['bytes']:,} B  sha256 {rec['sha256'][:16]}", flush=True)
    return rec


def fetch_files() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for key, url in FILES:
        fn = url.rsplit("/", 1)[-1]
        out = RAW / f"{key}__{fn}"
        if out.exists() and out.stat().st_size > 0:
            print(f"  have {out.name}")
            continue
        print(f"GET {url}", flush=True)
        try:
            body = http(url)
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED {key}: {e}", flush=True)
            with MANIFEST.open("a", encoding="utf-8") as fh:
                fh.write(json.dumps({"file": None, "key": key, "url": url, "retrieved_utc": now(),
                                     "error": str(e)}) + "\n")
            continue
        record(out, url, body, {"key": key})


def gql(query: str) -> dict:
    body = json.dumps({"query": query}).encode()
    raw = http(GNOMAD_API, data=body, headers={"Content-Type": "application/json"})
    return json.loads(raw)


POP = "populations { id ac an homozygote_count }"
Q_FULL = """
{ transcript(transcript_id: "%s", reference_genome: GRCh38) {
    transcript_id transcript_version gene_id gene { symbol }
    variants(dataset: gnomad_r4) {
      variant_id chrom pos ref alt rsids hgvsc hgvsp consequence lof lof_filter lof_flags flags
      exome { ac an homozygote_count filters %s }
      genome { ac an homozygote_count filters %s }
      joint { ac an homozygote_count filters %s }
    } } }"""

LOF_CSQ = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
           "start_lost", "stop_lost", "transcript_ablation"}


def fetch_gnomad(specs: list[str]) -> None:
    out_dir = RAW / "gnomad"
    out_dir.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        sym, enst = spec.split(":")
        out = out_dir / f"{sym}_{enst}_gnomad_r4_pop.json"
        if out.exists() and out.stat().st_size > 0:
            print(f"  have {out.name}")
            continue
        print(f"gnomAD {sym} {enst}", flush=True)
        q = Q_FULL % (enst, POP, POP, POP)
        d = gql(q)
        if d.get("errors"):
            print(f"  GraphQL errors (full query): {json.dumps(d['errors'])[:800]}", flush=True)
            raise SystemExit(2)
        tx = d["data"]["transcript"]
        keep = [v for v in tx["variants"] if v.get("lof") or v.get("consequence") in LOF_CSQ]
        payload = {"transcript_id": tx["transcript_id"], "transcript_version": tx["transcript_version"],
                   "gene_id": tx["gene_id"], "symbol": (tx.get("gene") or {}).get("symbol"),
                   "dataset": "gnomad_r4", "n_variants_total": len(tx["variants"]),
                   "n_lof_rows": len(keep), "variants": keep}
        record(out, GNOMAD_API, json.dumps(payload).encode(),
               {"query": "transcript variants gnomad_r4 with exome/genome/joint populations",
                "gene": sym, "transcript": enst})
        time.sleep(4)  # gnomAD API hiz siniri


Q_ALL = """
{ transcript(transcript_id: "%s", reference_genome: GRCh38) {
    transcript_id transcript_version
    variants(dataset: gnomad_r4) {
      variant_id ref alt hgvsc hgvsp consequence lof lof_filter
      joint { ac an homozygote_count filters populations { id ac an homozygote_count } }
    } } }"""


def fetch_gnomad_all(specs: list[str]) -> None:
    """Butun varyantlar (sinonim dahil), homozigot beklentisinin notr kalibrasyonu icin."""
    out_dir = RAW / "gnomad"
    out_dir.mkdir(parents=True, exist_ok=True)
    for spec in specs:
        sym, enst = spec.split(":")
        out = out_dir / f"{sym}_{enst}_gnomad_r4_allvariants_pop.json"
        if out.exists() and out.stat().st_size > 0:
            print(f"  have {out.name}")
            continue
        print(f"gnomAD all variants {sym} {enst}", flush=True)
        d = gql(Q_ALL % enst)
        if d.get("errors"):
            print(f"  GraphQL errors: {json.dumps(d['errors'])[:800]}", flush=True)
            raise SystemExit(2)
        tx = d["data"]["transcript"]
        payload = {"transcript_id": tx["transcript_id"], "transcript_version": tx["transcript_version"],
                   "symbol": sym, "dataset": "gnomad_r4", "n_variants_total": len(tx["variants"]),
                   "variants": tx["variants"]}
        record(out, GNOMAD_API, json.dumps(payload).encode(),
               {"query": "transcript variants gnomad_r4, all consequences, joint populations",
                "gene": sym, "transcript": enst})
        time.sleep(4)


def fetch_ensembl(ensts: list[str]) -> None:
    out_dir = RAW / "ensembl"
    out_dir.mkdir(parents=True, exist_ok=True)
    for enst in ensts:
        out = out_dir / f"{enst}_structure.json"
        if out.exists() and out.stat().st_size > 0:
            print(f"  have {out.name}")
            continue
        url = f"{ENSEMBL}/lookup/id/{enst}?expand=1;utr=1"
        body = http(url, headers={"Content-Type": "application/json"}, tries=9)
        record(out, url, body, {"transcript": enst})
        time.sleep(0.5)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "files":
        fetch_files()
    elif mode == "gnomad":
        fetch_gnomad(sys.argv[2:])
    elif mode == "gnomad_all":
        fetch_gnomad_all(sys.argv[2:])
    elif mode == "ensembl":
        fetch_ensembl(sys.argv[2:])
    else:
        raise SystemExit(f"unknown mode {mode}")
