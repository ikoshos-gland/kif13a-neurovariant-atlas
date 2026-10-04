"""Brain imaging genetics axis: does the KIF13A locus carry a white matter signal in UK Biobank?

Why this axis exists
--------------------
Two independent cell-type analyses agree that KIF13A is an oligodendrocyte-lineage gene, not a
neuronal one (HANDOFF.md 3.6). Oligodendrocytes build myelin, so if KIF13A variation has any
human brain consequence the phenotype should be *white matter*, not the neuronal axonal-transport
disease KIF1A causes. That prediction had never been tested against human data. This script tests
it against the Oxford BIG40 resource: GWAS of 3,935 brain imaging-derived phenotypes (IDPs) in up
to 33,224 UK Biobank subjects, summary statistics public and unauthenticated.

Discipline the script enforces
------------------------------
1. Phenotypes are classified into white matter microstructure / white matter macrostructure /
   grey matter and everything else BEFORE any p value is read (stage ``classify`` writes the
   table, stage ``scan`` reads it back and never re-derives it).
2. Both thresholds are applied and reported: the conventional 5e-8 and BIG40's own
   multiple-testing threshold across IDPs, whose published value is fetched live (stage
   ``thresholds``), never assumed.
3. Controls are mandatory. KIF1A gets the identical extraction; three established white matter
   loci are extracted to prove the pipeline can find a signal; negative-control loci and a set of
   random matched windows give the empirical null for "best p across 675 white matter IDPs in a
   500 kb window".
4. A null is only reported together with its detectable bound (stage ``report``).

Build discipline
----------------
BIG40 is **GRCh37/hg19** (stated on the BIG40 index page and on the PheWeb ``/about`` page, both
saved with provenance). Every locus coordinate is therefore fetched from the **GRCh37** Ensembl
REST endpoint, never converted by hand. The project has been bitten by a build mixup before
(gnomAD GraphQL returns v2.1.1 constraint for GRCh37 and v4 for GRCh38, HANDOFF.md section 7a).

Data routes, cheapest first
---------------------------
* ``pheweb33k`` region API: one HTTP call per (IDP, locus), returns every variant in the window
  with MAF >= 0.01 and **p <= 0.01** (measured, not assumed: stage ``probe`` prints the largest p
  and smallest MAF each IDP returns). Because everything below 0.01 is present, the minimum p per
  IDP inside a window is exactly recoverable, and an IDP that returns nothing has min p > 0.01.
  About 10 calls/s, so a full 3,935-IDP sweep costs ~7 min and ~1 MB. This is the primary route.
* ``pheweb33k`` gene page: one call per gene, returns the precomputed best phenotype per gene
  region (gene +- 100 kb) for phenotypes with p < 1e-6 (at least the top 4). Used as an
  independent cross-check of the sweep.
* ``release2/Table1.txt``: BIG40's own table of local peak associations with -log10(p) > 7.5
  across all IDPs, 22k discovery sample with 11k replication. One 0.9 MB download covers every
  IDP at once.
* ``release2/stats33k/NNNN.txt.gz``: the full per-IDP summary statistics, 344 MB each, single
  member gzip (checked: no BGZF sync points, so no random access) served at 0.1-0.45 MB/s. Used
  for exactly three IDPs: the two where MAF 0.001-0.01 coverage is load-bearing (the cingulum OD
  phenotype the AstraZeneca PheWAS lead named, and the locus's only 5e-8 hit) and the one whose
  region query answers HTTP 500. Streams the prefix and stops once the window is passed (~100-120
  MB each). Nothing is written to disk except the extracted window.

Usage (system python 3.13 from the repo root; the project venv is blocked by Windows, see
CLAUDE.md)::

    python -u scripts\\imaging_gwas_axis.py --stage sources
    python -u scripts\\imaging_gwas_axis.py --stage classify
    python -u scripts\\imaging_gwas_axis.py --stage sumstats --phenos 2012,1415
    python -u scripts\\imaging_gwas_axis.py --stage scan
    python -u scripts\\imaging_gwas_axis.py --stage thresholds
    python -u scripts\\imaging_gwas_axis.py --stage followup
    python -u scripts\\imaging_gwas_axis.py --stage report
    python -u scripts\\imaging_gwas_axis.py --stage figure

Every download goes through ``atlas_mcp.core.save_*``, which writes a ``.prov.json`` sidecar
(source URL, UTC timestamp, sha256, bytes) and appends to ``data/raw/PROVENANCE.jsonl``.
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import gzip
import io
import json
import math
import random
import re
import sys
import time
import urllib.parse
import zlib
from pathlib import Path

import httpx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

BIG40 = "https://open.win.ox.ac.uk/ukbiobank/big40"
PHEWEB = f"{BIG40}/pheweb33k"
ENSEMBL37 = "https://grch37.rest.ensembl.org"
ENSEMBL38 = "https://rest.ensembl.org"

RAW = ROOT / "data" / "raw" / "imaging_gwas"
INTERIM = ROOT / "data" / "interim" / "imaging_gwas"
DOCS = ROOT / "docs"
FIGS = ROOT / "figures"
UA = {"User-Agent": core.USER_AGENT}

# ---------------------------------------------------------------------------------------------
# Phenotype classification. Written down before any association is read; stage `scan` only reads
# the file this produces. The three classes are exactly the ones the analysis plan names.
# ---------------------------------------------------------------------------------------------

# Diffusion MRI: FA, MD/L1/L2/L3 ("diffusivity"), MO, OD, ICVF, ISOVF on tracts and on the FA
# skeleton. These are the phenotypes the oligodendrocyte hypothesis predicts.
WM_MICRO_CATEGORIES = {
    "WM tract FA",
    "WM tract diffusivity",
    "WM tract MO",
    "WM tract OD",
    "WM tract ICVF",
    "WM tract ISOVF",
}

# White matter volume and lesion load. Named explicitly because they sit inside broader BIG40
# categories ("regional and tissue volume") that are otherwise negative controls.
WM_MACRO_PHENOS = {
    "IDP_T2_FLAIR_BIANCA_WMH_volume",          # total white matter hyperintensity volume
    "aseg_global_volume_WM-hypointensities",   # T1 white matter hypointensities (lesion-like)
    "IDP_T1_SIENAX_white_normalised_volume",
    "IDP_T1_SIENAX_white_unnormalised_volume",
    "aseg_lh_volume_CerebralWhiteMatter",
    "aseg_rh_volume_CerebralWhiteMatter",
    "aseg_lh_volume_Cerebellum-White-Matter",
    "aseg_rh_volume_Cerebellum-White-Matter",
}

# Everything else is the negative-control class: regional grey matter volumes, cortical thickness
# and area, subcortical volumes, functional connectivity, task activation, T2*, and cortical
# grey-white intensity contrast. The last one is flagged: the literature reads grey-white contrast
# as partly reflecting intracortical myelination, so a critic can argue it belongs with myelin
# phenotypes. It is NOT diffusion and NOT a white matter volume, so under the plan's own
# definitions it stays in the control class, and the report states what changes if it moves.
GREY_WHITE_CONTRAST = "cortical grey-white contrast"
QC_CATEGORY = "QC"

CLASS_WM_MICRO = "ak_madde_mikroyapi"
CLASS_WM_MACRO = "ak_madde_makroyapi"
CLASS_OTHER = "gri_madde_ve_diger"
CLASS_QC = "QC"


def classify(phenocode: str, category: str) -> tuple[str, str]:
    """Return (class, rule) for one IDP. Deterministic, uses only name and category."""
    if category in WM_MICRO_CATEGORIES:
        return CLASS_WM_MICRO, f"kategori '{category}' difuzyon MRG"
    if phenocode in WM_MACRO_PHENOS:
        return CLASS_WM_MACRO, "ak madde hacmi / hiperintensite yuku (ad listesi)"
    if category == QC_CATEGORY:
        return CLASS_QC, "goruntuleme kalite olcutu, fenotip degil"
    if category == GREY_WHITE_CONTRAST:
        return CLASS_OTHER, "gri-ak kontrasti: difuzyon degil, hacim degil -> kontrol sinifi (TARTISMALI)"
    return CLASS_OTHER, f"kategori '{category}' -> negatif kontrol sinifi"


# ---------------------------------------------------------------------------------------------
# Loci. Coordinates are fetched live from Ensembl GRCh37 (BIG40 is GRCh37).
# ---------------------------------------------------------------------------------------------

LOCI: dict[str, dict] = {
    # target and its two protein-coding neighbours (a locus signal is not a gene signal)
    "KIF13A": {"gene": "KIF13A", "role": "hedef", "sweep": "all"},
    "NUP153": {"gene": "NUP153", "role": "komsu (daha kisitli: LOEUF 0.277)", "sweep": "none"},
    "NHLRC1": {"gene": "NHLRC1", "role": "komsu", "sweep": "none"},
    # the neuronal kinesin contrast
    "KIF1A": {"gene": "KIF1A", "role": "noronal kinezin karsilastirmasi", "sweep": "all"},
    # positive controls: loci where a white matter signal is established in the literature.
    # Chosen live, see docs/ak_madde_goruntuleme.md for the justification and the citations.
    "VCAN": {"gene": "VCAN", "role": "pozitif kontrol: dMRI mikroyapi", "sweep": "wm"},
    "TRIM47": {"gene": "TRIM47", "role": "pozitif kontrol: ak madde hiperintensitesi", "sweep": "wm"},
    "EFEMP1": {"gene": "EFEMP1", "role": "pozitif kontrol: ak madde hiperintensitesi", "sweep": "wm"},
    # negative control loci: no expected brain phenotype
    "LCT": {"gene": "LCT", "role": "negatif kontrol (laktaz, beyin fenotipi beklenmez)", "sweep": "wm"},
    "MYH7": {"gene": "MYH7", "role": "negatif kontrol (kalp miyozini)", "sweep": "wm"},
}

PADS = (100_000, 250_000)
PVAL_PRESENT_MAX = 0.01  # measured in stage `probe`: PheWeb holds only p <= 0.01, MAF >= 0.01

# random matched windows: empirical null for "best p across the 675 white matter IDPs in a window"
N_RANDOM_WINDOWS = 12
RANDOM_SEED = 20260917
# GRCh37 chromosome lengths, used only to place random windows inside gene-rich autosomes
GRCH37_LEN = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566,
}


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------

def _get(url: str, *, timeout: float = 180.0, headers: dict | None = None) -> httpx.Response:
    last: Exception | None = None
    for attempt in range(4):
        try:
            r = httpx.get(url, timeout=timeout, follow_redirects=True, headers={**UA, **(headers or {})})
            if r.status_code in (429, 500, 502, 503, 504):
                raise RuntimeError(f"HTTP {r.status_code}")
            r.raise_for_status()
            return r
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"{url}: {last}")


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    cols: list[str] = []
    for row in rows:
        for k in row:
            if k not in cols:
                cols.append(k)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"  -> {core.rel(path)}  ({len(rows)} satir)")


def gene_grch37(symbol: str) -> dict:
    """Live GRCh37 coordinates for a gene symbol. Cached under data/raw/imaging_gwas."""
    cache = RAW / f"ensembl37_{symbol}.json"
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    url = f"{ENSEMBL37}/lookup/symbol/homo_sapiens/{symbol}?expand=0"
    r = _get(url, headers={"Accept": "application/json"})
    core.save_bytes("imaging_gwas", cache.name, r.content, source=url,
                    note=f"{symbol} GRCh37 koordinatlari (BIG40 GRCh37)")
    return r.json()


# ---------------------------------------------------------------------------------------------
# stage: sources
# ---------------------------------------------------------------------------------------------

def stage_sources() -> None:
    """Pull and provenance-stamp every BIG40 resource this axis relies on."""
    print("== stage sources ==")
    for url, name, note in [
        (f"{BIG40}/", "BIG40_index.html", "BIG40 ana sayfa: surum gecmisi, GRCh37 beyani, dosya yollari"),
        (f"{PHEWEB}/about", "pheweb33k_about.html", "PheWeb 1.1.21, 'All positions are on GRCh37'"),
        (f"{BIG40}/release2/Table1.txt", "BIG40_Table1_local_peaks.txt",
         "BIG40 yerel tepe iliskileri tablosu, -log10(p) > 7.5, 22k kesif + 11k tekrar"),
        (f"{PHEWEB}/api/phenotypes.json", "pheweb33k_phenotypes.json",
         "3935 IDP: phenocode, kategori, genom capindaki en iyi p"),
    ]:
        r = _get(url)
        print(core.save_bytes("imaging_gwas", name, r.content, source=url, note=note))
    for spec in LOCI.values():
        g = gene_grch37(spec["gene"])
        print(f"  {spec['gene']}: chr{g['seq_region_name']}:{g['start']}-{g['end']} "
              f"({g['assembly_name']}, strand {g['strand']})")


# ---------------------------------------------------------------------------------------------
# stage: probe  (verify what the PheWeb copy actually contains, before trusting it)
# ---------------------------------------------------------------------------------------------

def stage_probe() -> None:
    print("== stage probe: PheWeb kopyasinin gercekten ne icerdigi ==")
    g = gene_grch37("KIF13A")
    lo, hi = g["start"] - PADS[1], g["end"] + PADS[1]
    rows = []
    for ph in ("IDP_dMRI_TBSS_L3_Anterior_corona_radiata_L", "IDP_T1_SIENAX_white_normalised_volume",
               "rfMRI_connectivity_ICA100_edge_198", "IDP_dMRI_TBSS_OD_Cingulum_cingulate_gyrus_L"):
        d = region_query_sync(ph, g["seq_region_name"], lo, hi)
        if not d:
            rows.append({"pheno": ph, "n": 0})
            continue
        p = [float(x) for x in d["pvalue"]]
        rows.append({"pheno": ph, "n": len(p), "min_p": min(p), "max_p": max(p),
                     "min_maf": min(float(x) for x in d["maf"])})
    for r in rows:
        print("  ", r)
    print("  Beklenti: her yanitta max_p <= 0.01 ve min_maf >= 0.01 -> PheWeb kopyasi esiklenmis.")


def region_query_sync(phenocode: str, chrom: str, lo: int, hi: int) -> dict | None:
    filt = f"analysis in 3 and chromosome in  '{chrom}' and position ge {lo} and position le {hi}"
    url = f"{PHEWEB}/api/region/{urllib.parse.quote(phenocode)}/lz-results/?filter={urllib.parse.quote(filt)}"
    d = _get(url).json()["data"]
    return d if d.get("position") else None


# ---------------------------------------------------------------------------------------------
# stage: classify
# ---------------------------------------------------------------------------------------------

def stage_classify() -> None:
    """Build the IDP classification table. Runs before any association is read."""
    print("== stage classify ==")
    import pandas as pd

    html = (RAW / "BIG40_IDPs_v4.html").read_text(encoding="utf-8")
    tab = pd.read_html(io.StringIO(html))[0]
    tab.columns = list(tab.iloc[0])
    tab = tab.iloc[1:].reset_index(drop=True)
    phenos = json.loads((RAW / "pheweb33k_phenotypes.json").read_text(encoding="utf-8"))
    # The BIG40 IDP table and PheWeb spell the same IDP differently: PheWeb sanitises spaces,
    # "+", brackets and parentheses. Normalising both sides to [A-Za-z0-9_] matches 3933/3935; the
    # two thalamic "MV(Re)" nuclei lose their brackets entirely in PheWeb, so they are hard-mapped.
    def norm(s: str) -> str:
        return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")

    by_norm = {norm(p["phenocode"]): p["phenocode"] for p in phenos}
    hard = {"ThalamNuclei_lh_volume_MV(Re)": "ThalamNuclei_lh_volume_MVRe",
            "ThalamNuclei_rh_volume_MV(Re)": "ThalamNuclei_rh_volume_MVRe"}

    rows = []
    for _, t in tab.iterrows():
        idp_name = str(t["IDP short name"]).strip()
        cat = str(t["Category name"]).strip()
        code = hard.get(idp_name) or by_norm.get(norm(idp_name))
        pw = code is not None
        code = code or idp_name
        cls, rule = classify(code, cat)
        rows.append({
            "phenocode": code,
            "idp_short_name": idp_name,
            "pheno_id": str(t["Pheno"]).strip(),
            "ukb_field": str(t["UKB ID"]).strip(),
            "big40_category": cat,
            "sinif": cls,
            "sinif_kurali": rule,
            "n_all": str(t["N(all)"]).strip(),
            "n_disc": str(t["N(disc)"]).strip(),
            "h2": str(t["Heritability"]).strip(),
            "h2_se": str(t["Heritability(SE)"]).strip(),
            "in_pheweb33k": "yes" if pw else "no",
            "big40_sumstats_url": f"{BIG40}/release2/stats33k/{str(t['Pheno']).strip()}.txt.gz",
            "aciklama": str(t["IDP description"]).strip(),
        })
    missing = [r["phenocode"] for r in rows if r["in_pheweb33k"] == "no"]
    write_tsv(INTERIM / "idp_classification.tsv", rows)
    from collections import Counter
    c = Counter(r["sinif"] for r in rows)
    for k, v in c.most_common():
        print(f"  {k}: {v}")
    print(f"  PheWeb33k'de bulunmayan IDP: {len(missing)}")
    if missing[:5]:
        print("   ", missing[:5])


# ---------------------------------------------------------------------------------------------
# stage: scan  (the extraction)
# ---------------------------------------------------------------------------------------------

async def _region_batch(client: httpx.AsyncClient, sem: asyncio.Semaphore, phenocode: str,
                        chrom: str, lo: int, hi: int) -> tuple[str, dict | None]:
    filt = f"analysis in 3 and chromosome in  '{chrom}' and position ge {lo} and position le {hi}"
    url = (f"{PHEWEB}/api/region/{urllib.parse.quote(phenocode)}/lz-results/"
           f"?filter={urllib.parse.quote(filt)}")
    async with sem:
        for attempt in range(6):
            try:
                r = await client.get(url)
                if r.status_code in (429, 500, 502, 503, 504):
                    raise RuntimeError(f"HTTP {r.status_code}")
                r.raise_for_status()
                d = r.json()["data"]
                return phenocode, (d if d.get("position") else None)
            except Exception:  # noqa: BLE001
                await asyncio.sleep(2.0 * (attempt + 1))
        # Never silently drop a phenotype: the caller collects failures and retries them serially,
        # and stage `report` refuses to summarise a locus whose sweep is incomplete.
        return phenocode, {"__failed__": True}


async def sweep(phenocodes: list[str], chrom: str, lo: int, hi: int, *, conc: int = 6,
                label: str = "") -> list[dict]:
    """One region query per phenotype. Returns every (pheno, variant) row with p < 0.05."""
    out: list[dict] = []
    t0 = time.time()
    sem = asyncio.Semaphore(conc)

    def absorb(ph: str, d: dict | None, failed: list[str]) -> None:
        if not d:
            return
        if d.get("__failed__"):
            failed.append(ph)
            return
        for i in range(len(d["position"])):
            out.append({
                "phenocode": ph, "chr": d["chr"][i], "pos": d["position"][i],
                "rsid": d["rsid"][i], "ref": d["ref"][i], "alt": d["alt"][i],
                "af": d["af"][i], "maf": d["maf"][i], "beta": d["beta"][i],
                "se": d["sebeta"][i], "p": d["pvalue"][i],
                "nearest_genes": d["nearest_genes"][i],
            })

    failed: list[str] = []
    async with httpx.AsyncClient(timeout=180, headers=UA, follow_redirects=True,
                                 limits=httpx.Limits(max_connections=conc + 2)) as cl:
        done = 0
        for chunk_start in range(0, len(phenocodes), 240):
            chunk = phenocodes[chunk_start:chunk_start + 240]
            res = await asyncio.gather(*[_region_batch(cl, sem, p, chrom, lo, hi) for p in chunk])
            for ph, d in res:
                absorb(ph, d, failed)
            done += len(chunk)
            print(f"    {label} {done}/{len(phenocodes)} IDP, {len(out)} satir, "
                  f"{len(failed)} basarisiz, {time.time()-t0:.0f}s", flush=True)
        if failed:
            print(f"    {label}: {len(failed)} IDP yeniden deneniyor (seri)", flush=True)
            retry, failed = failed, []
            sem2 = asyncio.Semaphore(1)
            for ph in retry:
                p2, d2 = await _region_batch(cl, sem2, ph, chrom, lo, hi)
                absorb(p2, d2, failed)
                await asyncio.sleep(0.2)
    if failed:
        # BIG40's PheWeb answers HTTP 500 for a small number of IDPs (one, as of 2026-09-17:
        # IDP_dMRI_TBSS_FA_Posterior_corona_radiata_L, pheno 1479, a white matter FA phenotype).
        # A silent hole in a white matter sweep would be exactly the kind of broken extraction
        # that makes a null unusable, so the names are written next to the sweep and stage
        # `report` fills them from the full summary statistics instead.
        (INTERIM / f"scan_{label}_failed.txt").write_text("\n".join(failed) + "\n", encoding="utf-8")
        print(f"    !! {label}: {len(failed)} IDP PheWeb'de HTTP 500 veriyor -> "
              f"scan_{label}_failed.txt, sumstats ile kapatilacak: {failed[:5]}", flush=True)
    return out


def stage_scan(only: str | None = None) -> None:
    print("== stage scan ==")
    cls_rows = read_tsv(INTERIM / "idp_classification.tsv")
    all_codes = [r["phenocode"] for r in cls_rows if r["in_pheweb33k"] == "yes"]
    wm_codes = [r["phenocode"] for r in cls_rows
                if r["in_pheweb33k"] == "yes" and r["sinif"] in (CLASS_WM_MICRO, CLASS_WM_MACRO)]
    print(f"  taranacak IDP: tum {len(all_codes)}, ak madde {len(wm_codes)}")

    targets: list[tuple[str, str, int, int, list[str], str]] = []
    for name, spec in LOCI.items():
        if spec["sweep"] == "none":
            continue
        if only and name != only:
            continue
        g = gene_grch37(spec["gene"])
        codes = all_codes if spec["sweep"] == "all" else wm_codes
        targets.append((name, g["seq_region_name"], int(g["start"]), int(g["end"]), codes, spec["role"]))

    # random matched windows, same size as the KIF13A +-250 kb window, for the empirical null
    if not only:
        rng = random.Random(RANDOM_SEED)
        kg = gene_grch37("KIF13A")
        span = (int(kg["end"]) + PADS[1]) - (int(kg["start"]) - PADS[1])
        chroms = list(GRCH37_LEN)
        for i in range(N_RANDOM_WINDOWS):
            c = rng.choice(chroms)
            start = rng.randrange(5_000_000, GRCH37_LEN[c] - 5_000_000 - span)
            targets.append((f"RANDOM{i+1:02d}", c, start + PADS[1], start + span - PADS[1],
                            wm_codes, "rastgele eslesmis pencere (bos hipotez)"))

    for name, chrom, gstart, gend, codes, role in targets:
        path = INTERIM / f"scan_{name}.tsv.gz"
        if path.exists():
            print(f"  {name}: onbellek var, atlaniyor ({core.rel(path)})")
            continue
        lo, hi = gstart - PADS[1], gend + PADS[1]
        print(f"  {name} [{role}] chr{chrom}:{lo}-{hi} ({len(codes)} IDP)")
        rows = asyncio.run(sweep(codes, chrom, lo, hi, label=name))
        buf = io.StringIO()
        cols = ["phenocode", "chr", "pos", "rsid", "ref", "alt", "af", "maf", "beta", "se", "p",
                "nearest_genes"]
        w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
        content = gzip.compress(buf.getvalue().encode("utf-8"))
        path.write_bytes(content)
        core.save_bytes("imaging_gwas", path.name, content, base=INTERIM.parent,
                        source=f"{PHEWEB}/api/region/<IDP>/lz-results/?filter=chr{chrom}:{lo}-{hi}",
                        note=f"BIG40 33k, {name} +-250kb, {len(codes)} IDP, p<0.05 satirlari",
                        n_idp=len(codes), window=f"{chrom}:{lo}-{hi}", build="GRCh37")
        print(f"  -> {core.rel(path)} ({len(rows)} satir)")


# ---------------------------------------------------------------------------------------------
# stage: sumstats  (the only bulk download, two IDPs, streamed prefix)
# ---------------------------------------------------------------------------------------------

def stage_sumstats(pheno_ids: list[str], chrom: str = "6", lo: int = 0, hi: int = 0) -> None:
    """Stream release2/stats33k/<pheno>.txt.gz and keep only the rows inside the window.

    The file is sorted by chromosome then position and is a single-member gzip, so the prefix has
    to be streamed; the loop stops as soon as the window has been passed. Nothing but the
    extracted window is written to disk.
    """
    print("== stage sumstats ==")
    if not hi:
        g = gene_grch37("KIF13A")
        chrom = g["seq_region_name"]
        lo, hi = int(g["start"]) - PADS[1], int(g["end"]) + PADS[1]
    ctag = f"{int(chrom):02d}" if chrom.isdigit() else chrom
    for pid in pheno_ids:
        out = INTERIM / f"sumstats_{pid}_chr{chrom}_{lo}_{hi}.tsv"
        if out.exists():
            print(f"  {pid}: onbellek var, atlaniyor")
            continue
        url = f"{BIG40}/release2/stats33k/{pid}.txt.gz"
        print(f"  {pid}: {url} -> chr{chrom}:{lo}-{hi}", flush=True)
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
        tail = b""
        kept: list[dict] = []
        header: list[str] = []
        got = 0
        t0 = time.time()
        passed = False
        with httpx.stream("GET", url, timeout=600, follow_redirects=True, headers=UA) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            for chunk in r.iter_bytes(1 << 20):
                got += len(chunk)
                data = tail + dec.decompress(chunk)
                *lines, tail = data.split(b"\n")
                for raw in lines:
                    line = raw.decode()
                    if not header:
                        header = line.split()
                        continue
                    f = line.split()
                    if f[0] != ctag:
                        if kept:
                            passed = True
                        continue
                    pos = int(f[2])
                    if pos < lo:
                        continue
                    if pos > hi:
                        passed = True
                        break
                    kept.append({"chr": f[0], "rsid": f[1], "pos": pos, "a1": f[3], "a2": f[4],
                                 "beta": f[5], "se": f[6], "neglog10p": f[7]})
                if passed:
                    break
                if got % (20 << 20) < (1 << 20):
                    print(f"    {got/1048576:.0f}/{total/1048576:.0f} MB, {len(kept)} satir, "
                          f"{time.time()-t0:.0f}s", flush=True)
        write_tsv(out, kept)
        core.save_bytes("imaging_gwas", out.name, out.read_bytes(), base=INTERIM.parent, source=url,
                        note=f"BIG40 33k sumstats IDP {pid}, chr{chrom}:{lo}-{hi} penceresi, "
                             f"MAF>=0.001 tam kapsam, akis onekinden ({got/1048576:.0f} MB indirildi)",
                        window=f"{chrom}:{lo}-{hi}", build="GRCh37", n_records=len(kept))


# ---------------------------------------------------------------------------------------------
# stage: followup  (a locus is not a gene: PheWAS profile, eQTL, LD, L2G, published hits)
# ---------------------------------------------------------------------------------------------

# GRCh37 positions, every one of them resolved live from Ensembl GRCh37 before use
FOLLOWUP_VARIANTS = [
    ("rs143422697", "6-17828055-T-G", "BIG40'in tek 5e-8 gecen hiti (gri-ak kontrasti)"),
    ("rs4712326", "6-18052156-C-T", "rs143422697 ile r2=0,83"),
    ("rs4716216", "6-18074064-G-A", "rs143422697 ile r2=0,66"),
    ("rs75428660", "6-17886045-A-G", "yayimlanmis kortikal kalinlik / sulkus derinligi hiti, KIF13A introniki"),
    ("rs76091407", "6-18020433-G-A", "yayimlanmis ak madde mikroyapi hiti (UDIP-FA), KIF13A'dan 32,6 kb"),
    ("rs144963738", "6-18012583-G-C", "BIG40'in en iyi ak madde mikroyapi hiti"),
]


def _hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for the hypergeometric distribution, exact, in log space."""
    from math import lgamma, exp
    def lc(a: int, b: int) -> float:
        if b < 0 or b > a:
            return -math.inf
        return lgamma(a + 1) - lgamma(b + 1) - lgamma(a - b + 1)
    denom = lc(N, n)
    tot = 0.0
    for i in range(k, min(K, n) + 1):
        tot += exp(lc(K, i) + lc(N - K, n - i) - denom)
    return min(1.0, tot)


def stage_followup() -> None:
    print("== stage followup ==")
    cls = {c["phenocode"]: c for c in read_tsv(INTERIM / "idp_classification.tsv")}
    n_by_class: dict[str, int] = {}
    for c in cls.values():
        n_by_class[c["sinif"]] = n_by_class.get(c["sinif"], 0) + 1
    rows = []
    for rsid, vid, note in FOLLOWUP_VARIANTS:
        r = _get(f"{PHEWEB}/api/variant/{vid}")
        j = r.json()
        if not j:
            print(f"  {rsid}: BIG40'ta yok")
            continue
        core.save_bytes("imaging_gwas", f"pheweb33k_variant_{rsid}.json", r.content,
                        source=f"{PHEWEB}/api/variant/{vid}",
                        note=f"BIG40 33k tum IDP PheWAS, {rsid} ({note})")
        ph = sorted(j["phenos"], key=lambda z: z["pval"])
        got = [x for x in ph if x["phenocode"] in cls]
        for topn in (20, len(got)):
            sub = got[:topn]
            for klass, K in n_by_class.items():
                k = sum(1 for x in sub if cls[x["phenocode"]]["sinif"] == klass)
                if k == 0:
                    continue
                rows.append({
                    "rsid": rsid, "grch37": vid, "not": note,
                    "kume": f"en iyi {topn}" if topn == 20 else f"p<=0.01 olan {topn}",
                    "sinif": klass, "gozlenen": k, "beklenen": round(topn * K / len(cls), 2),
                    "hipergeometrik_p_tek_yonlu": f"{_hypergeom_sf(k, len(cls), K, topn):.2g}",
                    "en_iyi_p": f"{min((x['pval'] for x in sub if cls[x['phenocode']]['sinif'] == klass), default=float('nan')):.3g}",
                })
        best = got[0]
        print(f"  {rsid} {vid}: p<=0.01 olan {len(got)} IDP, en iyi {best['phenocode']} "
              f"p={best['pval']:.3g} ({cls[best['phenocode']]['sinif']})")
    write_tsv(INTERIM / "variant_phewas_class_enrichment.tsv", rows)
    print("  UYARI: IDP'ler birbiriyle yuksek korelasyonlu (ayni traktin 9 difuzyon olcutu), "
          "dolayisiyla hipergeometrik p degeri gercek bir anlamlilik testi degil, yalniz "
          "listenin tepesinin ne kadar tek bir sinifa toplandiginin olcusudur.")
    _gwascat_window()
    _locus_to_gene()


def _gwascat_window() -> None:
    """Every published GWAS association in the locus, from the GWAS Catalog region endpoint.

    The Catalog reports GRCh38, BIG40 is GRCh37; at this locus the two builds differ by only
    231 bp, so a position from the wrong build silently looks valid. The window is therefore
    stated in GRCh38 and converted from the GRCh37 gene coordinates via the Ensembl GRCh38
    lookup, not by adding an offset.
    """
    print("  GWAS Catalog bolge taramasi (yayimlanmis her iliski)")
    g38 = _get(f"{ENSEMBL38}/lookup/symbol/homo_sapiens/KIF13A?expand=0",
               headers={"Accept": "application/json"}).json()
    lo, hi = int(g38["start"]) - PADS[1], int(g38["end"]) + PADS[1]
    base = "https://www.ebi.ac.uk/gwas/rest/api"
    url = (f"{base}/singleNucleotidePolymorphisms/search/findByChromBpLocationRange"
           f"?chrom={g38['seq_region_name']}&bpStart={lo}&bpEnd={hi}&size=500")
    snps = _get(url).json().get("_embedded", {}).get("singleNucleotidePolymorphisms", [])
    rows = []
    for s in snps:
        rs = s["rsId"]
        pos = next((l["chromosomePosition"] for l in s.get("locations", [])), None)
        ctx = {c["gene"]["geneName"]: c.get("distance") for c in s.get("genomicContexts", [])}
        try:
            assoc = _get(f"{base}/singleNucleotidePolymorphisms/{rs}/associations"
                         f"?projection=associationBySnp").json()
        except Exception as exc:  # noqa: BLE001
            print(f"    {rs}: iliski cekilemedi ({exc})")
            continue
        for a in assoc.get("_embedded", {}).get("associations", []):
            if not a.get("pvalue"):
                continue
            rows.append({"rsid": rs, "pos_grch38": pos, "p": a["pvalue"],
                         "beta": a.get("betaNum"), "beta_unit": a.get("betaUnit"),
                         "beta_yon": a.get("betaDirection"), "or": a.get("orPerCopyNum"),
                         "efo_traits": ";".join(t["trait"] for t in a.get("efoTraits", [])),
                         "kif13a_uzaklik": ctx.get("KIF13A"), "nhlrc1_uzaklik": ctx.get("NHLRC1"),
                         "nup153_uzaklik": ctx.get("NUP153")})
        time.sleep(0.15)
    rows.sort(key=lambda r: r["p"])
    write_tsv(INTERIM / "gwascat_window_associations.tsv", rows)
    core.save_bytes("imaging_gwas", "gwascat_window_associations.tsv",
                    (INTERIM / "gwascat_window_associations.tsv").read_bytes(),
                    base=INTERIM.parent, source=url,
                    note=f"GWAS Catalog, KIF13A +-250kb (GRCh38 chr6:{lo}-{hi}), {len(snps)} SNP, "
                         f"{len(rows)} iliski")
    wm = [r for r in rows if any(k in r["efo_traits"].lower() for k in
                                 ("white matter", "hyperintens", "anisotrop", "diffusi", "myelin",
                                  "microstructure"))]
    print(f"    {len(snps)} SNP, {len(rows)} iliski; ak madde fenotipi olan: {len(wm)}")
    for r in wm:
        print(f"      p={r['p']:.1e} {r['rsid']} chr6:{r['pos_grch38']} | {r['efo_traits']} "
              f"| KIF13A uzaklik {r['kif13a_uzaklik']}")


def _locus_to_gene() -> None:
    """eQTL, sQTL, LD and L2G evidence for the locus leads, all saved with provenance."""
    print("  lokus-gen kanitlari")
    gtex = "https://gtexportal.org/api/v2"
    meta = "https://eqtlgen.molgeniscloud.org/api/v2/metabrain_cis_eqtl_v2"
    ot = "https://api.platform.opentargets.org/api/v4/graphql"
    # GTEx: the variant must first be resolved to a GTEx id, which also reveals whether it was
    # eligible for eQTL testing at all (maf01 = false means "not tested", not "not an eQTL").
    for rsid, _vid, _note in FOLLOWUP_VARIANTS:
        r = _get(f"{gtex}/dataset/variant?snpId={rsid}", headers={"Accept": "application/json"})
        d = r.json().get("data", [])
        core.save_bytes("imaging_gwas", f"gtex_variant_{rsid}.json", r.content,
                        source=f"{gtex}/dataset/variant?snpId={rsid}",
                        note=f"GTEx v8 varyant kaydi {rsid}; maf01=false ise eQTL testine girmemis")
        if not d:
            print(f"    {rsid}: GTEx'te yok")
            continue
        v = d[0]
        if not v.get("maf01"):
            print(f"    {rsid} ({v['variantId']}): GTEx maf01=false -> eQTL testine girmemis")
            continue
        for ep in ("singleTissueEqtl", "singleTissueSqtl"):
            rr = _get(f"{gtex}/association/{ep}?variantId={v['variantId']}",
                      headers={"Accept": "application/json"})
            core.save_bytes("imaging_gwas", f"gtex_{ep}_{rsid}.json", rr.content,
                            source=f"{gtex}/association/{ep}?variantId={v['variantId']}",
                            note=f"GTEx v8 {ep}, {rsid}")
            hits = rr.json().get("data", [])
            brain = [h for h in hits if str(h.get("tissueSiteDetailId", "")).startswith("Brain")]
            print(f"    {rsid} GTEx {ep}: {len(hits)} anlamli satir, {len(brain)} beyin dokusu; "
                  f"genler {sorted({h['geneSymbol'] for h in hits})}")
    # MetaBrain: the only brain cis-eQTL query interface still answering (eQTL Catalogue's REST
    # API returns HTTP 410, PsychENCODE is bulk download only).
    for gene in ("KIF13A", "NUP153", "NHLRC1", "TPMT"):
        r = _get(f"{meta}?q=GeneSymbol=={gene}&num=50")
        core.save_bytes("imaging_gwas", f"metabrain_cis_eqtl_{gene}.json", r.content,
                        source=f"{meta}?q=GeneSymbol=={gene}",
                        note=f"MetaBrain cis-eQTL en iyi etkiler, {gene}, GRCh38")
        for x in r.json().get("items", []):
            print(f"    MetaBrain {gene:8s} {x.get('Brainregion'):16s} sira {x.get('QTLRank')} "
                  f"{x.get('SNP')} p={x.get('MetaP'):.2e} beta={x.get('MetaBeta')} N={x.get('MetaPN')}")
    # Open Targets credible sets and L2G (constraint numbers are never taken from Open Targets)
    q = """query V($v:String!){ variant(variantId:$v){ id rsIds
      credibleSets(page:{index:0,size:30}){ count rows{ pValueMantissa pValueExponent
        study{ id studyType traitFromSource projectId nSamples pubmedId
               biosample{ biosampleName } }
        l2GPredictions(page:{index:0,size:5}){ rows{ score target{ approvedSymbol } } } } } } }"""
    for rsid, vid38 in (("rs75428660", "6_17885814_A_G"), ("rs143422697", "6_17827824_T_G"),
                        ("rs76091407", "6_18020202_G_A")):
        body = json.dumps({"query": q, "variables": {"v": vid38}}).encode()
        resp = httpx.post(ot, content=body, timeout=180,
                          headers={**UA, "Content-Type": "application/json"})
        core.save_bytes("imaging_gwas", f"opentargets_credsets_{rsid}.json", resp.content,
                        source=f"{ot} (POST, variantId {vid38})",
                        note=f"Open Targets kredibl kumeler + L2G, {rsid}")
        var = (resp.json().get("data") or {}).get("variant")
        if not var:
            print(f"    {rsid}: Open Targets varyant dizininde yok")
            continue
        for row in var["credibleSets"]["rows"]:
            s = row["study"]
            l2g = [(round(x["score"], 3), x["target"]["approvedSymbol"])
                   for x in row["l2GPredictions"]["rows"]]
            print(f"    OT {rsid} {s['studyType']:5s} {str(s.get('traitFromSource'))[:34]:34s} "
                  f"p={row['pValueMantissa']}e{row['pValueExponent']} N={s['nSamples']} L2G {l2g}")
    # LD blocks, 1000 Genomes EUR, GRCh37 (the build BIG40 uses)
    ld_rows = []
    for lead in ("rs143422697", "rs75428660", "rs76091407", "rs9297056"):
        r = _get(f"{ENSEMBL37}/ld/human/{lead}/1000GENOMES:phase_3:EUR?window_size=500&r2=0.2",
                 headers={"Accept": "application/json"}, timeout=300)
        core.save_bytes("imaging_gwas", f"ld_block_{lead}.json", r.content,
                        source=f"{ENSEMBL37}/ld/human/{lead}/1000GENOMES:phase_3:EUR"
                               f"?window_size=500&r2=0.2",
                        note=f"{lead} LD blogu, 1000G faz 3 EUR, r2>=0.2, GRCh37")
        partners = r.json()
        ids = [x["variation2"] for x in partners]
        pos: dict[str, int] = {}
        for i in range(0, len(ids), 180):
            rr = httpx.post(f"{ENSEMBL37}/variation/homo_sapiens",
                            json={"ids": ids[i:i + 180]}, timeout=300,
                            headers={**UA, "Content-Type": "application/json",
                                     "Accept": "application/json"})
            for k, v in rr.json().items():
                for m in v.get("mappings", []):
                    if m["seq_region_name"] == "6":
                        pos[k] = m["start"]
        ps = sorted(pos.values())
        ld_rows.append({
            "oncu": lead, "r2_ge_0.2_ortak": len(ids),
            "blok_grch37": f"chr6:{ps[0]}-{ps[-1]}" if ps else "",
            "blok_kb": round((ps[-1] - ps[0]) / 1000) if ps else 0,
            "KIF13A_govdesinde": sum(1 for p in ps if 17_759_414 <= p <= 17_987_854),
            "NUP153_govdesinde": sum(1 for p in ps if 17_615_269 <= p <= 17_706_656),
            "NHLRC1_govdesinde": sum(1 for p in ps if 18_120_718 <= p <= 18_122_851),
            "TPMT_govdesinde": sum(1 for p in ps if 18_128_542 <= p <= 18_155_305),
            "NHLRC1_konum_olarak_kapsaniyor": "yes" if ps and ps[0] <= 18_120_718
                                              and ps[-1] >= 18_122_851 else "no",
        })
        print("    LD", ld_rows[-1])
    write_tsv(INTERIM / "ld_blocks.tsv", ld_rows)


# ---------------------------------------------------------------------------------------------
# stage: external  (published GWAS summary statistics from the GWAS Catalog FTP)
# ---------------------------------------------------------------------------------------------

# Accession -> (FTP file, what it is). Both are GRCh37, 1-based, GWAS-SSF v1.0, stated in the
# accompanying -meta.yaml which is pulled and saved next to the extract. Both are under 170 MB,
# so unlike BIG40's per-IDP files they can be streamed whole in a few minutes; neither is sorted,
# so the whole file has to pass through.
EXTERNAL_SUMSTATS = {
    "GCST90860937": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90860001-GCST90861000/GCST90860937/GCST90860937.tsv.gz",
        "UDIP-FA, ak madde mikroyapisi, N=11.736, PMID 42236496; rs76091407 p=4e-09 iddiasinin kaynagi"),
    "GCST90095129": (
        "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
        "GCST90095001-GCST90096000/GCST90095129/GCST90095129_buildGRCh37.tsv.gz",
        "Vertex-wise sulcal depth, N=33.748, PMID 34910505; rs75428660 p=1e-09 iddiasinin kaynagi"),
}


def stage_external(accessions: list[str], chrom: str = "6", lo: int = 0, hi: int = 0) -> None:
    print("== stage external ==")
    if not hi:
        g = gene_grch37("KIF13A")
        chrom, lo, hi = g["seq_region_name"], int(g["start"]) - PADS[1], int(g["end"]) + PADS[1]
    for acc in accessions:
        url, note = EXTERNAL_SUMSTATS[acc]
        out = INTERIM / f"external_{acc}_chr{chrom}_{lo}_{hi}.tsv"
        if out.exists():
            print(f"  {acc}: onbellek var, atlaniyor")
            continue
        meta = _get(url + "-meta.yaml")
        core.save_bytes("imaging_gwas", f"{acc}_meta.yaml", meta.content, source=url + "-meta.yaml",
                        note=f"{acc} ozet istatistik ust verisi (yapi, koordinat sistemi, N)")
        print(f"  {acc}: {url}", flush=True)
        dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
        tail = b""
        header: list[str] = []
        kept: list[dict] = []
        got = 0
        t0 = time.time()
        with httpx.stream("GET", url, timeout=1800, follow_redirects=True, headers=UA) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            for chunk in r.iter_bytes(1 << 20):
                got += len(chunk)
                data = tail + dec.decompress(chunk)
                *lines, tail = data.split(b"\n")
                for raw in lines:
                    line = raw.decode()
                    if not header:
                        header = line.rstrip("\r").split("\t")
                        continue
                    f = line.rstrip("\r").split("\t")
                    rec = dict(zip(header, f))
                    if rec.get("chromosome") != chrom:
                        continue
                    try:
                        pos = int(rec["base_pair_location"])
                    except (KeyError, ValueError):
                        continue
                    if lo <= pos <= hi:
                        kept.append(rec)
                if got % (30 << 20) < (1 << 20):
                    print(f"    {got/1048576:.0f}/{total/1048576:.0f} MB, {len(kept)} satir, "
                          f"{time.time()-t0:.0f}s", flush=True)
        write_tsv(out, kept)
        core.save_bytes("imaging_gwas", out.name, out.read_bytes(), base=INTERIM.parent, source=url,
                        note=f"{note}; chr{chrom}:{lo}-{hi} penceresi (GRCh37)",
                        window=f"{chrom}:{lo}-{hi}", build="GRCh37", n_records=len(kept))
        if kept and "p_value" in header:
            best = min(kept, key=lambda r: float(r["p_value"]))
            n5 = sum(1 for r in kept if float(r["p_value"]) < 5e-8)
            print(f"    en iyi p = {best['p_value']} @ chr{chrom}:{best['base_pair_location']} "
                  f"({best.get('rsid', best.get('variant_id', '?'))}); 5e-8 gecen {n5} varyant")


# ---------------------------------------------------------------------------------------------
# stage: thresholds  (never assume BIG40's own threshold, read it out of the paper)
# ---------------------------------------------------------------------------------------------

def stage_thresholds() -> dict[str, float]:
    print("== stage thresholds ==")
    xml = RAW / "Smith2021_NatNeuro_PMC7610742.xml"
    if not xml.exists():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               "?db=pmc&id=PMC7610742&retmode=xml")
        r = _get(url)
        core.save_bytes("imaging_gwas", xml.name, r.content, source=url,
                        note="Smith 2021 Nat Neurosci (BIG40 yayini) tam metin, anlamlilik esigi icin")
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml.read_text(encoding="utf-8")))
    # the paper writes the two thresholds as -Log10(P); numbers are extracted, not assumed
    single = re.search(r"standard GWAS threshold of .{0,12}Log10\(\s*P\s*\)\s*\D{0,5}(\d{1,2})\D{1,3}(\d)",
                       text)
    bonf = re.search(r"[Bb]onferroni corrected level of .{0,12}Log10\(\s*P\s*\)\s*\D{0,5}(\d{1,2})\.(\d)",
                     text)
    t_single = float(f"{single.group(1)}.{single.group(2)}") if single else None
    t_bonf = float(f"{bonf.group(1)}.{bonf.group(2)}") if bonf else None
    n_idp = len(read_tsv(INTERIM / "idp_classification.tsv"))
    derived = (t_single or 7.5) + math.log10(n_idp)
    out = {
        "conventional_5e-8": 5e-8,
        "big40_single_pheno_neglog10": t_single,
        "big40_single_pheno_p": 10 ** (-t_single) if t_single else None,
        "big40_bonferroni_neglog10": t_bonf,
        "big40_bonferroni_p": 10 ** (-t_bonf) if t_bonf else None,
        "bonferroni_recomputed_neglog10": round(derived, 3),
        "bonferroni_recomputed_p": 10 ** (-derived),
        "n_idp": n_idp,
        "source": "PMID 33875891 / PMC7610742 / doi 10.1038/s41593-021-00826-4",
    }
    for k, v in out.items():
        print(f"  {k}: {v}")
    (INTERIM / "thresholds.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out


# ---------------------------------------------------------------------------------------------
# stage: report
# ---------------------------------------------------------------------------------------------

def _load_scan(name: str) -> list[dict]:
    path = INTERIM / f"scan_{name}.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _best_per_pheno(rows: list[dict], lo: int, hi: int) -> dict[str, dict]:
    best: dict[str, dict] = {}
    for r in rows:
        pos = int(r["pos"])
        if pos < lo or pos > hi:
            continue
        p = float(r["p"])
        cur = best.get(r["phenocode"])
        if cur is None or p < cur["p"]:
            best[r["phenocode"]] = {"p": p, "pos": pos, "rsid": r["rsid"], "ref": r["ref"],
                                    "alt": r["alt"], "af": r["af"], "beta": r["beta"],
                                    "se": r["se"], "nearest_genes": r["nearest_genes"]}
    return best


def _z(alpha: float) -> float:
    """Two-sided normal critical value, via the inverse error function (no scipy dependency)."""
    from math import erf, sqrt
    lo_, hi_ = 0.0, 40.0
    target = 1.0 - alpha / 2.0
    for _ in range(200):
        mid = (lo_ + hi_) / 2
        if 0.5 * (1 + erf(mid / sqrt(2))) < target:
            lo_ = mid
        else:
            hi_ = mid
    return (lo_ + hi_) / 2


def stage_report() -> None:
    print("== stage report ==")
    thr = json.loads((INTERIM / "thresholds.json").read_text(encoding="utf-8"))
    cls = read_tsv(INTERIM / "idp_classification.tsv")
    by_code = {c["phenocode"]: c for c in cls}
    genes = {g: gene_grch37(g) for g in {v["gene"] for v in LOCI.values()}}

    t_conv, t_big40, t_bonf = 5e-8, thr["big40_single_pheno_p"], thr["big40_bonferroni_p"]

    summary: list[dict] = []
    tophits: list[dict] = []
    scans = sorted(p.name[5:-7] for p in INTERIM.glob("scan_*.tsv.gz"))
    for name in scans:
        spec = LOCI.get(name)
        if spec:
            g = genes[spec["gene"]]
            gstart, gend, role = int(g["start"]), int(g["end"]), spec["role"]
            sweep_kind = spec["sweep"]
        else:  # random window: the stored window already is gene-equivalent +- 250 kb
            role, sweep_kind = "rastgele eslesmis pencere", "wm"
            rows0 = _load_scan(name)
            gstart, gend = (min(int(r["pos"]) for r in rows0), max(int(r["pos"]) for r in rows0)) \
                if rows0 else (0, 0)
        rows = _load_scan(name)
        # the one IDP whose PheWeb region query answers HTTP 500 is filled from full sumstats
        failed = INTERIM / f"scan_{name}_failed.txt"
        patched: list[str] = []
        if failed.exists():
            for ph in failed.read_text(encoding="utf-8").split():
                pid = by_code[ph]["pheno_id"]
                cand = sorted(INTERIM.glob(f"sumstats_{pid}_chr*.tsv"))
                if not cand:
                    print(f"  !! {name}: {ph} icin sumstats yok, sinif minimumu eksik kalir")
                    continue
                for r in read_tsv(cand[0]):
                    rows.append({"phenocode": ph, "chr": r["chr"], "pos": r["pos"], "rsid": r["rsid"],
                                 "ref": r["a1"], "alt": r["a2"], "af": "", "maf": "",
                                 "beta": r["beta"], "se": r["se"],
                                 "p": str(10 ** (-float(r["neglog10p"]))), "nearest_genes": ""})
                patched.append(ph)
        for pad in PADS:
            lo, hi = gstart - pad, gend + pad
            best = _best_per_pheno(rows, lo, hi)
            for klass in (CLASS_WM_MICRO, CLASS_WM_MACRO, CLASS_OTHER, CLASS_QC):
                codes = [c["phenocode"] for c in cls if c["sinif"] == klass]
                if sweep_kind == "wm" and klass in (CLASS_OTHER, CLASS_QC):
                    continue
                hits = {c: best[c] for c in codes if c in best}
                if hits:
                    top = min(hits.items(), key=lambda kv: kv[1]["p"])
                else:
                    top = (None, None)
                summary.append({
                    "lokus": name, "rol": role, "pencere_kb": pad // 1000, "sinif": klass,
                    "n_idp": len(codes), "n_idp_p_le_0.01": len(hits),
                    "min_p": f"{top[1]['p']:.3g}" if top[1] else ">0.01",
                    "min_p_idp": top[0] or "", "lead_rsid": top[1]["rsid"] if top[1] else "",
                    "lead_pos": top[1]["pos"] if top[1] else "",
                    "n_p_lt_5e-8": sum(1 for v in hits.values() if v["p"] < t_conv),
                    "n_p_lt_big40_7.5": sum(1 for v in hits.values() if v["p"] < t_big40),
                    "n_p_lt_bonferroni_11.1": sum(1 for v in hits.values() if v["p"] < t_bonf),
                    "sumstats_ile_tamamlanan": ",".join(patched),
                })
                if pad == PADS[1]:
                    for code, v in sorted(hits.items(), key=lambda kv: kv[1]["p"])[:6]:
                        tophits.append({
                            "lokus": name, "sinif": klass, "phenocode": code,
                            "aciklama": by_code.get(code, {}).get("aciklama", "")[:110],
                            "ukb_field": by_code.get(code, {}).get("ukb_field", ""),
                            "n": by_code.get(code, {}).get("n_all", ""),
                            "h2": by_code.get(code, {}).get("h2", ""),
                            "rsid": v["rsid"], "chr_pos_grch37": f"{v.get('nearest_genes','')}"
                                                                 f" chr6:{v['pos']}" if False else v["pos"],
                            "ref": v["ref"], "alt": v["alt"], "af": v["af"], "beta": v["beta"],
                            "se": v["se"], "p": f"{v['p']:.3g}", "nearest_genes": v["nearest_genes"],
                        })
    write_tsv(INTERIM / "summary_by_class.tsv", summary)
    write_tsv(INTERIM / "top_hits.tsv", tophits)

    # ---- empirical null over the random windows, white matter class only
    nulls = []
    for row in summary:
        if row["lokus"].startswith("RANDOM") and row["pencere_kb"] == PADS[1] // 1000 \
                and row["sinif"] == CLASS_WM_MICRO:
            nulls.append(float(row["min_p"]) if row["min_p"] != ">0.01" else 0.01)
    kif = [r for r in summary if r["lokus"] == "KIF13A" and r["pencere_kb"] == PADS[1] // 1000
           and r["sinif"] == CLASS_WM_MICRO]
    if nulls and kif:
        obs = float(kif[0]["min_p"])
        rank = sum(1 for n in nulls if n <= obs)
        print(f"  ampirik bos hipotez ({len(nulls)} rastgele pencere, ak madde mikroyapi):")
        print(f"    rastgele pencerelerin en iyi p degerleri: "
              f"{', '.join(f'{n:.2g}' for n in sorted(nulls))}")
        print(f"    KIF13A gozlenen: {obs:.3g} -> {rank}/{len(nulls)} pencere bu kadar veya daha "
              f"iyi, ampirik p = {(rank+1)/(len(nulls)+1):.2f}")

    # ---- internally calibrated null for the class profile of a single variant.
    # The white matter enrichment seen at the KIF13A haplotype (stage followup) is only
    # interesting if it is not what every variant in BIG40 looks like: white matter IDPs are
    # more heritable than average, so any true signal might show it. The KIF13A and KIF1A
    # windows were both swept over all 3,935 IDPs, which gives a matched comparison at no extra
    # cost: the distribution of the white-matter-microstructure share among the phenotypes that
    # each variant hits at p <= 0.01.
    share_rows = []
    for name in ("KIF13A", "KIF1A"):
        path = INTERIM / f"scan_{name}.tsv.gz"
        if not path.exists():
            continue
        per_var: dict[tuple, dict[str, int]] = {}
        for r in _load_scan(name):
            c = by_code.get(r["phenocode"])
            if not c:
                continue
            key = (r["pos"], r["rsid"])
            per_var.setdefault(key, {})
            per_var[key][c["sinif"]] = per_var[key].get(c["sinif"], 0) + 1
        shares = []
        for (pos, rsid), counts in per_var.items():
            tot = sum(counts.values())
            if tot < 30:
                continue
            shares.append((counts.get(CLASS_WM_MICRO, 0) / tot, tot, pos, rsid))
        shares.sort(reverse=True)
        base = len([c for c in cls if c["sinif"] == CLASS_WM_MICRO]) / len(cls)
        print(f"  {name}: >=30 IDP'ye p<=0.01 ile dokunan {len(shares)} varyant, "
              f"ak madde mikroyapi payi: ortanca {shares[len(shares)//2][0]:.2f}, "
              f"maksimum {shares[0][0]:.2f} ({shares[0][3]}), beklenen taban {base:.3f}")
        for sh, tot, pos, rsid in shares[:8]:
            share_rows.append({"lokus": name, "rsid": rsid, "pos_grch37": pos, "n_idp_p_le_0.01": tot,
                               "ak_madde_mikroyapi_payi": f"{sh:.3f}", "taban_pay": f"{base:.3f}"})
    write_tsv(INTERIM / "variant_class_share.tsv", share_rows)

    # ---- coherence, comparable across every locus (all of them were swept over the 683 white
    # matter IDPs). For the variant with the best white matter p at each locus: how many white
    # matter IDPs does it hit at p <= 0.01, and are those IDPs one anatomical structure or
    # scattered? A real tract signal shows many diffusion measures of the same tract; noise does
    # not. Reported as the count plus the three best phenotype names.
    coh_rows = []
    for name in sorted(p.name[5:-7] for p in INTERIM.glob("scan_*.tsv.gz")):
        rows = _load_scan(name)
        wm = [r for r in rows if by_code.get(r["phenocode"], {}).get("sinif") == CLASS_WM_MICRO]
        if not wm:
            continue
        lead = min(wm, key=lambda r: float(r["p"]))
        same = [r for r in wm if r["pos"] == lead["pos"]]
        same.sort(key=lambda r: float(r["p"]))
        coh_rows.append({
            "lokus": name, "rol": (LOCI.get(name) or {}).get("role", "rastgele pencere"),
            "lead_rsid": lead["rsid"], "pos_grch37": lead["pos"], "af": lead["af"],
            "en_iyi_ak_madde_p": f"{float(lead['p']):.3g}",
            "ayni_varyantta_p<=0.01_ak_madde_IDP": len(same),
            "ilk_uc": " | ".join(r["phenocode"].replace("IDP_dMRI_TBSS_", "") for r in same[:3]),
        })
    write_tsv(INTERIM / "wm_coherence_by_locus.tsv", coh_rows)

    # ---- power: detectable effect from the real standard errors in the window
    print("  guc (gercek standart hatalardan, 80 % guc):")
    zb = 0.8416
    power_rows = []
    for pid, label in (("2012", "OD cingulum cingulate gyrus L (UKB 25427, N=31.355)"),
                       ("1415", "wg rh intensity-contrast medialorbitofrontal (UKB 27037, N=31.968)"),
                       ("1479", "FA posterior corona radiata L (UKB 25083, N=31.356)")):
        cand = sorted(INTERIM.glob(f"sumstats_{pid}_chr*.tsv"))
        if not cand:
            continue
        rows = read_tsv(cand[0])
        # MAF is not in the sumstats file; the PheWeb copy supplies af for p<=0.01 variants, so the
        # SE itself is binned by its own value, which is a monotone function of MAF at fixed N.
        ses = sorted(float(r["se"]) for r in rows)
        for q, qlabel in ((0.05, "en bilgili %5"), (0.25, "alt ceyrek"), (0.5, "ortanca")):
            se = ses[int(q * len(ses))]
            for alpha, alabel in ((t_conv, "5e-8"), (t_big40, f"BIG40 {thr['big40_single_pheno_neglog10']}"),
                                  (t_bonf, f"BIG40 Bonferroni {thr['big40_bonferroni_neglog10']}")):
                beta_min = (_z(alpha) + zb) * se
                power_rows.append({"idp": label, "se_kuantil": qlabel, "se": f"{se:.4g}",
                                   "esik": alabel, "saptanabilir_beta_sd": f"{beta_min:.3g}"})
        print(f"    {label}: SE ortanca {ses[len(ses)//2]:.4g}, en bilgili %5 {ses[int(0.05*len(ses))]:.4g}")
    # analytic bound: for a phenotype scaled to unit variance and additive coding,
    # chi-square ~= N * R^2, so the detectable variance explained is (z_alpha + z_beta)^2 / N and
    # does not depend on allele frequency. The frequency only decides which beta that R^2 means.
    for N, nlabel in ((31355, "N=31.355 (dMRI IDP'leri)"), (33224, "N=33.224 (T1 hacim IDP'leri)")):
        for alpha, alabel in ((t_conv, "5e-8"), (t_big40, "BIG40 7,5"), (t_bonf, "BIG40 Bonferroni 11,1")):
            r2 = (_z(alpha) + zb) ** 2 / N
            row = {"idp": nlabel, "se_kuantil": "analitik", "se": "", "esik": alabel,
                   "saptanabilir_beta_sd": "", "saptanabilir_R2_yuzde": f"{100*r2:.3f}"}
            for maf in (0.01, 0.05, 0.25):
                row[f"beta_maf_{maf}"] = f"{math.sqrt(r2 / (2 * maf * (1 - maf))):.3f}"
            power_rows.append(row)
            print(f"    {nlabel} {alabel}: saptanabilir R2 = {100*r2:.3f} %, "
                  f"beta(MAF 0,05) = {math.sqrt(r2/(2*0.05*0.95)):.3f} SD, "
                  f"beta(MAF 0,25) = {math.sqrt(r2/(2*0.25*0.75)):.3f} SD")
    write_tsv(INTERIM / "power.tsv", power_rows)
    print("  esikler:", {k: v for k, v in thr.items() if k.endswith("_p") or k == "conventional_5e-8"})


# ---------------------------------------------------------------------------------------------
# stage: figure
# ---------------------------------------------------------------------------------------------

def stage_figure() -> None:
    print("== stage figure ==")
    thr = json.loads((INTERIM / "thresholds.json").read_text(encoding="utf-8"))
    cls = {c["phenocode"]: c for c in read_tsv(INTERIM / "idp_classification.tsv")}
    g = gene_grch37("KIF13A")
    gstart, gend = int(g["start"]), int(g["end"])
    lo, hi = gstart - PADS[1], gend + PADS[1]
    rows = _load_scan("KIF13A")
    W, H = 1180, 560
    ml, mr, mt, mb = 88, 24, 58, 110
    pw, ph = W - ml - mr, H - mt - mb
    ymax = 9.0

    def x(pos: int) -> float:
        return ml + (pos - lo) / (hi - lo) * pw

    def y(p: float) -> float:
        return mt + ph - min(-math.log10(p), ymax) / ymax * ph

    colours = {CLASS_WM_MICRO: "#1f77b4", CLASS_WM_MACRO: "#17becf", CLASS_OTHER: "#bbbbbb",
               CLASS_QC: "#dddddd"}
    order = [CLASS_OTHER, CLASS_QC, CLASS_WM_MACRO, CLASS_WM_MICRO]
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'viewBox="0 0 {W} {H}" font-family="Arial, Helvetica, sans-serif">',
             f'<rect width="{W}" height="{H}" fill="white"/>']
    # gene bodies
    for sym, colour in (("KIF13A", "#333333"), ("NUP153", "#8c564b"), ("NHLRC1", "#8c564b")):
        gg = gene_grch37(sym)
        x0, x1 = x(int(gg["start"])), x(int(gg["end"]))
        parts.append(f'<rect x="{max(ml,x0):.1f}" y="{mt+ph+16}" width="{max(2.0,min(x1,ml+pw)-max(ml,x0)):.1f}" '
                     f'height="9" fill="{colour}" opacity="0.85"/>')
        parts.append(f'<text x="{(max(ml,x0)+min(x1,ml+pw))/2:.1f}" y="{mt+ph+40}" font-size="11" '
                     f'text-anchor="middle" fill="{colour}">{sym}</text>')
    # thresholds
    for val, lab, dash in ((-math.log10(5e-8), "5e-8", "6 4"),
                           (thr["big40_single_pheno_neglog10"], f"BIG40 {thr['big40_single_pheno_neglog10']}", "2 3")):
        yy = mt + ph - val / ymax * ph
        parts.append(f'<line x1="{ml}" y1="{yy:.1f}" x2="{ml+pw}" y2="{yy:.1f}" stroke="#d62728" '
                     f'stroke-width="1" stroke-dasharray="{dash}"/>')
        parts.append(f'<text x="{ml+pw-4}" y="{yy-4:.1f}" font-size="11" text-anchor="end" '
                     f'fill="#d62728">{lab}</text>')
    # One point per (variant, class): the best p that variant reaches inside that phenotype class.
    # Plotting all 164k (IDP, variant) pairs produces a 6 MB unreadable SVG; the per-class minimum
    # is what the analysis actually reports and it keeps the file small.
    best: dict[tuple, float] = {}
    for r in rows:
        c = cls.get(r["phenocode"])
        if not c:
            continue
        key = (int(r["pos"]), c["sinif"])
        p = float(r["p"])
        if p < best.get(key, 2.0):
            best[key] = p
    # The control class has 3,236 IDPs, so at a loose cut it paints the whole panel grey and hides
    # the two white matter classes; it is drawn at a stricter cut on purpose (stated in the caption).
    cut = {CLASS_WM_MICRO: 0.01, CLASS_WM_MACRO: 0.01, CLASS_OTHER: 5e-4, CLASS_QC: 5e-4}
    for klass in order:
        pts = [f'<circle cx="{x(pos):.0f}" cy="{y(p):.0f}" r="2.4"/>'
               for (pos, kl), p in best.items() if kl == klass and p <= cut[klass]]
        if pts:
            parts.append(f'<g fill="{colours[klass]}" opacity="0.8">' + "".join(pts) + "</g>")
    # label the variants the report argues about
    for rsid, pos, klass, dy in (("rs143422697", 17828055, CLASS_OTHER, -10),
                                 ("rs75428660", 17886045, CLASS_WM_MICRO, -10),
                                 ("rs144963738", 18012583, CLASS_WM_MICRO, -30),
                                 ("rs76091407", 18020433, CLASS_WM_MICRO, 16),
                                 ("rs9297056 (TPMT)", 18145411, CLASS_WM_MICRO, -10)):
        p = best.get((pos, klass))
        if not p:
            continue
        parts.append(f'<circle cx="{x(pos):.1f}" cy="{y(p):.1f}" r="4.2" fill="none" '
                     f'stroke="#222" stroke-width="1.2"/>')
        parts.append(f'<text x="{x(pos):.1f}" y="{y(p)+dy:.1f}" font-size="10.5" '
                     f'text-anchor="middle" fill="#222">{rsid}</text>')
    # +-100 kb window edges
    for edge in (gstart - PADS[0], gend + PADS[0]):
        parts.append(f'<line x1="{x(edge):.1f}" y1="{mt}" x2="{x(edge):.1f}" y2="{mt+ph}" '
                     f'stroke="#999" stroke-width="1" stroke-dasharray="3 4"/>')
    parts.append(f'<text x="{x(gend + PADS[0])-4:.1f}" y="{mt+13}" font-size="10" '
                 f'text-anchor="end" fill="#999">+-100 kb</text>')
    # axes
    parts.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph}" stroke="#333"/>')
    parts.append(f'<line x1="{ml}" y1="{mt+ph}" x2="{ml+pw}" y2="{mt+ph}" stroke="#333"/>')
    for v in range(0, int(ymax) + 1, 2):
        yy = mt + ph - v / ymax * ph
        parts.append(f'<line x1="{ml-5}" y1="{yy:.1f}" x2="{ml}" y2="{yy:.1f}" stroke="#333"/>')
        parts.append(f'<text x="{ml-9}" y="{yy+4:.1f}" font-size="11" text-anchor="end">{v}</text>')
    for pos in range(17_600_000, 18_300_000, 100_000):
        if lo <= pos <= hi:
            parts.append(f'<line x1="{x(pos):.1f}" y1="{mt+ph}" x2="{x(pos):.1f}" y2="{mt+ph+5}" stroke="#333"/>')
            parts.append(f'<text x="{x(pos):.1f}" y="{mt+ph+64}" font-size="10" text-anchor="middle">'
                         f'{pos/1e6:.1f} Mb</text>')
    parts.append(f'<text x="{18}" y="{mt+ph/2}" font-size="12" text-anchor="middle" '
                 f'transform="rotate(-90 18 {mt+ph/2})">-log10(p)</text>')
    parts.append(f'<text x="{ml}" y="24" font-size="15" font-weight="bold">'
                 f'KIF13A lokusu, UK Biobank beyin goruntuleme GWAS (BIG40 33k, GRCh37)</text>')
    parts.append(f'<text x="{ml}" y="42" font-size="12" fill="#555">3.935 IDP, KIF13A +-250 kb; '
                 f'nokta = bir varyantin o sinifta ulastigi en iyi p (ak madde p&lt;=0,01; kontrol sinifi p&lt;=5e-4)</text>')
    lx = ml
    for klass, lab in ((CLASS_WM_MICRO, "ak madde mikroyapi (675 IDP)"),
                       (CLASS_WM_MACRO, "ak madde makroyapi (8)"),
                       (CLASS_OTHER, "gri madde ve diger (3.236)")):
        parts.append(f'<circle cx="{lx+6}" cy="{H-22}" r="4" fill="{colours[klass]}"/>')
        parts.append(f'<text x="{lx+16}" y="{H-18}" font-size="11">{lab}</text>')
        lx += 20 + 7.0 * len(lab)
    parts.append("</svg>")
    FIGS.mkdir(parents=True, exist_ok=True)
    out = FIGS / "KIF13A_ak_madde_goruntuleme.svg"
    out.write_text("\n".join(parts), encoding="utf-8")
    print(f"  -> {core.rel(out)}")


# ---------------------------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["sources", "probe", "classify", "scan", "sumstats", "external",
                             "thresholds", "followup", "report", "figure"])
    ap.add_argument("--phenos", default="2012,1415")
    ap.add_argument("--only", default=None)
    ap.add_argument("--accessions", default="GCST90860937,GCST90095129")
    ap.add_argument("--chrom", default="6")
    ap.add_argument("--lo", type=int, default=0)
    ap.add_argument("--hi", type=int, default=0)
    a = ap.parse_args()
    RAW.mkdir(parents=True, exist_ok=True)
    INTERIM.mkdir(parents=True, exist_ok=True)
    if a.stage == "sources":
        stage_sources()
    elif a.stage == "probe":
        stage_probe()
    elif a.stage == "classify":
        stage_classify()
    elif a.stage == "scan":
        stage_scan(a.only)
    elif a.stage == "sumstats":
        stage_sumstats([p.strip() for p in a.phenos.split(",") if p.strip()], a.chrom, a.lo, a.hi)
    elif a.stage == "external":
        stage_external([x.strip() for x in a.accessions.split(",") if x.strip()],
                       a.chrom, a.lo, a.hi)
    elif a.stage == "thresholds":
        stage_thresholds()
    elif a.stage == "followup":
        stage_followup()
    elif a.stage == "report":
        stage_report()
    elif a.stage == "figure":
        stage_figure()
    else:
        raise SystemExit(f"bilinmeyen stage {a.stage}")


if __name__ == "__main__":
    main()
