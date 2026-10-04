"""The loss-of-function axis: does haploinsufficiency explain KIF13A better than missense does?

The project was built on motor-domain missense, transferred from KIF1A. Every human observation so
far points the other way, so this script tests the loss-of-function hypothesis on public data and is
written to be able to reject it.

  1. every gnomAD v4.1 LoF variant on the MANE transcript of KIF13A, KIF1A and KIF13B, broken down by
     consequence, allele count, homozygotes and position, with the ones that are predicted to escape
     nonsense-mediated decay separated out. NMD escape is computed from the Ensembl MANE exon
     structure (last exon, or within 50 nt of the last exon-exon junction), not guessed, and is
     cross-checked against LOFTEE's own END_TRUNC filter that gnomAD ships
  2. gnomAD constraint for a panel of genes whose haploinsufficiency ClinGen has actually curated
     (dosage score read live from the ClinGen curation list), so LOEUF 0.449 can be read in context
  3. every human loss-of-function observation: ClinVar SNV LoF, ClinVar CNVs with the protein-coding
     gene content of each one resolved on Ensembl, DECIPHER gene-level dosage metrics, Europe PMC
  4. IMPC mouse knockout phenotypes, asking specifically what the heterozygotes show

Nothing is re-downloaded if it is already in data/raw (see data/raw/PROVENANCE.jsonl).

    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\lof_axis.py
    mcp\\atlas_server\\.venv\\Scripts\\python.exe -u scripts\\lof_axis.py --refresh --no-figure
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core  # noqa: E402

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL = "https://rest.ensembl.org"
ENSEMBL37 = "https://grch37.rest.ensembl.org"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
IMPC_SOLR = "https://www.ebi.ac.uk/mi/impc/solr"
CLINGEN_DOSAGE = "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
DECIPHER_GENE = "https://www.deciphergenomics.org/gene/%s/overview/clinical-info"
BROWSER_UA = ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
              "(KHTML, like Gecko) Chrome/124 Safari/537.36")

GENES = ["KIF13A", "KIF1A", "KIF13B"]

# ClinGen haploinsufficiency score 3 ("sufficient evidence"), neurodevelopmental, autosomal.
# The score itself is re-read live from CLINGEN_DOSAGE, never hard-coded into the output.
PANEL = ["ARID1B", "SYNGAP1", "ANKRD11", "CHD8", "DYRK1A", "STXBP1", "SCN1A", "TCF4",
         "ADNP", "SETD5", "NSD1", "POGZ"]

# gnomAD/VEP consequences that LOFTEE considers; everything else is not a LoF candidate
LOF_CONSEQUENCES = {"stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant"}

# KIF13A domains, MANE numbering (same table as scripts/figure_candidates.py)
DOMAINS = [(5, 352, "motor"), (353, 380, "boyun"), (381, 468, "CC1"), (469, 519, "FHA"), (520, 1805, "kuyruk")]

NMD_WINDOW_NT = 50            # the textbook 50-nucleotide rule
ENSEMBL_MAX_SPAN = 4_500_000  # the overlap endpoint refuses regions of 5 Mb and above
CNV_ENUMERATE_MAX = 30_000_000  # above this a CNV is a chromosome arm event; counting its genes is pointless

# ClinVar records that carry a real neurological phenotype, verified one by one from the VCV XML
VCV_OF_INTEREST = {"442594": "KIF13A-only deletion, affected child",
                   "2570672": "p.Lys1122Glu tail missense VUS with a phenotype"}


# --------------------------------------------------------------------------- small helpers


def table(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


def cached(subdir: str, pattern: str, base: Path = core.RAW) -> Path | None:
    """Newest cached file matching the pattern, never a .prov.json sidecar."""
    folder = base / subdir
    if not folder.exists():
        return None
    hits = sorted(p for p in folder.glob(pattern) if not p.name.endswith(".prov.json"))
    return hits[-1] if hits else None


def as_int(value) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


async def gnomad(query: str) -> dict:
    data = await core.post_json(GNOMAD_API, {"query": query}, timeout=300)
    if data.get("errors"):
        raise RuntimeError(f"gnomAD GraphQL error: {data['errors'][:2]}")
    return data["data"]


# --------------------------------------------------------------------------- 1. exon structure


async def exon_structure(enst: str, refresh: bool) -> dict:
    """MANE exon structure from Ensembl, turned into CDS coordinates and an NMD boundary."""
    hit = None if refresh else cached("ensembl", f"{enst}_structure_*.json")
    if hit:
        raw = json.loads(hit.read_text(encoding="utf-8"))
        print(f"  using cached {core.rel(hit)}")
    else:
        url = f"{ENSEMBL}/lookup/id/{enst}?expand=1;utr=1"
        raw = await core.get_json(url, headers={"Content-Type": "application/json"})
        rec = core.save_json("ensembl", f"{enst}_structure_{core.today()}.json", raw,
                             source="Ensembl REST lookup (expand, utr)", urls=[url],
                             tool="scripts/lof_axis.py", params={"transcript": enst})
        print(f"  downloaded {rec['path']}")

    tr = raw["Translation"]
    cds_lo, cds_hi = min(tr["start"], tr["end"]), max(tr["start"], tr["end"])
    exons, cum = [], 0
    for i, ex in enumerate(raw["Exon"], start=1):
        lo, hi = min(ex["start"], ex["end"]), max(ex["start"], ex["end"])
        coding = max(0, min(hi, cds_hi) - max(lo, cds_lo) + 1)
        exons.append({"exon_number": i, "exon_id": ex["id"], "genomic_start": lo, "genomic_end": hi,
                      "length": hi - lo + 1, "coding_length": coding,
                      "cds_start": cum + 1 if coding else None, "cds_end": cum + coding if coding else None})
        cum += coding
    n_exons = len(exons)
    # the last exon-exon junction, expressed in CDS coordinates: everything after it escapes NMD
    last_junction_cds = sum(e["coding_length"] for e in exons[:-1])
    coding_exons = [e for e in exons if e["coding_length"]]
    return {"transcript": f"{raw['id']}.{raw['version']}", "strand": raw["strand"], "n_exons": n_exons,
            "cds_length": cum, "protein_length": tr["length"], "exons": exons,
            "last_coding_exon": coding_exons[-1]["exon_number"],
            "last_junction_cds": last_junction_cds,
            "nmd_escape_from_cds": last_junction_cds - NMD_WINDOW_NT + 1}


CDS_POS = re.compile(r"^c\.(?:\*|-)?(\d+)")


def cds_position(hgvsc: str | None) -> int | None:
    """First CDS coordinate in an HGVS c. string. c.1234+1G>A and c.1235-1del both resolve."""
    if not hgvsc:
        return None
    m = CDS_POS.match(hgvsc.split(":")[-1].strip())
    return int(m.group(1)) if m else None


def locate(cds_pos: int | None, struct: dict) -> dict:
    """Exon number and NMD verdict for a CDS position."""
    if cds_pos is None:
        return {"exon_number": None, "in_last_exon": None, "nmd_escape_rule": None}
    exon_number = None
    for e in struct["exons"]:
        if e["cds_start"] and e["cds_start"] <= cds_pos <= e["cds_end"]:
            exon_number = e["exon_number"]
            break
    escape = cds_pos >= struct["nmd_escape_from_cds"]
    return {"exon_number": exon_number,
            "in_last_exon": bool(exon_number == struct["n_exons"]),
            "nmd_escape_rule": escape}


def domain_of(pos: int | None) -> str:
    if not pos:
        return ""
    for lo, hi, name in DOMAINS:
        if lo <= pos <= hi:
            return name
    return "N-uc" if pos < DOMAINS[0][0] else "C-uc"


# --------------------------------------------------------------------------- 2. gnomAD LoF


async def gnomad_lof_variants(symbol: str, enst: str, refresh: bool) -> list[dict]:
    """Every variant on the MANE transcript that carries a LoF annotation, with LOFTEE's verdict.

    The project's existing pull (variants.gnomad_variants) drops lof_filter and lof_flags, and those
    two fields are exactly what separates a real haploinsufficiency allele from a 3' truncation.
    """
    hit = None if refresh else cached("gnomad", f"{symbol}_{enst}_gnomad_r4_lofaxis_*.tsv")
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return table(hit)

    query = """
    { transcript(transcript_id: "%s", reference_genome: GRCh38) {
        transcript_id transcript_version
        variants(dataset: gnomad_r4) {
          variant_id chrom pos ref alt rsids hgvsc hgvsp consequence lof lof_filter lof_flags flags
          exome { ac an homozygote_count filters }
          genome { ac an homozygote_count filters }
          joint { ac an homozygote_count filters }
        } } }""" % enst
    tx = (await gnomad(query))["transcript"]
    rows = []
    for v in tx["variants"]:
        if not v.get("lof") and v.get("consequence") not in LOF_CONSEQUENCES:
            continue
        joint = v.get("joint") or {}
        ac, an = joint.get("ac"), joint.get("an")
        rows.append({
            "variant_id": v["variant_id"], "chrom": v["chrom"], "pos": v["pos"],
            "rsids": v.get("rsids"), "transcript": f"{tx['transcript_id']}.{tx['transcript_version']}",
            "consequence": v.get("consequence"), "hgvsc": v.get("hgvsc"), "hgvsp": v.get("hgvsp"),
            "loftee": v.get("lof"), "lof_filter": v.get("lof_filter"), "lof_flags": v.get("lof_flags"),
            "flags": v.get("flags"),
            "joint_ac": ac, "joint_an": an, "joint_af": (ac / an) if ac is not None and an else None,
            "joint_hom": joint.get("homozygote_count"), "joint_filters": joint.get("filters"),
            "exome_ac": (v.get("exome") or {}).get("ac"), "genome_ac": (v.get("genome") or {}).get("ac"),
        })
    rec = core.save_tsv("gnomad", f"{symbol}_{enst}_gnomad_r4_lofaxis_{core.today()}.tsv", rows,
                        source="gnomAD GraphQL (transcript query, LoF fields)", urls=[GNOMAD_API],
                        tool="scripts/lof_axis.py",
                        params={"gene": symbol, "transcript": enst, "fields": "lof,lof_filter,lof_flags"})
    print(f"  downloaded {rec['path']} ({len(rows)} LoF rows)")
    return table(cached("gnomad", f"{symbol}_{enst}_gnomad_r4_lofaxis_*.tsv"))


def annotate_lof(symbol: str, rows: list[dict], struct: dict) -> list[dict]:
    out = []
    for r in rows:
        cds = cds_position(r.get("hgvsc"))
        where = locate(cds, struct)
        ac, hom = as_int(r.get("joint_ac")), as_int(r.get("joint_hom"))
        filters = r.get("joint_filters") or ""
        confident = ac >= 1 and "AC0" not in filters
        m = re.search(r"p\.[A-Za-z]{3}(\d+)", r.get("hgvsp") or "")
        protein_pos = int(m.group(1)) if m else None
        loftee = r.get("loftee") or ""
        lof_filter = r.get("lof_filter") or ""
        lof_flags = r.get("lof_flags") or ""
        last50 = bool(protein_pos and protein_pos > struct["protein_length"] - 50)
        out.append({
            "gene": symbol, "variant_id": r["variant_id"], "consequence": r["consequence"],
            "hgvsc": r.get("hgvsc"), "hgvsp": r.get("hgvsp"),
            "cds_pos": cds, "exon_number": where["exon_number"], "n_exons": struct["n_exons"],
            "protein_pos": protein_pos,
            "protein_fraction": round(protein_pos / struct["protein_length"], 4) if protein_pos else None,
            "domain": domain_of(protein_pos) if symbol == "KIF13A" else "",
            "in_last_exon": where["in_last_exon"], "nmd_escape_rule": where["nmd_escape_rule"],
            "in_last_50_codons": last50,
            "loftee": loftee, "lof_filter": lof_filter, "lof_flags": lof_flags, "flags": r.get("flags"),
            "joint_ac": ac, "joint_an": as_int(r.get("joint_an")), "joint_af": r.get("joint_af"),
            "joint_hom": hom, "joint_filters": filters, "confident_carrier": confident,
            # the class the haploinsufficiency argument actually rests on
            "true_lof_candidate": bool(confident and loftee == "HC" and not lof_filter
                                       and where["nmd_escape_rule"] is False),
        })
    out.sort(key=lambda r: (r["cds_pos"] or 0))
    return out


def lof_summary(symbol: str, rows: list[dict], struct: dict, constraint: dict) -> dict:
    conf = [r for r in rows if r["confident_carrier"]]
    true_lof = [r for r in conf if r["true_lof_candidate"]]
    escape = [r for r in conf if r["nmd_escape_rule"]]
    return {
        "gene": symbol,
        "transcript": struct["transcript"],
        "protein_length": struct["protein_length"],
        "n_exons": struct["n_exons"],
        "cds_length": struct["cds_length"],
        "nmd_escape_from_cds": struct["nmd_escape_from_cds"],
        "nmd_escape_from_codon": (struct["nmd_escape_from_cds"] + 2) // 3,
        "lof_rows_total": len(rows),
        "lof_rows_confident": len(conf),
        "lof_rows_ac0_or_zero": len(rows) - len(conf),
        "confident_by_consequence": dict(Counter(r["consequence"] for r in conf).most_common()),
        "confident_loftee": dict(Counter(r["loftee"] or "(none)" for r in conf).most_common()),
        "confident_lof_filter": dict(Counter(r["lof_filter"] or "(none)" for r in conf).most_common()),
        "confident_lof_flags": dict(Counter(r["lof_flags"] or "(none)" for r in conf).most_common()),
        "confident_allele_count_sum": sum(r["joint_ac"] for r in conf),
        "confident_max_allele_count": max((r["joint_ac"] for r in conf), default=0),
        "confident_singletons": sum(1 for r in conf if r["joint_ac"] == 1),
        "confident_homozygotes": sum(r["joint_hom"] for r in conf),
        "confident_in_last_exon": sum(1 for r in conf if r["in_last_exon"]),
        "confident_nmd_escape": len(escape),
        "confident_in_last_50_codons": sum(1 for r in conf if r["in_last_50_codons"]),
        "confident_end_trunc": sum(1 for r in conf if r["lof_filter"] == "END_TRUNC"),
        "true_lof_candidates": len(true_lof),
        "true_lof_allele_count_sum": sum(r["joint_ac"] for r in true_lof),
        "true_lof_max_allele_count": max((r["joint_ac"] for r in true_lof), default=0),
        "true_lof_carrier_frequency_per_100k": round(
            1e5 * sum(r["joint_ac"] for r in true_lof)
            / max(1, max((r["joint_an"] for r in true_lof), default=1)), 2),
        "loeuf": constraint.get("oe_lof_upper"),
        "loeuf_percentile": constraint.get("oe_lof_percentile"),
        "pli": constraint.get("pli"),
        "obs_lof": constraint.get("obs_lof"), "exp_lof": constraint.get("exp_lof"),
        "mis_z": constraint.get("mis_z"), "oe_mis": constraint.get("oe_mis"),
    }


# --------------------------------------------------------------------------- 3. constraint panel


CONSTRAINT_QUERY = """{ gene(gene_symbol: "%s", reference_genome: GRCh38) {
  gene_id symbol
  mane_select_transcript { ensembl_id ensembl_version refseq_id }
  gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper oe_lof_percentile
    exp_mis obs_mis oe_mis oe_mis_lower oe_mis_upper exp_syn obs_syn oe_syn syn_z
    lof_z mis_z pli flags } } }"""


async def constraint_panel(symbols: list[str], refresh: bool) -> dict[str, dict]:
    hit = None if refresh else cached("gnomad", "constraint_panel_gnomad_r4_*.json")
    if hit:
        out = json.loads(hit.read_text(encoding="utf-8"))
        if all(s in out for s in symbols):
            print(f"  using cached {core.rel(hit)}")
            return out
    out = {}
    for s in symbols:
        gene = (await gnomad(CONSTRAINT_QUERY % s))["gene"]
        c = dict(gene.get("gnomad_constraint") or {})
        c["gene_id"] = gene["gene_id"]
        c["mane_select_transcript"] = (gene.get("mane_select_transcript") or {}).get("ensembl_id")
        out[s] = c
    rec = core.save_json("gnomad", f"constraint_panel_gnomad_r4_{core.today()}.json", out,
                         source="gnomAD GraphQL (gene_symbol constraint)", urls=[GNOMAD_API],
                         tool="scripts/lof_axis.py", params={"genes": symbols})
    print(f"  downloaded {rec['path']}")
    return out


async def clingen_dosage(refresh: bool) -> tuple[dict[str, dict], str]:
    """ClinGen dosage sensitivity curation, read live. Absence from this file is itself a finding."""
    hit = None if refresh else cached("clingen", "ClinGen_gene_curation_list_GRCh38_*.tsv")
    if hit:
        text = hit.read_text(encoding="utf-8")
        print(f"  using cached {core.rel(hit)}")
    else:
        raw = await core.get_bytes(CLINGEN_DOSAGE)
        rec = core.save_bytes("clingen", f"ClinGen_gene_curation_list_GRCh38_{core.today()}.tsv", raw,
                              source="ClinGen dosage sensitivity curation list", urls=[CLINGEN_DOSAGE],
                              tool="scripts/lof_axis.py", params={})
        print(f"  downloaded {rec['path']}")
        text = raw.decode("utf-8", "replace")

    lines = text.splitlines()
    stamp = next((l.lstrip("#").strip() for l in lines[:5] if re.match(r"#\d{1,2} \w+,\d{4}", l)), "UNVERIFIED")
    header = next(l.lstrip("#").split("\t") for l in lines if l.startswith("#Gene Symbol"))
    idx = {h: i for i, h in enumerate(header)}
    out = {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        out[f[0]] = {
            "hi_score": f[idx["Haploinsufficiency Score"]],
            "hi_description": f[idx["Haploinsufficiency Description"]],
            "hi_pmids": [f[idx[f"Haploinsufficiency PMID{i}"]] for i in range(1, 7)
                         if f[idx[f"Haploinsufficiency PMID{i}"]]],
            "ts_score": f[idx["Triplosensitivity Score"]],
            "last_evaluated": f[idx["Date Last Evaluated"]],
            "hi_disease_id": f[idx["Haploinsufficiency Disease ID"]],
        }
    return out, stamp


# --------------------------------------------------------------------------- 4. human observations


async def clinvar_cnvs(gene: str, refresh: bool) -> list[dict]:
    """Every non-SNV ClinVar record on the gene, live. The single-gene filter hides all of them."""
    hit = None if refresh else cached("clinvar", f"{gene}_clinvar_structural_*.json")
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return json.loads(hit.read_text(encoding="utf-8"))
    key = core.secret("NCBI_API_KEY")
    common = {"db": "clinvar", "retmode": "json", **({"api_key": key} if key else {})}
    term = f'{gene}[gene] NOT "single nucleotide variant"[Type of variation]'
    search = await core.get_json(f"{EUTILS}/esearch.fcgi", params={**common, "term": term, "retmax": "500"})
    ids = search["esearchresult"]["idlist"]
    records = []
    for block in core.chunks(ids, 100):
        summ = await core.get_json(f"{EUTILS}/esummary.fcgi", params={**common, "id": ",".join(block)})
        res = summ["result"]
        for uid in res["uids"]:
            r = res[uid]
            vset = (r.get("variation_set") or [{}])[0]
            locs = {loc.get("assembly_name"): loc for loc in vset.get("variation_loc", [])}
            germ = r.get("germline_classification") or {}
            records.append({
                "accession": r.get("accession"), "title": r.get("title"),
                "variant_type": vset.get("variant_type") or r.get("obj_type"),
                "classification": germ.get("description") or "(none)",
                "review_status": germ.get("review_status") or "",
                "last_evaluated": germ.get("last_evaluated") or "",
                "traits": ";".join(t.get("trait_name", "") for t in germ.get("trait_set") or []),
                "assembly": "GRCh38" if "GRCh38" in locs else ("GRCh37" if "GRCh37" in locs else ""),
                "chr": (locs.get("GRCh38") or locs.get("GRCh37") or {}).get("chr"),
                "start": as_int((locs.get("GRCh38") or locs.get("GRCh37") or {}).get("start")),
                "stop": as_int((locs.get("GRCh38") or locs.get("GRCh37") or {}).get("stop")),
                "grch37_start": as_int((locs.get("GRCh37") or {}).get("start")),
                "grch37_stop": as_int((locs.get("GRCh37") or {}).get("stop")),
                "n_scv": len((r.get("supporting_submissions") or {}).get("scv") or []),
            })
    rec = core.save_json("clinvar", f"{gene}_clinvar_structural_{core.today()}.json", records,
                         source="NCBI E-utilities clinvar esearch+esummary", urls=[f"{EUTILS}/esearch.fcgi"],
                         tool="scripts/lof_axis.py", params={"term": term})
    print(f"  downloaded {rec['path']} ({len(records)} records)")
    return records


async def protein_coding_genes(chrom: str, start: int, stop: int, assembly: str) -> list[str]:
    """Protein-coding genes inside an interval. Ensembl refuses >=5 Mb, so ask in windows."""
    base = ENSEMBL if assembly == "GRCh38" else ENSEMBL37
    names: set[str] = set()
    lo = start
    while lo <= stop:
        hi = min(stop, lo + ENSEMBL_MAX_SPAN - 1)
        feats = await core.get_json(f"{base}/overlap/region/human/{chrom}:{lo}-{hi}",
                                    params={"feature": "gene", "content-type": "application/json"})
        for g in feats:
            if g.get("biotype") == "protein_coding":
                names.add(g.get("external_name") or g.get("gene_id") or g.get("id"))
        lo = hi + 1
    return sorted(names)


async def cnv_gene_content(gene: str, records: list[dict], refresh: bool) -> list[dict]:
    hit = None if refresh else cached("lof", f"{gene}_cnv_gene_content_*.tsv", base=core.INTERIM)
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return table(hit)
    rows = []
    for r in records:
        if not r["chr"] or not r["start"] or r["stop"] - r["start"] < 1000:
            continue  # small indels are handled as sequence variants, not CNVs
        span = r["stop"] - r["start"] + 1
        if span > CNV_ENUMERATE_MAX:
            rows.append({**r, "span_bp": span, "n_protein_coding_genes": "SAYILMADI",
                         "contains_target": "",
                         "protein_coding_genes": f"kromozom kolu olcegi ({span / 1e6:.0f} Mb), sayilmadi"})
            continue
        genes = await protein_coding_genes(r["chr"], r["start"], r["stop"], r["assembly"])
        rows.append({**r, "span_bp": span,
                     "n_protein_coding_genes": len(genes),
                     "contains_target": gene in genes,
                     "protein_coding_genes": ";".join(genes) if len(genes) <= 40 else
                     ";".join(genes[:40]) + f";...(+{len(genes) - 40})"})
    rows.sort(key=lambda r: r["span_bp"])
    rec = core.save_tsv("lof", f"{gene}_cnv_gene_content_{core.today()}.tsv", rows, base=core.INTERIM,
                        source="derived: ClinVar esummary + Ensembl REST overlap",
                        urls=[f"{ENSEMBL}/overlap/region/human", f"{ENSEMBL37}/overlap/region/human"],
                        tool="scripts/lof_axis.py", params={"gene": gene})
    print(f"  wrote {rec['path']} ({len(rows)} CNVs)")
    return rows


def clinvar_sequence_lof(gene: str) -> list[dict]:
    """LoF sequence variants out of the already-pulled single-gene ClinVar table."""
    path = cached("clinvar", f"{gene}_clinvar_*.tsv")
    if not path:
        return []
    keep = []
    for r in table(path):
        mc = (r.get("molecular_consequence") or "").lower()
        if any(k in mc for k in ("nonsense", "frameshift", "splice acceptor", "splice donor")):
            keep.append({"gene": gene, "accession": r["accession"], "title": r["title"],
                         "molecular_consequence": r["molecular_consequence"],
                         "classification": r["classification"], "review_status": r["review_status"],
                         "stars": r["stars"], "traits": r["traits"], "hgvsc": r.get("hgvsc"),
                         "hgvsp": r.get("hgvsp"), "protein_pos": r.get("protein_pos"),
                         "last_evaluated": r.get("last_evaluated")})
    return keep


async def vcv_details(refresh: bool) -> list[dict]:
    """Read the phenotype out of the VCV XML itself. esummary only says "See cases"."""
    hit = None if refresh else cached("clinvar", "KIF13A_vcv_phenotyped_*.json")
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return json.loads(hit.read_text(encoding="utf-8"))
    key = core.secret("NCBI_API_KEY")
    out = []
    for uid, note in VCV_OF_INTEREST.items():
        raw = await core.get_bytes(f"{EUTILS}/efetch.fcgi",
                                   params={"db": "clinvar", "rettype": "vcv", "is_variationid": "true",
                                           "id": uid, **({"api_key": key} if key else {})})
        core.save_bytes("clinvar/vcv_xml", f"VCV{int(uid):09d}_{core.today()}.xml", raw,
                        source="NCBI E-utilities clinvar efetch (rettype=vcv)",
                        urls=[f"{EUTILS}/efetch.fcgi"], tool="scripts/lof_axis.py", params={"id": uid})
        x = raw.decode("utf-8", "replace")
        traits = sorted(set(re.findall(r'<ElementValue Type="Preferred">([^<]*)</ElementValue>', x)))
        out.append({
            "variation_id": uid, "accession": f"VCV{int(uid):09d}", "note": note,
            "title": (re.search(r"<VariationName>([^<]*)</VariationName>", x)
                      or re.search(r'VariationName="([^"]*)"', x) or [None, ""])[1],
            "traits": ";".join(t for t in traits if t != "See cases"),
            "affected_status": ";".join(sorted(set(re.findall(r'AffectedStatus="([^"]*)"', x)))),
            "collection_method": ";".join(sorted(set(re.findall(r"<MethodType>([^<]*)</MethodType>", x)))),
            "origin": ";".join(sorted(set(re.findall(r"<Origin>([^<]*)</Origin>", x)))),
            "submitters": ";".join(sorted(set(re.findall(r'SubmitterName="([^"]*)"', x)))),
            "cited_pmids": ";".join(sorted(set(re.findall(r'<ID Source="PubMed">(\d+)</ID>', x)))),
        })
    rec = core.save_json("clinvar", f"KIF13A_vcv_phenotyped_{core.today()}.json", out,
                         source="NCBI E-utilities clinvar efetch (rettype=vcv)", urls=[f"{EUTILS}/efetch.fcgi"],
                         tool="scripts/lof_axis.py", params={"ids": list(VCV_OF_INTEREST)})
    print(f"  downloaded {rec['path']}")
    return out


async def decipher_gene_metrics(gene: str, refresh: bool) -> dict:
    """DECIPHER publishes gene-level dosage metrics in the server-rendered page, no account needed.

    The JSON API is key-gated (HTTP 403), so the patient CNV browser cannot be pulled into a table;
    only these gene-level numbers are machine-readable today.
    """
    hit = None if refresh else cached("decipher", f"{gene}_gene_ssr_*.json")
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return json.loads(hit.read_text(encoding="utf-8"))
    url = DECIPHER_GENE % gene
    html = (await core.get_bytes(url, headers={"User-Agent": BROWSER_UA, "Accept": "text/html"})).decode("utf-8", "replace")
    m = re.search(r"window\.__NUXT__=(.*?)</script>", html, re.S)
    payload = m.group(1) if m else ""
    fields = {}
    for key in ("hi_score", "hi_index", "s_het", "p_li", "loeuf", "missense_z_score", "p_hi", "p_ts",
                "obs_missense", "exp_missense", "pdn", "pgof", "plof", "lof_pctl", "dn_pctl", "gof_pctl",
                "rank_max", "public_patients", "public_snv_count"):
        mm = re.search(rf'\b{key}:\s*"?([\-\w.]+)"?', payload)
        fields[key] = mm.group(1) if mm else None
    fields["clingen_dosage_record"] = "absent" if re.search(r"ClingenDSGene:\s*a\b", payload) else "present_or_unknown"
    fields["source_url"] = url
    fields["retrieved"] = core.today()
    rec = core.save_json("decipher", f"{gene}_gene_ssr_{core.today()}.json", fields,
                         source="DECIPHER gene page, server-rendered payload (no account)", urls=[url],
                         tool="scripts/lof_axis.py", params={"gene": gene})
    print(f"  downloaded {rec['path']}")
    return fields


async def europepmc(queries: dict[str, str], refresh: bool) -> dict:
    hit = None if refresh else cached("literature", "lof_axis_epmc_*.json")
    if hit:
        out = json.loads(hit.read_text(encoding="utf-8"))
        if set(queries) <= set(out):
            print(f"  using cached {core.rel(hit)}")
            return out
    out = {}
    for label, q in queries.items():
        d = await core.get_json(EPMC, params={"query": q, "format": "json", "pageSize": "25",
                                              "resultType": "core" if label == "lam2021" else "lite"})
        out[label] = {"query": q, "hitCount": d["hitCount"], "results": [
            {"pmid": r.get("pmid"), "year": r.get("pubYear"), "title": r.get("title"),
             "journal": ((r.get("journalInfo") or {}).get("journal") or {}).get("title"),
             "doi": r.get("doi"), "isOpenAccess": r.get("isOpenAccess"),
             "abstract": r.get("abstractText") if label == "lam2021" else None}
            for r in d["resultList"]["result"][:10]]}
    rec = core.save_json("literature", f"lof_axis_epmc_{core.today()}.json", out,
                         source="Europe PMC REST search", urls=[EPMC], tool="scripts/lof_axis.py",
                         params={"queries": list(queries)})
    print(f"  downloaded {rec['path']}")
    return out


# --------------------------------------------------------------------------- 5. mouse


async def impc(gene: str, refresh: bool) -> dict:
    """IMPC: the statistically called phenotypes, and what was tested in heterozygotes."""
    hit = None if refresh else cached("impc", f"{gene}_impc_*.json")
    if hit:
        print(f"  using cached {core.rel(hit)}")
        return json.loads(hit.read_text(encoding="utf-8"))
    mouse = gene[0] + gene[1:].lower()
    calls = await core.get_json(f"{IMPC_SOLR}/genotype-phenotype/select",
                                params={"q": f"marker_symbol:{mouse}", "rows": "500", "wt": "json"})
    facets = await core.get_json(f"{IMPC_SOLR}/statistical-result/select",
                                 params={"q": f"marker_symbol:{mouse}", "rows": "0", "wt": "json",
                                         "facet": "true", "facet.limit": "60",
                                         "facet.field": ["zygosity", "significant", "status", "procedure_name"]})
    het = await core.get_json(f"{IMPC_SOLR}/statistical-result/select",
                              params={"q": f"marker_symbol:{mouse} AND zygosity:heterozygote",
                                      "rows": "100", "wt": "json"})
    via = await core.get_json(f"{IMPC_SOLR}/experiment/select",
                              params={"q": f'gene_symbol:{mouse} AND procedure_name:"Viability Primary Screen"',
                                      "rows": "60", "wt": "json",
                                      "fl": "parameter_name,data_point,zygosity,category"})
    gene_doc = await core.get_json(f"{IMPC_SOLR}/gene/select",
                                   params={"q": f"marker_symbol:{mouse}", "rows": "1", "wt": "json"})

    def ff(name):
        vals = facets["facet_counts"]["facet_fields"].get(name, [])
        return {vals[i]: vals[i + 1] for i in range(0, len(vals), 2) if vals[i + 1]}

    g = (gene_doc["response"]["docs"] or [{}])[0]
    out = {
        "marker_symbol": mouse,
        "mgi_accession_id": g.get("mgi_accession_id"),
        "alleles": g.get("allele_name"),
        "significant_top_level_mp_terms": g.get("significant_top_level_mp_terms"),
        "not_significant_top_level_mp_terms": g.get("not_significant_top_level_mp_terms"),
        "n_significant_calls": calls["response"]["numFound"],
        "statistical_results_by_zygosity": ff("zygosity"),
        "statistical_results_by_significance": ff("significant"),
        "statistical_results_by_status": ff("status"),
        "procedures_tested": ff("procedure_name"),
        "calls": [{k: d.get(k) for k in ("mp_term_id", "mp_term_name", "parameter_name", "procedure_name",
                                         "zygosity", "sex", "p_value", "effect_size", "allele_symbol",
                                         "phenotyping_center")}
                  for d in calls["response"]["docs"]],
        "heterozygote_tests": [{k: d.get(k) for k in ("parameter_name", "procedure_name", "p_value",
                                                      "significant", "status", "allele_symbol")}
                               for d in het["response"]["docs"]],
        "viability_counts": {d.get("parameter_name"): d.get("data_point")
                             for d in via["response"]["docs"] if d.get("data_point") is not None},
        "retrieved": core.today(),
    }
    rec = core.save_json("impc", f"{gene}_impc_{core.today()}.json", out,
                         source="IMPC Solr (genotype-phenotype, statistical-result, experiment, gene)",
                         urls=[f"{IMPC_SOLR}/genotype-phenotype/select"], tool="scripts/lof_axis.py",
                         params={"marker_symbol": mouse})
    print(f"  downloaded {rec['path']}")
    return out


# --------------------------------------------------------------------------- 6. figure


FIG_W, FIG_H = 1180, 620
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
C_TRUE, C_ESCAPE, C_LC, C_AC0 = "#e34948", "#7fb2e5", "#eb9a34", "#c9c8c0"


def build_figure(per_gene: dict, summaries: dict, panel_rows: list[dict], path: Path) -> None:
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{FIG_W}" height="{FIG_H}" '
         f'viewBox="0 0 {FIG_W} {FIG_H}" font-family="Arial, Helvetica, sans-serif">',
         f'<rect width="{FIG_W}" height="{FIG_H}" fill="#ffffff"/>',
         f'<text x="40" y="30" font-size="14" fill="{INK}">KIF13A fonksiyon kaybi ekseni: '
         f'gnomAD v4.1 LoF varyantlari, NMD kacisi ve kisitlama baglami</text>']

    # ---- panel A: one track per gene, LoF positions along the protein
    L, R = 70, 40
    trackw = FIG_W - L - R
    top = 70
    o.append(f'<text x="{L}" y="{top - 14}" font-size="12" fill="{INK}">A. MANE proteini boyunca '
             f'guvenilir tasiyicisi olan LoF varyantlari (nokta capi allel sayisiyla artar)</text>')
    for gi, gene in enumerate(GENES):
        rows = [r for r in per_gene[gene] if r["confident_carrier"]]
        s = summaries[gene]
        length = s["protein_length"]
        y0 = top + gi * 86
        esc_codon = s["nmd_escape_from_codon"]

        def gx(p):
            return L + trackw * (p - 1) / (length - 1)

        o.append(f'<rect x="{L}" y="{y0 + 22}" width="{trackw}" height="14" fill="#f2f1ec" stroke="{GRID}"/>')
        o.append(f'<rect x="{gx(esc_codon):.1f}" y="{y0 + 22}" width="{trackw - (gx(esc_codon) - L):.1f}" '
                 f'height="14" fill="{C_ESCAPE}" opacity="0.35"/>')
        if gene == "KIF13A":
            for lo, hi, name in DOMAINS[:2]:
                o.append(f'<rect x="{gx(lo):.1f}" y="{y0 + 22}" width="{gx(hi) - gx(lo):.1f}" height="14" '
                         f'fill="#cfe3f7"/>')
            o.append(f'<text x="{gx(178):.1f}" y="{y0 + 33}" font-size="9" fill="{MUTED}" '
                     f'text-anchor="middle">motor</text>')
        o.append(f'<text x="{L}" y="{y0 + 16}" font-size="11" fill="{INK}">{gene}  '
                 f'({length} aa, LOEUF {s["loeuf"]:.3f}, {s["lof_rows_confident"]} guvenilir LoF, '
                 f'{s["true_lof_candidates"]} NMD-duyarli)</text>')
        for r in rows:
            pos = r["protein_pos"] or (((r["cds_pos"] or 1) + 2) // 3)
            pos = min(max(pos, 1), length)
            ac = max(1, r["joint_ac"])
            rad = 2.2 + min(4.4, 1.35 * (ac ** 0.5) / 2)
            colour = C_TRUE if r["true_lof_candidate"] else (C_ESCAPE if r["nmd_escape_rule"] else C_LC)
            o.append(f'<circle cx="{gx(pos):.1f}" cy="{y0 + 29}" r="{rad:.1f}" fill="{colour}" '
                     f'fill-opacity="0.75" stroke="{colour}" stroke-width="0.6"/>')
        for tick in (1, length):
            o.append(f'<text x="{gx(tick):.1f}" y="{y0 + 50}" font-size="9" fill="{MUTED}" '
                     f'text-anchor="{"start" if tick == 1 else "end"}">{tick}</text>')
        o.append(f'<text x="{gx(esc_codon):.1f}" y="{y0 + 50}" font-size="9" fill="#2f6ea8" '
                 f'text-anchor="middle">NMD kacisi &gt;= {esc_codon}</text>')

    # ---- panel B: stacked bars of confident LoF composition
    by, bh, bx0 = top + 3 * 86 + 34, 150, L
    bw = 330
    o.append(f'<text x="{bx0}" y="{by - 10}" font-size="12" fill="{INK}">B. Guvenilir LoF bilesimi</text>')
    maxn = max(summaries[g]["lof_rows_confident"] for g in GENES)
    slot = bw / len(GENES)
    for gi, gene in enumerate(GENES):
        s = summaries[gene]
        parts = [("NMD-duyarli (HC, filtresiz)", s["true_lof_candidates"], C_TRUE),
                 ("NMD kacisi", s["confident_nmd_escape"], C_ESCAPE),
                 ("diger (LC / LOFTEE filtresi)", s["lof_rows_confident"] - s["true_lof_candidates"]
                  - s["confident_nmd_escape"], C_LC)]
        x = bx0 + gi * slot + slot * 0.22
        w = slot * 0.56
        yy = by + bh
        for _, n, col in parts:
            h = bh * n / maxn
            yy -= h
            o.append(f'<rect x="{x:.1f}" y="{yy:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{col}"/>')
        o.append(f'<text x="{x + w / 2:.1f}" y="{by + bh + 14}" font-size="10" fill="{INK}" '
                 f'text-anchor="middle">{gene}</text>')
        o.append(f'<text x="{x + w / 2:.1f}" y="{yy - 5:.1f}" font-size="10" fill="{MUTED}" '
                 f'text-anchor="middle">{s["lof_rows_confident"]}</text>')
    ly = by + bh + 34
    for label, col in (("NMD-duyarli (HC, filtresiz)", C_TRUE), ("NMD kacisi", C_ESCAPE),
                       ("LC / LOFTEE filtresi", C_LC)):
        o.append(f'<rect x="{bx0}" y="{ly - 8}" width="9" height="9" fill="{col}"/>')
        o.append(f'<text x="{bx0 + 14}" y="{ly}" font-size="10" fill="{MUTED}">{label}</text>')
        ly += 14

    # ---- panel C: LOEUF of the ClinGen haploinsufficiency panel vs our three genes
    cx0 = bx0 + bw + 90
    cw = FIG_W - cx0 - R
    o.append(f'<text x="{cx0}" y="{by - 10}" font-size="12" fill="{INK}">C. LOEUF: ClinGen '
             f'haployetmezlik paneli (HI=3) ve bizim genlerimiz</text>')
    values = [(r["gene"], float(r["loeuf"]), r["group"]) for r in panel_rows if r["loeuf"] not in (None, "")]
    lo_v, hi_v = 0.0, max(v for _, v, _ in values) * 1.15

    def cxp(v):
        return cx0 + cw * (v - lo_v) / (hi_v - lo_v)

    axis_y = by + bh + 6
    o.append(f'<line x1="{cx0}" y1="{axis_y}" x2="{cx0 + cw}" y2="{axis_y}" stroke="{MUTED}"/>')
    for t in (0.0, 0.2, 0.4, 0.6, 0.8):
        if t > hi_v:
            continue
        o.append(f'<line x1="{cxp(t):.1f}" y1="{by}" x2="{cxp(t):.1f}" y2="{axis_y}" stroke="{GRID}"/>')
        o.append(f'<text x="{cxp(t):.1f}" y="{axis_y + 14}" font-size="10" fill="{MUTED}" '
                 f'text-anchor="middle">{t:.1f}</text>')
    step = bh / (len(values) + 1)
    for i, (g, v, group) in enumerate(sorted(values, key=lambda x: x[1])):
        yy = by + step * (i + 1)
        col = C_TRUE if g == "KIF13A" else ("#2a78d6" if group == "target" else MUTED)
        o.append(f'<circle cx="{cxp(v):.1f}" cy="{yy:.1f}" r="{4.5 if group == "target" else 3}" fill="{col}"/>')
        o.append(f'<text x="{cxp(v) + 8:.1f}" y="{yy + 3.5:.1f}" font-size="10" '
                 f'fill="{INK if group == "target" else MUTED}">{g} {v:.2f}</text>')
    o.append(f'<text x="{cx0}" y="{axis_y + 30}" font-size="10" fill="{MUTED}">LOEUF (kucuk = LoF '
             f'intoleransi daha guclu)</text>')
    o.append(f'<text x="40" y="{FIG_H - 12}" font-size="10" fill="{MUTED}">Kaynak: gnomAD v4.1 '
             f'(GraphQL, transcript sorgusu), Ensembl REST MANE ekzon yapisi, ClinGen doz duyarliligi '
             f'kurasyon listesi. Olcum tarihi {core.today()}.</text>')
    o.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(o), encoding="utf-8")
    print(f"  wrote {core.rel(path)} ({path.stat().st_size} bytes)")


# --------------------------------------------------------------------------- main


async def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes", default=",".join(GENES), help="genes for the LoF breakdown")
    ap.add_argument("--panel", default=",".join(PANEL), help="ClinGen haploinsufficiency comparison panel")
    ap.add_argument("--refresh", action="store_true", help="ignore every cached pull")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    core.load_dotenv()
    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    panel = [g.strip() for g in args.panel.split(",") if g.strip()]
    stamp = core.today()

    print("== gnomAD constraint")
    con = await constraint_panel(sorted(set(panel + genes)), args.refresh)

    print("\n== exon structure and gnomAD LoF variants")
    per_gene, summaries, structs = {}, {}, {}
    for gene in genes:
        t = core.target(gene)
        enst = core.unversioned(t["mane_transcript"])
        print(f"{gene} ({t['mane_transcript']}, {t['protein_length']} aa)")
        struct = await exon_structure(enst, args.refresh)
        if struct["protein_length"] != t["protein_length"]:
            raise RuntimeError(f"{gene}: Ensembl says {struct['protein_length']} aa, "
                               f"config/targets.yaml says {t['protein_length']}")
        raw = await gnomad_lof_variants(gene, enst, args.refresh)
        rows = annotate_lof(gene, raw, struct)
        per_gene[gene], structs[gene] = rows, struct
        summaries[gene] = lof_summary(gene, rows, struct, con.get(gene, {}))
        s = summaries[gene]
        print(f"  {s['lof_rows_total']} LoF rows, {s['lof_rows_confident']} with a confident carrier, "
              f"{s['confident_nmd_escape']} escape NMD, {s['true_lof_candidates']} NMD-competent HC")

    all_rows = [r for g in genes for r in per_gene[g]]
    core.save_tsv("lof", f"gnomad_lof_variants_{stamp}.tsv", all_rows, base=core.INTERIM,
                  source="derived: gnomAD v4.1 LoF pull + Ensembl MANE exon structure",
                  urls=[GNOMAD_API, f"{ENSEMBL}/lookup/id"], tool="scripts/lof_axis.py",
                  params={"genes": genes, "nmd_window_nt": NMD_WINDOW_NT})
    core.save_tsv("lof", f"gnomad_lof_summary_{stamp}.tsv",
                  [{k: v for k, v in summaries[g].items()} for g in genes], base=core.INTERIM,
                  source="derived", urls=[GNOMAD_API], tool="scripts/lof_axis.py", params={"genes": genes})

    # per-exon table for KIF13A, so the reader can see where the surviving LoF sits
    exon_rows = []
    for e in structs["KIF13A"]["exons"]:
        hits = [r for r in per_gene["KIF13A"]
                if r["confident_carrier"] and r["exon_number"] == e["exon_number"]]
        exon_rows.append({**e, "confident_lof": len(hits),
                          "true_lof_candidates": sum(1 for r in hits if r["true_lof_candidate"]),
                          "allele_count_sum": sum(r["joint_ac"] for r in hits),
                          "codon_start": ((e["cds_start"] or 0) + 2) // 3 if e["cds_start"] else None,
                          "codon_end": ((e["cds_end"] or 0) + 2) // 3 if e["cds_end"] else None})
    core.save_tsv("lof", f"KIF13A_exon_lof_{stamp}.tsv", exon_rows, base=core.INTERIM,
                  source="derived: Ensembl MANE exon structure + gnomAD LoF",
                  urls=[f"{ENSEMBL}/lookup/id"], tool="scripts/lof_axis.py", params={"gene": "KIF13A"})

    print("\n== ClinGen dosage sensitivity")
    dosage, dosage_date = await clingen_dosage(args.refresh)
    print(f"  curation list dated {dosage_date}; KIF13A "
          f"{'present' if 'KIF13A' in dosage else 'ABSENT (never curated)'}")
    panel_rows = []
    for gene in panel + genes:
        c = con.get(gene, {})
        d = dosage.get(gene)
        panel_rows.append({
            "gene": gene, "group": "target" if gene in genes else "clingen_hi3",
            "clingen_hi_score": d["hi_score"] if d else "not curated",
            "clingen_hi_description": d["hi_description"] if d else "gen ClinGen doz listesinde yok",
            "clingen_last_evaluated": d["last_evaluated"] if d else "",
            "loeuf": c.get("oe_lof_upper"), "loeuf_percentile": c.get("oe_lof_percentile"),
            "pli": c.get("pli"), "lof_z": c.get("lof_z"),
            "obs_lof": c.get("obs_lof"), "exp_lof": c.get("exp_lof"),
            "mis_z": c.get("mis_z"), "oe_mis": c.get("oe_mis"), "syn_z": c.get("syn_z"),
            "flags": ";".join(c.get("flags") or []), "mane_transcript": c.get("mane_select_transcript"),
        })
    core.save_tsv("lof", f"constraint_panel_{stamp}.tsv", panel_rows, base=core.INTERIM,
                  source=f"gnomAD v4.1 constraint + ClinGen dosage curation list ({dosage_date})",
                  urls=[GNOMAD_API, CLINGEN_DOSAGE], tool="scripts/lof_axis.py",
                  params={"panel": panel, "targets": genes})
    ranked = sorted((r for r in panel_rows if r["loeuf"] is not None), key=lambda r: r["loeuf"])
    below = sum(1 for r in ranked if r["group"] == "clingen_hi3" and r["loeuf"] < (con["KIF13A"]["oe_lof_upper"]))
    n_hi3 = sum(1 for r in ranked if r["group"] == "clingen_hi3")
    print(f"  KIF13A LOEUF {con['KIF13A']['oe_lof_upper']:.3f}; "
          f"{below}/{n_hi3} panel genes are more LoF-constrained")

    print("\n== human loss-of-function observations")
    seq_lof = clinvar_sequence_lof("KIF13A")
    print(f"  ClinVar single-gene LoF sequence variants: {len(seq_lof)}")
    cnv_records = await clinvar_cnvs("KIF13A", args.refresh)
    cnvs = await cnv_gene_content("KIF13A", cnv_records, args.refresh)
    solo = [c for c in cnvs if as_int(c["n_protein_coding_genes"]) == 1 and str(c["contains_target"]) == "True"]
    print(f"  ClinVar non-SNV records: {len(cnv_records)}; CNVs >= 1 kb: {len(cnvs)}; "
          f"KIF13A-only: {len(solo)}")
    core.save_tsv("lof", f"KIF13A_clinvar_lof_sequence_{stamp}.tsv", seq_lof, base=core.INTERIM,
                  source="derived from data/raw/clinvar KIF13A single-gene pull",
                  urls=[], tool="scripts/lof_axis.py", params={"gene": "KIF13A"})

    phenotyped = await vcv_details(args.refresh)
    for p in phenotyped:
        print(f"  {p['accession']}: {p['traits']} | affected {p['affected_status']} | "
              f"{p['collection_method']} | {p['submitters'][:60]}")
    core.save_tsv("lof", f"KIF13A_phenotyped_records_{stamp}.tsv", phenotyped, base=core.INTERIM,
                  source="ClinVar VCV XML (efetch)", urls=[f"{EUTILS}/efetch.fcgi"],
                  tool="scripts/lof_axis.py", params={"ids": list(VCV_OF_INTEREST)})

    dec = await decipher_gene_metrics("KIF13A", args.refresh)
    print(f"  DECIPHER pHaplo {dec.get('p_hi')}, pTriplo {dec.get('p_ts')}, "
          f"s_het {dec.get('s_het')}, mechanism rank {dec.get('rank_max')}, "
          f"ClinGen dosage record {dec['clingen_dosage_record']}")

    lit = await europepmc({
        "lam2021": "EXT_ID:33526817",
        "kif13a_lof": 'KIF13A AND (haploinsufficiency OR "loss of function" OR frameshift OR nonsense)',
        "kif13a_deletion_case": 'KIF13A AND ("6p22.3" OR microdeletion) AND (patient OR proband OR case)',
        "kif13a_mouse": '"Kif13a" AND (knockout OR "null mice" OR heterozygous)',
    }, args.refresh)
    lam = lit["lam2021"]["results"][0] if lit["lam2021"]["results"] else {}
    print(f"  Lam 2021 verified: PMID {lam.get('pmid')} {lam.get('journal')} {lam.get('year')} "
          f"open access {lam.get('isOpenAccess')}")

    print("\n== mouse")
    mouse = await impc("KIF13A", args.refresh)
    het_sig = [h for h in mouse["heterozygote_tests"] if h.get("significant")]
    print(f"  {mouse['n_significant_calls']} significant calls, all zygosities "
          f"{Counter(c['zygosity'] for c in mouse['calls'])}")
    print(f"  heterozygote statistical results: {len(mouse['heterozygote_tests'])} "
          f"({len(het_sig)} significant)")
    mouse_kif13b = await impc("KIF13B", args.refresh)
    mouse_kif1a = await impc("KIF1A", args.refresh)

    verdict = {
        "measured": stamp,
        "gnomad_lof": {g: summaries[g] for g in genes},
        "constraint_panel_rank": [{"gene": r["gene"], "loeuf": r["loeuf"], "hi": r["clingen_hi_score"]}
                                  for r in ranked],
        "clingen_dosage_list_date": dosage_date,
        "kif13a_in_clingen_dosage": "KIF13A" in dosage,
        "clinvar_sequence_lof": len(seq_lof),
        "clinvar_cnv_records": len(cnv_records),
        "clinvar_cnv_kif13a_only": [c["accession"] for c in solo],
        "clinvar_phenotyped_records": phenotyped,
        "decipher": dec,
        "literature": {k: {"query": v["query"], "hitCount": v["hitCount"]} for k, v in lit.items()},
        "impc": {g: {"n_significant_calls": m["n_significant_calls"],
                     "zygosity_of_calls": dict(Counter(c["zygosity"] for c in m["calls"])),
                     "heterozygote_tests": len(m["heterozygote_tests"]),
                     "heterozygote_significant": sum(1 for h in m["heterozygote_tests"] if h.get("significant")),
                     "viability_counts": m["viability_counts"]}
                 for g, m in (("KIF13A", mouse), ("KIF13B", mouse_kif13b), ("KIF1A", mouse_kif1a))},
    }
    rec = core.save_json("lof", f"lof_axis_summary_{stamp}.json", verdict, base=core.INTERIM,
                         source="derived", urls=[], tool="scripts/lof_axis.py", params={"genes": genes})
    print(f"\n  wrote {rec['path']}")

    if not args.no_figure:
        print("\n== figure")
        build_figure(per_gene, summaries, panel_rows, ROOT / "figures" / "KIF13A_lof_axis.svg")
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
