"""Out-of-sample validation of the paralog-calibrated prioritisation method on NSD1 (teacher) / NSD3 (student).

The multi-pair screen (`scripts/paralog_pairs_scan.py`, `docs/olcekleme_fizibilite.md`) claimed that cheap
metrics predict in advance which paralog pairs transfer. NSD1/NSD3 was its best neurodevelopmental
candidate. This script runs the same pipeline that was run on KIF1A/KIF13A, blind to the prediction, so the
measured transfer signal can be compared against the predicted one.

Stages (each cached under data/interim/nsd_pair/; rerun one with --stage):
  identity  MANE transcript, RefSeq, UniProt accession and exact length for both genes, live, plus the
            isoform check (UniProt canonical sequence vs Ensembl MANE translation, residue by residue)
  domains   UniProt features + InterPro entries -> domain table, per-domain paralog identity
  universe  gnomAD v4.1 transcript query (never the gene query: trap 4) and ClinVar for both genes
  map       NSD1 -> NSD3 residue map, same aligner parameters as atlas_mcp.structure.paralog_map
  curate    ClinVar VCV XML for NSD1 -> positives with review status, PP3-dependent hold-out, and
            gnomAD negatives drawn from the SAME structural region as the positives (trap 10)
  score     one Ensembl VEP run over positives and negatives together (trap 11), MANE-indexed dbNSFP
  measure   AUROC with residue-grouped bootstrap (trap 13), sens@90%spec, thresholds, mandatory baselines
  constraint observed/expected synonymous and missense per region, the calibration check
  lift      the measured counterpart of the screen's transfer_lift, same definition
  transfer  calibrated threshold + paralog evidence on NSD3 gnomAD missense -> tiers A/B/C/D

    python -u scripts/nsd_pair_calibration.py --stage all
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "mcp" / "atlas_server" / "src"))

from atlas_mcp import core, variants  # noqa: E402

from Bio import Align  # noqa: E402
from Bio.Align import substitution_matrices  # noqa: E402

core.load_dotenv()

TEACHER, STUDENT = "NSD1", "NSD3"
OUT_SUB = "nsd_pair"                       # data/interim/nsd_pair/
ENSEMBL = "https://rest.ensembl.org"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
INTERPRO = "https://www.ebi.ac.uk/interpro/api"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
VEP_REGION = "https://rest.ensembl.org/vep/human/region"

RNG = np.random.default_rng(20260917)
N_BOOT = 2000
POPULATION_MIN_AC = 3     # main negative set; AC >= 2 kept as a sensitivity set
AF_DOMINANT_MAX = 1e-5    # a fully penetrant dominant de novo allele should be rarer than this

P_LP = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
B_LB = {"benign", "likely benign", "benign/likely benign"}

# model key in the VEP-derived row -> (label, +1 if higher means more damaging else -1)
MODELS = {
    "am_pathogenicity": ("AlphaMissense", 1), "cadd_phred": ("CADD v1.7", 1),
    "esm1b": ("ESM1b (dbNSFP)", -1), "primateai": ("PrimateAI", 1),
    "revel": ("REVEL", 1), "metarnn": ("MetaRNN", 1),
}
LABEL_TRAINED = {"revel", "metarnn"}   # trained on clinical labels: reference rows only (trap 12)

TEXT_FLAGS = {
    "in_silico": r"\bPP3\b|\bBP4\b|in[- ]silico|computational|REVEL|CADD|SIFT|PolyPhen|AlphaMissense|algorithm|predict(ed|ion) (tool|to be)",
    "functional": r"\bPS3\b|\bBS3\b|functional (stud|assay|evidence|analys)|in vitro|experimental stud|enzymatic activ|methyltransferase activ",
    "frequency": r"\bBA1\b|\bBS1\b|\bBS2\b|\bPM2\b|gnomAD|population database|allele frequency|ExAC|1000 Genomes",
    "de_novo": r"\bPS2\b|\bPM6\b|de novo",
    "same_residue": r"\bPS1\b|\bPM5\b|same amino acid change|different (amino acid|missense) change",
    "segregation": r"\bPP1\b|\bBS4\b|segregat",
}


def log(*parts: Any) -> None:
    """Provenance and progress go to stderr so stdout stays a clean report."""
    print(*parts, file=sys.stderr, flush=True)


def out_dir(args: argparse.Namespace) -> Path:
    return Path(args.out).resolve() if args.out else (core.INTERIM / OUT_SUB)


def save_tsv(args: argparse.Namespace, name: str, rows: list[dict], **prov: Any) -> dict:
    base = out_dir(args)
    rec = core.save_tsv(".", name, rows, base=base, tool="scripts/nsd_pair_calibration.py", **prov)
    log(f"  wrote {rec['path']}  ({rec['bytes']} B, {len(rows)} rows)")
    return rec


def save_json(args: argparse.Namespace, name: str, obj: Any, **prov: Any) -> dict:
    base = out_dir(args)
    folder = base / "."
    folder.mkdir(parents=True, exist_ok=True)
    rec = core.save_bytes(".", name, json.dumps(obj, indent=1, ensure_ascii=False, default=str).encode("utf-8"),
                          base=base, tool="scripts/nsd_pair_calibration.py", **prov)
    log(f"  wrote {rec['path']}  ({rec['bytes']} B)")
    return rec


def load_json(args: argparse.Namespace, name: str) -> Any:
    path = out_dir(args) / name
    if not path.exists():
        raise SystemExit(f"missing {core.rel(path)} - run the earlier stage first")
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(args: argparse.Namespace, name: str) -> list[dict]:
    path = out_dir(args) / name
    if not path.exists():
        raise SystemExit(f"missing {core.rel(path)} - run the earlier stage first")
    return list(csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"))


async def ens_get(path: str, **params: Any) -> Any:
    """Ensembl REST answered 500/503 intermittently on 2026-09-17; give it more attempts than core's default."""
    return await core.get_json(f"{ENSEMBL}/{path.lstrip('/')}", params=params or None, retries=9, timeout=180)


async def gnomad(query: str) -> dict[str, Any]:
    data = await core.post_json(GNOMAD_API, {"query": query}, timeout=300)
    if data.get("errors"):
        raise RuntimeError(f"gnomAD GraphQL error: {data['errors'][:2]}")
    return data["data"]


# =========================================================================== stage: identity

CONSTRAINT_FIELDS = ("oe_lof_upper pli mis_z oe_mis oe_mis_lower oe_mis_upper obs_mis exp_mis "
                     "obs_syn exp_syn oe_syn syn_z obs_lof exp_lof oe_lof flags")


async def stage_identity(args: argparse.Namespace) -> dict[str, Any]:
    report: dict[str, Any] = {"genes": {}}
    q = " ".join(
        f'{tag}: gene(gene_symbol: "{sym}", reference_genome: GRCh38) {{ gene_id symbol chrom start stop strand '
        f"canonical_transcript_id mane_select_transcript {{ ensembl_id ensembl_version refseq_id refseq_version }} "
        f"gnomad_constraint {{ {CONSTRAINT_FIELDS} }} }}"
        for tag, sym in (("teacher", TEACHER), ("student", STUDENT))
    )
    data = await gnomad("{ " + q + " }")
    save_json(args, f"{TEACHER}_{STUDENT}_gnomad_gene_{core.today()}.json", data,
              source="gnomAD v4.1 GraphQL gene query (MANE Select + constraint)", urls=[GNOMAD_API])

    for tag, sym in (("teacher", TEACHER), ("student", STUDENT)):
        g = data[tag]
        mane = g["mane_select_transcript"]
        enst = f"{mane['ensembl_id']}.{mane['ensembl_version']}"
        refseq = f"{mane['refseq_id']}.{mane['refseq_version']}"

        gene = await ens_get(f"lookup/symbol/homo_sapiens/{sym}")
        tx = await ens_get(f"lookup/id/{mane['ensembl_id']}", expand=1)
        seq_doc = await ens_get(f"sequence/id/{mane['ensembl_id']}", type="protein", object_type="transcript")
        ens_prot = seq_doc["seq"].rstrip("*")
        core.save_bytes("sequences", f"{sym}_{mane['ensembl_id']}.fasta",
                        f">{sym}|{enst}|{refseq}\n{ens_prot}\n".encode("utf-8"),
                        source="Ensembl REST protein sequence of the MANE Select transcript",
                        urls=[f"{ENSEMBL}/sequence/id/{mane['ensembl_id']}"],
                        tool="scripts/nsd_pair_calibration.py", params={"gene": sym, "type": "protein"})

        # UniProt: reviewed human entry for this gene, then the isoform whose sequence equals the MANE one
        hits = await core.get_json(f"{UNIPROT}/search", params={
            "query": f"gene_exact:{sym} AND organism_id:9606 AND reviewed:true", "format": "json",
            "fields": "accession,id,protein_name,length,sequence,xref_ensembl,xref_refseq"})
        entry = hits["results"][0]
        acc = entry["primaryAccession"]
        full = await core.get_json(f"{UNIPROT}/{acc}.json")
        core.save_json("uniprot", f"{acc}_{core.today()}.json", full, source="UniProt REST entry",
                       urls=[f"{UNIPROT}/{acc}.json"], tool="scripts/nsd_pair_calibration.py",
                       params={"accession": acc})
        canon = full["sequence"]["value"]

        matching_isoform, offsets = acc, []
        if canon != ens_prot:
            # isoform trap: find the isoform that does match, and record the offset table
            iso_accs = [alt for iso in (full.get("alternativeProducts") or {}).get("isoforms", [])
                        for alt in iso.get("isoformIds", []) if alt != acc]
            matching_isoform = None
            for iso in iso_accs:
                fa = (await core.request("GET", f"{UNIPROT}/{iso}.fasta", headers={"Accept": "text/plain"})).text
                iso_seq = "".join(fa.splitlines()[1:])
                if iso_seq == ens_prot:
                    matching_isoform = iso
                    break
            offsets = _offset_table(canon, ens_prot)

        report["genes"][sym] = {
            "role": "paralog_reference (clinical labels)" if tag == "teacher" else "target (clinically empty sibling)",
            "ensembl_gene": g["gene_id"],
            "locus": f"chr{g['chrom']}:{g['start']}-{g['stop']}",
            "strand": g["strand"],
            "ensembl_gene_version": gene.get("version"),
            "mane_transcript": enst,
            "refseq_transcript": refseq,
            "ensembl_canonical_transcript": g["canonical_transcript_id"],
            "mane_equals_canonical": g["canonical_transcript_id"] == mane["ensembl_id"],
            "translation_id": (tx.get("Translation") or {}).get("id"),
            "translation_length": (tx.get("Translation") or {}).get("length"),
            "ensembl_protein_length": len(ens_prot),
            "uniprot": acc,
            "uniprot_entry_name": full.get("uniProtkbId"),
            "uniprot_canonical_length": len(canon),
            "uniprot_isoform_matching_mane": matching_isoform,
            "uniprot_canonical_equals_mane": canon == ens_prot,
            "isoform_offsets": offsets,
            "sequence_sha256_prefix": hashlib.sha256(ens_prot.encode()).hexdigest()[:16],
            "constraint": g["gnomad_constraint"],
        }
        log(f"  {sym}: {enst} / {refseq} / {acc}  {len(ens_prot)} aa  "
            f"MANE==canonical:{g['canonical_transcript_id'] == mane['ensembl_id']}  "
            f"UniProt canonical==MANE:{canon == ens_prot}")

    save_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json", report,
              source="derived: gnomAD gene query + Ensembl REST lookup/sequence + UniProt REST",
              urls=[GNOMAD_API, ENSEMBL, UNIPROT])
    return report


def _offset_table(canonical: str, mane: str) -> list[dict[str, int]]:
    """Blocks of (canonical position -> MANE position) offsets, for the isoform trap."""
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score, aligner.end_gap_score = -10.0, -0.5, 0.0
    aln = aligner.align(canonical, mane)[0]
    blocks = []
    for (c0, c1), (m0, _m1) in zip(*aln.aligned.tolist()):
        blocks.append({"canonical_from": c0 + 1, "canonical_to": c1, "offset": m0 - c0})
    return blocks


def sequences(args: argparse.Namespace) -> dict[str, str]:
    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    out = {}
    for sym, info in ident.items():
        enst = info["mane_transcript"].split(".")[0]
        path = core.RAW / "sequences" / f"{sym}_{enst}.fasta"
        out[sym] = "".join(path.read_text(encoding="utf-8").splitlines()[1:])
    return out


# =========================================================================== stage: domains

async def stage_domains(args: argparse.Namespace) -> dict[str, Any]:
    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    rows: list[dict[str, Any]] = []
    for sym, info in ident.items():
        acc = info["uniprot"]
        entry = json.loads((core.RAW / "uniprot" / f"{acc}_{core.today()}.json").read_text(encoding="utf-8"))
        for f in entry.get("features", []):
            if f["type"] not in ("Domain", "Region", "Zinc finger", "Motif", "Binding site", "Compositional bias"):
                continue
            loc = f.get("location", {})
            rows.append({"gene": sym, "uniprot": acc, "source": "UniProt", "type": f["type"],
                         "name": f.get("description") or (f.get("ligand") or {}).get("name"),
                         "start": loc.get("start", {}).get("value"), "end": loc.get("end", {}).get("value")})
        ip = await core.get_json(f"{INTERPRO}/entry/all/protein/uniprot/{acc}/",
                                 params={"page_size": 200, "type": "domain"})
        core.save_json("interpro", f"{acc}_domains_{core.today()}.json", ip, source="InterPro API entries for protein",
                       urls=[f"{INTERPRO}/entry/all/protein/uniprot/{acc}/"],
                       tool="scripts/nsd_pair_calibration.py", params={"accession": acc, "type": "domain"})
        for item in ip.get("results", []):
            meta = item["metadata"]
            for prot in item.get("proteins", []):
                for frag in prot.get("entry_protein_locations", [{}]):
                    for f in frag.get("fragments", [{}]):
                        rows.append({"gene": sym, "uniprot": acc, "source": f"InterPro/{meta.get('source_database')}",
                                     "type": meta.get("type"), "name": f"{meta.get('accession')} {meta.get('name')}",
                                     "start": f.get("start"), "end": f.get("end")})
    save_tsv(args, f"{TEACHER}_{STUDENT}_domains_{core.today()}.tsv", rows,
             source="UniProt features + InterPro domain entries", urls=[UNIPROT, INTERPRO])
    for sym in ident:
        sub = [r for r in rows if r["gene"] == sym and r["source"] == "UniProt"]
        log(f"  {sym}: {len(sub)} UniProt features, "
            f"{sum(1 for r in rows if r['gene'] == sym) - len(sub)} InterPro fragments")
    return {"n_rows": len(rows)}


# =========================================================================== stage: universe

GNOMAD_TX_QUERY = """
{ transcript(transcript_id: "%s", reference_genome: GRCh38) {
    transcript_id transcript_version
    variants(dataset: gnomad_r4) {
      variant_id chrom pos ref alt rsids hgvsc hgvsp consequence lof lof_filter flags
      exome { ac an af homozygote_count filters }
      genome { ac an af homozygote_count filters }
      joint { ac an homozygote_count filters }
      in_silico_predictors { id value }
    } } }"""


async def stage_universe(args: argparse.Namespace) -> dict[str, Any]:
    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    summary: dict[str, Any] = {}
    for sym, info in ident.items():
        enst = info["mane_transcript"].split(".")[0]
        tx = (await gnomad(GNOMAD_TX_QUERY % enst))["transcript"]
        rows = []
        for v in tx["variants"]:
            joint = v.get("joint") or {}
            ref_aa, prot_pos, alt_aa = variants._protein_change(v.get("hgvsp"))
            preds = {p["id"]: p["value"] for p in (v.get("in_silico_predictors") or [])}
            filters = list(joint.get("filters") or [])
            rows.append({
                "variant_id": v["variant_id"], "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
                "transcript": f"{tx['transcript_id']}.{tx['transcript_version']}", "consequence": v.get("consequence"),
                "hgvsc": v.get("hgvsc"), "hgvsp": v.get("hgvsp"), "protein_pos": prot_pos,
                "aa_ref": ref_aa, "aa_alt": alt_aa, "lof": v.get("lof"),
                "joint_ac": joint.get("ac"), "joint_an": joint.get("an"),
                "joint_af": (joint["ac"] / joint["an"]) if joint.get("an") else None,
                "joint_hom": joint.get("homozygote_count"), "joint_filters": filters,
                "exome_ac": (v.get("exome") or {}).get("ac"), "genome_ac": (v.get("genome") or {}).get("ac"),
                "flags": v.get("flags"), "cadd_phred": preds.get("cadd"), "revel": preds.get("revel_max"),
                "spliceai_ds_max": preds.get("spliceai_ds_max"), "phylop": preds.get("phylop"),
            })
        save_tsv(args, f"{sym}_{enst}_gnomad_r4_all_{core.today()}.tsv", rows,
                 source="gnomAD v4.1 GraphQL transcript query (MANE Select; never the gene query)",
                 urls=[GNOMAD_API], params={"transcript": enst})
        mis = [r for r in rows if r["consequence"] == "missense_variant"]
        conf = [r for r in mis if (r["joint_ac"] or 0) >= 1 and "AC0" not in (r["joint_filters"] or [])]
        summary[sym] = {
            "gnomad": {
                "transcript": f"{tx['transcript_id']}.{tx['transcript_version']}",
                "n_variants_on_transcript": len(rows), "n_missense": len(mis),
                "n_missense_confidently_observed_ac_ge1_no_AC0": len(conf),
                "n_missense_AC0_filtered": sum(1 for r in mis if "AC0" in (r["joint_filters"] or [])),
                "n_missense_ac0_value": sum(1 for r in mis if (r["joint_ac"] or 0) == 0),
                "consequence_counts": dict(Counter(r["consequence"] for r in rows).most_common(10)),
                "distinct_missense_positions": len({r["protein_pos"] for r in conf if r["protein_pos"]}),
            }
        }
        log(f"  {sym} gnomAD: {len(rows)} variants on {enst}, {len(mis)} missense, {len(conf)} confidently observed")

        cv_rows, cv_counts = await _clinvar_records(sym, info)
        save_tsv(args, f"{sym}_clinvar_{core.today()}.tsv", cv_rows,
                 source="NCBI E-utilities ClinVar esearch+esummary (single_gene)",
                 urls=[f"{EUTILS}/esummary.fcgi"], params={"gene": sym})
        summary[sym]["clinvar"] = cv_counts
        log(f"  {sym} ClinVar: {cv_counts['n_records']} records, {cv_counts['n_missense']} missense, "
            f"{cv_counts['n_missense_P_LP']} P/LP missense (esearch Properties count "
            f"{cv_counts['esearch_plp_missense_properties']}, inflated field form "
            f"{cv_counts['esearch_plp_missense_wrong_field']})")
    save_json(args, f"{TEACHER}_{STUDENT}_universe_{core.today()}.json", summary,
              source="derived: gnomAD transcript query + ClinVar esummary", urls=[GNOMAD_API, EUTILS])
    return summary


async def _esearch_count(term: str) -> int:
    key = core.secret("NCBI_API_KEY")
    res = await core.get_json(f"{EUTILS}/esearch.fcgi", params={
        "db": "clinvar", "retmode": "json", "term": term, "retmax": 0, **({"api_key": key} if key else {})})
    return int(res["esearchresult"]["count"])


async def _clinvar_records(sym: str, info: dict) -> tuple[list[dict], dict]:
    key = core.secret("NCBI_API_KEY")
    base = {"db": "clinvar", "retmode": "json", **({"api_key": key} if key else {})}
    term = f"{sym}[gene] AND single_gene[prop]"
    search = await core.get_json(f"{EUTILS}/esearch.fcgi",
                                 params={**base, "term": term, "usehistory": "y", "retmax": 0})
    res = search["esearchresult"]
    total = int(res["count"])
    mane_nm = info["refseq_transcript"].split(".")[0]
    rows: list[dict] = []
    for start in range(0, total, 400):
        page = await core.get_json(f"{EUTILS}/esummary.fcgi", params={
            **base, "query_key": res["querykey"], "WebEnv": res["webenv"], "retstart": start, "retmax": 400},
            timeout=300)
        result = page["result"]
        for uid in result.get("uids", []):
            r = result[uid]
            germ = r.get("germline_classification") or {}
            vs = (r.get("variation_set") or [{}])[0]
            loc = next((l for l in vs.get("variation_loc", []) if l.get("assembly_name") == "GRCh38"), {})
            m = variants.TITLE_RE.match(r.get("title", ""))
            title_nm, cdna, prot = (m.group(1), m.group(3), m.group(4)) if m else (None, None, None)
            ref_aa, prot_pos, alt_aa = variants._protein_change(prot)
            spdi = vs.get("canonical_spdi") or ""
            parts = spdi.split(":")
            review = germ.get("review_status", "")
            rows.append({
                "uid": uid, "accession": r.get("accession"), "title": r.get("title"),
                "variant_type": vs.get("variant_type"), "classification": germ.get("description"),
                "review_status": review, "stars": variants.STARS.get(review, 0),
                "last_evaluated": (germ.get("last_evaluated") or "")[:10],
                "traits": [t.get("trait_name") for t in germ.get("trait_set", [])],
                "molecular_consequence": r.get("molecular_consequence_list"),
                "title_transcript": title_nm, "title_is_mane": title_nm == mane_nm,
                "hgvsc": cdna, "hgvsp": prot, "protein_pos": prot_pos if title_nm == mane_nm else None,
                "aa_ref": ref_aa, "aa_alt": alt_aa,
                "chrom": loc.get("chr"), "start": loc.get("start"), "stop": loc.get("stop"), "spdi": spdi,
                "ref": parts[2] if len(parts) == 4 else None, "alt": parts[3] if len(parts) == 4 else None,
            })
    mis = [r for r in rows if r["molecular_consequence"] and "missense variant" in r["molecular_consequence"]]
    plp = [r for r in mis if (r["classification"] or "").lower() in P_LP]
    blb = [r for r in mis if (r["classification"] or "").lower() in B_LB]
    # trap 14: the only Entrez form that is a real field is "clinsig ..."[Properties]
    right = (f'{sym}[gene] AND single_gene[prop] AND ("clinsig pathogenic"[Properties] OR '
             f'"clinsig likely pathogenic"[Properties]) AND "missense variant"[molecular consequence]')
    wrong = f'{sym}[gene] AND single_gene[prop] AND pathogenic[clinsig] AND "missense variant"[molecular consequence]'
    counts = {
        "n_records": total, "n_parsed": len(rows), "query": term,
        "classification_counts": dict(Counter(r["classification"] for r in rows).most_common(12)),
        "n_missense": len(mis), "n_missense_P_LP": len(plp), "n_missense_B_LB": len(blb),
        "n_missense_VUS": sum(1 for r in mis if (r["classification"] or "").lower() == "uncertain significance"),
        "n_missense_P_LP_on_mane_title": sum(1 for r in plp if r["protein_pos"]),
        "records_not_on_mane_title": sum(1 for r in rows if r["title_transcript"] and not r["title_is_mane"]),
        "esearch_plp_missense_properties": await _esearch_count(right),
        "esearch_plp_missense_wrong_field": await _esearch_count(wrong),
    }
    return rows, counts


# =========================================================================== stage: map

def align_pair(src_seq: str, tgt_seq: str) -> dict[int, int]:
    """Exactly atlas_mcp.structure.paralog_map's parameters, so the two case studies are comparable."""
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


def stage_map(args: argparse.Namespace) -> dict[str, Any]:
    seqs = sequences(args)
    src, tgt = seqs[TEACHER], seqs[STUDENT]
    mapping = align_pair(src, tgt)
    rows = [{f"{TEACHER}_pos": s, f"{TEACHER}_aa": src[s - 1], f"{STUDENT}_pos": t, f"{STUDENT}_aa": tgt[t - 1],
             "identical": src[s - 1] == tgt[t - 1]} for s, t in sorted(mapping.items())]
    core.save_tsv("paralog_maps", f"{TEACHER}_to_{STUDENT}_pairwise_MANE.tsv", rows, base=core.INTERIM,
                  source="derived: Biopython PairwiseAligner (BLOSUM62, gap -10/-0.5, free end gaps) on Ensembl MANE proteins",
                  urls=[f"{ENSEMBL}/sequence/id"], tool="scripts/nsd_pair_calibration.py",
                  params={"source": TEACHER, "target": STUDENT})
    ident = sum(1 for r in rows if r["identical"])
    log(f"  aligned pairs {len(mapping)}, identical {ident} ({ident / len(mapping):.3f} over aligned pairs, "
        f"{ident / min(len(src), len(tgt)):.3f} over the shorter protein)")

    # identity per domain, never one global number: a single figure hides the structure of the relationship
    per_domain = []
    for a, b, name in DOMAINS[TEACHER]:
        span = [pos for pos in range(a, b + 1)]
        aligned = [pos for pos in span if pos in mapping]
        same = [pos for pos in aligned if src[pos - 1] == tgt[mapping[pos] - 1]]
        per_domain.append({
            "domain": name, f"{TEACHER}_start": a, f"{TEACHER}_end": b, "n_aa": len(span),
            f"{STUDENT}_start": mapping.get(aligned[0]) if aligned else None,
            f"{STUDENT}_end": mapping.get(aligned[-1]) if aligned else None,
            "n_aligned": len(aligned), "gap_fraction": round(1 - len(aligned) / len(span), 3),
            "n_identical": len(same),
            "identity_over_aligned": round(len(same) / len(aligned), 3) if aligned else None,
            "identity_over_domain": round(len(same) / len(span), 3),
            "is_shared_module_part": name in MODULE_NAMES,
        })
        log(f"    {name:<12} {a:>5}-{b:<5} aligned {len(aligned):>4}/{len(span):<4} "
            f"identity {per_domain[-1]['identity_over_aligned']}")
    save_tsv(args, f"{TEACHER}_{STUDENT}_domain_identity_{core.today()}.tsv", per_domain,
             source="derived: paralog residue map, identity computed per annotated domain", urls=[])

    # landmark check: the map is only trustworthy if UniProt's own ligand-binding residues line up
    landmarks = []
    dom_path = out_dir(args) / f"{TEACHER}_{STUDENT}_domains_{core.today()}.tsv"
    if dom_path.exists():
        for r in csv.DictReader(open(dom_path, encoding="utf-8"), delimiter="	"):
            if r["gene"] != TEACHER or r["type"] not in ("Binding site", "Site"):
                continue
            for pos in range(int(r["start"]), int(r["end"]) + 1):
                t = mapping.get(pos)
                landmarks.append({"feature": f'{r["type"]}: {r["name"]}', f"{TEACHER}_pos": pos,
                                  f"{TEACHER}_aa": src[pos - 1], f"{STUDENT}_pos": t,
                                  f"{STUDENT}_aa": tgt[t - 1] if t else None,
                                  "identical": bool(t) and src[pos - 1] == tgt[t - 1]})
        n_id = sum(1 for r in landmarks if r["identical"])
        log(f"  landmark check: {n_id}/{len(landmarks)} UniProt ligand-binding residues of {TEACHER} "
            f"map to an identical {STUDENT} residue")
        save_tsv(args, f"{TEACHER}_{STUDENT}_landmark_check_{core.today()}.tsv", landmarks,
                 source="derived: UniProt binding-site features mapped through the paralog map", urls=[])
    lo, hi = SHARED_MODULE[TEACHER]
    mod = [pos for pos in range(lo, hi + 1) if pos in mapping]
    mod_ident = sum(1 for pos in mod if src[pos - 1] == tgt[mapping[pos] - 1])
    return {"aligned_pairs": len(mapping), "identical": ident,
            "identity_over_aligned": round(ident / len(mapping), 3),
            "identity_over_shorter": round(ident / min(len(src), len(tgt)), 3),
            f"{TEACHER}_residues_without_partner": len(src) - len(mapping),
            "landmark_check": {"n": len(landmarks), "n_identical": sum(1 for r in landmarks if r["identical"])},
            "shared_module": {"span": SHARED_MODULE[TEACHER], "n_aa": hi - lo + 1, "n_aligned": len(mod),
                              "n_identical": mod_ident,
                              "identity_over_aligned": round(mod_ident / len(mod), 3) if mod else None},
            "per_domain": per_domain}


def paralog_map_table(args: argparse.Namespace) -> dict[int, dict]:
    path = core.INTERIM / "paralog_maps" / f"{TEACHER}_to_{STUDENT}_pairwise_MANE.tsv"
    if not path.exists():
        raise SystemExit("run --stage map first")
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter="\t"):
        out[int(r[f"{TEACHER}_pos"])] = {
            "student_pos": int(r[f"{STUDENT}_pos"]), "teacher_aa": r[f"{TEACHER}_aa"],
            "student_aa": r[f"{STUDENT}_aa"], "identical": r["identical"] == "True"}
    return out


# =========================================================================== domain tables

# NSD1 modules in MANE numbering (UniProt Q96L73 features + InterPro, both fetched live 2026-09-17).
# Non-overlapping intervals: where UniProt and InterPro disagree by a few residues the union is used.
# InterPro/CDD resolves five PHD fingers while UniProt annotates four zinc fingers, so the PHD naming
# follows InterPro (cd15648/15650/15653/15656/15659 for NSD1, cd15649/15652/15655/15658/15661 for NSD3).
NSD1_DOMAINS = [
    (1, 320, "N_ucu"), (321, 408, "PWWP1"), (409, 1538, "ara_bolge"),
    (1539, 1589, "PHD1"), (1590, 1638, "PHD2"), (1639, 1706, "PHD3"), (1707, 1751, "PHD4"),
    (1752, 1850, "PWWP2"), (1851, 1889, "AWS_onu"), (1890, 1941, "AWS"),
    (1942, 2065, "SET"), (2066, 2082, "PostSET"), (2083, 2117, "SET_sonrasi"),
    (2118, 2160, "PHD5"), (2161, 2210, "CysHis"), (2211, 2696, "C_ucu"),
]
NSD3_DOMAINS = [
    (1, 267, "N_ucu"), (268, 349, "PWWP1"), (350, 696, "ara_bolge"),
    (697, 747, "PHD1"), (748, 797, "PHD2"), (798, 851, "PHD3"), (852, 955, "PHD4"),
    (956, 1048, "PWWP2"), (1049, 1092, "AWS_onu"), (1093, 1144, "AWS"),
    (1145, 1268, "SET"), (1269, 1285, "PostSET"), (1286, 1320, "SET_sonrasi"),
    (1321, 1363, "PHD5"), (1364, 1413, "CysHis"), (1414, 1437, "C_ucu"),
]
# The analogue of KIF1A's motor domain: the contiguous catalytic module from the first PHD finger to the
# Cys-His rich domain. Everything downstream of the curation happens inside it.
SHARED_MODULE = {"NSD1": (1539, 2210), "NSD3": (697, 1413)}
MODULE_NAMES = ("PWWP1", "PHD1", "PHD2", "PHD3", "PHD4", "PWWP2", "AWS", "SET", "PostSET", "PHD5", "CysHis")
DOMAINS = {"NSD1": NSD1_DOMAINS, "NSD3": NSD3_DOMAINS}


def domain_of(gene: str, pos: int | None) -> str | None:
    if not pos:
        return None
    return next((name for a, b, name in DOMAINS[gene] if a <= int(pos) <= b), None)


def in_module(gene: str, pos: int | None) -> bool:
    lo, hi = SHARED_MODULE[gene]
    return bool(pos) and lo <= int(pos) <= hi


# =========================================================================== stage: curate


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


async def fetch_vcv(gene: str, ids: list[str]) -> list[ET.Element]:
    key = core.secret("NCBI_API_KEY")
    archives: list[ET.Element] = []
    for i, batch in enumerate(core.chunks(ids, 50)):
        cache = core.RAW / "clinvar" / "vcv_xml" / f"{gene}_vcv_{core.today()}_{i:03d}.xml"
        if cache.exists():
            content = cache.read_bytes()
        else:
            params = {"db": "clinvar", "id": ",".join(batch), "rettype": "vcv", "is_variationid": "",
                      **({"api_key": key} if key else {})}
            content = (await core.request("GET", f"{EUTILS}/efetch.fcgi", params=params,
                                          headers={"Accept": "application/xml"}, timeout=300)).content
            core.save_bytes("clinvar/vcv_xml", cache.name, content,
                            source="NCBI E-utilities ClinVar VCV XML", urls=[f"{EUTILS}/efetch.fcgi"],
                            tool="scripts/nsd_pair_calibration.py", params={"gene": gene, "n_ids": len(batch)})
        archives += ET.fromstring(content).findall("VariationArchive")
    log(f"    VCV XML: {len(archives)} archives for {len(ids)} ids")
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


def assign_tier(r: dict) -> tuple[str, str]:
    """Same decision tree as scripts/kif1a_clinvar_curation.py, with the shared module in place of the
    motor domain. Nothing here is tuned on the outcome."""
    c, g, stars = r["class"], r["group"], int(r["stars"] or 0)
    af = r["gnomad_af_max"]
    if r["variant_type"] != "single nucleotide variant" and g in ("missense", "splice_canonical"):
        return "exclude", "not an SNV"
    if c == "P/LP":
        if g == "lof":
            return "lof_reference", "LoF: haploinsufficiency reference, not a missense calibration label"
        if g == "inframe":
            return "inframe_reference", "in-frame indel: not scorable by AlphaMissense/ESM"
        if stars == 0:
            return "exclude", "no assertion criteria"
        if g == "splice_canonical":
            return "positive_splice", "canonical splice P/LP: SpliceAI / Evo 2 positive"
        if g != "missense":
            return "exclude", f"consequence {g}"
        if r["protein_pos"] is None:
            return "exclude", "no MANE protein position"
        if af is not None and af >= AF_DOMINANT_MAX:
            return "positive_review", f"gnomAD AF {af:.2g} too common for a dominant de novo disorder"
        if (r["classification"] or "").lower() == "likely pathogenic" and int(r["n_submitters"] or 0) < 2 \
                and r["text_in_silico"] and not (r["text_de_novo"] or r["text_functional"]):
            return "positive_circular", "single-submitter LP whose comment cites in-silico evidence (PP3 may be decisive)"
        independent = stars >= 2 or r["text_de_novo"] or r["text_functional"]
        return ("positive_gold" if independent else "positive_silver",
                "independent evidence (>=2 stars, de novo or functional)" if independent else "P/LP, single submitter")
    if c == "B/LB":
        if g != "missense":
            return "negative_noncoding_reference", f"{g} B/LB"
        if stars == 0:
            return "exclude", "no assertion criteria"
        if (af is not None and af >= 1e-4) or (r["text_frequency"] and not r["text_in_silico"]):
            return "negative_frequency", "benign by population frequency"
        if r["text_in_silico"] and not r["text_frequency"]:
            return "negative_circular", "benign largely by BP4 / in-silico evidence"
        return "negative_other", "benign by other evidence"
    if c == "conflict" and g in ("missense", "splice_canonical"):
        e = (r["explanation"] or "").lower()
        has_p, has_b = "pathogenic" in e, "benign" in e
        if has_p and not has_b:
            return "conflict_P_vs_VUS", "sensitivity analysis only"
        if has_b and not has_p:
            return "conflict_B_vs_VUS", "sensitivity analysis only"
        return "exclude", "conflict spans P and B"
    return "exclude", f"{c} / {g}"


async def vep_batch(lines: list[str], tag: str) -> dict[str, dict]:
    """One VEP run for everything that needs scoring: the two classes of a comparison never come from
    two sources (trap 11). Input is sorted because VEP rejects unsorted VCF-style lines (trap 15)."""
    lines = sorted(lines, key=variants.vcf_sort_key)
    out: dict[str, dict] = {}
    raw_all = []
    for i, batch in enumerate(core.chunks(lines, 200)):
        raw = await core.post_json(VEP_REGION, {"variants": batch}, params=variants.VEP_PARAMS,
                                   headers={"Content-Type": "application/json"}, timeout=900)
        raw_all += raw
        for v in raw:
            chrom, pos, _dot, ref, alt = v["input"].split()[:5]
            # _flatten_vep indexes the dbNSFP per-transcript lists against Ensembl_transcriptid (trap 16)
            out[f"{chrom}-{pos}-{ref}-{alt}"] = variants._flatten_vep(v)
        log(f"    VEP batch {i + 1}: {len(out)}/{len(lines)}")
    core.save_json("vep", f"{tag}_{core.today()}.raw.json", raw_all, compress=True, source="Ensembl VEP REST",
                   urls=[VEP_REGION], tool="scripts/nsd_pair_calibration.py", params=variants.VEP_PARAMS)
    return out


def _paralog_columns(pmap: dict[int, dict], pos, ref, alt) -> dict[str, Any]:
    m = pmap.get(int(pos)) if pos else None
    return {"protein_change": f"{ref}{pos}{alt}" if pos else "", "protein_pos": pos,
            "cv_group": f"res{pos}" if pos else "",
            "student_pos": m["student_pos"] if m else "", "student_aa": m["student_aa"] if m else "",
            "residue_identical_in_student": bool(m) and m["student_aa"] == ref}


async def stage_curate(args: argparse.Namespace) -> dict[str, Any]:
    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    chrom = ident[TEACHER]["locus"].split(":")[0].removeprefix("chr")
    seq = sequences(args)[TEACHER]
    pmap = paralog_map_table(args)

    rows = load_tsv(args, f"{TEACHER}_clinvar_{core.today()}.tsv")
    for r in rows:
        r["class"], r["group"] = cls(r["classification"]), csq_group(r["molecular_consequence"] or "")
    cand = [r for r in rows if r["class"] == "P/LP"
            or (r["class"] == "B/LB" and r["group"] in ("missense", "splice_canonical", "inframe"))
            or (r["class"] == "conflict" and r["group"] in ("missense", "splice_canonical"))]
    log(f"  ClinVar {len(rows)} records -> {len(cand)} curation candidates")

    parsed = {p["variation_id"]: p for p in map(parse_vcv, await fetch_vcv(TEACHER, [r["uid"] for r in cand]))}
    lines, key_of = [], {}
    for r in cand:
        if r["spdi"] and r["ref"] and r["alt"] and len(r["ref"]) == 1 and len(r["alt"]) == 1:
            _acc, p0, ref, alt = r["spdi"].split(":")
            lines.append(f"{chrom} {int(p0) + 1} . {ref} {alt} . . .")
            key_of[r["uid"]] = f"{chrom}-{int(p0) + 1}-{ref}-{alt}"
    vep = await vep_batch(lines, f"{TEACHER}_clinvar_candidates_vep")

    curated = []
    for r in cand:
        p = parsed.get(r["uid"], {})
        v = vep.get(key_of.get(r["uid"], ""), {})
        afs = [x for x in (v.get("gnomad_exomes_af"), v.get("gnomad_genomes_af")) if isinstance(x, (int, float))]
        pos = v.get("protein_pos") if v.get("mane_select") else None
        row = {
            "variation_id": r["uid"], "accession": r["accession"], "hgvsc": r["hgvsc"], "hgvsp_vep": v.get("hgvsp"),
            "protein_pos": pos, "aa_ref": v.get("aa_ref"), "aa_alt": v.get("aa_alt"),
            "domain": domain_of(TEACHER, pos), "in_shared_module": in_module(TEACHER, pos),
            "classification": r["classification"], "class": r["class"], "group": r["group"],
            "stars": r["stars"], "review_status": r["review_status"], "variant_type": r["variant_type"],
            "variant_key": key_of.get(r["uid"]), "gnomad_af_max": max(afs) if afs else None,
            "am_pathogenicity": v.get("am_pathogenicity"), "cadd_phred": v.get("cadd_phred"),
            "revel": v.get("revel"), "spliceai_ds_max": v.get("spliceai_ds_max"),
            "last_evaluated": r["last_evaluated"], "traits": r["traits"],
            **{k: p.get(k) for k in ("date_created", "n_submitters", "explanation", "submitters", "methods",
                                     "origins", "n_pmids", "literature_only",
                                     *[f"text_{t}" for t in TEXT_FLAGS])},
        }
        row["tier"], row["reason"] = assign_tier(row)
        curated.append(row)
    save_tsv(args, f"{TEACHER}_clinvar_curated_{core.today()}.tsv", curated,
             source="derived: ClinVar VCV XML + one Ensembl VEP run", urls=[EUTILS, VEP_REGION])

    log("  tiers: " + ", ".join(f"{t}={n}" for t, n in Counter(r["tier"] for r in curated).most_common()))
    pos_core = [r for r in curated if r["tier"] in ("positive_gold", "positive_silver")]
    log("  positives by domain: " + str(Counter(r["domain"] for r in pos_core).most_common()))
    log(f"  positives inside the shared module {SHARED_MODULE[TEACHER]}: "
        f"{sum(1 for r in pos_core if r['in_shared_module'])}/{len(pos_core)}")

    # ------------------------------------------------------------------ negatives, same region as positives
    labelled = {(r["protein_pos"], r["aa_alt"]) for r in curated if r["class"] in ("P/LP", "conflict")}
    enst = ident[TEACHER]["mane_transcript"].split(".")[0]
    gn = [r for r in load_tsv(args, f"{TEACHER}_{enst}_gnomad_r4_all_{core.today()}.tsv")
          if r["consequence"] == "missense_variant"]
    controls = []
    for r in gn:
        if not r["protein_pos"] or not in_module(TEACHER, int(r["protein_pos"])):
            continue
        if (int(r["protein_pos"]), r["aa_alt"]) in labelled:
            continue
        if "AC0" in (r["joint_filters"] or ""):
            continue
        if seq[int(r["protein_pos"]) - 1] != r["aa_ref"]:
            raise RuntimeError(f"MANE sequence mismatch at {r['hgvsp']} ({TEACHER})")
        controls.append(r)
    ac_bins = Counter("AC1" if int(r["joint_ac"] or 0) == 1 else "AC2" if int(r["joint_ac"] or 0) == 2
                      else "AC3-4" if int(r["joint_ac"] or 0) < 5 else "AC>=5" for r in controls)
    log(f"  gnomAD module missense not labelled in ClinVar: {len(controls)}  {ac_bins.most_common()}")

    sets: list[dict] = []
    for r in curated:
        set_name = {
            "positive_gold": "missense_positive_core", "positive_silver": "missense_positive_core",
            "negative_other": "missense_negative_clinvar_outside_module",
            "negative_frequency": "missense_negative_clinvar_outside_module",
            "positive_splice": "splice_positive", "positive_circular": "heldout_pp3_dependent_positive",
            "conflict_P_vs_VUS": "sensitivity_conflict_P_vs_VUS", "lof_reference": "reference_lof",
        }.get(r["tier"])
        if not set_name:
            continue
        sets.append({
            "set": set_name,
            "label": 1 if set_name in ("missense_positive_core", "splice_positive", "heldout_pp3_dependent_positive")
            else 0 if set_name.startswith("missense_negative") else "",
            "source": "ClinVar", "id": r["accession"], "hgvsc": r["hgvsc"], "variant_key": r["variant_key"],
            "domain": r["domain"], "in_shared_module": r["in_shared_module"],
            **_paralog_columns(pmap, r["protein_pos"], r["aa_ref"], r["aa_alt"]),
            "classification": r["classification"], "stars": r["stars"], "review_status": r["review_status"],
            "tier": r["tier"], "gnomad_af_max": r["gnomad_af_max"],
            "text_in_silico": r["text_in_silico"], "text_functional": r["text_functional"],
            "text_de_novo": r["text_de_novo"],
        })
    for c in controls:
        ac = int(c["joint_ac"] or 0)
        set_name = ("missense_negative_module_population" if ac >= args.min_ac else
                    "missense_negative_module_population_ac2" if ac == 2 else None)
        if not set_name:
            continue
        sets.append({
            "set": set_name, "label": 0, "source": "gnomAD v4.1", "id": c["variant_id"], "hgvsc": c["hgvsc"],
            "variant_key": c["variant_id"], "domain": domain_of(TEACHER, int(c["protein_pos"])),
            "in_shared_module": True,
            **_paralog_columns(pmap, c["protein_pos"], c["aa_ref"], c["aa_alt"]),
            "classification": f"population AC={ac} hom={c['joint_hom']}", "stars": "", "review_status": "",
            "tier": "", "gnomad_af_max": c["joint_af"],
            "text_in_silico": "", "text_functional": "", "text_de_novo": "",
        })
    save_tsv(args, f"{TEACHER}_calibration_sets_{core.today()}.tsv", sets,
             source="derived: curated ClinVar + gnomAD module controls + paralog map", urls=[],
             params={"population_min_ac": args.min_ac, "shared_module": str(SHARED_MODULE[TEACHER])})

    report: dict[str, Any] = {"shared_module": SHARED_MODULE[TEACHER], "sets": {},
                              "gnomad_control_allele_count_bins": dict(ac_bins.most_common())}
    for s, n in Counter(r["set"] for r in sets).most_common():
        sub = [r for r in sets if r["set"] == s]
        report["sets"][s] = {
            "n": n, "n_residues": len({r["protein_pos"] for r in sub if r["protein_pos"]}),
            "n_in_module": sum(1 for r in sub if r["in_shared_module"] is True),
            "n_identical_residue_in_student": sum(1 for r in sub if r["residue_identical_in_student"] is True),
            "by_domain": dict(Counter(r["domain"] for r in sub).most_common()),
        }
        log(f"  {n:>4}  {s:<46} residues={report['sets'][s]['n_residues']:<4} "
            f"module={report['sets'][s]['n_in_module']:<4} "
            f"identical={report['sets'][s]['n_identical_residue_in_student']}")
    report["tiers"] = dict(Counter(r["tier"] for r in curated).most_common())
    report["exclusions"] = dict(Counter(r["reason"] for r in curated if r["tier"] == "exclude").most_common(10))
    report["positives_by_domain"] = dict(Counter(r["domain"] for r in pos_core).most_common())
    report["n_gnomad_module_controls"] = len(controls)
    report["evidence_flags_missense_P_LP"] = {
        k: sum(1 for r in curated if r["class"] == "P/LP" and r["group"] == "missense" and r[k])
        for k in ("text_in_silico", "text_functional", "text_de_novo", "text_frequency", "text_same_residue")}
    report["n_missense_P_LP_candidates"] = sum(1 for r in curated if r["class"] == "P/LP" and r["group"] == "missense")
    save_json(args, f"{TEACHER}_curation_summary_{core.today()}.json", report,
              source="derived: curation summary", urls=[])
    return report


# =========================================================================== stage: score


async def stage_score(args: argparse.Namespace) -> dict[str, Any]:
    sets = load_tsv(args, f"{TEACHER}_calibration_sets_{core.today()}.tsv")
    todo = [r for r in sets if r["variant_key"] and r["set"] not in ("splice_positive", "reference_lof")]
    keys = sorted({r["variant_key"] for r in todo})
    lines = [f"{k.split('-')[0]} {k.split('-')[1]} . {k.split('-')[2]} {k.split('-')[3]} . . ." for k in keys]
    vep = await vep_batch(lines, f"{TEACHER}_calibration_scores_vep")
    rows = []
    for r in todo:
        v = vep.get(r["variant_key"], {})
        row = {k: r[k] for k in ("set", "label", "source", "id", "protein_change", "protein_pos", "cv_group",
                                 "domain", "in_shared_module", "residue_identical_in_student",
                                 "student_pos", "student_aa", "gnomad_af_max", "tier", "stars",
                                 "text_in_silico", "text_functional", "text_de_novo")}
        for key in MODELS:
            try:
                row[key] = float(v.get(key))
            except (TypeError, ValueError):
                row[key] = None
        row["vep_consequence"] = v.get("most_severe_consequence")
        row["vep_hgvsp"] = v.get("hgvsp")
        rows.append(row)
    save_tsv(args, f"{TEACHER}_calibration_scores_{core.today()}.tsv", rows,
             source="derived: one Ensembl VEP run over positives and negatives together",
             urls=[VEP_REGION], params=variants.VEP_PARAMS)
    cov = {MODELS[k][0]: f"{sum(1 for r in rows if r[k] is not None)}/{len(rows)}" for k in MODELS}
    log(f"  score coverage: {cov}")
    return {"n_scored": len(rows), "coverage": cov}


# =========================================================================== stage: measure


def auroc(pos: np.ndarray, neg: np.ndarray) -> float | None:
    if pos.size == 0 or neg.size == 0:
        return None
    values = np.concatenate([pos, neg])
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(order.size, dtype=float)
    ranks[order] = np.arange(1, order.size + 1)
    sorted_vals = values[order]
    i = 0
    while i < sorted_vals.size:
        j = i
        while j + 1 < sorted_vals.size and sorted_vals[j + 1] == sorted_vals[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return float((ranks[:pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size))


def grouped_bootstrap_ci(pos_rows: list[tuple[str, float]], neg_rows: list[tuple[str, float]]):
    """Residues, not variants, are resampled: several substitutions of one residue are not
    independent observations (trap 13)."""
    if not pos_rows or not neg_rows:
        return None
    pg: dict[str, list[float]] = defaultdict(list)
    ng: dict[str, list[float]] = defaultdict(list)
    for g, s in pos_rows:
        pg[g].append(s)
    for g, s in neg_rows:
        ng[g].append(s)
    pk, nk = list(pg), list(ng)
    draws = []
    for _ in range(N_BOOT):
        p = np.concatenate([pg[pk[i]] for i in RNG.integers(0, len(pk), len(pk))])
        n = np.concatenate([ng[nk[i]] for i in RNG.integers(0, len(nk), len(nk))])
        a = auroc(np.asarray(p, float), np.asarray(n, float))
        if a is not None:
            draws.append(a)
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def sens_at_spec(pos: np.ndarray, neg: np.ndarray, spec: float = 0.90):
    """Sensitivity at a fixed specificity and the score threshold that achieves it: the number that
    actually gets transferred to the sibling."""
    if pos.size == 0 or neg.size == 0:
        return None
    threshold = float(np.quantile(neg, spec))
    return float((pos >= threshold).mean()), threshold


def stage_measure(args: argparse.Namespace) -> dict[str, Any]:
    rows = load_tsv(args, f"{TEACHER}_calibration_scores_{core.today()}.tsv")
    for r in rows:
        for k in MODELS:
            r[k] = float(r[k]) if r[k] not in ("", None) else None
        r["module"] = r["in_shared_module"] == "True"
        r["identical"] = r["residue_identical_in_student"] == "True"

    def pick(name: str, module_only: bool = False) -> list[dict]:
        out = [r for r in rows if r["set"] == name]
        return [r for r in out if r["module"]] if module_only else out

    pos_core_mod = pick("missense_positive_core", True)
    neg_pop = pick("missense_negative_module_population")
    neg_pop_ac2 = pick("missense_negative_module_population_ac2")
    comparisons = {
        f"PRIMARY module: core positives vs gnomAD AC>={args.min_ac}": (pos_core_mod, neg_pop),
        "module: core positives vs gnomAD AC>=2": (pos_core_mod, neg_pop + neg_pop_ac2),
        f"module: held-out PP3-dependent positives vs gnomAD AC>={args.min_ac}":
            (pick("heldout_pp3_dependent_positive", True), neg_pop),
        f"module: positives with functional evidence vs gnomAD AC>={args.min_ac}":
            ([r for r in pos_core_mod if r["text_functional"] == "True"], neg_pop),
        "WHOLE PROTEIN (confounded): core positives vs ClinVar B/LB":
            (pick("missense_positive_core"), pick("missense_negative_clinvar_outside_module")),
    }
    # boundary sensitivity: the module was fixed from UniProt/InterPro annotation before any label was
    # placed in it, but a reader will still ask whether the result rides on that choice
    for lo, hi, label in ((1968, 2219, "screen hotspot window"), (1942, 2082, "SET + Post-SET only")):
        def inside(rs, lo=lo, hi=hi):
            return [r for r in rs if r["protein_pos"] and lo <= int(r["protein_pos"]) <= hi]
        comparisons[f"SENSITIVITY {label} ({lo}-{hi})"] = (inside(pos_core_mod), inside(neg_pop))
    report: dict[str, Any] = {}
    for name, (p_rows, n_rows) in comparisons.items():
        log(f"\n  {name}  n+={len(p_rows)} n-={len(n_rows)}")
        entry: dict[str, Any] = {"n_pos": len(p_rows), "n_neg": len(n_rows),
                                 "n_pos_residues": len({r["protein_pos"] for r in p_rows}),
                                 "n_neg_residues": len({r["protein_pos"] for r in n_rows}), "models": {}}
        for key, (label, sign) in MODELS.items():
            p = [(r["cv_group"], sign * r[key]) for r in p_rows if r[key] is not None]
            n = [(r["cv_group"], sign * r[key]) for r in n_rows if r[key] is not None]
            a = auroc(np.array([s for _, s in p]), np.array([s for _, s in n]))
            if a is None:
                continue
            ci = grouped_bootstrap_ci(p, n)
            sens = sens_at_spec(np.array([s for _, s in p]), np.array([s for _, s in n]))
            entry["models"][label] = {"auroc": round(a, 3), "ci95": [round(ci[0], 3), round(ci[1], 3)],
                                      "sens_at_90_spec": round(sens[0], 3), "threshold": round(sign * sens[1], 4),
                                      "n_pos": len(p), "n_neg": len(n), "label_trained": key in LABEL_TRAINED}
            log(f"    {label:<16}{' (label-trained)' if key in LABEL_TRAINED else '':<16} AUROC {a:.3f} "
                f"95% CI {ci[0]:.3f}-{ci[1]:.3f}  sens@90%spec {sens[0]:.2f} (esik {sign * sens[1]:.4g})"
                f"  n+={len(p)} n-={len(n)}")
        for base_name, getter in (
            (f"paralog residue identical in {STUDENT}", lambda r: 1.0 if r["identical"] else 0.0),
            ("is inside the shared module", lambda r: 1.0 if r["module"] else 0.0),
        ):
            a = auroc(np.array([getter(r) for r in p_rows]), np.array([getter(r) for r in n_rows]))
            if a is not None:
                entry["models"][f"baseline: {base_name}"] = {"auroc": round(a, 3)}
                log(f"    baseline {base_name:<42} AUROC {a:.3f}")
        report[name] = entry

    # the empirical analogue of the screen's transfer lift: identity at pathogenic residues minus
    # identity at population-tolerated residues in the same region
    p_ident = [r for r in pos_core_mod]
    n_ident = [r for r in neg_pop]
    report["paralog_identity"] = {
        "positives_identical_fraction": round(sum(1 for r in p_ident if r["identical"]) / len(p_ident), 3)
        if p_ident else None,
        "negatives_identical_fraction": round(sum(1 for r in n_ident if r["identical"]) / len(n_ident), 3)
        if n_ident else None,
        "n_pos": len(p_ident), "n_neg": len(n_ident),
    }
    save_json(args, f"{TEACHER}_calibration_metrics_{core.today()}.json", report,
              source="derived: AUROC with residue-grouped bootstrap on one VEP score table", urls=[])
    return report


# =========================================================================== stage: constraint


async def stage_constraint(args: argparse.Namespace) -> dict[str, Any]:
    """Observed/expected synonymous and missense per region, the calibration check the method requires.

    Reuses scripts/regional_constraint.py wholesale (Samocha 2014 trinucleotide rates, PMID 25086666,
    expectation calibrated on synonymous sites), so the numbers are computed the same way as the
    KIF13A regional constraint axis. The calibration here is fitted on the pooled synonymous sites of
    NSD1 and NSD3 only, which makes the two-gene average synonymous o/e 1.0 by construction; the
    informative quantity is therefore the synonymous o/e of a single region.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    import regional_constraint as rc  # noqa: PLC0415

    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    rates = await rc.mutation_rates(refresh=False)
    per_gene: dict[str, list[dict]] = {}
    structs: dict[str, dict] = {}
    for gene in (TEACHER, STUDENT):
        enst = ident[gene]["mane_transcript"].split(".")[0]
        struct = await rc.transcript_structure(enst, False)
        pad_lo, pad_hi = struct["gene_start"] - 2, struct["gene_end"] + 2
        plus = await rc.genomic_sequence(struct["chrom"], pad_lo, pad_hi, f"{gene}_{enst}", False)
        cds, gpos, contexts = rc.build_cds(struct, plus, pad_lo)
        if len(cds) % 3:
            raise RuntimeError(f"{gene}: CDS length {len(cds)} is not a multiple of 3")
        possible = rc.enumerate_snvs(gene, cds, gpos, contexts, rates, struct["protein_length"])
        obs = await rc.gnomad_variants(gene, enst, False)
        rows, qc = rc.attach_observed(gene, possible, obs, cds, rc.MAX_AF)
        per_gene[gene], structs[gene] = rows, struct
        log(f"  {gene}: CDS {len(cds)} nt, {len(rows)} possible SNVs, gnomAD rows {len(obs)}, qc {qc}")
    cal = rc.fit_calibration(per_gene[TEACHER] + per_gene[STUDENT], rc.N_RATE_BINS)
    log(f"  calibration fitted on {cal['n_syn_possible']} possible synonymous sites "
        f"({cal['n_syn_observed']} observed) in {len(cal['bins'])} rate bins")
    stats = []
    for gene, rows in per_gene.items():
        rc.apply_calibration(rows, cal)
        length = structs[gene]["protein_length"]
        lo, hi = SHARED_MODULE[gene]
        regions = [(1, length, "tum protein"), (lo, hi, "paylasilan modul"),
                   *[(a, b, name) for a, b, name in DOMAINS[gene] if name in MODULE_NAMES]]
        for a, b, name in regions:
            stats.append(rc.region_stats(rows, a, b, name, gene))
        outside = [r for r in rows if not (lo <= r["codon"] <= hi)]
        st = rc.region_stats(outside, 1, length, "modul disi", gene)
        st["n_aa"] = length - (hi - lo + 1)
        stats.append(st)
    save_tsv(args, f"{TEACHER}_{STUDENT}_regional_constraint_{core.today()}.tsv", stats,
             source="derived: Samocha 2014 trinucleotide rates + gnomAD v4.1 exomes, synonymous-calibrated "
                    "expectation (scripts/regional_constraint.py helpers)",
             urls=[rc.RATES_URL, GNOMAD_API], params={"max_af": rc.MAX_AF, "n_rate_bins": rc.N_RATE_BINS})
    for st in stats:
        log(f"  {st['gene']:<5} {st['region']:<18} {st['n_aa']:>5} aa  "
            f"o/e mis {st['oe_mis']} ({st['oe_mis_lower']}-{st['oe_mis_upper']}, obs {st['obs_mis']})  "
            f"o/e syn {st['oe_syn']} (obs {st['obs_syn']})")
    return {"n_regions": len(stats),
            "calibration": {k: cal[k] for k in ("n_syn_possible", "n_syn_observed")},
            "rows": stats}


# =========================================================================== stage: transfer


def _tier(over: bool, identical: bool, in_mod: bool, teacher_path: list[str]) -> str:
    """The project's tier definitions, unchanged from scripts/build_kif13a_candidates.py."""
    if over and identical and teacher_path:
        return "A"
    if over and identical and in_mod:
        return "B"
    if over:
        return "C"
    if teacher_path and identical:
        return "D"
    return "-"


async def stage_transfer(args: argparse.Namespace) -> dict[str, Any]:
    ident = load_json(args, f"{TEACHER}_{STUDENT}_identity_{core.today()}.json")["genes"]
    metrics = load_json(args, f"{TEACHER}_calibration_metrics_{core.today()}.json")
    primary = next(k for k in metrics if k.startswith("PRIMARY"))
    am = metrics[primary]["models"]["AlphaMissense"]
    cutoff = args.am_cutoff if args.am_cutoff is not None else am["threshold"]
    generic = 0.564   # AlphaMissense's own likely-pathogenic cutoff, kept for the comparison
    log(f"  calibrated AlphaMissense threshold {cutoff} (90% specificity in {TEACHER}, "
        f"sens {am['sens_at_90_spec']}), generic {generic}")

    seq = sequences(args)[STUDENT]
    inv: dict[int, dict] = {}
    for t_pos, m in paralog_map_table(args).items():
        inv[m["student_pos"]] = {"teacher_pos": t_pos, "teacher_aa": m["teacher_aa"],
                                 "student_aa": m["student_aa"], "identical": m["identical"]}
    sets = load_tsv(args, f"{TEACHER}_calibration_sets_{core.today()}.tsv")
    teacher_path: dict[int, list[str]] = defaultdict(list)
    teacher_tol: dict[int, list[str]] = defaultdict(list)
    for r in sets:
        if not r["protein_pos"]:
            continue
        pos = int(r["protein_pos"])
        if r["set"] == "missense_positive_core":
            teacher_path[pos].append(r["protein_change"])
        elif r["set"].startswith("missense_negative_module_population"):
            teacher_tol[pos].append(r["protein_change"])

    enst = ident[STUDENT]["mane_transcript"].split(".")[0]
    gn = [r for r in load_tsv(args, f"{STUDENT}_{enst}_gnomad_r4_all_{core.today()}.tsv")
          if r["consequence"] == "missense_variant" and r["protein_pos"]]
    keys = sorted({r["variant_id"] for r in gn})
    cache = out_dir(args) / f"{STUDENT}_missense_scores_{core.today()}.tsv"
    if cache.exists():
        log(f"  using cached {core.rel(cache)}")
        scores = {r["variant_key"]: r for r in load_tsv(args, cache.name)}
    else:
        vep = await vep_batch([f"{k.split('-')[0]} {k.split('-')[1]} . {k.split('-')[2]} {k.split('-')[3]} . . ."
                               for k in keys], f"{STUDENT}_missense_vep")
        scores = {k: {"variant_key": k, **{m: v.get(m) for m in MODELS},
                      "protein_pos": v.get("protein_pos"), "aa_ref": v.get("aa_ref"), "aa_alt": v.get("aa_alt"),
                      "hgvsp": v.get("hgvsp"), "am_class": v.get("am_class"),
                      "consequence": v.get("most_severe_consequence")}
                  for k, v in vep.items()}
        save_tsv(args, cache.name, list(scores.values()),
                 source="Ensembl VEP REST on every gnomAD missense variant of the student gene",
                 urls=[VEP_REGION], params=variants.VEP_PARAMS)

    rows = []
    for g in gn:
        s = scores.get(g["variant_id"], {})
        pos = int(s.get("protein_pos") or g["protein_pos"] or 0) or None
        if not pos:
            continue
        if seq[pos - 1] != g["aa_ref"]:
            raise RuntimeError(f"MANE sequence mismatch at {g['hgvsp']} ({STUDENT})")
        try:
            am_score = float(s.get("am_pathogenicity"))
        except (TypeError, ValueError):
            am_score = None
        p = inv.get(pos, {})
        identical = bool(p) and p.get("identical")
        path_here = teacher_path.get(p.get("teacher_pos", -1), [])
        tol_here = teacher_tol.get(p.get("teacher_pos", -1), [])
        ac = int(g["joint_ac"] or 0)
        confident = ac >= 1 and "AC0" not in (g["joint_filters"] or "")
        over = am_score is not None and am_score >= cutoff
        rows.append({
            "variant_id": g["variant_id"], "hgvsc": g["hgvsc"], "hgvsp": g["hgvsp"],
            "protein_pos": pos, "aa_ref": g["aa_ref"], "aa_alt": g["aa_alt"],
            "protein_change": f"{g['aa_ref']}{pos}{g['aa_alt']}",
            "domain": domain_of(STUDENT, pos), "in_shared_module": in_module(STUDENT, pos),
            "joint_ac": ac, "joint_af": g["joint_af"], "joint_hom": g["joint_hom"],
            "confidently_observed": confident,
            "am_pathogenicity": am_score, "am_class": s.get("am_class"),
            "cadd_phred": s.get("cadd_phred"), "esm1b": s.get("esm1b"), "primateai": s.get("primateai"),
            "over_calibrated_cutoff": over,
            "over_generic_cutoff": am_score is not None and am_score >= generic,
            "teacher_pos": p.get("teacher_pos"), "teacher_aa": p.get("teacher_aa"),
            "residue_identical": identical,
            "teacher_pathogenic_here": ";".join(path_here),
            "teacher_tolerated_here": ";".join(tol_here),
            "tier": _tier(over, bool(identical), in_module(STUDENT, pos), path_here),
        })
    prov = dict(source=f"derived: gnomAD v4.1 {STUDENT} missense + VEP scores + {TEACHER}-calibrated threshold "
                       f"+ {TEACHER}->{STUDENT} paralog map",
                urls=[GNOMAD_API, VEP_REGION], tool="scripts/nsd_pair_calibration.py",
                params={"am_cutoff": cutoff, "am_generic": generic})
    rec = core.save_tsv("processed", f"{STUDENT}_candidates_{core.today()}.tsv", rows,
                        base=ROOT / "data", **prov)
    log(f"  wrote {rec['path']} ({len(rows)} rows)")

    conf = [r for r in rows if r["confidently_observed"]]
    tiers = {t: {"all": sum(1 for r in rows if r["tier"] == t),
                 "confident": sum(1 for r in conf if r["tier"] == t)} for t in ("A", "B", "C", "D", "-")}
    report = {
        "am_cutoff_calibrated": cutoff, "am_cutoff_generic": generic,
        "n_missense": len(rows), "n_confidently_observed": len(conf),
        "n_over_calibrated": sum(1 for r in rows if r["over_calibrated_cutoff"]),
        "n_over_generic": sum(1 for r in rows if r["over_generic_cutoff"]),
        "n_over_calibrated_confident": sum(1 for r in conf if r["over_calibrated_cutoff"]),
        "n_over_generic_confident": sum(1 for r in conf if r["over_generic_cutoff"]),
        "tiers": tiers,
        "tier_A_examples": [r["protein_change"] + f" (AM {r['am_pathogenicity']}, {TEACHER} "
                            f"{r['teacher_pathogenic_here']})" for r in
                            sorted((r for r in rows if r["tier"] == "A"),
                                   key=lambda r: -(r["am_pathogenicity"] or 0))[:12]],
        "tier_A_by_domain": dict(Counter(r["domain"] for r in rows if r["tier"] == "A").most_common()),
        "n_identical_residue": sum(1 for r in rows if r["residue_identical"]),
        "n_in_shared_module": sum(1 for r in rows if r["in_shared_module"]),
        "saved": rec["path"],
    }
    log(f"  tiers (all / confidently observed): " + ", ".join(
        f"{t} {v['all']} ({v['confident']})" for t, v in tiers.items()))
    log(f"  generic threshold {generic} would call {report['n_over_generic']} variants, "
        f"calibrated {cutoff} calls {report['n_over_calibrated']}")
    save_json(args, f"{STUDENT}_transfer_summary_{core.today()}.json", report,
              source="derived: transfer summary", urls=[])
    return report


# =========================================================================== stage: lift


def hotspot_window(positions: list[int], length: int) -> tuple[int, int]:
    """The narrowest window holding half of the unique pathogenic residues, exactly as
    scripts/paralog_pairs_scan.py::hotspot defines it, so the measured and predicted lift are
    the same quantity."""
    uniq = sorted(set(positions))
    half = max(1, (len(uniq) + 1) // 2)
    best = min(((uniq[i + half - 1] - uniq[i], uniq[i], uniq[i + half - 1])
                for i in range(len(uniq) - half + 1)), key=lambda t: t[0])
    return best[1], best[2]


def stage_lift(args: argparse.Namespace) -> dict[str, Any]:
    """The measured counterpart of the screen's `transfer_lift`, computed two ways:
    on the curated positive set, and on the raw ClinVar P/LP residues the screen itself used."""
    seqs = sequences(args)
    src, tgt = seqs[TEACHER], seqs[STUDENT]
    mapping = {p: m["student_pos"] for p, m in paralog_map_table(args).items()}

    def signal(positions: list[int]) -> dict[str, Any]:
        uniq = sorted(set(positions))
        aligned = [p for p in uniq if p in mapping]
        ident = [p for p in aligned if src[p - 1] == tgt[mapping[p] - 1]]
        w0, w1 = hotspot_window(uniq, len(src))
        win = [p for p in range(w0, w1 + 1) if p in mapping]
        bg = sum(1 for p in win if src[p - 1] == tgt[mapping[p] - 1]) / len(win) if win else None
        sig = len(ident) / len(uniq)
        return {"n_unique_residues": len(uniq), "n_aligned": len(aligned),
                "gap_fraction_at_plp_residues": round(1 - len(aligned) / len(uniq), 3),
                "transfer_signal_residues": round(sig, 3),
                "hotspot_window": [w0, w1], "hotspot_width_fraction": round((w1 - w0 + 1) / len(src), 3),
                "background_identity_in_window": round(bg, 3) if bg is not None else None,
                "transfer_lift": round(sig - bg, 3) if bg is not None else None}

    sets = load_tsv(args, f"{TEACHER}_calibration_sets_{core.today()}.tsv")
    curated = [int(r["protein_pos"]) for r in sets
               if r["set"] == "missense_positive_core" and r["protein_pos"]]
    curated_mod = [p for p in curated if in_module(TEACHER, p)]
    raw = load_tsv(args, f"{TEACHER}_clinvar_{core.today()}.tsv")
    screen_like = [int(r["protein_pos"]) for r in raw
                   if r["protein_pos"] and (r["classification"] or "").lower() in P_LP
                   and r["molecular_consequence"] and "missense variant" in r["molecular_consequence"]]
    neg = [r for r in sets if r["set"] == "missense_negative_module_population"]
    out = {
        "curated_positives": signal(curated),
        "curated_positives_inside_module": signal(curated_mod),
        "screen_definition_raw_clinvar_plp": signal(screen_like),
        "shared_module_background_identity": None,
        "population_control_identity_in_module": {
            "n": len(neg),
            "identical_fraction": round(sum(1 for r in neg if r["residue_identical_in_student"] == "True")
                                        / len(neg), 3) if neg else None},
    }
    lo, hi = SHARED_MODULE[TEACHER]
    win = [p for p in range(lo, hi + 1) if p in mapping]
    out["shared_module_background_identity"] = round(
        sum(1 for p in win if src[p - 1] == tgt[mapping[p] - 1]) / len(win), 3)
    out["lift_vs_module_background"] = round(
        out["curated_positives_inside_module"]["transfer_signal_residues"]
        - out["shared_module_background_identity"], 3)
    out["lift_vs_population_controls"] = round(
        out["curated_positives_inside_module"]["transfer_signal_residues"]
        - out["population_control_identity_in_module"]["identical_fraction"], 3)
    for k, v in out.items():
        log(f"  {k}: {v}")
    save_json(args, f"{TEACHER}_{STUDENT}_measured_lift_{core.today()}.json", out,
              source="derived: paralog map + curated positives, screen-identical lift definition", urls=[])
    return out


# =========================================================================== entry point

STAGES = {
    "identity": ("async", lambda a: stage_identity(a)),
    "domains": ("async", lambda a: stage_domains(a)),
    "universe": ("async", lambda a: stage_universe(a)),
    "map": ("sync", lambda a: stage_map(a)),
    "curate": ("async", lambda a: stage_curate(a)),
    "score": ("async", lambda a: stage_score(a)),
    "measure": ("sync", lambda a: stage_measure(a)),
    "constraint": ("async", lambda a: stage_constraint(a)),
    "lift": ("sync", lambda a: stage_lift(a)),
    "transfer": ("async", lambda a: stage_transfer(a)),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", default="all",
                   help="comma-separated: identity, domains, universe, map, curate, score, measure, "
                        "constraint, lift, transfer, or all")
    p.add_argument("--out", default=None, help="output directory (default data/interim/nsd_pair)")
    p.add_argument("--min-ac", type=int, default=POPULATION_MIN_AC, help="allele count floor for the main negative set")
    p.add_argument("--am-cutoff", type=float, default=None, help="override the calibrated AlphaMissense threshold")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    order = list(STAGES) if args.stage == "all" else [s.strip() for s in args.stage.split(",")]
    out_dir(args).mkdir(parents=True, exist_ok=True)
    log(f"# nsd_pair_calibration {TEACHER} -> {STUDENT}  stages={order}  out={core.rel(out_dir(args))}  {core.today()}")
    results = {}
    for name in order:
        if name not in STAGES:
            raise SystemExit(f"unknown stage {name!r}; known: {', '.join(STAGES)}")
        kind, fn = STAGES[name]
        log(f"\n== stage {name}")
        results[name] = asyncio.run(fn(args)) if kind == "async" else fn(args)
    print(json.dumps(results, indent=1, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
