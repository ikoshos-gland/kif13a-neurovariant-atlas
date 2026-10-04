"""Variant-level sources: gnomAD v4.1, ClinVar, Ensembl VEP (+plugins), AlphaMissense, CADD.

Numbering rule for the whole project: protein positions are MANE Select.
gnomAD's gene-level query reports HGVSp on the most severe transcript (not MANE),
so variants are pulled with the transcript-level query on the MANE transcript.
"""

from __future__ import annotations

import csv
import io
import re
from collections import Counter
from typing import Any

from . import core

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
VEP = "https://rest.ensembl.org/vep/human"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
CADD_API = "https://cadd.gs.washington.edu/api/v1.0"

HGVSP_RE = re.compile(r"p\.\(?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|=|\*|Ter|fs|del|dup|ins)?")
CODING = {
    "missense_variant", "synonymous_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "inframe_insertion", "inframe_deletion", "splice_donor_variant",
    "splice_acceptor_variant", "splice_region_variant", "protein_altering_variant",
}


def _protein_change(hgvsp: str | None) -> tuple[str | None, int | None, str | None]:
    if not hgvsp:
        return None, None, None
    m = HGVSP_RE.search(hgvsp)
    if not m:
        return None, None, None
    ref3, pos, alt3 = m.groups()
    alt = core.AA3.get(alt3, alt3) if alt3 else None
    return core.AA3.get(ref3, ref3), int(pos), alt


async def _gnomad(query: str) -> dict[str, Any]:
    data = await core.post_json(GNOMAD_API, {"query": query}, timeout=300)
    if data.get("errors"):
        raise RuntimeError(f"gnomAD GraphQL error: {data['errors'][:2]}")
    return data["data"]


# --------------------------------------------------------------------------- gnomAD


async def gnomad_constraint(gene: str, save: bool = True) -> dict[str, Any]:
    t = core.target(gene)
    query = """
    { gene(gene_id: "%s", reference_genome: GRCh38) {
        gene_id symbol canonical_transcript_id
        mane_select_transcript { ensembl_id ensembl_version refseq_id refseq_version }
        gnomad_constraint { exp_lof obs_lof oe_lof oe_lof_lower oe_lof_upper oe_lof_percentile
          exp_mis obs_mis oe_mis oe_mis_lower oe_mis_upper exp_syn obs_syn oe_syn syn_z
          lof_z mis_z pli flags }
    } }""" % t["ensembl_gene"]
    data = (await _gnomad(query))["gene"]
    c = data["gnomad_constraint"] or {}
    out: dict[str, Any] = {
        "gene": t["symbol"],
        "dataset": "gnomAD v4.1 constraint",
        "mane_select_transcript": data.get("mane_select_transcript"),
        "pLI": c.get("pli"),
        "LOEUF (oe_lof_upper)": c.get("oe_lof_upper"),
        "oe_lof": c.get("oe_lof"),
        "missense_z": c.get("mis_z"),
        "oe_mis": c.get("oe_mis"),
        "oe_mis_upper": c.get("oe_mis_upper"),
        "obs_mis/exp_mis": [c.get("obs_mis"), c.get("exp_mis")],
        "synonymous_z": c.get("syn_z"),
        "flags": c.get("flags"),
        "interpretation_hint": "LOEUF < 0.6 is the usual loss-of-function intolerance cutoff.",
    }
    if save:
        out["saved"] = core.save_json(
            "gnomad", f"{t['symbol']}_constraint_gnomad_r4_{core.today()}.json", data,
            source="gnomAD GraphQL", urls=[GNOMAD_API], tool="gnomad_constraint", params={"gene": gene},
        )
    return out


async def gnomad_variants(gene: str, consequence: str = "missense", save: bool = True) -> dict[str, Any]:
    t = core.target(gene)
    enst = core.unversioned(t["mane_transcript"])
    query = """
    { transcript(transcript_id: "%s", reference_genome: GRCh38) {
        transcript_id transcript_version
        variants(dataset: gnomad_r4) {
          variant_id chrom pos ref alt rsids hgvsc hgvsp consequence lof lof_filter flags
          exome { ac an af homozygote_count filters }
          genome { ac an af homozygote_count filters }
          joint { ac an homozygote_count filters }
          in_silico_predictors { id value }
        } } }""" % enst
    tx = (await _gnomad(query))["transcript"]
    rows: list[dict[str, Any]] = []
    for v in tx["variants"]:
        joint = v.get("joint") or {}
        ref_aa, prot_pos, alt_aa = _protein_change(v.get("hgvsp"))
        preds = {p["id"]: p["value"] for p in (v.get("in_silico_predictors") or [])}
        rows.append({
            "variant_id": v["variant_id"], "chrom": v["chrom"], "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
            "rsids": v.get("rsids"), "transcript": f"{tx['transcript_id']}.{tx['transcript_version']}",
            "consequence": v.get("consequence"), "hgvsc": v.get("hgvsc"), "hgvsp": v.get("hgvsp"),
            "protein_pos": prot_pos, "aa_ref": ref_aa, "aa_alt": alt_aa, "lof": v.get("lof"),
            "joint_ac": joint.get("ac"), "joint_an": joint.get("an"),
            "joint_af": (joint["ac"] / joint["an"]) if joint.get("an") else None,
            "joint_hom": joint.get("homozygote_count"), "joint_filters": joint.get("filters"),
            "exome_af": (v.get("exome") or {}).get("af"), "genome_af": (v.get("genome") or {}).get("af"),
            "flags": v.get("flags"), "cadd_phred": preds.get("cadd"), "revel": preds.get("revel_max"),
            "spliceai_ds_max": preds.get("spliceai_ds_max"), "pangolin_largest_ds": preds.get("pangolin_largest_ds"),
            "phylop": preds.get("phylop"),
        })

    if consequence == "missense":
        keep = [r for r in rows if r["consequence"] == "missense_variant"]
    elif consequence == "coding":
        keep = [r for r in rows if r["consequence"] in CODING]
    else:
        keep = rows

    afs = [r["joint_af"] for r in keep if r["joint_af"] is not None]
    out: dict[str, Any] = {
        "gene": t["symbol"],
        "transcript": f"{tx['transcript_id']}.{tx['transcript_version']} (MANE Select)",
        "dataset": "gnomad_r4 (v4.1, exomes+genomes joint)",
        "n_variants_on_transcript": len(rows),
        "consequence_filter": consequence,
        "n_selected": len(keep),
        "consequence_counts": dict(Counter(r["consequence"] for r in rows).most_common(12)),
        "selected_af_bins": {
            "singleton_ac1": sum(1 for r in keep if r["joint_ac"] == 1),
            "af_lt_1e-4": sum(1 for a in afs if a < 1e-4),
            "af_1e-4_to_1e-2": sum(1 for a in afs if 1e-4 <= a < 1e-2),
            "af_ge_1e-2": sum(1 for a in afs if a >= 1e-2),
        },
        "distinct_protein_positions": len({r["protein_pos"] for r in keep if r["protein_pos"]}),
        "top_by_af": [
            {k: r[k] for k in ("variant_id", "hgvsp", "joint_af", "joint_hom")}
            for r in sorted(keep, key=lambda r: r["joint_af"] or 0, reverse=True)[:8]
        ],
    }
    if save:
        stem = f"{t['symbol']}_{enst}_gnomad_r4_{consequence}_{core.today()}"
        out["saved"] = core.save_tsv(
            "gnomad", stem + ".tsv", keep, source="gnomAD GraphQL (transcript query)", urls=[GNOMAD_API],
            tool="gnomad_variants", params={"gene": gene, "transcript": enst, "consequence": consequence},
        )
    return out


# --------------------------------------------------------------------------- ClinVar

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
}
TITLE_RE = re.compile(r"^(NM_\d+)\.\d+\(([^)]+)\):(c\.[^ ]+)(?: \((p\.[^)]+)\))?")


async def clinvar_variants(gene: str, save: bool = True) -> dict[str, Any]:
    t = core.target(gene)
    key = core.secret("NCBI_API_KEY")
    base = {"db": "clinvar", "retmode": "json", **({"api_key": key} if key else {})}
    term = f"{t['symbol']}[gene] AND single_gene[prop]"
    search = await core.get_json(f"{EUTILS}/esearch.fcgi", params={**base, "term": term, "usehistory": "y", "retmax": 0})
    res = search["esearchresult"]
    total = int(res["count"])
    mane_nm = core.unversioned(t["refseq_transcript"])

    rows: list[dict[str, Any]] = []
    for start in range(0, total, 400):
        page = await core.get_json(
            f"{EUTILS}/esummary.fcgi",
            params={**base, "query_key": res["querykey"], "WebEnv": res["webenv"], "retstart": start, "retmax": 400},
            timeout=300,
        )
        result = page["result"]
        for uid in result.get("uids", []):
            r = result[uid]
            germ = r.get("germline_classification") or {}
            vs = (r.get("variation_set") or [{}])[0]
            loc = next((l for l in vs.get("variation_loc", []) if l.get("assembly_name") == "GRCh38"), {})
            m = TITLE_RE.match(r.get("title", ""))
            title_nm, cdna, prot = (m.group(1), m.group(3), m.group(4)) if m else (None, None, None)
            ref_aa, prot_pos, alt_aa = _protein_change(prot)
            spdi = vs.get("canonical_spdi") or ""
            spdi_parts = spdi.split(":")
            review = germ.get("review_status", "")
            rows.append({
                "uid": uid, "accession": r.get("accession"), "title": r.get("title"),
                "variant_type": vs.get("variant_type"),
                "classification": germ.get("description"), "review_status": review,
                "stars": STARS.get(review, 0), "last_evaluated": (germ.get("last_evaluated") or "")[:10],
                "traits": [tr.get("trait_name") for tr in germ.get("trait_set", [])],
                "molecular_consequence": r.get("molecular_consequence_list"),
                "title_transcript": title_nm, "title_is_mane": title_nm == mane_nm,
                "hgvsc": cdna, "hgvsp": prot, "protein_pos": prot_pos if title_nm == mane_nm else None,
                "aa_ref": ref_aa, "aa_alt": alt_aa,
                "chrom": loc.get("chr"), "start": loc.get("start"), "stop": loc.get("stop"),
                "spdi": spdi, "ref": spdi_parts[2] if len(spdi_parts) == 4 else None,
                "alt": spdi_parts[3] if len(spdi_parts) == 4 else None,
            })

    missense = [r for r in rows if r["molecular_consequence"] and "missense variant" in r["molecular_consequence"]]
    plp = [r for r in missense if (r["classification"] or "").lower() in
           {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}]
    out: dict[str, Any] = {
        "gene": t["symbol"], "query": term, "n_records": total, "n_parsed": len(rows),
        "mane_refseq": t["refseq_transcript"],
        "classification_counts": dict(Counter(r["classification"] for r in rows).most_common(10)),
        "n_missense": len(missense), "n_missense_P_LP": len(plp),
        "n_missense_P_LP_with_mane_position": sum(1 for r in plp if r["protein_pos"]),
        "records_not_on_mane_title": sum(1 for r in rows if r["title_transcript"] and not r["title_is_mane"]),
        "example_P_LP": [{k: r[k] for k in ("accession", "hgvsp", "classification", "stars")} for r in plp[:8]],
    }
    if save:
        out["saved"] = core.save_tsv(
            "clinvar", f"{t['symbol']}_clinvar_{core.today()}.tsv", rows, source="NCBI E-utilities ClinVar esummary",
            urls=[f"{EUTILS}/esearch.fcgi", f"{EUTILS}/esummary.fcgi"], tool="clinvar_variants", params={"term": term},
        )
    return out


# --------------------------------------------------------------------------- VEP

VEP_PARAMS = {
    "AlphaMissense": 1, "CADD": 1, "SpliceAI": 1, "REVEL": 1,
    "dbNSFP": "Ensembl_transcriptid,ESM1b_score,PrimateAI_score,MetaRNN_score",
    "mane": 1, "canonical": 1, "hgvs": 1, "protein": 1, "numbers": 1, "variant_class": 1,
    "af_gnomade": 1, "af_gnomadg": 1,
}
REGION_RE = re.compile(r"^(?:chr)?([0-9XYM]{1,2})[-:_ ](\d+)[-:_ ]([ACGTN]+)[-:_>/ ]([ACGTN]+)$", re.I)


def vcf_sort_key(line: str) -> tuple[int, int]:
    chrom, pos = line.split()[:2]
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (int(chrom) if chrom.isdigit() else order.get(chrom.upper(), 26), int(pos))


def _dbnsfp_for(tc: dict[str, Any], field: str, enst: str) -> str | None:
    tids = str(tc.get("ensembl_transcriptid", "")).split(",")
    values = str(tc.get(field.lower(), "")).split(",")
    if enst in tids and len(values) == len(tids):
        value = values[tids.index(enst)]
        return None if value in (".", "") else value
    return None


def _pick_transcript(consequences: list[dict[str, Any]]) -> dict[str, Any] | None:
    for key in ("mane_select", "canonical"):
        hits = [tc for tc in consequences if tc.get(key)]
        if hits:
            return hits[0]
    return consequences[0] if consequences else None


def _flatten_vep(v: dict[str, Any]) -> dict[str, Any]:
    tc = _pick_transcript(v.get("transcript_consequences", [])) or {}
    enst = tc.get("transcript_id", "")
    am = tc.get("alphamissense") or {}
    splice = tc.get("spliceai") or {}
    ds = [splice.get(k) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL") if isinstance(splice.get(k), (int, float))]
    alt = (v.get("allele_string") or "/").split("/")[-1]
    gnomade = gnomadg = None
    clin_sig: list[str] = []
    for cv in v.get("colocated_variants", []):
        freqs = (cv.get("frequencies") or {}).get(alt) or {}
        gnomade = freqs.get("gnomade", gnomade)
        gnomadg = freqs.get("gnomadg", gnomadg)
        clin_sig += cv.get("clin_sig") or []
    ref_aa, prot_pos, alt_aa = _protein_change(tc.get("hgvsp"))
    return {
        "input": v.get("input"), "location": f"{v.get('seq_region_name')}:{v.get('start')}",
        "allele_string": v.get("allele_string"), "most_severe_consequence": v.get("most_severe_consequence"),
        "gene": tc.get("gene_symbol"), "transcript": enst, "mane_select": tc.get("mane_select"),
        "consequence_terms": tc.get("consequence_terms"), "hgvsc": tc.get("hgvsc"), "hgvsp": tc.get("hgvsp"),
        "protein_pos": prot_pos, "aa_ref": ref_aa, "aa_alt": alt_aa, "exon": tc.get("exon"),
        "am_pathogenicity": am.get("am_pathogenicity"), "am_class": am.get("am_class"),
        "cadd_phred": tc.get("cadd_phred"), "revel": tc.get("revel"),
        "spliceai_ds_max": max(ds) if ds else None,
        "esm1b": _dbnsfp_for(tc, "ESM1b_score", enst), "primateai": tc.get("primateai_score"),
        "metarnn": _dbnsfp_for(tc, "MetaRNN_score", enst),
        "sift": tc.get("sift_prediction"), "polyphen": tc.get("polyphen_prediction"),
        "gnomad_exomes_af": gnomade, "gnomad_genomes_af": gnomadg,
        "clinvar_clin_sig": sorted(set(clin_sig)) or None,
    }


async def vep_annotate(variants: list[str], label: str = "batch", save: bool = True) -> dict[str, Any]:
    hgvs, region, ids, bad = [], [], [], []
    for raw in variants:
        s = raw.strip()
        if not s:
            continue
        if re.match(r"^rs\d+$", s, re.I):
            ids.append(s)
        elif (m := REGION_RE.match(s)):
            chrom, pos, ref, alt = m.groups()
            region.append(f"{chrom} {pos} . {ref.upper()} {alt.upper()} . . .")
        elif ":" in s and re.search(r":[cgnpm]\.", s):
            hgvs.append(s)
        else:
            bad.append(s)

    # VEP rejects VCF-style region input that is not sorted by chromosome and position ("appears to be unsorted")
    region.sort(key=vcf_sort_key)
    raw_results: list[dict[str, Any]] = []
    for items, endpoint, key in ((hgvs, "hgvs", "hgvs_notations"), (region, "region", "variants"), (ids, "id", "ids")):
        for batch in core.chunks(items, 200):
            raw_results += await core.post_json(
                f"{VEP}/{endpoint}", {key: batch}, params=VEP_PARAMS,
                headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=300,
            )
    rows = [_flatten_vep(v) for v in raw_results]
    out: dict[str, Any] = {
        "n_input": len(variants), "n_annotated": len(rows), "unparsed_inputs": bad,
        "scores": "AlphaMissense, CADD v1.7, SpliceAI (max DS), REVEL, dbNSFP ESM1b/PrimateAI/MetaRNN (MANE-indexed), gnomAD AF",
        "rows": rows if len(rows) <= 12 else rows[:12],
    }
    if len(rows) > 12:
        out["note"] = f"showing 12 of {len(rows)} rows; full table saved"
    if save and rows:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", label)[:60]
        prov = dict(source="Ensembl VEP REST", urls=[VEP], tool="vep_annotate", params=VEP_PARAMS)
        out["saved"] = core.save_tsv("vep", f"vep_{safe}_{core.today()}.tsv", rows, **prov)
        out["saved_raw"] = core.save_json("vep", f"vep_{safe}_{core.today()}.raw.json", raw_results, compress=True, **prov)
    return out


# --------------------------------------------------------------------------- AlphaMissense


async def alphamissense(
    uniprot: str, positions: str | None = None, variants: list[str] | None = None, save: bool = True
) -> dict[str, Any]:
    entries = await core.get_json(f"{AFDB_API}/{uniprot}")
    entry = entries[0]
    url = entry.get("amAnnotationsUrl")
    if not url:
        raise ValueError(
            f"AlphaFold DB has no AlphaMissense table for {uniprot} (AlphaMissense covers UniProt canonical "
            "isoforms only). For a non-canonical MANE isoform such as KIF1A Q12756-3, use vep_annotate, "
            "whose AlphaMissense plugin is transcript-aware."
        )
    content = await core.get_bytes(url, timeout=300)
    table = list(csv.DictReader(io.StringIO(content.decode("utf-8"))))

    by_pos: dict[int, list[float]] = {}
    lpath_by_pos: Counter[int] = Counter()
    lookup: dict[str, dict[str, Any]] = {}
    for row in table:
        pv = row["protein_variant"]
        pos = int(pv[1:-1])
        score = float(row["am_pathogenicity"])
        by_pos.setdefault(pos, []).append(score)
        if row["am_class"] == "LPath":
            lpath_by_pos[pos] += 1
        lookup[pv] = {"variant": pv, "am_pathogenicity": score, "am_class": row["am_class"]}

    wanted = core.parse_positions(positions)
    ref_at = {int(k[1:-1]): k[0] for k in lookup}
    per_pos = [
        {"pos": p, "ref": ref_at.get(p), "mean_am": round(sum(s) / len(s), 4), "n_lpath": lpath_by_pos[p]}
        for p, s in sorted(by_pos.items())
        if wanted is None or p in wanted
    ]

    out: dict[str, Any] = {
        "uniprot": uniprot, "source": url, "protein_length": entry.get("uniprotEnd"),
        "n_substitutions": len(table),
        "class_counts": dict(Counter(r["am_class"] for r in table)),
    }
    if variants:
        out["variant_scores"] = [lookup.get(v.strip(), {"variant": v, "error": "not in table"}) for v in variants]
    if wanted is not None:
        out["positions"] = per_pos[:120]
    else:
        out["most_sensitive_positions"] = sorted(per_pos, key=lambda r: r["mean_am"], reverse=True)[:15]
    if save:
        out["saved"] = core.save_bytes(
            "alphamissense", f"AF-{uniprot}-F1-aa-substitutions_{core.today()}.csv", content,
            source="AlphaFold DB AlphaMissense per-protein table", urls=[url], tool="alphamissense",
            params={"uniprot": uniprot}, n_records=len(table),
        )
    return out


# --------------------------------------------------------------------------- CADD


async def cadd_region(chrom: str, start: int, end: int, version: str = "GRCh38-v1.7", save: bool = True) -> dict[str, Any]:
    chrom = chrom.removeprefix("chr")
    if end < start:
        start, end = end, start
    if end - start > 2000:
        raise ValueError("Keep CADD region queries to 2 kb or less; loop over exons for larger spans.")
    rows: list[dict[str, Any]] = []
    for s in range(start, end + 1, 100):
        e = min(s + 99, end)
        data = await core.get_json(f"{CADD_API}/{version}/{chrom}:{s}-{e}", timeout=180)
        if data and isinstance(data[0], list):  # range queries: header row + value rows
            header, *values = data
            rows += [dict(zip(header, v)) for v in values]
        else:  # single-position queries: list of dicts
            rows += data
    out: dict[str, Any] = {
        "region": f"{chrom}:{start}-{end}", "version": version, "n_scores": len(rows),
        "max_phred": max((float(r["PHRED"]) for r in rows), default=None),
        "n_phred_ge_20": sum(1 for r in rows if float(r["PHRED"]) >= 20),
        "rows": rows[:10],
    }
    if save and rows:
        out["saved"] = core.save_tsv(
            "cadd", f"cadd_{version}_{chrom}_{start}_{end}_{core.today()}.tsv", rows, source="CADD API",
            urls=[f"{CADD_API}/{version}"], tool="cadd_region", params={"chrom": chrom, "start": start, "end": end},
        )
    return out
