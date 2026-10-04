"""Protein and structure sources: AlphaFold DB files, UniProt features, PDB coverage,
and a quick MANE-to-MANE paralog residue map."""

from __future__ import annotations

from collections import defaultdict
from typing import Any

# Imported at module load on purpose: on Windows, importing numpy (via Biopython) for the first
# time inside a tool call deadlocks against the stdio reader thread of the MCP server.
from Bio import Align
from Bio.Align import substitution_matrices

from . import core

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
PDBE_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
RCSB_FILES = "https://files.rcsb.org/download"
ENSEMBL_SEQ = "https://rest.ensembl.org/sequence/id"


# --------------------------------------------------------------------------- AlphaFold DB


def _confident_segments(residues: list[int], scores: list[float], cutoff: float = 70.0, min_len: int = 20) -> list[str]:
    segments, start, prev = [], None, None
    for r, s in zip(residues, scores):
        if s >= cutoff and start is None:
            start = r
        if s < cutoff and start is not None:
            if prev - start + 1 >= min_len:
                segments.append(f"{start}-{prev}")
            start = None
        prev = r
    if start is not None and prev - start + 1 >= min_len:
        segments.append(f"{start}-{prev}")
    return segments


async def alphafold_model(uniprot: str, include_pae: bool = True, save: bool = True) -> dict[str, Any]:
    entry = (await core.get_json(f"{AFDB_API}/{uniprot}"))[0]
    version = entry.get("latestVersion")
    folder = f"alphafold/{entry['entryId']}_v{version}"
    conf = await core.get_json(entry["plddtDocUrl"])
    residues, scores = conf["residueNumber"], conf["confidenceScore"]
    out: dict[str, Any] = {
        "uniprot": uniprot, "entry": entry["entryId"], "model_version": version,
        "model_created": entry.get("modelCreatedDate"), "length": entry.get("uniprotEnd"),
        "mean_plddt": round(sum(scores) / len(scores), 2),
        "fraction_very_high_gt90": entry.get("fractionPlddtVeryHigh"),
        "fraction_confident_70_90": entry.get("fractionPlddtConfident"),
        "fraction_low_50_70": entry.get("fractionPlddtLow"),
        "fraction_very_low_lt50": entry.get("fractionPlddtVeryLow"),
        "confident_segments_plddt_ge70_len_ge20": _confident_segments(residues, scores),
        "has_alphamissense": bool(entry.get("amAnnotationsUrl")),
    }
    if save:
        prov = dict(source="AlphaFold DB", tool="alphafold_model", params={"uniprot": uniprot})
        saved = [
            core.save_bytes(folder, entry["cifUrl"].rsplit("/", 1)[-1], await core.get_bytes(entry["cifUrl"], timeout=300),
                            urls=[entry["cifUrl"]], **prov),
            core.save_json(folder, entry["plddtDocUrl"].rsplit("/", 1)[-1], conf, urls=[entry["plddtDocUrl"]], **prov),
        ]
        if include_pae:
            pae = await core.get_bytes(entry["paeDocUrl"], timeout=600)
            saved.append(core.save_bytes(folder, entry["paeDocUrl"].rsplit("/", 1)[-1], pae, urls=[entry["paeDocUrl"]], **prov))
        saved.append(core.save_json(folder, "api_entry.json", entry, urls=[f"{AFDB_API}/{uniprot}"], **prov))
        out["saved"] = saved
    return out


# --------------------------------------------------------------------------- UniProt

FEATURE_TYPES = [
    "Domain", "Region", "Coiled coil", "Motif", "Compositional bias", "Binding site", "Site",
    "Modified residue", "Natural variant", "Mutagenesis", "Alternative sequence", "Sequence conflict",
]


async def uniprot_features(accession: str, save: bool = True) -> dict[str, Any]:
    data = await core.get_json(f"{UNIPROT}/{accession}.json")
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for f in data.get("features", []):
        loc = f.get("location", {})
        grouped[f["type"]].append({
            "start": loc.get("start", {}).get("value"), "end": loc.get("end", {}).get("value"),
            "description": f.get("description") or (f.get("ligand") or {}).get("name"),
            **({"change": f"{f['alternativeSequence'].get('originalSequence')}>{','.join(f['alternativeSequence'].get('alternativeSequences', []))}"}
               if f.get("alternativeSequence") else {}),
        })
    verbose = {"Domain", "Region", "Coiled coil", "Motif", "Binding site", "Site", "Compositional bias"}
    out: dict[str, Any] = {
        "accession": data.get("primaryAccession"), "entry_name": data.get("uniProtkbId"),
        "protein": ((data.get("proteinDescription") or {}).get("recommendedName") or {}).get("fullName", {}).get("value"),
        "length": (data.get("sequence") or {}).get("length"),
        "sequence_version": (data.get("entryAudit") or {}).get("sequenceVersion"),
        "entry_version": (data.get("entryAudit") or {}).get("entryVersion"),
        "features": {t: grouped[t] for t in FEATURE_TYPES if t in grouped and t in verbose},
        "feature_counts": {t: len(grouped[t]) for t in FEATURE_TYPES if t in grouped},
        "mutagenesis": grouped.get("Mutagenesis", [])[:25],
    }
    if save:
        out["saved"] = core.save_json(
            "uniprot", f"{accession}_{core.today()}.json", data, source="UniProt REST",
            urls=[f"{UNIPROT}/{accession}.json"], tool="uniprot_features", params={"accession": accession},
        )
    return out


# --------------------------------------------------------------------------- PDB


async def pdb_structures(uniprot: str, download: list[str] | None = None, save: bool = True) -> dict[str, Any]:
    try:
        best = (await core.get_json(f"{PDBE_BEST}/{uniprot}")).get(uniprot, [])
    except Exception as exc:  # PDBe answers 404 when there is no mapped structure
        if "404" in str(exc):
            best = []
        else:
            raise
    ids = sorted({b["pdb_id"].upper() for b in best})
    meta: dict[str, dict[str, Any]] = {}
    if ids:
        query = """query($ids:[String!]!){ entries(entry_ids:$ids){ rcsb_id struct{title} exptl{method}
                   rcsb_accession_info{initial_release_date} rcsb_entry_info{resolution_combined} } }"""
        data = await core.post_json(RCSB_GRAPHQL, {"query": query, "variables": {"ids": ids}})
        for e in data["data"]["entries"]:
            meta[e["rcsb_id"]] = e
    rows = []
    for b in best:
        m = meta.get(b["pdb_id"].upper(), {})
        rows.append({
            "pdb_id": b["pdb_id"].upper(), "chain": b["chain_id"], "unp_start": b["unp_start"], "unp_end": b["unp_end"],
            "coverage": b.get("coverage"), "resolution": b.get("resolution"), "method": b.get("experimental_method"),
            "released": ((m.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10],
            "title": (m.get("struct") or {}).get("title"),
        })
    unique = {}
    for r in rows:
        unique.setdefault(r["pdb_id"], r)
    out: dict[str, Any] = {
        "uniprot": uniprot, "n_chains": len(rows), "n_entries": len(unique),
        "entries": sorted(unique.values(), key=lambda r: (r["resolution"] or 99))[:25],
    }
    if not rows:
        out["note"] = "No experimental structure mapped to this UniProt accession; use alphafold_model."
    if save and rows:
        out["saved"] = core.save_tsv(
            "pdb", f"{uniprot}_pdb_coverage_{core.today()}.tsv", rows, source="PDBe SIFTS best_structures + RCSB",
            urls=[f"{PDBE_BEST}/{uniprot}", RCSB_GRAPHQL], tool="pdb_structures", params={"uniprot": uniprot},
        )
    if download:
        out["downloaded"] = []
        for pdb_id in download:
            url = f"{RCSB_FILES}/{pdb_id.upper()}.cif"
            out["downloaded"].append(core.save_bytes(
                "pdb", f"{pdb_id.upper()}.cif", await core.get_bytes(url, timeout=600),
                source="RCSB PDB", urls=[url], tool="pdb_structures", params={"pdb_id": pdb_id},
            ))
    return out


# --------------------------------------------------------------------------- paralog map


async def mane_protein(gene: str) -> tuple[dict[str, Any], str]:
    t = core.target(gene)
    enst = core.unversioned(t["mane_transcript"])
    cache = core.RAW / "sequences" / f"{t['symbol']}_{enst}.fasta"
    if cache.exists():
        seq = "".join(cache.read_text(encoding="utf-8").splitlines()[1:])
    else:
        url = f"{ENSEMBL_SEQ}/{enst}"
        seq = (await core.request("GET", url, params={"type": "protein"}, headers={"Accept": "text/plain"})).text.strip()
        core.save_bytes("sequences", cache.name, f">{t['symbol']}|{t['mane_transcript']}|MANE\n{seq}\n".encode(),
                        source="Ensembl REST", urls=[url], tool="paralog_map", params={"transcript": enst})
    if len(seq) != t["protein_length"]:
        raise RuntimeError(f"{gene}: MANE protein is {len(seq)} aa but targets.yaml says {t['protein_length']}")
    return t, seq


async def paralog_map(source_gene: str, target_gene: str, positions: str, save: bool = True) -> dict[str, Any]:
    src, src_seq = await mane_protein(source_gene)
    tgt, tgt_seq = await mane_protein(target_gene)
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score, aligner.extend_gap_score = -10.0, -0.5
    aligner.end_gap_score = 0.0
    aln = aligner.align(src_seq, tgt_seq)[0]

    mapping: dict[int, int] = {}
    for (s0, s1), (t0, _t1) in zip(*aln.aligned.tolist()):
        for i in range(s1 - s0):
            mapping[s0 + i + 1] = t0 + i + 1
    identical = sum(1 for s, t in mapping.items() if src_seq[s - 1] == tgt_seq[t - 1])

    def window_identity(pos: int, half: int = 7) -> float | None:
        pairs = [(p, mapping[p]) for p in range(pos - half, pos + half + 1) if p in mapping]
        if not pairs:
            return None
        return round(sum(src_seq[s - 1] == tgt_seq[t - 1] for s, t in pairs) / (2 * half + 1), 2)

    rows = []
    for pos in sorted(core.parse_positions(positions) or []):
        if not 1 <= pos <= len(src_seq):
            continue
        tpos = mapping.get(pos)
        rows.append({
            f"{src['symbol']}_pos": pos, f"{src['symbol']}_aa": src_seq[pos - 1],
            f"{tgt['symbol']}_pos": tpos, f"{tgt['symbol']}_aa": tgt_seq[tpos - 1] if tpos else None,
            "identical": bool(tpos) and src_seq[pos - 1] == tgt_seq[tpos - 1],
            "window15_identity": window_identity(pos),
        })
    out: dict[str, Any] = {
        "method": "pairwise global Needleman-Wunsch-Gotoh, BLOSUM62, gap -10/-0.5, free end gaps (preliminary; "
                  "the definitive WP1 map comes from the multi-paralog MSA plus structural alignment)",
        "source": f"{src['symbol']} {src['mane_transcript']} ({len(src_seq)} aa)",
        "target": f"{tgt['symbol']} {tgt['mane_transcript']} ({len(tgt_seq)} aa)",
        "aligned_pairs": len(mapping),
        "identity_over_aligned_pairs": round(identical / len(mapping), 3) if mapping else None,
        "positions": rows,
    }
    if save:
        full = [{f"{src['symbol']}_pos": s, f"{src['symbol']}_aa": src_seq[s - 1], f"{tgt['symbol']}_pos": t,
                 f"{tgt['symbol']}_aa": tgt_seq[t - 1]} for s, t in sorted(mapping.items())]
        out["saved"] = core.save_tsv(
            "paralog_maps", f"{src['symbol']}_to_{tgt['symbol']}_pairwise_MANE.tsv", full, base=core.INTERIM,
            source="derived: Biopython PairwiseAligner on Ensembl MANE proteins", urls=[ENSEMBL_SEQ],
            tool="paralog_map", params={"source": source_gene, "target": target_gene},
        )
    return out
